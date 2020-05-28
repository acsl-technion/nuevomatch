	/*
 * MIT License
 * Copyright (c) 2019 Alon Rashelbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <list>
#include <vector>
#include <set>
#include <iostream>
#include <algorithm> // Set intersection

#include <argument_handler.h>
#include <logging.h>
#include <generic_classifier.h>
#include <nuevomatch.h>
#include <algorithms.h>
#include <rule_db.h>
#include <object_io.h>

using namespace std;

trace_packet gen_packet_in_rule(const list<openflow_rule>& rule_db, uint32_t rule_idx, uint32_t tries);
vector<trace_packet> generate_unique_packets(const list<openflow_rule>& rule_db, uint32_t num, set<uint32_t>& malformed_packets);
vector<trace_packet> generate_random_packets(const list<openflow_rule>& rule_db, uint32_t num, set<uint32_t>& malformed_packets);

// Holds arguments information
static argument_t my_arguments[] = {
		// Name,			Required,	IsBoolean,	Default,	Help
		{"-f",				1,			0,			NULL,		"Filter filename. Supported formats: [Classbench, Classbnech-ng, Binary]"},
		{"--indices",		0,			0,			NULL,		"Load an indices file to apply a filter on the ruleset"},
		{"--reverse",		0,			1,			NULL,		"Reverse the effect of the indices file"},
		{"-n",				1,			0,			"0",		"Number of packets to generate"},
		{"-p",				0,			1,			NULL,		"Prints rule-set to stdout"},
		{"-o",				0,			0,			NULL,		"Output trace filename"},
		{"-s",				0,			1,			NULL,		"Shuffle trace"},
		{"-r",				0,			1,			NULL,		"Randomize trace"},
		{NULL,				0,			0,			NULL,		"Reads a ruleset file / NuevoMatch classifier and generates accurate packet trace file"} /* Sentinel */
};

/**
 * @brief Represents the group of all integer intervals that match a specific rule
 */

class IntervalList {
private:

	list<range> _intervals;

	/**
	 * @brief Initiate empty list
	 */
	IntervalList() {};

public:

	/**
	 * @brief Initiate new interval (with field boundaries)
	 */
	IntervalList(uint32_t low, uint32_t high) {
		_intervals.push_back(range(low, high));
	}

	/**
	 * @brief Subtract a rule region from this
	 * @param low The rule field low value
	 * @param high The rule field high value
	 * @returns The union IntervalList between this and the rule
	 */
	IntervalList apply_rule(uint32_t low, uint32_t high) {

		IntervalList output;
		uint32_t left_cursor = low;
		uint32_t right_cursor = high;
		uint32_t maximum = high;

		for (list<range>::iterator it = _intervals.begin(); it != _intervals.end(); ++it) {
			range& r = *it;

			// Skip intervals in case they do not intersect the rule
			if (r.high < left_cursor) continue;
			if (right_cursor < r.low) break;

			// At this point, the interval r intersect the rule

			// Get the valid lower bound between cursor and the interval
			uint32_t min = std::max(left_cursor, r.low);
			// Get the valid upper bound between the cursor and the interval
			uint32_t max = std::min(right_cursor, r.high);
			if (max > maximum) max = maximum;
			// Add the current interval as valid interval to the output
			output._intervals.push_back(range(min, max));

			// Update additional interval in case required
			if (max < r.high) {
				auto position = it;
				_intervals.insert(++position, range(max+1, r.high));
			}

			// Update the current interval in case required
			r.high = min - 1;
			if (min == 0 || r.high < r.low) {
				it = _intervals.erase(it);
				--it;
			}

			// Update the cursor
			left_cursor = max + 1;
			if (right_cursor > maximum) break;
		}

		return output;
	}

	/**
	 * @brief Returns a valid random value inside this
	 */
	uint32_t random_value() const {
		// TODO
		// In case the rule covers nothing, return 0
		if (_intervals.size() == 0) return 0;
		uint32_t x = gen_uniform_random_uint32(0, _intervals.size()-1);
		// Get any interval within this
		auto it = _intervals.begin();
		for (uint32_t i=0; i<x; ++i, ++it);
		// Get any value within the interval
		const range* intvl = &(*it);
		return gen_uniform_random_uint32(intvl->low, intvl->high);
	}

	/**
	 * @brief Returns the number of intervals in this
	 */
	uint32_t size() const { return _intervals.size(); }

	/**
	 * @brief Used for debugging. Print this
	 */
	void print() const {
		for (auto it : _intervals) {
			fprintf(stderr, "[%u, %u] ", it.low, it.high);
		}
		fprintf(stderr,"\n");
	}
};

// Holds the rule-set attributes
static uint32_t field_num;
static ruleset_type_t ruleset_type;

/**
 * @brief Application entry point
 */
int main(int argc, char** argv) {

	// Parse arguments
	parse_arguments(argc, argv, my_arguments);
	const char* ruleset_filename = ARG("-f")->value;
	const char* indices_filename = ARG("--indices")->value;

	// Rule-set size & rules
	list<openflow_rule> rule_db;

	// Parse rule-set file
	messagef("Reading ruleset file %s...", ruleset_filename);
	ruleset_type = classify_ruleset_file(ruleset_filename);
	if (ruleset_type == UNKNOWN) {
		throw error("Cannot parse rule-set file: file format not recognized");
	} else if (ruleset_type == CLASSBENCH) {
		messagef("Recognized rule-set as Classbench text file");
		rule_db = read_classbench_file(ruleset_filename);
	} else if (ruleset_type == CLASSBENCHNG) {
		messagef("Recognized rule-set as Classbench-ng text file");
		rule_db = read_classbench_ng_file(ruleset_filename);
	} else if (ruleset_type == BINARY) {
		messagef("Recognized rule-set as binary format");
		ObjectReader classifier_file(ruleset_filename);
		rule_db = load_rule_database(classifier_file);
	}

	// Read ruleset indices file
	if (indices_filename != NULL) {
		ObjectReader reader(indices_filename);
		bool reverse = ARG("--reverse")->available;
		std::set<uint32_t> indices = load_indices_database(reader);
		messagef("Rule set had %lu rules before truncating by indices", rule_db.size());
		messagef("Found %u indices in file. Applying on ruleset...", indices.size());
		rule_db = apply_indices_on_ruleset(rule_db, indices, reverse);
	}

	messagef("Rule set has %lu rules", rule_db.size());

	// Check size of rule-set
	if (rule_db.size() == 0) {
		throw error("Rule-set has zero rules");
	}
	field_num = rule_db.front().fields.size();

	// Print rule-table to stderr
	int print_table = ARG("-p")->available;
	if (print_table) {
		print_rule_db(rule_db);
		return 0;
	}

	// Open the output filename for write
	const char* output_filename = ARG("-o")->value;
	FILE* out_file_ptr = fopen(output_filename, "w");
	if (!out_file_ptr) {
		throw error("cannot open output filename for writing");
	}

	// Allocate memory
	uint32_t packets_num = atoi( ARG("-n")->value );
	set<uint32_t> malformed_packets;
	vector<trace_packet> packets;

	// Generate rules
	if (ARG("-r")->available) {
		fprintf(stderr, "Generating random packets...\n");
		packets = generate_random_packets(rule_db, packets_num, malformed_packets);
	} else {
		fprintf(stderr, "Generating unique packets...\n");
		packets = generate_unique_packets(rule_db, packets_num, malformed_packets);
		// Generate additional packets to fill gaps in malformed packets
		fprintf(stderr, "Creating additional %lu random packets\n", malformed_packets.size());
		vector<trace_packet> additional_packets = generate_random_packets(rule_db, malformed_packets.size(), malformed_packets);
		std::move(additional_packets.begin(), additional_packets.end(), std::back_inserter(packets));
	}

	// Initiate checkpoint value
	uint32_t checkpoint = packets.size() < 100 ? 1 : packets.size()/100;

	// Shuffle packets
	bool shuffle = ARG("-s")->available;
	uint32_t* perm = new uint32_t[packets.size()];
	random_permutation(perm, packets.size());

	// Print packets
	for (uint32_t i=0; i<packets.size(); ++i) {

		uint32_t packet_idx = shuffle ? perm[i] : i;
		if (malformed_packets.find(packet_idx) != malformed_packets.end()) {
			continue;
		}

		// Update status
		if (i%checkpoint==0) {
			fprintf(stderr, "\rPrinting rules (%u%%)...", i/checkpoint);
		}

		// Print packet fields
		for (uint32_t j=0; j<field_num; ++j) {
			fprintf(out_file_ptr, "%u\t", packets[packet_idx].header[j]);
		}

		// Print rule number
		fprintf(out_file_ptr, "%u\n", packets[packet_idx].match_priority);
	}

	delete perm;
	fprintf(stderr, "\nDone\n");
	return 0;
}

/**
 * @brief Generate at least a single packet per rule
 * @param rule_db The rule database
 * @param packet_num The number of packets
 * @param A set of invalid packet indices
 * return A list of packets
 */
vector<trace_packet> generate_unique_packets(const list<openflow_rule>& rule_db, uint32_t packets_num, set<uint32_t>& malformed_packets) {

	// Count how many non-unique rules are there
	set<uint32_t> non_unique;

	// The output packets
	vector<trace_packet> output;
	output.resize(packets_num);
	for (uint32_t i=0; i<packets_num; ++i) {
		output[i].header.resize(field_num);
	}

	// Regular upper bounds for 5-tuple classifier
	uint32_t checkpoint, counter;

	// For each dimension...
	for (uint32_t f=0; f<field_num; ++f) {

		// Initiate checkpoint value
		checkpoint = rule_db.size() < 100 ? 1 : rule_db.size()/100;
		counter = 0;

		// Build the IntervalList
		IntervalList interval(0, get_field_bound(f, ruleset_type));
		vector<IntervalList> rule_intervals;
		vector<uint32_t> priorities;

		for(auto rule : rule_db) {
			// Update status
			if (counter%checkpoint==0) {
				fprintf(stderr, "\rBuilding IntervalList for field %u... (%u%%)...", f, counter/checkpoint);
			}
			rule_intervals.push_back( interval.apply_rule(rule.fields[f].low, rule.fields[f].high) );
			priorities.push_back(rule.priority);
			++counter;
		}
		fprintf(stderr, "\rBuilding IntervalList for field %u - done           \n", f);

		// Initiate checkpoint value
		checkpoint = packets_num < 100 ? 1 : packets_num/100;

		// Count how many non-unique rules are there in the current field
		set<uint32_t> current_non_unique;

		// Rule iterator
		auto rule_it = rule_db.begin();

		// Build packets
		for (uint32_t i=0; i<packets_num; ++i) {
			// Update status
			if (i%checkpoint==0) {
				fprintf(stderr, "\rGenerating packets for field %u... (%u%%)...", f, i/checkpoint);
			}

			uint32_t rule_idx = i % rule_db.size();

			// In case the current rule is non_unique
			if (rule_intervals[rule_idx].size() == 0) {
				current_non_unique.insert(rule_idx);

				// Randomize the value of the current packet
				uint32_t low = rule_it->fields[f].low;
				uint32_t high = rule_it->fields[f].high;
				output[i].header[f] = gen_uniform_random_uint32(low, high);
				// Match index is not known yet
			}
			// Otherwise, randomize packet value in field
			else {
				output[i].header[f] = rule_intervals[rule_idx].random_value();
				output[i].match_priority = priorities[rule_idx];
			}

			// Update iterator
			if (++rule_it == rule_db.end()) {
				rule_it = rule_db.begin();
			}
		}

		// Update the non_unique rule set
		if (f == 0) {
			non_unique = current_non_unique;
		} else {
			set<uint32_t> intersect;
			set_intersection(
					non_unique.begin(), non_unique.end(),
					current_non_unique.begin(), current_non_unique.end(),
					std::inserter(intersect,intersect.begin()));
			non_unique = intersect;
		}

		fprintf(stderr, "\rGenerating packets for field %u - done.      \n", f);
	}

	fprintf(stderr, "Non-unique rules: %u\n", (uint32_t)non_unique.size());

	// Count how many non-reachable rules are there
	set<uint32_t> unreachable_rules;

	// Initiate checkpoint value
	checkpoint = non_unique.size() < 100 ? 1 : non_unique.size()/100;
	counter=0;

	// Handle non-unique rules...
	for (auto idx : non_unique) {

		// Update status
		if (counter%checkpoint==0) {
			fprintf(stderr, "\rHandling non-unique rules (%u%%)...", counter/checkpoint);
		}
		++counter;

		// Skip unreachable rules
		if (unreachable_rules.find(idx) != unreachable_rules.end()) {
			continue;
		}

		// For each packet of this rule..
		for (uint32_t i=idx; i<packets_num; i+= rule_db.size()) {
			trace_packet packet = gen_packet_in_rule(rule_db, idx, 5);
			output[i] = packet;
			if (packet.match_priority == 0xffffffff) {
				malformed_packets.insert(i);
				unreachable_rules.insert(idx);
			} else {
				// The current rule is reachable
				unreachable_rules.erase(idx);
			}
		}
	}

	fprintf(stderr, "\rHandling non-unique rules - done. Unreachable rules: %u. Malformed packets: %u \n",
			(uint32_t)unreachable_rules.size(), (uint32_t)malformed_packets.size());

	return output;
}


/**
 * @brief Generate at least a single packet per rule
 * @param rule_db The rule database
 * @param packet_num The number of packets
 * @param malformed_packets A set of invalid packet indices
 */
vector<trace_packet> generate_random_packets(const list<openflow_rule>& rule_db, uint32_t packets_num, set<uint32_t>& malformed_packets) {

	// Store pointers to rules
	vector<const openflow_rule*> rule_ptr_list;

	// The output packets
	vector<trace_packet> output;
	output.resize(packets_num);

	// Randomize packets
	infof("Randomizing %u packets...", packets_num);
	for (uint32_t i=0; i<packets_num; ++i) {
		// Choose rule by random
		uint32_t rule_idx = gen_uniform_random_uint32(0, rule_db.size()-1);

		// Get the rule
		auto rule_it = rule_db.begin();
		for (uint32_t i=0; i<rule_idx; ++i) {
			++rule_it;
		}
		const openflow_rule& rule = *rule_it;
		rule_ptr_list.push_back(&rule);

		// Choose field values by random
		output[i].header.resize(field_num);
		for (uint32_t j=0; j<field_num; ++j) {
			uint32_t field_start = rule_it->fields[j].low;
			uint32_t field_end = rule_it->fields[j].high;
			// Check for errors
			if (field_start > field_end) {
				throw errorf("field start (%u) > field end (%u) for rule %u",
						field_start, field_end, rule_idx);
			}
			output[i].header[j] = gen_uniform_random_uint32(field_start, field_end);
			output[i].match_priority = rule_it->priority; // Temporary match index, will be verified later
		}
	}

	fprintf(stderr, "Done. Starting to perform rule-matching...");

	uint32_t checkpoint = packets_num/100;
	if (checkpoint==0) checkpoint=1;

	// Search exact rule-match for each packet
	for (uint32_t i=0; i<packets_num; ++i) {
		// Update status
		if (i%checkpoint==0) {
			fprintf(stderr, "\rPerforming rule-matching for packets (%u%%)...", i/checkpoint);
		}

		int found=0;

		// Search for the first match rule
		uint32_t rule_idx = 0;
		for(auto rule : rule_db) {
			int match=1;
			// For each filed
			for (uint32_t j=0; j<field_num; ++j) {
				// Get rule boundaries
				uint32_t field_start = rule.fields[j].low;
				uint32_t field_end = rule.fields[j].high;
				// Check collision
				if ( (output[i].header[j] < field_start) ||
					 (output[i].header[j] > field_end ) )
				{
					match=0;
					break;
				}
			}
			// In case of match
			if (match) {
				output[i].match_priority = rule.priority;
				found=1;
				break;
			}
			++rule_idx;
		}

		// Check for errors
		if (found==0) {
			// Print packet
			cerr << endl <<
					"Warning: Packet " << i << "[" << output[i].to_string() << "]" <<
					"did not match any rule. Packet origin rule: [" << rule_ptr_list[i]->to_string() << "]" <<
					"Skipping packet" << endl;
			malformed_packets.insert(i);
		}
	}
	return output;
}

/**
 * @brief Generates a packet within the required rule index
 * @return 1 On success, 0 after maximum tries
 */
trace_packet gen_packet_in_rule(const list<openflow_rule>& rule_db, uint32_t rule_idx, uint32_t tries) {
	uint32_t r, counter = 0;
	uint32_t priority = 0xffffffff;
	// Initiate an invalid
	trace_packet packet;
	packet.match_priority = 0xffffffff;
	packet.header.resize(field_num);

	// Get the rule
	auto rule_it = rule_db.begin();
	for (uint32_t i=0; i<rule_idx; ++i) {
		++rule_it;
	}
	const openflow_rule& rule = *rule_it;

	do {

		// Choose field values by random
		for (uint32_t j=0; j<field_num; ++j) {
			packet.header[j] = gen_uniform_random_uint32(rule.fields[j].low, rule.fields[j].high);
		}

		// Validate the rule is indeed the requested index
		rule_it = rule_db.begin();
		for (r=0; r<=rule_idx; ++r) {

			int match=1;

			// For each filed
			for (uint32_t j=0; j<field_num; ++j) {
				// Get rule boundaries
				uint32_t field_start = rule_it->fields[j].low;
				uint32_t field_end = rule_it->fields[j].high;
				// Check collision
				if ( (packet.header[j] < field_start) || (packet.header[j] > field_end ) ) {
					match=0;
					break;
				}
			}

			if (match){
				priority = rule_it->priority;
				break;
			}
			++rule_it;
		}

		// Returns an invalid apcket
		if (++counter == tries) {
			return packet;
		}

		// Break only if the matching rule is with the desired index
	} while (r != rule_idx);

	// Found match!
	// Set the packet index
	packet.match_priority = priority;
	return packet;
}


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
#include <map>
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

// Holds arguments informatio
static argument_t my_arguments[] = {
		// Name,			Required,	IsBoolean,	Default,	Help
		{"-f",				1,			0,			NULL,		"Filter filename. Supported formats: [Classbench, Classbnech-ng, Binary]"},
		{"--indices",		0,			0,			NULL,		"Load an indices file to apply a filter on the ruleset"},
		{"--reverse",		0,			1,			NULL,		"Reverse the effect of the indices file"},
		{"-p",				0,			1,			NULL,		"Prints rule-set to stdout"},
		{"-o",				0,			0,			NULL,		"Output trace filename"},
		{"-s",				0,			1,			NULL,		"Shuffle trace"},
		{"-n",				0,			0,			"0",		"Number of packets to generate (used only when no locality file is specified)."},
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


typedef vector<trace_packet> header_options;
typedef map<int, header_options> rule_mapping_t;

rule_mapping_t generate_mapping( const list<openflow_rule>& rule_db);
vector<int> generate_uniform_locality(size_t num_of_packets, size_t num_of_rules);
void generate_trace(rule_mapping_t mapping, vector<int>& trace_indices, FILE* file);
trace_packet gen_packet_in_rule(const list<openflow_rule>& rule_db, int rule_idx, int tries);
vector<int> load_custom_locality(const char* filename);

// How many different trace options to generate for each
// rule in the ruleset?
#define OPTIONS 1

// Holds the rule-set attributes
static uint32_t field_num;
static ruleset_type_t ruleset_type;

/**
 * @brief Prints progres to the screen
 */
void print_progress(int counter, const char* message, size_t size) {
	if ( (size ==0) || (counter < 0) ) {
		fprintf(stderr, "\r%s... Done   \n", message);
	} else {
		int checkpoint = size < 100 ? 1 : size/100;
		if (counter%checkpoint==0) {
			fprintf(stderr, "\r%s... (%u%%)", message, counter/checkpoint);
		}
	}
}

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
		std::set<uint32_t> indices = read_indices_file(reader);
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

	// Generate a unique mapping between the rules and packet header
 	auto mapping = generate_mapping(rule_db);

 	// Randomize trace indices
 	uint32_t num_of_packets = atoi( ARG("-n")->value );	
 	vector<int> trace_indices = generate_uniform_locality(num_of_packets, rule_db.size());

	// Shuffle trace indices, if necessary
	if (ARG("-s")->available) {
		messagef("Shuffeling trace...");
		uint32_t* perm = new uint32_t[trace_indices.size()];
		random_permutation(perm, trace_indices.size());
		vector<int> trace_indices_new(trace_indices.size());
		for (size_t i=0; i<trace_indices.size(); ++i) {
			trace_indices_new[i] = trace_indices[perm[i]];
		}
		trace_indices = trace_indices_new;
		delete[] perm;
	}

	generate_trace(mapping, trace_indices, out_file_ptr);

	return 0;
}


/**
 * @brief Generates a mapping between a rule index to a
 * random header packet that matches the rule.
 * @param rule_db The ruleset
 * @param A set of unreachable rule indices
 * @returns A mapping rule_idx->(priority, header)
 */
rule_mapping_t generate_mapping(const list<openflow_rule>& rule_db) {
	rule_mapping_t output;

	// Count how many non-unique rules are there
	set<int> non_unique;

	auto rule = rule_db.begin();
	for (size_t i=0; i<rule_db.size(); ++i) {
		output[i] = vector<trace_packet>();
		for (int j=0; j<OPTIONS; ++j) {
			output[i].push_back(trace_packet());
			output[i].back().header.resize(field_num);
		}	
		++rule;
	}

	// For each field
	for (uint32_t f=0; f<field_num; ++f) {
		// Build an interval-list for the current field
		IntervalList interval(0, get_field_bound(f, ruleset_type));
		auto rule = rule_db.begin();

		// Count how many non-unique rules are there in the current field
		set<int> current_non_unique;

		char message[256];
		snprintf(message, 256, "Calculating interval-list for field %d", f);

		for (size_t i=0; i<rule_db.size(); ++i) {
			// Print progress to screen
			print_progress(i, message, rule_db.size());
			// Calculate the interval for the current rule
			auto sub_interval = interval.apply_rule(rule->fields[f].low, rule->fields[f].high);
			// In case we can guarantee unique value for the current rule
			if (sub_interval.size() > 0) {
				// Fill values for the current field
				for (int j=0; j<OPTIONS; ++j) {
					output[i][j].header[f] = sub_interval.random_value();
				}
			}
			// We cannot guarantee a unique mapping
			else {
				uint32_t low = rule->fields[f].low;
				uint32_t high = rule->fields[f].high;
				current_non_unique.insert(i);
				for (int j=0; j<OPTIONS; ++j) {
					output[i][j].header[f] = gen_uniform_random_uint32(low, high);
				}
			}
			// Update match priority
			for (int j=0; j<OPTIONS; ++j) {
				output[i][j].match_priority = rule->priority;
			}
			++rule;
		}

		// Update the non_unique rule set
		if (f == 0) {
			non_unique = current_non_unique;
		} else {
			set<int> intersect;
			set_intersection(
					non_unique.begin(), non_unique.end(),
					current_non_unique.begin(), current_non_unique.end(),
					std::inserter(intersect, intersect.begin()));
			non_unique = intersect;
		}
		
		// Print progress to screen
		print_progress(-1, message, 1);
	}

	// Update mapping for non-unique rules
	fprintf(stderr, "Non-unique rules: %lu\n", non_unique.size());

	set<int> unreachable_rules;

	// Handle non-unique rules...
	int counter = 0;
	for (auto idx : non_unique) {
		print_progress(counter++, "Handling non-unique rules", non_unique.size());
		// Skip unreachable rules
		if (unreachable_rules.find(idx) != unreachable_rules.end()) {
			continue;
		}
		// Try to generate a header that matches the rule. No more than 5 times.
		uint32_t desired_priority = output[idx][0].match_priority;
		trace_packet packet = gen_packet_in_rule(rule_db, idx, 5);
		// We don't care for duplicates in this case
		for (int j=0; j<OPTIONS; ++j) {
			output[idx][j] = packet;
		}
		if (packet.match_priority != desired_priority) {
			unreachable_rules.insert(idx);
		} else {
			// The current rule is reachable
			unreachable_rules.erase(idx);
		}
	}
	print_progress(0, "Handling non-unique rules", 0);
	fprintf(stderr, "Unreachable rules: %lu\n", unreachable_rules.size());
	return output;
}

/**
 * @brief Generates a locality vector
 * @param num_of_packets Number of packets in trace
 * @param num_of_rules Number of rules in the ruleset
 */
vector<int> generate_uniform_locality(size_t num_of_packets, size_t num_of_rules) {
	vector<int> output(num_of_packets);
	messagef("Generating uniform locality...");
	for (size_t i=0; i<num_of_packets; ++i) {
		output[i] = gen_uniform_random_uint32(0, num_of_rules-1);
	}
	return output;
}

/**
 * @brief Generate trace from mapping and indices. Prints to output file
 * @param mapping A mapping from rule-index to packet header
 * @param trace_indices The indices of the packets in trace
 * @param file The output file
 */
void generate_trace(rule_mapping_t mapping, vector<int>& trace_indices, FILE* file) {
	for (size_t i=0; i<trace_indices.size(); ++i) {
		// Print progress
		print_progress(i, "Generating trace", trace_indices.size());

		// Which rule are we writing?
		int rule_idx = trace_indices[i] % mapping.size();

		// Select a packet from the rule options
		vector<trace_packet>& packet_options = mapping[rule_idx];
		int option_idx = gen_uniform_random_uint32(0, packet_options.size()-1);
		trace_packet& packet = packet_options[option_idx];

		// Print packet fields
		for (uint32_t j=0; j<field_num; ++j) {
			fprintf(file, "%u\t", packet.header[j]);
		}

		// Print rule number
		fprintf(file, "%u\n", packet.match_priority);
	}
	print_progress(0, "Generating trace", 0);
}

/**
 * @brief Generates a packet within the required rule index. Does not always succeed.
 * @param rule_db The ruleset
 * @param rule_idx The required rule index
 * @param tries Number of tries
 * @returns A trace packet
 */
trace_packet gen_packet_in_rule(const list<openflow_rule>& rule_db, int rule_idx, int tries) {

	trace_packet output;
	output.header.resize(field_num);

	// Get a reference to the desired rule
	auto rule_it = rule_db.begin();
	for (int i=0; i<rule_idx; ++i) {
		++rule_it;
	}
	const openflow_rule& rule = *rule_it;

	while (tries > 0) {
		// Choose field values by random
		for (uint32_t j=0; j<field_num; ++j) {
			output.header[j] = gen_uniform_random_uint32(rule.fields[j].low, rule.fields[j].high);
		}

		// Validate the rule is indeed the requested index
		rule_it = rule_db.begin();
		for (int r=0; r<=rule_idx; ++r) {

			int match=1;

			// For each filed
			for (uint32_t j=0; j<field_num; ++j) {
				// Get rule boundaries
				uint32_t field_start = rule_it->fields[j].low;
				uint32_t field_end = rule_it->fields[j].high;
				// Check collision
				if ( (output.header[j] < field_start) || (output.header[j] > field_end ) ) {
					match=0;
					break;
				}
			}

			if (match){
				output.match_priority = rule_it->priority;
				break;
			}
			++rule_it;
		}

		// Can break - found a desired header!
		if ((int)output.match_priority == rule_idx) {
			break;
		}

		--tries;
	}
	return output;
}


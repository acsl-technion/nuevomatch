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

#include <tuple_merge.h>
#include <math.h>
#include <vector>
#include <unordered_map>
#include <time.h>
#include <sstream>
#include <ElementaryClasses.h>
#include <TupleMerge/TupleMergeOnline.h>

using namespace std;

// Macro for addressing classifier as TupleMerge
#define CLASSIFIER reinterpret_cast<TupleMergeOnline*>(this->tm_classifier)

// Macro for addressing rules 
#define RULES reinterpret_cast<std::vector<Rule>*>(this->my_rules)

TupleMerge::TupleMerge(int limit) : build_time(0), limit(limit) {
	std::unordered_map<string,string> m;
	stringstream ss;
	ss << limit;
	m.insert(std::pair<string, string>("TM.Limit.Collide", ss.str()));
	tm_classifier = new TupleMergeOnline(m);
	my_rules = new std::vector<Rule>();
}

TupleMerge::~TupleMerge() {
	delete CLASSIFIER;
	delete RULES;
};

/**
 * @brief Build the classifier data structure
 * @returns 1 On success, 0 on fail
 */
int TupleMerge::build(const std::list<openflow_rule>& rule_db) {
	for (auto r : rule_db) {
		if (r.fields.size() != 5) {
			throw errorf("Cannot build TupleMerge classifier: number of rule fields != 5");
		}
		std::vector<uint32_t> prefix_length;
		// Prefix from range:
		// perform XOR: finds the bits that differ (not part of prefix)
		// perform floor(log2)+1: find the largest differ-bit number
		// perform 32-: the prefix is the remaining number of bits
		for (int i=0; i<2; ++i) {
			int differ_bits=0;
			if (r.fields[i].high != r.fields[i].low) {
				differ_bits = floor(log2((r.fields[i].high^r.fields[i].low)))+1;
			}
			prefix_length.push_back(32-differ_bits);
		}
		// Port prefixes are easy. high == low ? 32 : 16
		prefix_length.push_back(r.fields[2].high == r.fields[2].low ? 32 : 16);
		prefix_length.push_back(r.fields[3].high == r.fields[3].low ? 32 : 16);
		// Protocol prefixes are also easy.
		prefix_length.push_back(r.fields[4].high == r.fields[4].low ? 32 : 24);
		// Check whether the rule prefixes were calculated
		// If yes, validate prefixes
		if (r.prefixes.size() != 0) {
			for (uint32_t i=0; i<5; i++) {
				if (prefix_length[i] != r.prefixes[i]) {
					throw errorf("Prefix calculation error");
				}
			}
		}
		Rule tm_rule;
		// Note: as TM chooses priority by MAX and not MIN,
		// we flip them (small prios become large)
		tm_rule.priority = 0x7fffffff - r.priority;
		tm_rule.id = r.priority;
		for (uint32_t i=0; i<5; i++) {
			tm_rule.range[i][0] = r.fields[i].low;
			tm_rule.range[i][1] = r.fields[i].high;
			tm_rule.prefix_length[i] = prefix_length[i];
		}
		RULES->push_back(tm_rule);
	}
	struct timespec start_time, end_time;
	clock_gettime(CLOCK_MONOTONIC, &start_time);
	CLASSIFIER->ConstructClassifier(*RULES);
	clock_gettime(CLOCK_MONOTONIC, &end_time);
	this->build_time = (double)((end_time.tv_sec * 1e9 + end_time.tv_nsec) -
			  (start_time.tv_sec * 1e9 + start_time.tv_nsec)) / 1e6;
	return 1;
}

/**
 * @brief Packs this to byte array
 * @returns An object-packer with the binary data
 */
ObjectPacker TupleMerge::pack() const {
	ObjectPacker out;
	out << limit;
	out << (uint32_t)RULES->size();
	// TODO - this is a simplistic packning without changing TM source code
	for (auto r : *RULES) {
		out << r.id << r.dim << r.tag
			<< r.priority << r.markedDelete;
		for (int i=0; i<r.dim; ++i) {
			out << r.prefix_length[i]
			    << r.range[i][0] << r.range[i][1];
		}
	}
	return out;
}

/**
 * @brief Creates this from a memory location
 * @param object An object-reader instance
 */
void TupleMerge::load(ObjectReader& object) {
	RULES->clear();

	// Load the limit value
	limit = object.read<uint32_t>();
	std::unordered_map<string,string> m;
	stringstream ss;
	ss << limit;
	m.insert(std::pair<string, string>("TM.Limit.Collide", ss.str()));
	delete CLASSIFIER;
	tm_classifier = new TupleMergeOnline(m);

	uint32_t num_of_rules = object.read<uint32_t>();
	for (uint32_t i =0; i<num_of_rules; ++i) {
		Rule r;
		object >> r.id >> r.dim >> r.tag >> r.priority >> r.markedDelete;
		for (uint32_t j=0; j<r.dim; ++j) {
			object >> r.prefix_length[j] >> r.range[j][0] >> r.range[j][1];
		}
		RULES->push_back(r);
	}
	struct timespec start_time, end_time;
	clock_gettime(CLOCK_MONOTONIC, &start_time);
	CLASSIFIER->ConstructClassifier(*RULES);
	clock_gettime(CLOCK_MONOTONIC, &end_time);
	this->build_time = (double)((end_time.tv_sec * 1e9 + end_time.tv_nsec) -
			  (start_time.tv_sec * 1e9 + start_time.tv_nsec)) / 1e6;
}

/**
 * @brief Returns the number of rules
 */
unsigned int TupleMerge::get_num_of_rules() const {
	return RULES->size();
}

/**
 * @brief Returns the memory size of this in bytes
 */
unsigned int TupleMerge::get_size() const {
	return CLASSIFIER->MemSizeBytes();
}

/**
 * @brief Returns the building time of this in milliseconds
 */
unsigned int TupleMerge::get_build_time() const {
	return this->build_time;
}

/**
 * @brief Returns the maximum supported number of fields this can classify
 */
const unsigned int TupleMerge::get_supported_number_of_fields() const {
	return 5;
}

/**
 * @brief Starts the performance measurement of this
 */
void TupleMerge::start_performance_measurement() {
	clock_gettime(CLOCK_MONOTONIC, &start_time);
}

	/**
	 * @brief Stops the performance measurement of this
	 */
void TupleMerge::stop_performance_measurement() {
	clock_gettime(CLOCK_MONOTONIC, &end_time);
}

/**
 * @brief clones this to another instance
 */
GenericClassifier* TupleMerge::clone() {
	TupleMerge* new_classifier = new TupleMerge();
	ObjectReader reader = ObjectReader(pack());
	new_classifier->load(reader);
	return new_classifier;
}

/**
 * @brief Start an asynchronous process of classification for an input packet.
 * @param header An array of 32bit integers according to the number of supported fields.
 * @returns A unique id for the packet
 */
unsigned int TupleMerge::classify_async(const unsigned int* header, int priority) {
	uint32_t match_id = -1;
	uint32_t packet_id = 0xffffffff;

	// Lookup only valid packets
	if (header != nullptr) {
		// Perform lookup
		match_id = classify_sync(header, priority);
		packet_id = _packet_counter++;
	}

	// Notify all listeners
	for (auto it : _listeners) {
		it->on_new_result(packet_id, match_id, match_id, _additional_args);
	}
	return packet_id;
}

/**
 * @brief Start a synchronous process of classification an input packet.
 * @param header An array of 32bit integers according to the number of supported fields.
 * @returns The matching rule action/priority (or 0xffffffff if not found)
 */
unsigned int TupleMerge::classify_sync(const unsigned int* header, int priority) {
	// Packet not found
	if (header == nullptr) return -1;
	// Must use adapter for packet
	Packet p(5);
	for (int i=0; i<5; ++i) {
		p[i] = header[i];
	}
	// Note: as TM chooses priority by MAX and not MIN,
	// we flip them (small prio become large)
	uint32_t output = CLASSIFIER->ClassifyAPacket(p, 0x7fffffff - priority);
	return output == -1 ? -1 : 0x7fffffff - output;
}

/**
 * @brief Prints debug information
 * @param verbose Set the verbosity level of printing
 */
void TupleMerge::print(uint32_t verbose) const {
	messagef("Tuple Merge Classifier");
	if (verbose > 1) {
		messagef("There are %d tables:", CLASSIFIER->NumTables());
		uint32_t rule_num=0;
		for (auto i=0; i<CLASSIFIER->NumTables(); ++i) {
			messagef("Table %d with %d rules",i ,CLASSIFIER->RulesInTable(i));
			rule_num+=CLASSIFIER->RulesInTable(i);
		}
		messagef("Collision Limit: %u", limit);
		messagef("Total rules: %u", rule_num);
		messagef("Build time %.3f ms", this->build_time);
	}

	// Measure performance
	uint32_t total_usec = (end_time.tv_sec * 1e6 + end_time.tv_nsec / 1e3) -
						  (start_time.tv_sec * 1e6 + start_time.tv_nsec /1e3);
	messagef("Performance: total time %u usec. Average time: %.3f usec per packet.",
			total_usec, (double)total_usec / _packet_counter);
}

/**
 * @brief Returns a string representation of this
 */
const std::string TupleMerge::to_string() const {
	return "OnlineTupleMerge";
}



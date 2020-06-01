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

#pragma once

#include <object_io.h>
#include <logging.h>

#include <vector>
#include <set>
#include <list>
#include <string>


#define FIVE_TUPLE_FIELDS 5

/**
 * @brief Simulates an input packet when using packet-traces
 */
struct trace_packet{
    std::vector<uint32_t> header;
    uint32_t match_priority;

    /**
     * @brief Returns a pointer to the header values
     */
    const uint32_t* get() const { return &header[0]; }

    /**
     * @brief Returns a string representation of this
     */
    std::string to_string() const;
};


/**
 * @brief A 32-bit range [low, high] (inclusive)
 */
struct range {
	uint32_t low;
	uint32_t high;
	uint32_t calc() { return high - low; }

	range() : low(0), high(0) {};
	range(uint32_t low, uint32_t high) : low(low), high(high) {};

	/**
	 * @brief Build range from a vector of values
	 */
	range(std::vector<uint32_t> v) {
		if (v.size() != 2) {
			throw error("Cannot set range from vector with size != 2");
		}
		low = v[0];
		high = v[1];
	}

	/**
	 * @brief returns true iff low <= high
	 */
	bool is_valid() const {
		return low <= high;
	}

    /**
     * @brief Returns a string representation of this
     */
    std::string to_string() const;

    /**
     * @brief Forces order between ranges. Ranges are ordered
     *        by their start value, then by their length.
     *        e.g., (0, 3) < (0, 10) < (1, 2)
     * @note Required for using range within sets
     *
     */
    bool operator<(const range& rhs) const {
    	return (this->low < rhs.low) || (this->high < rhs.high);
    }
};

/**
 * @brief A 5-tuple matching rule (used with Classbench)
 */
struct matching_rule {
	unsigned int priority;
	range field[5];
};

/**
 * @brief An OpenFlow rule with arbitrary number of fields
 */
struct openflow_rule {
	unsigned int priority;
	std::vector<range> fields;
	std::vector<uint32_t> prefixes; // Used for hash-based algos, such as TM

	openflow_rule() : priority(0) {};

	/**
	 * @brief Create new OpenFlow rule from a classic 5-tuple rule
	 */
	openflow_rule(const matching_rule& rule);

	/**
	 * @brief Converts this to a 5-tuple rule
	 * @throws In case the rule cannot be converted
	 */
	matching_rule convert_to_five_tuple() const;

    /**
     * @brief Returns a string representation of this
     */
    std::string to_string() const;

    /**
     * @brief Used for sorting rules in STL containers
     */
    bool operator<(const openflow_rule& rhs) const;
};

/**
 * @brief Used to verify the rule-set origin
 */
typedef enum { UNKNOWN = 0, CLASSBENCH, CLASSBENCHNG, BINARY } ruleset_type_t;

/**
 * @brief Loads rules database (as generated using NuevoMatch) from memory
 * @param reader An object-reader with binary data
 * @return A list of open-flow rules
 */
std::list<openflow_rule> load_rule_database(ObjectReader& reader);

/**
 * @brief Reads Classbench file into rule_table_t data structure
 * @param filename Path to Classbench file
 * @return A list of open-flow rules
 * @throws In case of an error
 */
std::list<openflow_rule> read_classbench_file(const char* filename);

/**
 * @brief Reads Classbench-ng file with OpenFlow 1.0 rules
 * @param filename Path to Classbench-ng file
 * @return A list of open-flow rules
 * @note The resulting rule-table has at most 9 fields
 */
std::list<openflow_rule> read_classbench_ng_file(const char* filename);

/**
 * @brief Read a textual rule-set file and identify its origins.
 * @param filename Path to the rule-set file
 * @return The type of the rule-set
 */
ruleset_type_t classify_ruleset_file(const char* filename);

/**
 * @brief Returns the upper bound for a field of a rule-set
 * @param field_idx The index of the field within the rule-db
 * @param type The rule-set type
 */
uint64_t get_field_bound(int field_idx, ruleset_type_t type);

/**
 * @brief Truncate the ruleset according to a set of indices
 * @param ruleset Original ruleset
 * @param indices A set of indices
 * @param reverse If true, select the rules with indices not in the indices set
 */
std::list<openflow_rule> apply_indices_on_ruleset(std::list<openflow_rule>& ruleset, std::set<uint32_t>& indices, bool reverse);

/**
 * @brief Reads a ruleset indices file and return a set with all indices available in file
 * @param reader An object-reader with binary data
 */
std::set<uint32_t> read_indices_file(ObjectReader& reader);

/**
 * @brief Reads a textual trace file into memory
 * @param[in] trace_filename The textual trace filename
 * @param[in] indices A vector of custom fields to look at
 * @param[out] num_of_packets The number of packet in trace
 * @returns An array of trace packet headers, or NULL in case of an error
 */
trace_packet* read_trace_file(const char* trace_filename, const std::vector<uint32_t>& indices, uint32_t* num_of_packets);

/**
 * @brief Prints a rule-db to stdout
 * @param rule_db The rule database
 */
void print_rule_db(std::list<openflow_rule> rules);


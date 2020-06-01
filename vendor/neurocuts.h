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

#include <time.h>

#include <generic_classifier.h>
#include <rule_db.h>

/**
 * @brief A NeuroCuts classifier implementation in C++
 */
class NeuroCuts : public GenericClassifier {
private:

	struct node {
		bool is_partition;
		uint32_t id;
		uint32_t depth;
		uint32_t num_of_children;
		uint32_t num_of_rules;
		int max_priority;
		node** children;
		node* parent;
		matching_rule** rules;
		range boundaries[FIVE_TUPLE_FIELDS];
	};

	uint32_t _num_of_rules;
	uint32_t _size;
	uint32_t _build_time;
	uint32_t _max_binth;

	node* _nodes;
	matching_rule* _rules;
	node* _root;

	// Clones should not delete nodes and rules
	bool _is_clone;

	// Performance
	struct timespec perf_start_time, perf_end_time;

	/**
	 * @brief Helper recursive method for finding matching node
	 * @param current The node to check
	 * @param header The input packet header
	 * @returns A pointer to a matching rule
	 */
	matching_rule* match(node* current, const uint32_t* header, int priority);

	/**
	 * @brief Used for debugging. Prints a NeuroCuts node
	 */
	void print_node(node* n);

	/**
	 * @brief Computes the memory access of this
	 * @param n The node from which start calculating. At normal cases, should be root.
	 */
	uint32_t compute_mem_access(node* n) const;

public:

	NeuroCuts();
	~NeuroCuts();

	/**
	 * @brief Build the classifier data structure
	 * @returns 1 On success, 0 on fail
	 */
	int build(const std::list<openflow_rule>& rule_db) {
		throw std::runtime_error("NeuroCuts build should be done using NeuroCuts library");
	}

	/**
	 * @brief Packs this to byte array
	 * @returns An object-packer with the binary data
	 */
	ObjectPacker pack() const {
		throw std::runtime_error("NeuroCuts pack should be done using Python library");
	}

	/**
	 * @brief Creates this from a memory location
	 * @param object An object-reader instance
	 */
	void load(ObjectReader& object);

	/**
	 * @brief Returns the number of rules
	 */
	uint32_t get_num_of_rules() const { return _num_of_rules; }

	/**
	 * @brief Returns the memory size of this in bytes
	 */
	uint32_t get_size() const { return _size; }

	/**
	 * @brief Returns the building time of this in milliseconds
	 */
	uint32_t get_build_time() const { return _build_time; }

	/**
	 * @brief Returns the maximum supported number of fields this can classify
	 */
	virtual const unsigned int get_supported_number_of_fields() const { return FIVE_TUPLE_FIELDS; }

	/**
	 * @brief Starts the performance measurement of this
	 */
	void start_performance_measurement();

	/**
	 * @brief Stops the performance measurement of this
	 */
	void stop_performance_measurement();

	/**
	 * @brief clones this to another instance
	 */
	virtual GenericClassifier* clone() {
		NeuroCuts* clone = new NeuroCuts(*this);
		clone->_is_clone = true;
		return clone;
	}

	/**
	 * @brief Start an asynchronous process of classification for an input packet.
	 * @param header An array of 32bit integers according to the number of supported fields.
	 * @param priority The priority of a previous matching rule.
	 * Stops classifiying when there is no potential better priority
	 * @returns A unique id for the packet
	 */
	virtual unsigned int classify_async(const unsigned int* header, int priority);

	/**
	 * @brief Start a synchronous process of classification an input packet.
	 * @param header An array of 32bit integers according to the number of supported fields.
	 * @param priority The priority of a previous matching rule.
	 * Stops classifiying when there is no potential better priority
	 * @returns The matching rule action/priority (or 0xffffffff if not found)
	 */
	virtual unsigned int classify_sync(const unsigned int* header, int priority);

	/**
	 * @brief Prints statistical information
	 * @param verbose Set the verbosity level of printing
	 */
	virtual void print(uint32_t verbose=1) const;

	/**
	 * @brief Returns a string representation of this
	 */
	virtual const std::string to_string() const { return "NueroCuts"; }
};

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

#include <vector>
#include <queue>
#include <list>
#include <map>
#include <algorithm>

#include <generic_classifier.h>
#include <object_io.h>
#include <rule_db.h>


struct node {
	int depth;
	int problematic;
	int node_has_rule;
	int row;
	int column;
	int index;
	bool is_compressed;
	int count;

	matching_rule boundary;

	std::list <matching_rule*> classifier;
	std::list <node *> children;
	std::list <node *> actual_children;
	std::vector<int> cuts;

	node() : depth(0), problematic(0), node_has_rule(0),
		  row(0), column(0), index(0), is_compressed(0), count(0)
	{
	  cuts.resize(CLASSIFIER_FIELDS);
	  for (int i = 0; i < CLASSIFIER_FIELDS; i++) {
		  cuts[i] = 0;
	  }
	}
};

struct TreeDetails {
	node* root;
	std::vector<bool> wideFields;
	TreeDetails() : root(nullptr) {
		wideFields.resize(CLASSIFIER_FIELDS);
	}
};

/**
 * @brief A CutSplit object matching the vector interface
 */
class EffiCuts : public GenericClassifier {
private:

	uint32_t _num_of_rules;
	uint32_t _binth;
	uint32_t _size;
	uint32_t _build_time;

	matching_rule* rule_db;

	// Performance
	struct timespec perf_start_time, perf_end_time;

	// Trees
	std::list<TreeDetails> _trees;

public:

	EffiCuts(uint32_t binth);
	~EffiCuts();

	/**
	 * @brief Build the classifier data structure
	 * @returns 1 On success, 0 on fail
	 */
	int build(const std::list<openflow_rule>& rule_db);

	/**
	 * @brief Packs this to byte array
	 * @returns An object-packer with the binary data
	 */
	ObjectPacker pack() const;

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
		return new EffiCuts(*this);
	}

	/**
	 * @brief Start an asynchronous process of classification for an input packet.
	 * @param header An array of 32bit integers according to the number of supported fields.
	 * @returns A unique id for the packet
	 */
	uint32_t classify_async(const uint32_t* header, int priority);

	/**
	 * @brief Start a synchronous process of classification an input packet.
	 * @param header An array of 32bit integers according to the number of supported fields.
	 * @returns The matching rule action/priority (or 0xffffffff if not found)
	 */
	uint32_t classify_sync(const uint32_t* header, int priority);

	/**
	 * @brief Prints statistical information
	 * @param verbose Set the verbosity level of printing
	 */
	virtual void print(uint32_t verbose=1) const;

	/**
	 * @brief Returns a string representation of this
	 */
	virtual const std::string to_string() const { return "EffiCuts"; }

	/**
	 * @brief Returns the maximum supported number of fields this can classify
	 */
	virtual const unsigned int get_supported_number_of_fields() const { return FIVE_TUPLE_FIELDS; }
};

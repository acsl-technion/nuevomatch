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

#include <generic_classifier.h>

/**
 * @brief Adapter for TupleMergeOnline Classifier
 */
class TupleMerge : public GenericClassifier {
public:


	TupleMerge(int limit = 10);
	~TupleMerge();

	/**
	 * @brief Build the classifier data structure
	 * @returns 1 On success, 0 on fail
	 */
	virtual int build(const std::list<openflow_rule>& rule_db);

	/**
	 * @brief Packs this to byte array
	 * @returns An object-packer with the binary data
	 */
	virtual ObjectPacker pack() const;

	/**
	 * @brief Creates this from a memory location
	 * @param object An object-reader instance
	 */
	virtual void load(ObjectReader& object);

	/**
	 * @brief Returns the number of rules
	 */
	virtual unsigned int get_num_of_rules() const;

	/**
	 * @brief Returns the memory size of this in bytes
	 */
	virtual unsigned int get_size() const;

	/**
	 * @brief Returns the building time of this in milliseconds
	 */
	virtual unsigned int get_build_time() const;

	/**
	 * @brief Returns the maximum supported number of fields this can classify
	 */
	virtual const unsigned int get_supported_number_of_fields() const;

	/**
	 * @brief Starts the performance measurement of this
	 */
	virtual void start_performance_measurement();

	/**
	 * @brief Stops the performance measurement of this
	 */
	virtual void stop_performance_measurement();

	/**
	 * @brief clones this to another instance
	 */
	virtual GenericClassifier* clone();

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
	 * @brief Prints debug information
	 * @param verbose Set the verbosity level of printing
	 */
	virtual void print(uint32_t verbose=1) const;

	/**
	 * @brief Returns a string representation of this
	 */
	virtual const std::string to_string() const;

protected:
	struct timespec start_time, end_time;
	void* my_rules;
	void* tm_classifier;
	uint32_t build_time;
	uint32_t limit;
};




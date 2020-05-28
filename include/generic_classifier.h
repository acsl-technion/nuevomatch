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
#include <list>

#include <object_io.h>
#include <rule_db.h>


#define CLASSIFIER_FIELDS 5

/**
 * @brief Stores output pair for any classifier
 */
typedef struct {
	int priority;
	int action;
} classifier_output_t;

class GenericClassifierListener {
public:

	/**
	 * @brief Is invoked by the classifier when new result is available
	 * @param id A unique packet id
	 * @param priority The priority of the rule
	 * @param action The action to take on the packet
	 * @param args Any additional information
	 */
	virtual void on_new_result(unsigned int id, int priority, int action, void* args) = 0;
	virtual ~GenericClassifierListener() {}
};

/**
 * @brief Abstract class, used to specify a common interface
 * for different packet classification vectors
 */
class GenericClassifier {
protected:

	// Listeners of results of this
	std::vector<GenericClassifierListener*> _listeners;

	// Counts packets
	volatile unsigned int _packet_counter;

	// Additional information to pass on new results
	void* _additional_args;

public:

	GenericClassifier() : _packet_counter(0), _additional_args(nullptr) {}

	/**
	 * @brief Build the classifier data structure
	 * @returns 1 On success, 0 on fail
	 */
	virtual int build(const std::list<openflow_rule>& rule_db) = 0;

	/**
	 * @brief Packs this to byte array
	 * @returns An object-packer with the binary data
	 */
	virtual ObjectPacker pack() const = 0;

	/**
	 * @brief Creates this from a memory location
	 * @param object An object-reader instance
	 */
	virtual void load(ObjectReader& object) = 0;

	/**
	 * @brief Returns the number of rules
	 */
	virtual unsigned int get_num_of_rules() const = 0;

	/**
	 * @brief Returns the memory size of this in bytes
	 */
	virtual unsigned int get_size() const = 0;

	/**
	 * @brief Returns the building time of this in milliseconds
	 */
	virtual unsigned int get_build_time() const = 0;

	/**
	 * @brief Returns the maximum supported number of fields this can classify
	 */
	virtual const unsigned int get_supported_number_of_fields() const = 0;

	/**
	 * @brief Starts the performance measurement of this
	 */
	virtual void start_performance_measurement() = 0;

	/**
	 * @brief Stops the performance measurement of this
	 */
	virtual void stop_performance_measurement() = 0;

	/**
	 * @brief clones this to another instance
	 */
	virtual GenericClassifier* clone() = 0;

	/**
	 * @brief Start an asynchronous process of classification for an input packet.
	 * @param header An array of 32bit integers according to the number of supported fields.
	 * @returns A unique id for the packet
	 */
	virtual unsigned int classify_async(const unsigned int* header) = 0;

	/**
	 * @brief Start a synchronous process of classification an input packet.
	 * @param header An array of 32bit integers according to the number of supported fields.
	 * @returns The matching rule action/priority (or 0xffffffff if not found)
	 */
	virtual unsigned int classify_sync(const unsigned int* header) = 0;

	/**
	 * @brief Prints debug information
	 * @param verbose Set the verbosity level of printing
	 */
	virtual void print(uint32_t verbose=1) const = 0;

	virtual ~GenericClassifier() {}

	/**
	 * @brief Returns a string representation of this
	 */
	virtual const std::string to_string() const = 0;

	/**
	 * @brief Adds new listener to this
	 */
	void add_listener(GenericClassifierListener& listener) {
		_listeners.push_back(&listener);
	}

	/**
	 * @brief Resets the all classifier counters
	 */
	virtual void reset_counters() { _packet_counter = 0; }

	/**
	 * @brief Sets any additional arguments to pass to listeners on new result
	 */
	virtual void set_additional_args (void* args) {
		_additional_args = args;
	}
};

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

#include <set>
#include <time.h>
#include <bits/stdc++.h> // UINT_MAX

#include <basic_types.h>
#include <pipeline_thread.h>
#include <cpu_core_tools.h>
#include <generic_classifier.h>

template <uint32_t N>
class ParallelClassifier : public GenericClassifier, public GenericClassifierListener {
private:

	// Each worker process job
	typedef struct {
		WorkBatch<const uint32_t*, N> packets;
		uint32_t worker_id;
		uint32_t first_packet_counter;
	} worker_job_t;

	// Worker threads
	PipelineThread<worker_job_t>** _workers;

	// The classifier objects
	uint32_t _size;
	GenericClassifier** _classifiers;

	// Store next batch
	worker_job_t _next_batch;
	uint32_t _next_batch_size;

	// Hold information for each classifier
	typedef struct {
		uint32_t first_packet_counter;
		classifier_output_t results[N];
		uint32_t valid_results;
	} worker_info_t;

	// Hold lock for publishing new results
	volatile uint32_t _lock;

	worker_info_t* _workers_info;

	// Performance
	struct timespec start_time, end_time;
	uint32_t _total_serial_packets;

	/**
	 * @brief Is invoked by the classifiers when new result is available
	 * @param id A unique packet id
	 * @param priority The priority of the rule
	 * @param action The action to take on the packet
	 * @param args Any additional information
	 */
	virtual void on_new_result(uint32_t id, int priority, int action, void* args) {
		// Extract the information
		worker_info_t* info = static_cast<worker_info_t*>(args);
		// Update result
		if (id != 0xffffffff) {
			info->results[id] = {priority, action};
			++info->valid_results;
		}
	}

	/**
	 * @brief The worker thread main method
	 * @param job The next job to process
	 * @param args A pointer to ParallelClassifier instance
	 */
	static bool worker_method(worker_job_t& job, void* args) {

		// Get instance and classifier
		ParallelClassifier* instance = static_cast<ParallelClassifier*>(args);
		GenericClassifier* classifier = instance->_classifiers[job.worker_id];

		// Reset the classifier counter
		classifier->reset_counters();

		// Set additional information
		worker_info_t* info = &instance->_workers_info[job.worker_id];
		info->first_packet_counter = job.first_packet_counter;
		info->valid_results = 0;
		classifier->set_additional_args(info);

		// Process all packets in the batch
		for (uint32_t i=0; i<N; ++i) {
			// Perform classification
			classifier->classify_async(job.packets[i], -1);
		}

		// Request lock for publishing
		while(__sync_val_compare_and_swap(&instance->_lock, 0, 1));

		// Publish results
		for (uint32_t i=0; i<info->valid_results; ++i) {
			for (auto it : instance->_listeners) {
				it->on_new_result(info->first_packet_counter + i,
						info->results[i].priority,
						info->results[i].action,
						nullptr);
			}
		}

		// Release lock for publishing
		instance->_lock = 0;

		// Finished job
		return true;
	}

public:

	/**
	 * @brief Initialize this with classifiers
	 * @param queue_size The size of the working queue
	 * @param size The number of classifiers
	 * @param classifiers An array of initialized classifiers
	 */
	ParallelClassifier(uint32_t queue_size, uint32_t size, GenericClassifier** classifiers) :
		_size(size), _classifiers(classifiers), _next_batch_size(0), _lock(0), _total_serial_packets(0)
	{
		loggerf("Initializing ParallelClassifier with %u workers", size);

		// Initialize the worker threads
		_workers = new PipelineThread<worker_job_t>*[_size-1];
		_workers_info = new worker_info_t[_size];

		// Hold all utilized CPUs
		std::set<int> cpu_set;
		int core_idx = cpu_core_tools_get_index_of_current_thread();
		cpu_set.insert(core_idx);
		loggerf("Classifier 0 on CPU %d", core_idx);
		_classifiers[0]->add_listener(*this);

		// Initialize each worker on distinct core
		for (uint32_t i=0; i<size-1; ++i) {

			// Request core for the current iSet
			core_idx = cpu_core_tools_get_next_physical_core(core_idx);

			// Check for performance hazards
			if (cpu_set.find(core_idx) != cpu_set.end()) {
				warningf("Classifier %u has no available free CPU core. performance may degenerate", i);
			}
			cpu_set.insert(core_idx);

			// Initialize new worker
			_workers[i] = new PipelineThread<worker_job_t>(queue_size, core_idx, worker_method, this);
			loggerf("Classifier %u thread on CPU %d", i+1, core_idx);

			// Register this as listener to classifiers
			_classifiers[i+1]->add_listener(*this);
		}

		// Set the first packet counter to be zero
		_next_batch.first_packet_counter = 0;
	}

	~ParallelClassifier() {
		for (uint32_t i=0; i<_size-1; ++i) {
			delete _workers[i];
		}
		for (uint32_t i=0; i<_size; ++i) {
			delete _classifiers[i];
		}
		delete _workers;
		delete _classifiers;
		delete _workers_info;
	}

	/**
	 * @brief Build the classifier data structure
	 * @returns 1 On success, 0 on fail
	 */
	virtual int build(const std::list<openflow_rule>& rule_db) {
		throw std::runtime_error("ParallelClassifier does not support build");
	}

	/**
	 * @brief Packs this to byte array
	 * @returns An object-packer with the binary data
	 */
	virtual ObjectPacker pack() const {
		throw std::runtime_error("ParallelClassifier does not support pack");
	}

	/**
	 * @brief Creates this from a memory location
	 * @param object An object-reader instance
	 */
	virtual void load(ObjectReader& object) {
		throw std::runtime_error("ParallelClassifier does not support load");
	}

	/**
	 * @brief Returns the number of rules
	 */
	virtual unsigned int get_num_of_rules() const {
		return _classifiers[0]->get_num_of_rules();
	}

	/**
	 * @brief Returns the memory size of this in bytes
	 */
	virtual unsigned int get_size() const {
		return _classifiers[0]->get_size();
	}

	/**
	 * @brief Returns the building time of this in milliseconds
	 */
	virtual unsigned int get_build_time() const {
		return _classifiers[0]->get_build_time();
	}


	/**
	 * @brief Starts the performance measurement of this
	 */
	virtual void start_performance_measurement() {
		for (uint32_t i=0; i<_size-1; ++i) {
			_workers[i]->start_performance_measurements();
		}
		for (uint32_t i=0; i<_size; ++i) {
			_classifiers[i]->start_performance_measurement();
		}
		clock_gettime(CLOCK_MONOTONIC, &start_time);
	}

	/**
	 * @brief Stops the performance measurement of this
	 */
	virtual void stop_performance_measurement() {
		for (uint32_t i=0; i<_size-1; ++i) {
			_workers[i]->stop_performance_measurements();
		}
		for (uint32_t i=0; i<_size; ++i) {
			_classifiers[i]->stop_performance_measurement();
		}
		clock_gettime(CLOCK_MONOTONIC, &end_time);
	}

	/**
	 * @brief clones this to another instance
	 */
	virtual GenericClassifier* clone() {
		throw std::runtime_error("Clone is not supported for ParallelClassifier");
	}


	/**
	 * @brief Prints debug information
	 * @param verbose Set the verbosity level of printing
	 */
	virtual void print(uint32_t verbose=1) const {

		// Measure performance
		double total_usec = (double)((end_time.tv_sec * 1e9 + end_time.tv_nsec) -
							  (start_time.tv_sec * 1e9 + start_time.tv_nsec)) / 1e3;
		messagef("Performance: total time %.3lf usec. Average time: %.3lf usec per packet.", total_usec, total_usec / _packet_counter);

		messagef("Serial classifier avg time per batch: %.3lf us", total_usec / _total_serial_packets);

		for (uint32_t i=0; i<_size-1; ++i) {
			messagef("Classifier %u statistics: utilization: %.2lf%%, throughput: %.2lf rpus, "
					"backpressure: %.2lf rpus, avg time per batch: %.2lf us", i+1,
					_workers[i]->get_utilization(), _workers[i]->get_throughput(),
					_workers[i]->get_backpressure(), _workers[i]->get_average_work_time());
		}
	}

	/**
	 * @brief Resets the all classifier counters
	 */
	virtual void reset_counters() {
		_packet_counter = 0;
		_next_batch.first_packet_counter = 0;
		_total_serial_packets = 0;
	}

	/**
	 * @brief Start an asynchronous process of classification for an input packet.
	 * @param header An array of 32bit integers according to the number of supported fields.
	 * @returns A unique id for the packet
	 */
	virtual unsigned int classify_async(const unsigned int* header, int priority) {
		// Default packet counter is invalid packet
		uint32_t packet_counter = 0xffffffff;
		// Feed next batch
		if (header != NULL) {
			_next_batch.packets[_next_batch_size++] = header;
			packet_counter = _packet_counter++;
		}

		// Fire next batch
		if (_next_batch_size == N || header == NULL) {
			// Fill invalid packets with NULL headers
			for (uint32_t i=_next_batch_size; i<N; ++i) {
				_next_batch.packets[i] = NULL;
			}

			// Produce batch with next available parallel worker
			bool produced = false;
			for (uint32_t i=0; i<_size-1; ++i) {
				_next_batch.worker_id = i+1;
				if (_workers[i]->produce(_next_batch)) {
					produced = true;
					break;
				}
			}

			// In case no parallel worker available - do as serial work
			if (!produced) {
				_next_batch.worker_id = 0;
				++_total_serial_packets;
				ParallelClassifier<N>::worker_method(_next_batch, this);
			}

			_next_batch_size = 0;

			// Set the counter value for the first packet
			_next_batch.first_packet_counter = _packet_counter;
		}

		// Return packet counter
		return packet_counter;
	}

	/**
	 * @brief Start a synchronous process of classification an input packet.
	 * @param header An array of 32bit integers according to the number of supported fields.
	 * @returns The matching rule action/priority (or 0xffffffff if not found)
	 */
	uint32_t classify_sync(const uint32_t* header, int priority) {
		throw std::runtime_error("NuevoMatch parallel classifier can only perform asynchronous classification");
	}


	/**
	 * @brief Returns a string representation of this
	 */
	const std::string to_string() const { return "NuevoMatchParallelClassifier"; }

	/**
	 * @brief Returns the maximum supported number of fields this can classify
	 */
	const unsigned int get_supported_number_of_fields() const { return UINT_MAX; }

};

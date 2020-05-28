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

#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <exception>
#include <stdexcept>

#include <basic_types.h>
#include <pipeline_thread.h>
#include <cpu_core_tools.h>
#include <logging.h>

template<typename T, uint32_t N>
struct WorkBatch {
	T items[N];

	/**
	 * @brief Fast getters and setters to the items of this
	 */
	inline T& operator[] (uint32_t idx) { return items[idx]; }
	inline const T& operator[] (uint32_t idx) const { return items[idx]; }
};

/**
 * @brief A thread that takes place in a chin of pipeline events
 * @tparam T The consumed job type of this thread
 */
template<typename T>
class PipelineThread {
public:

	/**
	 * @brief The thread worker method type
	 * @param job A job to consume
	 * @param args Any additional arguments (initiated with the constructor of PipelineThread)
	 * @returns True iff the job was consumed
	 */
	typedef bool(*thread_worker_method_t)(T& job, void* args);

private:

	// The worker thread
	pthread_t _thread;

	// Running flag
	volatile enum { READY, RUNNING, STOP, ERROR } _state;

	// Thread input jobs
	T* _ringbuffer;
	uint32_t _size;
	volatile uint32_t _cursor_read, _cursor_write;

	// Used for measuring performance
	double _idle_time_nsec, _start_time_nsec, _stop_time_nsec;
	double _work_time_nsec;
	double _request_handeled, _requests_to_measure;
	double _declined_requests;

	bool _measure_performance;

	// The worker method
	thread_worker_method_t _worker_method;
	void* _worker_method_args;

	/**
	 * @brief The worker thread main loop
	 */
	static void* worker_start(void* args) {
		// Initiate state
		PipelineThread* instance = (PipelineThread*)args;
		if (instance->_state == READY) {
			instance->_state = RUNNING;
		}

		// Measure time
		struct timespec start_time, end_time;
		clock_gettime(CLOCK_MONOTONIC, &end_time);
		bool measure = false, last_iteration_measure = false;;

		try {
			// Main loop
			while (instance->_state == RUNNING) {

				// Measure performance for current iteration
				last_iteration_measure = measure;
				measure = instance->_measure_performance;

				// Measure time only when requred
				if (measure) {
					clock_gettime(CLOCK_MONOTONIC, &start_time);
				}

				// Measure the work time only if the last iteration was measured
				if (last_iteration_measure) {
					instance->_work_time_nsec +=
							(start_time.tv_sec * 1e9 + start_time.tv_nsec) -
							(end_time.tv_sec * 1e9 + end_time.tv_nsec);
				}

				// Wait for available job
				while ( instance->_state == RUNNING && (instance->_cursor_read == instance->_cursor_write) );

				// Measure time and update idle statistics only in case the current iteration was measured
				if (measure) {
					clock_gettime(CLOCK_MONOTONIC, &end_time);
					instance->_idle_time_nsec +=
							(end_time.tv_sec * 1e9 + end_time.tv_nsec) -
							(start_time.tv_sec * 1e9 + start_time.tv_nsec);

					// Statistics
					++instance->_request_handeled;
				}

				// Get next available job
				T& job = instance->_ringbuffer[instance->_cursor_read];

				// Wait until the job is consumed (backpressure)
				while ( instance->_state == RUNNING && !(instance->_worker_method(job, instance->_worker_method_args)) );

				// Update cursor (using fast modulo)
				instance->_cursor_read = (instance->_cursor_read+1) & (instance->_size - 1);
			}
		} catch (const std::exception& e) {
			std::cerr << "error in worker thread (" << instance->_thread << "): " << e.what() << std::endl;
		}

		// Stop this thread
		pthread_exit(NULL);
	}

public:

	/**
	 * @brief Initiate new pipeline thread
	 * @param size The size of the input ring-buffer. Must be a power of two!
	 * @param core_idx The core on which to run the thread, or -1 for auto management,
	 * @param worker_method The consumer method
	 * @param args Arguments to pass to the consumer method
	 */
	PipelineThread(uint32_t size, int core_idx, thread_worker_method_t worker_method, void* args) :
		_thread(0), _state(READY), _ringbuffer(nullptr),
		_size(size), _cursor_read(0), _cursor_write(0), _idle_time_nsec(0),
		_start_time_nsec(0), _stop_time_nsec(0), _work_time_nsec(0), _declined_requests(0),
		_measure_performance(false),
		_worker_method(worker_method), _worker_method_args(args)
	{

		// Validate size
		if (_size == 1) _size = 2;
		if (_size % 2 != 0) {
			throw std::runtime_error("queue size must be a power of two!");
		}

		// Allocate ring buffer
		_ringbuffer = new T[_size];

		// Initiate thread
		if(pthread_create(&_thread, NULL, worker_start, this)) {
			throw std::runtime_error(std::string("cannot initiate thread: ") + strerror(errno));
		}

		// Migrate thread to CPU
		if (core_idx >= 0) {
			try {
				cpu_core_tools_set_thread_affinity(_thread, (int)core_idx);
			} catch (const std::exception& e) {
				// Log error
				throw std::runtime_error(std::string("error while migrating thread to CPU ") + std::to_string(core_idx) + e.what());
			}
		}

	}

	~PipelineThread() {
		stop();
		delete[] _ringbuffer;
	}

	/**
	 * @brief Produces new job for this to consume
	 * @param job The job to produce
	 * @returns True iff the job was produced
	 */
	bool produce(T job) {
		// No slot available
		if (!available()) {
			return false;
		}
		// Update buffer
		_ringbuffer[_cursor_write] = job;
		_cursor_write = (_cursor_write+1)%_size;
		return true;
	}

	/**
	 * @brief Stops the worker thread
	 * @throws In case the thread stopped with error
	 */
	void stop() {
		if (_state == ERROR) {
			throw std::runtime_error("the thread stopped with internal error");
		}
		// Wait for worker thread
		if (_state == RUNNING) {
			_state = STOP;
			pthread_join(_thread, NULL);
		}
	}

	/**
	 * @brief Start performance measurements
	 */
	void start_performance_measurements() {
		struct timespec t;
		clock_gettime(CLOCK_MONOTONIC, &t);
		_request_handeled = 0;
		_start_time_nsec = (double)t.tv_sec * 1e9 + t.tv_nsec;
		_declined_requests = 0;
		_measure_performance = true;
		_work_time_nsec = 0;
	}

	/**
	 * @brief Stop performance measurements
	 */
	void stop_performance_measurements() {
		struct timespec t;
		clock_gettime(CLOCK_MONOTONIC, &t);
		_requests_to_measure = _request_handeled;
		_stop_time_nsec = (double)t.tv_sec * 1e9 + t.tv_nsec;
		_measure_performance = false;
	}

	/**
	 * @brief Returns the utilization percent
	 */
	double get_utilization() const { return 100.0 - (_idle_time_nsec / (_idle_time_nsec+_work_time_nsec) * 100); }

	/**
	 * @brief Returns the throughput of this (requests per us)
	 */
	double get_throughput() const { return _requests_to_measure / (_stop_time_nsec-_start_time_nsec - _idle_time_nsec) * 1e3; }

	/**
	 * @brief Returns the ratio of declined requests (requests per us). Used to measure backpressure.
	 */
	double get_backpressure() const { return _declined_requests / (_stop_time_nsec-_start_time_nsec) * 1e3; }

	/**
	 * @brief Returns the average work time per request (in us).
	 */
	double get_average_work_time() const { return (double)_work_time_nsec / _requests_to_measure / 1e3; }

	/**
	 * @brief Checks whether this has available slot
	 * @returns True iff a new job can be produced for this
	 * @note Updates backpressure statistics
	 */
	inline bool available() {
		bool result = ((_cursor_write+1) & (_size - 1)) != _cursor_read;
		if (!result) ++_declined_requests;
		return result;
	};
};

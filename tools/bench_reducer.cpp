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

#include <stdlib.h>
#include <time.h>

#include <basic_types.h>
#include <pipeline_thread.h>
#include <cpu_core_tools.h>
#include <argument_handler.h>
#include <logging.h>

// Holds arguments information
static argument_t my_arguments[] = {
		// Name,	Required,	IsBoolean,	Default,	Help
		{"--jobs",			0,			0,			"1024",		"Number of jobs to generate"},
		{"--queue-size",	0,			0,			"4",		"Queue size per thread"},
		{"--threads",		0,			0,			"2",		"Number of threads to test"},
		{NULL,				0,			0,			NULL,		"Benchmark for reducer thread. Used for testing ideal setting."} /* Sentinel */
};

// Used as flag type
typedef union {
	//uint8_t __align__[16];
	volatile uint32_t value;
} cache_line_t;


// Global variables
uint32_t num_of_threads, queue_size, num_of_jobs;

// Shared variables
cache_line_t *num_of_results;
cache_line_t *flags;

/**
 * @brief The method for the worker threads
 */
bool worker_func(uint32_t& job_id, void* args) {
	int thread_id = (long)args;
	++num_of_results[thread_id].value;
	// Update the flag to notify the reducer
	uint32_t idx = num_of_threads*(job_id % queue_size)+thread_id;
	flags[idx].value = 1;
	return true;
}

/**
 * @brief The method for the reducer thread
 */
bool reducer_func(uint32_t& job_id, void* args) {
	uint32_t counter = 0;
	// Wait for all workers to finish
	while (counter < num_of_threads) {
		// Go over all flags of the current job..
		for (uint32_t i=0; i<num_of_threads; ++i) {
			uint32_t idx = num_of_threads*(job_id % queue_size)+i;
			// Check whether the flag of the current job was set
			if (flags[idx].value == 1) {
				flags[idx].value = 0;
				++counter;
			}
		}
	}
	return true;
}

/**
 * @brief Entry point
 */
int main(int argc, char** argv) {

	// Parse arguments
	parse_arguments(argc, argv, my_arguments);

	// Set variables
	queue_size = atoi(ARG("--queue-size")->value);
	num_of_threads = atoi(ARG("--threads")->value);
	num_of_jobs = atoi(ARG("--jobs")->value);

	flags = new cache_line_t[num_of_threads*queue_size];
	num_of_results = new cache_line_t[num_of_threads];

	// Allocate and align
	/*
	if (posix_memalign((void**)&flags, sizeof(cache_line_t), sizeof(cache_line_t)*num_of_threads*queue_size)) {
		throw std::runtime_error("cannot allocate aligned memory");
	}

	if (posix_memalign((void**)&num_of_results, sizeof(cache_line_t), sizeof(cache_line_t)*num_of_threads)) {
		throw std::runtime_error("cannot allocate aligned memory");
	}*/


	messagef("Allocated flags on %p. Size of flag_t: %lu", flags, sizeof(cache_line_t));

	// Initiate all data
	messagef("Initiate all data");
	for (uint32_t i=0; i<queue_size; ++i) {
		for (uint32_t j=0; j<num_of_threads; ++j) {
			flags[i*num_of_threads+j].value=0;
		}
	}
	for (uint32_t j=0; j<num_of_threads; ++j) {
		num_of_results[j].value=0;
	}


	// Initiate pipeline-threads
	messagef("Initiate worker threads");
	PipelineThread<uint32_t> *threads[num_of_threads];
	for (uint32_t i=0; i<num_of_threads; ++i) {
		threads[i] = new PipelineThread<uint32_t>(queue_size, 2*(i+1), worker_func, (void*)(long)i);
	}

	// Starting reducer
	messagef("Initiate reducer thread");
	PipelineThread<uint32_t> reducer = PipelineThread<uint32_t>(queue_size, 0, reducer_func, nullptr);

	messagef("Starting benchmark");

	for (uint32_t i=0; i<num_of_threads; ++i) {
		threads[i]->start_performance_measurements();
	}
	reducer.start_performance_measurements();
	struct timespec start_time, end_time;
	clock_gettime(CLOCK_MONOTONIC, &start_time);

	// Produce all jobs
	for (uint32_t i=0; i<num_of_jobs; ++i) {

		// Produce item for the reducer
		while(!reducer.produce(i));

		// Produce all workers
		for (uint32_t j=0; j<num_of_threads; ++j) {
			while(!threads[j]->produce(i));
		}
	}

	// Wait for all jobs to complete
	for (uint32_t j=0; j<num_of_threads; ++j) {
		while(num_of_results[j].value<num_of_jobs);
	}
	clock_gettime(CLOCK_MONOTONIC, &end_time);
	for (uint32_t i=0; i<num_of_threads; ++i) {
			threads[i]->stop_performance_measurements();
	}
	reducer.stop_performance_measurements();

	// Print statistics
	double total_us = (end_time.tv_sec - start_time.tv_sec) * 1e6 + (double)(end_time.tv_nsec - start_time.tv_nsec) / 1e3;
	messagef("Total time: %.3lf us, average time: %.3lf us", total_us, total_us/num_of_jobs);

	for (uint32_t i=0; i<num_of_threads; ++i) {
		messagef("Thread %u statistics: utilization: %.2lf%%, throughput: %.2lf rpus, "
				"backpressure: %.2lf rpus, avg time per batch: %.2lf us", i,
				threads[i]->get_utilization(), threads[i]->get_throughput(),
				threads[i]->get_backpressure(), threads[i]->get_average_work_time());
	}

	messagef("Reducer statistics: utilization: %.2lf%%, throughput: %.2lf rpus, "
			"backpressure: %.2lf rpus, avg time per batch: %.2lf us",
			reducer.get_utilization(), reducer.get_throughput(),
			reducer.get_backpressure(), reducer.get_average_work_time());
}

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

/**
 * This micro-benchmark is used to check the RQRMI inference.
 * Use it to load a model from the file-system, and generate samples for inference
 */

#include <time.h>
#include <logging.h>
#include <argument_handler.h>
#include <basic_types.h>
#include <pipeline_thread.h>
#include <cpu_core_tools.h>

// The work batch size
#define BATCH 8
typedef WorkBatch<scalar_t, BATCH> batch_t;

bool* initiated;
bool worker(batch_t input, void* args);
void measure_clock(uint32_t samples);


// Holds arguments information
static argument_t my_arguments[] = {
		// Name,	Required,	IsBoolean,	Default,	Help
		{"-n",			1,			0,		NULL,		"Number of samples"},
		{"--single",	0,			1,		NULL,		"Run both threads on a single core"},
		{NULL,			0,			0,		NULL,		"Core-echo benchmark tool. "} /* Sentinel */
};

/**
 * @brief The thread worker method
 */
bool worker(batch_t& input, void* args) {
	for (int i=0; i<BATCH; ++i) {
		initiated[(int)input.items[i]] = true;
	}
	return true;
}

/**
 * @brief Entry point
 */
int main(int argc, char** argv) {

	// Initialize library with print to stderr
	SimpleLogger::get().set_sticky_force(true);

	// Parse arguments
	parse_arguments(argc, argv, my_arguments);
	uint32_t num_of_samples = atoi( my_arguments[0].value );

	// Measure time for clock syscall
	measure_clock(num_of_samples);

	// Get next available core
	int curent_core_idx = cpu_core_tools_get_index_of_current_thread();
	int worker_core_idx = cpu_core_tools_get_next_physical_core(curent_core_idx);

	// Run both threads on a single core?
	if ( my_arguments[1].available ) worker_core_idx = curent_core_idx;

	message_s("Compiled to use batched jobs of " << BATCH << ". To alter batch size, re-set #define");

	message_s("Main function on core " << curent_core_idx << ", worker thread on core " << worker_core_idx);
	PipelineThread<batch_t> worker_thread(1, worker_core_idx, worker, nullptr);

	// Initiate the array
	initiated = new bool[(uint32_t)num_of_samples];
	for (uint32_t i=0; i<num_of_samples; ++i) {
		initiated[i] = false;
	}

	message_s("Starting echo benchmark within two cores...");
	struct timespec start_time, end_time;

	worker_thread.start_performance_measurements();
	clock_gettime(CLOCK_MONOTONIC, &start_time);

	batch_t job;
	// First handle the remainder
	uint32_t rem = num_of_samples % BATCH;
	for (uint32_t i=0; i<rem; ++i) {
		job[i] = i;
	}
	for (uint32_t i=rem; i<BATCH; ++i) {
		job[i] = 0; // Dummy
	}
	// Produce batch
	while (!worker_thread.produce(job));

	// Handle all the rest
	for (uint32_t i=rem; i<num_of_samples; i+=BATCH) {
		// Build job
		for (uint32_t j=0; j<BATCH; ++j) {
			job[j] = rem+i+j;
		}
		// Produce batch
		while (!worker_thread.produce(job));
	}

	clock_gettime(CLOCK_MONOTONIC, &end_time);
	worker_thread.stop_performance_measurements();

	// Validate
	for (uint32_t i=0; i<num_of_samples; ++i) {
		if (initiated[i] == false) {
			throw error("not all items were produced");
		}
	}

	// Print statistics
	double total_time = (end_time.tv_sec - start_time.tv_sec) * 1e6 + (double)(end_time.tv_nsec - start_time.tv_nsec) / 1e3;
	double avg_time = total_time/ num_of_samples;

	message_s("Total time: " << total_time << " us. Average time: " << avg_time << " us.");
	message_s("Worker utilization: " << worker_thread.get_utilization() << "%, throughput: " <<
			worker_thread.get_throughput() << " rpus, backpressure: " <<
			worker_thread.get_backpressure() << "rpus, worker avg time for batch: " <<
			worker_thread.get_average_work_time() << " us");

	delete initiated;
}

/**
 * @brief Measure hot many micro-seconds the clock method requires
 */
void measure_clock(uint32_t samples) {
	struct timespec start_time, end_time, test;

	clock_gettime(CLOCK_MONOTONIC, &start_time);
	for (uint32_t i=0; i<samples; ++i) {
		clock_gettime(CLOCK_MONOTONIC, &test);
	}
	clock_gettime(CLOCK_MONOTONIC, &end_time);

	double total_time = (end_time.tv_sec - start_time.tv_sec) * 1e6 + (double)(end_time.tv_nsec - start_time.tv_nsec) / 1e3;
	double avg_time = total_time/ samples;
	message_s("Clock Measure: Total time: " << total_time << " us. Average time: " << avg_time << " us.");
}

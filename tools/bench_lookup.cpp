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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <time.h>

#include <logging.h>
#include <argument_handler.h>
#include <object_io.h>
#include <matrix_operations.h>
#include <algorithms.h>
#include <lookup.h>

using namespace std;

// Define batch size for lookups
#define BATCH 128

scalar_t *samples;
map<scalar_t, typename LookupCPU<BATCH>::result_t> expected_results;
volatile uint32_t num_of_results = 0;


LookupCPU<BATCH>::result_t *results;

void generate_dataset(uint32_t num_of_samples, LookupCPU<BATCH>* lookup);
void perform_lookup(uint32_t num_of_samples, Lookup<BATCH>* lookup);

// Holds arguments information
static argument_t my_arguments[] = {
		// Name,	Required,	IsBoolean,	Default,	Help
		{"-f",		1,			0,			NULL,		"Lookup data filename to load"},
		{"-s",		0,			0,			"256",		"Lookup queue size"},
		{"-n",		1,			0,			NULL,		"Number of random generated samples"},
		{"-r",		0,			0,			NULL,		"Number of repeats"},
		{NULL,		0,			0,			NULL,		"Lookup benchmark tool."} /* Sentinel */
};


/**
 * @brief Callback object for lookup results
 */
class MyListener : public LookupListener {
protected:
	void on_new_result(scalar_t input, uint32_t index, int found) {
		if (input < 0) return;
		results[num_of_results].found = found;
		results[num_of_results].input = input;
		results[num_of_results].index = index;
		++num_of_results;
	}
};


int main(int argc, char** argv) {

	// Print message buffer to stderr
	SimpleLogger::get().set_sticky_force(true);

	// Parse arguments
	parse_arguments(argc, argv, my_arguments);
	const char* filename=my_arguments[0].value;
	int queue_size=atoi(my_arguments[1].value);
	uint32_t num_of_samples=atoi(my_arguments[2].value);
	int repeats = my_arguments[3].available ? atoi(my_arguments[3].value) : 1;

	// Open the file
	ObjectReader handler(filename);

	// Create a new instance
	LookupCPU<BATCH>* lookup = new LookupCPU<BATCH>();
	lookup->set_queue_size(queue_size);
	lookup->load(handler);

	// Register callback
	MyListener listener;
	lookup->add_listener(listener);

	// Measure simulation
	struct timespec start_time, end_time;

	// Catch model evaluation errors
	try {

		// Generate the samples
		generate_dataset(num_of_samples, lookup);

		// Warm cache (is it necessary?)
		message_s("Warming cache...");
		for (int j=0; j<repeats*3; ++j) {
			perform_lookup(num_of_samples, lookup);
		}

		// Simulate until number of packets
		message_s("Starting simulation (batching: " << BATCH << ")...");
		lookup->start_performance_measurement();
		clock_gettime(CLOCK_MONOTONIC, &start_time);
		// Repeat R times
		for(int r=0; r<repeats; r++) {
			perform_lookup(num_of_samples, lookup);
		}
		clock_gettime(CLOCK_MONOTONIC, &end_time);
		lookup->stop_performance_measurement();

	} catch (const std::exception& e) {
		warning(e.what());
		throw error("Lookup failed. Compile with DEBUG flag for extended info. Exiting");
	}

	double total_clock_us = (end_time.tv_sec * 1e9 + end_time.tv_nsec - start_time.tv_sec * 1e9 - start_time.tv_nsec)/1000;

	long total_samples=num_of_samples;

	message_s("Total time to simulate " << num_of_samples << " samples: " << total_clock_us <<
			" us. Average time per sample: " << total_clock_us/total_samples/repeats << " us.");

	// Print worker statistics of lookup
	lookup->print();

	// Check correctness
	message_s("Validating correctness...");
	for (uint32_t i=0; i<num_of_results; ++i) {
		// Search for current result
		auto it = expected_results.find(results[i].input);
		if (it == expected_results.end() && results[i].found == 1) {
			message_s("Error with input: sample " << results[i].input << " was found in lookup but not in database");
		} else if (it != expected_results.end() && results[i].found == 0) {
			message_s("Error with input: sample " << results[i].input << " was not found in lookup but was in database");
		} else if (it != expected_results.end() && results[i].found) {
			if (results[i].index != it->second.index) {
				scalar_pair_t actual_result = lookup->get_record(results[i].index);
				scalar_pair_t expected_result = lookup->get_record(it->second.index);
				message_s("Error with input " << results[i].input << ". Actual index: " <<
						results[i].index << " [" <<  actual_result.first << "," << actual_result.second <<
						"). Expected index: %u [" << expected_result.first << "," << expected_result.second << ")");
			}
		}
	}

	// Free memory
	message_s("Done");
	delete lookup;
}

/**
 * @brief Generates database for benchmark
 * @param num_of_samples Number of required samples
 * @param lookup The lookup object
 */
void generate_dataset(uint32_t num_of_samples, LookupCPU<BATCH>* lookup) {

	message_s("Generating " << num_of_samples << " samples for database with " << lookup->get_size() << " pairs...");

	// Generate the samples DB
	samples = (scalar_t*)malloc(sizeof(scalar_t) * num_of_samples);
	results = (LookupCPU<BATCH>::result_t*)malloc(sizeof(LookupCPU<BATCH>::result_t) * num_of_samples);
	if (samples == NULL || results == NULL) {
		throw error("Cannot allocate memory for generated samples. Exiting.");
	}

	// Initialize randomization
	srand(time(NULL));

	uint32_t counter=0;
	for (uint32_t i=0; i<num_of_samples; ++i) {

		scalar_pair_t range = lookup->get_record(counter);
		scalar_t range_start = range.first;
		scalar_t range_end = SCALAR_PREV(range.second); //inclusive
		counter = (counter+1) % lookup->get_size();

		// Get random key between start and end
		scalar_t key = gen_uniform_random_scalar(range_start, range_end);
		samples[i] = key;
		expected_results[key].index = i;
	}
}

/**
 * @brief Performs the lookup
 */
void perform_lookup(uint32_t num_of_samples, Lookup<BATCH>* lookup) {
	Lookup<BATCH>::batch_t job;

	// First handle the remainder
	uint32_t rem = num_of_samples % BATCH;
	for (uint32_t i=0; i<rem; ++i) {
		job[i] = i;
	}
	for (uint32_t i=rem; i<BATCH; ++i) {
		job[i] = -1; // Invalid
	}
	while (!lookup->search(job));

	// Handle all the rest
	for (uint32_t i=rem; i<num_of_samples; i+=BATCH) {
		// Build job
		for (uint32_t j=0; j<BATCH; ++j) {
			job[j] = rem+i+j;
		}
		// Perform lookup
		while (!lookup->search(job));
	}

	// Wait for all results
	while (num_of_results < num_of_samples) {
		continue;
	}

	num_of_results=0;
}

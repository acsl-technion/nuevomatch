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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <algorithms.h>
#include <logging.h>
#include <object_io.h>
#include <matrix_operations.h>
#include <rqrmi_model.h>
#include <rqrmi_fast.h>
#include <argument_handler.h>

void generate_dataset(int num_of_samples);

// Holds information regarding the simulation
int num_of_samples = 0;
bool use_fast = false;
scalar_t* samples = NULL;
scalar_t* results = NULL;

// Holds arguments information
static argument_t my_arguments[] = {
		// Name,		Required,	IsBoolean,	Default,	Help
		{"-f",			1,			0,			NULL,		"RQRMI filename to load"},
		{"-n",			1,			0,			NULL,		"Number of random generated samples"},
		{"-p",			0,			1,			NULL,		"Store & print results at the end (alters performance measurements)"},

		/* Use RQRMIfast for performing calculation */
		{"--fast",		0,			1,			NULL,		"(Optimization) Set to true for using fast RQRMI evaluation. "
															"The SIMD engine should be set in compilation time. "
															"Works only for RQRMI models, not RMI."},
		{NULL,			0,			0,			NULL,		"RQRMI benchmark tool. Compile library with 'make release'. "} /* Sentinel */
};

/**
 * @brief Performs fast RQRMI inference using SIMD acceleration
 * @param rqrmi_fast An RQRMIFast instance
 */
void fast_evaluate(void* arg) {
	RQRMIFast& rqrmi_fast = *(RQRMIFast*)arg;
	wide_scalar_t inputs, outputs, status, error;
	// Perform fast evaluation using vector of inputs
	infof("Running batches from main samples...");
	int i;
	for (i=0; i<num_of_samples-(int)RQRMIFast::input_width(); i+=(int)RQRMIFast::input_width()) {
		for (uint32_t k=0; k<rqrmi_fast.input_width(); ++k) {
			inputs.scalars[k] = samples[i+k];
		}
		rqrmi_fast.evaluate(inputs, status, outputs, error);
		for (uint32_t k=0; k<rqrmi_fast.input_width(); ++k) {
			results[i+k] = outputs.scalars[k];
		}
	}
	// Handle the remainder
	infof("Running the remainder batches...");
	uint32_t counter = 0;
	int start = i;
	for (; i<num_of_samples; ++i) {
		inputs.scalars[counter++] = samples[i];
	}
	rqrmi_fast.evaluate(inputs, status, outputs, error);
	for (uint32_t k=0; k<counter; ++k) {
		results[start+k] = outputs.scalars[k];
	}
}

/**
 * @brief Performs slow RMI/RQRMI inference using generic engine
 * @param model A pointer to an RQRMI model
 */
void slow_evaluate(void* arg) {
	rqrmi_model_t* model = (rqrmi_model_t*)arg;
	for (int i=0; i<num_of_samples; ++i) {
		results[i] = rqrmi_evaluate_model(model, samples[i]);
	}
}

int main(int argc, char** argv) {

	// Parse arguments
	parse_arguments(argc, argv, my_arguments);
	const char* model_filename=ARG("-f")->value;
	num_of_samples=atoi( ARG("-n")->value );
	bool print_results = ARG("-p")->available;

	// Open the file
	ObjectReader reader(model_filename);

	// Load the model
	rqrmi_model_t* model = rqrmi_load_model(reader.buffer(), reader.size());
	if (model==RQRMI_MODEL_ERROR) {
		error("Cannot read model. Compile with DEBUG flag for extended info. Exiting");
		return 1;
	}

	// Measure simulation
	struct timespec start_time, end_time;

	// Do we use AVX, SSR or normal computation?
	use_fast = ARG("--fast")->available;
	RQRMIFast rqrmi_fast = RQRMIFast(model);

	// Catch model evaluation errors
	// Generate the samples
	loggerf("Generating %d samples...", num_of_samples);
	generate_dataset(num_of_samples);

	// Set inference method
	void(*inference_method)(void*);
	void* arg;
	if (use_fast) {
		inference_method = fast_evaluate;
		arg = &rqrmi_fast;
	} else {
		inference_method = slow_evaluate;
		arg = model;
	}

	// Warm cache (is it necessary?)
	logger("Warming cache...");
	for (int j=0; j<100000/num_of_samples; ++j) {
		inference_method(arg);
	}

	// Simulate until number of packets
	logger("Starting simulation...");
	clock_gettime(CLOCK_MONOTONIC, &start_time);
	inference_method(arg);
	clock_gettime(CLOCK_MONOTONIC, &end_time);


	double total_clock_us = (end_time.tv_sec * 1e9 + end_time.tv_nsec - start_time.tv_sec * 1e9 - start_time.tv_nsec)/1000;
	long total_packets=num_of_samples;

	loggerf("Total time to simulate %d samples: %.3lf us. "
		"Average time per sample: %.7lf us.",
		num_of_samples, total_clock_us, total_clock_us/total_packets);

	// Print results to stdout only if necessary
	if (print_results) {
		loggerf("Printing results to stdout");
		for (int i=0; i<num_of_samples; ++i) {
			printf("%f %f\n", samples[i], results[i]);
		}
	}

	// Free memory
	rqrmi_free_model(model);
}

void generate_dataset(int num_of_samples) {
	// Generate the samples DB
	samples = (scalar_t*)malloc(sizeof(scalar_t) * num_of_samples);
	results = (scalar_t*)malloc(sizeof(scalar_t) * num_of_samples);
	if (samples == NULL || results == NULL) {
		throw error("Cannot allocate memory for generated samples. Exiting.");
	}

	// Initialize randomization
	srand(time(NULL));

	for (int i=0; i<num_of_samples; ++i) {
		samples[i] = gen_uniform_random_uint32(0x0, 0xFFFFFFFF);
	}
}


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
#include <math.h>
#include <assert.h>
#include <string.h>

#include <basic_types.h>
#include <matrix_operations.h>
#include <rqrmi_tools.h>
#include <vector_list.h>

#include <algorithms.h>
#include <logging.h>
#include <rqrmi_model.h>

#define MAX_UINT 0xffffffff

// Tools information is often not required for debugging.
// Set dedicated debugging flag for models
#ifdef DEBUG_TOOLS
#  define tools_info(...) info(__VA_ARGS__)
#else
#  define tools_info(...)
#endif

// Used to store probing datasets for fast operations and cache responsibilities
struct rqrmi_probing {
	uint32_t  			num_of_stages;
	uint32_t  			num_of_records;
	uint32_t*  			stage_width;
	matrix_t*			records;
	vector_list_t***	responsibilities;
	vector_list_t**		transition_sets;
};

// Private methods
void print_responsibility(rqrmi_probing_t* probe, uint32_t stage_idx, uint32_t bucket_index);
vector_list_t* rqrmi_tools_get_records_in_responsibility(rqrmi_probing_t* probe, uint32_t stage_idx, uint32_t bucket_idx);

/**
 * @brief Generates an RQRMI probing data structure.
 * @param records The records the RQRMI index. Nx1 matrix. The key values / range start values. Should be sorted.
 * @param stage_num The number of stages in the RQRMI model
 * @param stages_width An array of integers with the width of each stage
 * @note On error returns RQRMI_PROBING_ERROR
 */
rqrmi_probing_t* rqrmi_tools_probe_new(matrix_t* records, uint32_t stage_num, uint32_t* stages_width) {

	// Check inputs
	if (stage_num < 1 || stages_width[0] != 1) {
		warning("bad inputs");
		return RQRMI_PROBING_ERROR;
	}

	// Allocate output
	rqrmi_probing_t* output = (rqrmi_probing_t*)malloc(sizeof(rqrmi_probing_t));
	if (output== NULL) {
		warning("Cannot allocate memory for rqrmi_probing_t");
		return RQRMI_PROBING_ERROR;
	}

	try {
		// Allocate responsibility caches
		output->num_of_stages = stage_num;
		output->num_of_records = records->rows;
		output->records = records;
		output->responsibilities=(vector_list_t***)malloc(sizeof(vector_list_t**)*stage_num);
		output->transition_sets=(vector_list_t**)malloc(sizeof(vector_list_t*)*stage_num);
		output->stage_width=(uint32_t*)malloc(sizeof(uint32_t)*stage_num);

		if (output->responsibilities == NULL || output->transition_sets == NULL || output->stage_width == NULL) {
			throw error("Cannot allocate memory for probe caches");
		}

		// Reset all points
		for (uint32_t i=0; i<stage_num; ++i) {
			output->responsibilities[i] = NULL;
			output->transition_sets[i] = VECTOR_LIST_ERROR;
			output->stage_width[i] = stages_width[i];
		}

		// Allocate responsibilities
		for (uint32_t i=0; i<stage_num; ++i) {
			output->responsibilities[i] =(vector_list_t**)malloc(sizeof(vector_list_t*)*stages_width[i]);
			if (!output->responsibilities[i]) {
				throw error("cannot allocate memory for responsibilities");
			}
			for (uint32_t j=0; j<stages_width[i]; ++j) {
				output->responsibilities[i][j]=VECTOR_LIST_ERROR;
			}
		}

		// Create responsibility of stage 0 (always true by definition)
		output->responsibilities[0][0] = vector_list_create(2);
		if (output->responsibilities[0][0] == VECTOR_LIST_ERROR) {
			throw error("Cannot allocate memory for responsibility of stage 0");
		}

		// Add the entire input domain to the responsibility
		scalar_t* interval = (scalar_t*)vector_list_push_back_and_get(output->responsibilities[0][0]);
		interval[0] = GET_SCALAR(records, 0 ,0);
		interval[1] = GET_SCALAR(records, output->num_of_records-1, 0);
		info("Responsibility R<0,0> is the input domain: [" << interval[0] << "," << interval[1] << "]");

	} catch (const std::exception& e) {
		warning(e.what());
		for (uint32_t i=0; i<stage_num; ++i) {
			free(output->responsibilities[i]);
		}
		free(output->responsibilities);
		free(output->transition_sets);
		free(output->stage_width);
		free(output);
		output = RQRMI_PROBING_ERROR;
	}

	return output;
}

/**
 * @brief Frees the allocated memory of rqrmi_probing_t structure
 */
void rqrmi_tools_probe_free(rqrmi_probing_t* probe) {
	if (probe == RQRMI_PROBING_ERROR) return;
	// Free responsibilities caches
	for (uint32_t i=0; i<probe->num_of_stages; ++i) {
		for (uint32_t j=0; j<probe->stage_width[i]; ++j) {
			vector_list_free(probe->responsibilities[i][j]);
		}
		free(probe->responsibilities[i]);
		vector_list_free(probe->transition_sets[i]);
	}
	// Free  pointers
	free(probe->stage_width);
	free(probe->responsibilities);
	free(probe->transition_sets);
	free(probe);
}

/**
 * @brief Calculate the transition set of a stage
 * @param model An RQRMI model
 * @param probe An RQRMI probing data structure
 * @param stage_idx The stage of the submodel
 * @returns A pointer to the transition set. VECTOR_LIST_ERROR on error.
 * @note  The user should not free the transition set
 * @note  The responsibility of stage stage_idx must be defined prior to this method
 * @note  The transition set is an ordered vector list and stored at the probe data structure
 * 		  Vector format: [ x,  B_i( M_i(x-epsilon) ),  B_i( M_i(x+epsilon) ) ]
 */
vector_list_t* rqrmi_tools_calculate_transition_set(rqrmi_model_t* model, rqrmi_probing_t* probe, uint32_t stage_idx) {

	vector_list_t* transition_inputs = VECTOR_LIST_ERROR;

	info("Calculating U_" << stage_idx << "...");

	try {

		// Clear previous allocated memory
		if (probe->transition_sets[stage_idx] != VECTOR_LIST_ERROR) {
			vector_list_free(probe->transition_sets[stage_idx]);
		}

		// Create the transition set vector list
		// Format: [ x,  B(M(x-eps)),  B(M(x+eps)) ]
		probe->transition_sets[stage_idx] = vector_list_create(3);
		if (probe->transition_sets[stage_idx] == VECTOR_LIST_ERROR) {
			throw error("Cannot create memory for transition set of stage " << stage_idx);
		}

		// The next stage width
		uint32_t next_width = stage_idx+1 == probe->num_of_stages ? probe->num_of_records : probe->stage_width[stage_idx+1];
		scalar_pair_t input_domain = rqrmi_get_input_domain(model);

		// Insert the input domain minimum to the transition set
		scalar_t S = rqrmi_evaluate_model(model, input_domain.first);
		uint32_t B = next_width * (S < 0 ? 0 : S >= 1 ? 1 - SCALAR_EPS : S);
		scalar_t* new_pt = (scalar_t*)vector_list_push_back_and_get(probe->transition_sets[stage_idx]);
		new_pt[0] = input_domain.first;
		new_pt[1] = B;
		new_pt[2] = B;

		// Store last bucket value for stage output discontinuity inputs
		uint32_t last_bucket = B;

		// Calculate transition inputs for each submodel
		for (uint32_t i=0; i<probe->stage_width[stage_idx]; ++i) {

			// Skip non compiled submodels
			if (!rqrmi_submodel_compiled(model, stage_idx, i)) {
				info("Skipping T<" << stage_idx << "," << i << "> since has empty responsibility");
				continue;
			}

			// Calculate the transition inputs of the current submodel
			transition_inputs = rqrmi_calculate_transition_inputs(model, stage_idx, i, next_width);
			if (transition_inputs == MATRIX_ERROR) {
				throw error("Cannot calculate transition inputs");
			}

		    // Debug Print
#ifndef NDEBUG
				// Print only T_i,j of internal stages
				if (stage_idx < probe->num_of_stages-1) {
					SimpleLogger::get().format("T<%u,%u>: {", stage_idx, i);
					scalar_t* cursor = (scalar_t*)vector_list_begin(transition_inputs);
					uint32_t position = 0;
					for(; cursor; cursor = (scalar_t*)vector_list_iterate(transition_inputs)) {
						SimpleLogger::get().format("%f (%.0f, %.0f)", cursor[0], cursor[1], cursor[2]);
						if (position++ < vector_list_get_size(transition_inputs)-1) SimpleLogger::get().add(", ");
					}
					SimpleLogger::get().add("}\n");
				}
#endif

			// Go over the responsibility intervals
			vector_list_t* responsibility = probe->responsibilities[stage_idx][i];
			scalar_t* responsibility_interval = (scalar_t*)vector_list_begin(responsibility);
			for (; responsibility_interval; responsibility_interval = (scalar_t*)vector_list_iterate(responsibility)) {

				// Add interval edges to points
				scalar_t*  transition_input;

				// Go over all transition inputs
				transition_input = (scalar_t*)vector_list_begin(transition_inputs);
				for (; transition_input; transition_input = (scalar_t*)vector_list_iterate(transition_inputs)) {

					bool condition_left  = transition_input[0] >= responsibility_interval[0];
					bool condition_right = transition_input[0] <= responsibility_interval[1];

					// In both conditions are met, copy value to transition set
					if (condition_left && condition_right) {

						// Debug information
						/*info("Transition set U_%u contains x=%f (x in T<%u,%u>) as x in S=[%f, %f] subsetof R<%u,%u>. f(x-eps)=%u, f(x+eps)=%u",
								stage_idx, transition_input[0], stage_idx, i,
								responsibility_interval[0], responsibility_interval[1],
								stage_idx, i, transition_input[1], transition_input[2]);*/

						new_pt = (scalar_t*)vector_list_push_back_and_get(probe->transition_sets[stage_idx]);
						new_pt[0] = transition_input[0]; // x
						new_pt[1] = transition_input[1];
						new_pt[2] = transition_input[2];

						// Update last bucket value for discontinuity inputs
						last_bucket = transition_input[2];
					}
				}

				// Add the responsibility end point to the transition set, as the stage output
				// has discontinuity point
				new_pt = (scalar_t*)vector_list_push_back_and_get(probe->transition_sets[stage_idx]);
				new_pt[0] = responsibility_interval[1]; // x
				new_pt[1] = last_bucket;
				// The right bucket value is still unknown.
				new_pt[2] = (scalar_t)(-1);
			}

			// Free used memory
			vector_list_free(transition_inputs);
		}

		// Sort the transition set by value
		vector_list_sort(probe->transition_sets[stage_idx], 0, scalar_compare_asc);

		// Fill blank bucket values
		scalar_t* current = (scalar_t*)vector_list_begin(probe->transition_sets[stage_idx]);
		scalar_t* next = (scalar_t*)vector_list_iterate(probe->transition_sets[stage_idx]);
		while (next) {
			if (current[2] < 0) current[2] = next[1];
			current = next;
			next = (scalar_t*)vector_list_iterate(probe->transition_sets[stage_idx]);
		}
		if (current[2] < 0) current[2] = current[1];


	    // Debug Print
#ifndef NDEBUG
			// Print only U_i of internal stages
			if (stage_idx < probe->num_of_stages-1) {
				SimpleLogger::get().format("U_%u: {", stage_idx);
				scalar_t* cursor = (scalar_t*)vector_list_begin(probe->transition_sets[stage_idx]);
				uint32_t position = 0;
				for(; cursor; cursor = (scalar_t*)vector_list_iterate(probe->transition_sets[stage_idx])) {
					SimpleLogger::get().format("%f (%.0f, %.0f)", cursor[0], cursor[1], cursor[2]);
					if (position++ < vector_list_get_size(probe->transition_sets[stage_idx])-1) SimpleLogger::get().add(", ");
				}
				SimpleLogger::get().add("}\n");
			}
#endif

	} catch (const std::exception& e) {
		warning(e.what());
		vector_list_free(transition_inputs);
		vector_list_free(probe->transition_sets[stage_idx]);
	}
	return probe->transition_sets[stage_idx];
}

/**
 * @brief Calculates the responsibility of stage stage_idx
 * @param model An RQRMI model
 * @param probe An RQRMI probing data structure
 * @param stage_idx The stage of to calculate
 * @return A set of responsibilities, each is a list of intervals, or VECTOR_LIST_ERROR on error.
 * @note  User should not free the returned pointer
 * @note  First stage has width of 1
 * @note  The responsibilities of stages {0,1,...,x-1} must be calculated before that of stage x
 */
vector_list_t** rqrmi_tools_calculate_responsibility(rqrmi_model_t* model, rqrmi_probing_t* probe, uint32_t stage_idx) {

	info("Calculating responsibilities for stage " << stage_idx << "...");

	// Check inputs for first stage
	if (stage_idx == 0) {
		error("Responsibility of stage 0 is already available");
		return 0;
	}

	// Free previous resources
	if (probe->responsibilities[stage_idx] != VECTOR_LIST_ERROR) {
		for (uint32_t i=0; i<probe->stage_width[stage_idx]; ++i) {
			vector_list_free(probe->responsibilities[stage_idx][i]);
		}
	}

	// Check that the responsibility of the previous stage is defined
	if (probe->responsibilities[stage_idx-1] == VECTOR_LIST_ERROR) {
		error("cannot calculate responsibility of stage " << stage_idx
			  << "before the responsibility of the previous stage!");
		return 0;
	}

	// Calculate the transition set of the previous stage
	if (probe->transition_sets[stage_idx-1] == VECTOR_LIST_ERROR) {
		error("the transition set of stage " << stage_idx-1 << " was not calculated");
		return 0;
	}

	uint32_t stage_width = probe->stage_width[stage_idx];

	// Set return status
	vector_list_t** output = VECTOR_LIST_ERROR;

	// Temporarily set the RQRMI stages
	uint32_t orig_num_of_stages = rqrmi_get_num_of_stages(model);
	rqrmi_set_num_of_stages(model, stage_idx);

	try {
		// Allocate memory for all responsibilities in stage
		for (uint32_t i=0; i<stage_width; ++i) {
			probe->responsibilities[stage_idx][i] = vector_list_create(2);
			if (probe->responsibilities[stage_idx][i] == VECTOR_LIST_ERROR) {
				throw error("Cannot allocate memory for responsibility <" << stage_idx << "," << i << ">");
			}
		}

		scalar_t last_value = rqrmi_get_input_domain(model).first;

		// Go over all values in the transition set
		vector_list_t* U_i = probe->transition_sets[stage_idx-1];
		scalar_t* vec = (scalar_t*)vector_list_begin(U_i);
		for (; vec; vec = (scalar_t*)vector_list_iterate(U_i)) {

			// Get next transition point
			scalar_t pt = vec[0];

			// Get the bucket just before current point B(S(pt-epsilon))
			uint32_t bucket = vec[1];

			// Add new responsibility interval to the submodel corresponds to the last bucket
			scalar_t* new_interval = (scalar_t*)vector_list_push_back_and_get(probe->responsibilities[stage_idx][bucket]);

			// Check whether the current value is relevant to the bucket as well
			scalar_t M = rqrmi_evaluate_model(model, pt);
			uint32_t B = probe->stage_width[stage_idx] * (M < 0 ? 0 : M >= 1 ? 1 - SCALAR_EPS : M);

			// In case the input is not part of the bucket
			if (B != bucket) {
				new_interval[0] = last_value;
				new_interval[1] = SCALAR_PREV(pt);
				last_value = pt;
			}
			// In case the input as a part of the bucket
			else {
				new_interval[0] = last_value;
				new_interval[1] = pt;
				last_value = SCALAR_NEXT(pt);
			}
		}

		// Compress intersecting intervals
		for (uint32_t i=0; i<stage_width; ++i) {
			vector_list_t* responsibility = probe->responsibilities[stage_idx][i];
			// On each iteration, work on two intervals. Fix current with next and delete next if necessary.
			scalar_t* cursor = (scalar_t*)vector_list_begin(responsibility);
			scalar_t* next = (scalar_t*)vector_list_iterate(responsibility);
			uint32_t position = 1;
			while (next) {
				// Check whether the responsibility in cursor_next intersects cursor
				if (SCALAR_PREV(next[0]) <= cursor[1]) {
					// Merge cursor with cursor next
					cursor[1] = next[1];
					// Delete cursor_next
					vector_list_remove_at(responsibility, position);
				} else {
					// Continue to next cursor
					cursor = next;
				}
				next = (scalar_t*)vector_list_iterate(responsibility);
				++position;
			}
		}

		// Print debug information
#		ifndef 	NDEBUG
		for (uint32_t i=0; i<stage_width; ++i) {
			print_responsibility(probe, stage_idx, i);
		}
#		endif

		output = probe->responsibilities[stage_idx];
	} catch (const std::exception& e) {
		warning(e.what());
		output = VECTOR_LIST_ERROR;
	}

	// Return RQRMI to original state
	rqrmi_set_num_of_stages(model, orig_num_of_stages);
	return output;
}

/**
 * @brief Returns the records that match responsibility
 * @param probe An RQRMI probing data structure
 * @param stage_idx The required stage
 * @param bucket_idx The required bucket in stage
 * @returns A vector list. Format: [start_input, end_input, record_idx, interval_idx]. On error returns VECTOR_LIST_ERROR.
 * @note The user should free the record list
 */
vector_list_t* rqrmi_tools_get_records_in_responsibility(rqrmi_probing_t* probe, uint32_t stage_idx, uint32_t bucket_idx) {

	vector_list_t* responsibility;
	vector_list_t* record_list = VECTOR_LIST_ERROR;

	try {
		// Check bucket is valid
		if (bucket_idx >= probe->stage_width[stage_idx]) {
			throw error("bucket " << bucket_idx << " is not valid for stage " << stage_idx);
		}

		// Check responsibility was calculated
		if (probe->responsibilities[stage_idx][bucket_idx] == VECTOR_LIST_ERROR) {
			throw error("responsibility <" << stage_idx << "," << bucket_idx << "> is not available");
		}

		// Allocate list
		record_list = vector_list_create(4);
		if (record_list == VECTOR_LIST_ERROR) {
			throw error("cannot allocate record list");
		}

		// Get the responsibility of the bucket
		responsibility = probe->responsibilities[stage_idx][bucket_idx];

		scalar_t* interval = (scalar_t*)vector_list_begin(responsibility);

		// Go over all responsibility intervals
		for (uint32_t j=0; j<vector_list_get_size(responsibility); ++j) {

			scalar_t interval_srt = interval[0]; // inclusive
			scalar_t interval_end = interval[1]; // inclusive

			// Go over all records, match
			for (uint32_t r=0; r<probe->num_of_records-1; ++r) {

				scalar_t record_start = GET_SCALAR(probe->records, r, 0);
				scalar_t record_end   = GET_SCALAR(probe->records, r+1, 0);

				// In case the current record does not fit in the interval
				if (record_end < interval_srt) continue;
				else if (record_start > interval_end) break;

				// Add a marker
				scalar_t* cursor = (scalar_t*)vector_list_push_back_and_get(record_list);
				cursor[0] = MAX(record_start, interval_srt);
				cursor[1] = MIN(record_end, interval_end);
				cursor[1] = MAX(cursor[0], cursor[1]); // This is necessary since uint32 -> float has precision error
				cursor[2] = r;
				cursor[3] = j;
			}

			// Get next interval
			interval = (scalar_t*)vector_list_iterate(responsibility);
		}

	} catch (const std::exception& e) {
		warning(e.what());
		// Free resources
		vector_list_free(record_list);
		record_list = VECTOR_LIST_ERROR;
	}

	return record_list;
}

/**
 * @brief Generate dataset from stage bucket
 * @param probe An RQRMI probing data structure
 * @param stage_idx The required stage
 * @param bucket_idx The required bucket in stage
 * @param num_of_samples Dataset maximum size
 * @param random Smart random sampling
 * @param shuffle Shuffle dataset
 * @returns A matrix. Each row is a sample with format [input, expected_output], both are scalars not normalized!
 * 		    On error returns MATRIX_ERROR.
 */
matrix_t* rqrmi_tools_generate_dataset(rqrmi_probing_t* probe,
		uint32_t stage_idx, uint32_t bucket_idx, uint32_t num_of_samples, bool random, bool shuffle)
{
	// Check bucket is valid
	if (bucket_idx >= probe->stage_width[stage_idx]) {
		error("bucket " << bucket_idx << " is not valid for stage " << stage_idx);
		return MATRIX_ERROR;
	}

	// Check responsibility was calculated
	if (probe->responsibilities[stage_idx][bucket_idx] == VECTOR_LIST_ERROR) {
		error("responsibility <" << stage_idx << "," << bucket_idx << "> is not available");
		return MATRIX_ERROR;
	}

	// Get the records that match the responsibility
	vector_list_t* matching_records = rqrmi_tools_get_records_in_responsibility(probe, stage_idx, bucket_idx);
	if (matching_records == VECTOR_LIST_ERROR) {
		error("cannot calculate matching records");
		return MATRIX_ERROR;
	}

	// Allocate markers
	// Format: [ start_input, end_input, record_idx, responsibility_idx, num_of_samples ] (scalars)
	vector_list_t* markers = vector_list_create(5);
	if (markers == VECTOR_LIST_ERROR) {
		error("cannot allocate marker list");
		return MATRIX_ERROR;
	}

	matrix_t* output = MATRIX_ERROR;

	try {

		vector_list_t* responsibility = probe->responsibilities[stage_idx][bucket_idx];
		info("Generating dataset according to R<" << stage_idx << "," << bucket_idx << "> (which is here: " << responsibility);
		print_responsibility(probe, stage_idx, bucket_idx);


		// Calculate each responsibility coverage of record indices
		scalar_t total_space = 0;
		scalar_t responsibility_space[vector_list_get_size(responsibility)];
		scalar_t responsibility_markers[vector_list_get_size(responsibility)];
		int responsibility_idx = -1;
		uint32_t min_record_idx = MAX_UINT, max_record_idx = 0;

		// Go over all records and calculate coverage
		scalar_t* record = (scalar_t*)vector_list_begin(matching_records);
		for (; record; record = (scalar_t*)vector_list_iterate(matching_records)) {

			// Handle new responsibility
			if (record[3] > responsibility_idx) {

				// Close current responsibility
				if (responsibility_idx >= 0 && responsibility_markers[responsibility_idx] > 0) {
					assert(max_record_idx >= min_record_idx);
					responsibility_space[responsibility_idx] = (max_record_idx-min_record_idx+1);
					total_space += responsibility_space[responsibility_idx];
				}

				// Initialize new responsibility counters
				responsibility_idx = record[3];
				min_record_idx = MAX_UINT;
				max_record_idx = 0;
				responsibility_markers[responsibility_idx] = 0;
				responsibility_space[responsibility_idx] = 0;
			}

			// Update responsibility statistics
			min_record_idx = MIN(min_record_idx, record[2]);
			max_record_idx = MAX(max_record_idx, record[2]);
			++responsibility_markers[responsibility_idx];

			// Add a marker
			scalar_t* marker = (scalar_t*)vector_list_push_back_and_get(markers);
			for (int i=0; i<4; ++i) marker[i] = record[i];
			marker[4] = 0;
		}

		// Close current responsibility
		if (responsibility_idx >= 0 && responsibility_markers[responsibility_idx] > 0) {
			assert(max_record_idx >= min_record_idx);
			responsibility_space[responsibility_idx] = (max_record_idx-min_record_idx+1);
			total_space += responsibility_space[responsibility_idx];
		}

		info("Number of markers: " << vector_list_get_size(markers) <<
				". Total space: " << total_space);

		scalar_t* marker = (scalar_t*)vector_list_begin(markers);
		uint32_t actual_samples = 0;

		// Sample from each responsibility according to its size
		for (uint32_t j=0; j<vector_list_get_size(responsibility); ++j) {

			// Advance marker until reaching current responsibility
			while (marker && marker[3] < j) marker = (scalar_t*)vector_list_iterate(markers);
			if (marker == VECTOR_LIST_ERROR) break;

			// In case the next marker is in the next responsibility
			if (marker[3] > j) continue;
			assert(marker[3] == j);

			// How many samples are available to the current responsibility?
			uint32_t samples_for_responsibility = num_of_samples * responsibility_space[j] / total_space;

			// Skip empty responsibilities
			if (samples_for_responsibility == 0) continue;
			if (responsibility_markers[j] == 0) continue;

			// In case there are less markers than samples
			if (responsibility_markers[j] <= samples_for_responsibility) {
				// Divide the samples between the markers
				uint32_t samples_per_marker = SCALAR_NEXT(samples_for_responsibility) / responsibility_markers[j];
				while (marker && marker[3] == j) {
					marker[4] = samples_per_marker;
					actual_samples += samples_per_marker;
					marker = (scalar_t*)vector_list_iterate(markers);
				}
			}
			// In case there are more markers than samples
			else {
				// Semi-equal sample the current responsibility
				scalar_t step = SCALAR_NEXT(responsibility_space[j]) / samples_for_responsibility;
				scalar_t current = marker[2];
				while (marker && marker[3] == j) {
					// Search for the next marker to sample
					while (marker && marker[2] <= current) marker = (scalar_t*)vector_list_iterate(markers);
					// In case we got out of markers for the current responsibility
					if (!marker || marker[3] > j) break;
					// Sample the current marker
					marker[4] = 1;
					++actual_samples;
					current += step;
				}
			}
		}

		// Allocate output
		info("Generating dataset of size: " << actual_samples);
		output = new_matrix(actual_samples, 2);
		if (output == MATRIX_ERROR) {
			throw error("cannot allocate memory for output");
		}

		// Populate output
		uint32_t cursor=0;
		marker = (scalar_t*)vector_list_begin(markers);
		for(; marker; marker = (scalar_t*)vector_list_iterate(markers)) {
			// Sample according to demand
			// Random sample
			if (random) {
				for (uint32_t i=0; i<marker[4]; ++i) {
					GET_SCALAR(output, cursor, 0) = gen_uniform_random_scalar(marker[0], marker[1]);
					GET_SCALAR(output, cursor, 1) = marker[2];
					++cursor;
				}
			}
			// Linspace sample
			else {
				// Special case of only one sample, sample the mraker mid point
				if (marker[4] == 1) {
					GET_SCALAR(output, cursor, 0) = (marker[0] + marker[1])/2;
					GET_SCALAR(output, cursor, 1) = marker[2];
					++cursor;
				}
				// All other cases
				else {
					scalar_t current = marker[0];
					scalar_t step = (marker[1] - marker[0]) / marker[4];
					for (uint32_t i=0; i<marker[4]; ++i) {
						GET_SCALAR(output, cursor, 0) = current;
						GET_SCALAR(output, cursor, 1) = marker[2];
						current += step;
						++cursor;
					}
				}
			}
		}

		// Shuffle dataset if necessary
		if (shuffle) {
			// Create random permutation
			uint32_t *permutation = (uint32_t*)malloc(sizeof(uint32_t)*actual_samples);
			random_permutation(permutation, actual_samples);

			// Create new output
			matrix_t* new_output = new_matrix(actual_samples, 2);

			// Populate new output
			for (uint32_t i=0; i<actual_samples; ++i) {
				*(uint32_t*)get_element(new_output, i, 0) = *(uint32_t*)get_element(output, permutation[i], 0);
				*(uint32_t*)get_element(new_output, i, 1) = *(uint32_t*)get_element(output, permutation[i], 1);
			}

			// Free previous dataset
			free_matrix(output);
			free(permutation);
			output=new_output;
		}

	} catch (const std::exception& e) {
		warning(e.what());
		free_matrix(output);
		output = MATRIX_ERROR;
	}

	// Free resources
	vector_list_free(markers);
	vector_list_free(matching_records);
	return output;
}

/**
 * @brief Calculates the maximum error and bucket coverage of a submodel
 * @param model An RQRMI model
 * @param probe An RQRMI probing data structure
 * @param stage_idx The required stage
 * @param bucket_idx The required submodel index
 * @returns A scalar pair of {maximum_error, coverage}. The error may be negative in case the submodels covers no records.
 * @note May throw exceptions
 */
scalar_pair_t rqrmi_tools_calculate_submodel_error(rqrmi_model_t* model, rqrmi_probing_t* probe, uint32_t stage_idx, uint32_t submodel_idx) {


	// Calculate the transition set of the previous stage
	if (probe->transition_sets[stage_idx] == VECTOR_LIST_ERROR) {
		throw error("the transition set of stage " << stage_idx << " was not yet calculated");
	}

	// Get U_i
	vector_list_t* U_i = probe->transition_sets[stage_idx];

	// Get the records in the submodel responsibility
	// Format: [start_input, end_input, record_idx, interval_idx]
	vector_list_t* matching_records = rqrmi_tools_get_records_in_responsibility(probe, stage_idx, submodel_idx);
	if (matching_records == VECTOR_LIST_ERROR) {
		throw error("cannot calculate submodel <" << stage_idx << "," << submodel_idx << "> matching records");
	}

	uint32_t next_width = stage_idx+1 == probe->num_of_stages ? probe->num_of_records : probe->stage_width[stage_idx+1];

	// Statistics
	int max_error = -1;
	scalar_t coverage = 0;
	uint8_t  bucket_is_active[next_width];
	uint8_t  bucket_should_be_active[next_width];

	// Reset the flags
	memset(bucket_is_active, 0, next_width*sizeof(uint8_t));
	memset(bucket_should_be_active, 0, next_width*sizeof(uint8_t));

	try {

		// In case of internal stage, calculate the the bucket for each record
		if (stage_idx < probe->num_of_stages - 1) {
			// Go over all matching records
			scalar_t* cursor = (scalar_t*)vector_list_begin(matching_records);
			for(; cursor; cursor = (scalar_t*)vector_list_iterate(matching_records)) {
				cursor[2] = floor(cursor[2] / probe->num_of_records * next_width);
				bucket_should_be_active[(int)cursor[2]] = 1;
			}
		}

		// Iterate through all records
		scalar_t* cursor = (scalar_t*)vector_list_begin(matching_records);

		// Go over all U_i in responsibility
		vector_list_t* responsibility = probe->responsibilities[stage_idx][submodel_idx];
		scalar_t* interval = (scalar_t*)vector_list_begin(responsibility);
		for(; interval; interval = (scalar_t*)vector_list_iterate(responsibility)) {

			scalar_t interval_srt = interval[0]; // inclusive
			scalar_t interval_end = interval[1]; // inclusive

			scalar_t* pt = (scalar_t*)vector_list_begin(U_i);
			for(; pt; pt = (scalar_t*)vector_list_iterate(U_i)) {

				// Pt format:
				// [ x,  B_i( M_i(x-epsilon) ),  B_i( M_i(x+epsilon) ) ]

				// Skip transition inputs that are not in the responsibility of this
				if (pt[0] < interval_srt) continue;
				else if (pt[0] > interval_end) break;

				// Go over all records smaller than current transition input
				for(; cursor; cursor = (scalar_t*)vector_list_iterate(matching_records)) {
					// The record starts after the transition input
					if (cursor[0] > pt[0]) {
						break;
					}
					// At this point, either the transition input is inside the record,
					// or the entire record is smaller than the transition input.
					// In either case, the bucket is activated and the error should be calculated.
					int required_bucket = cursor[2];

					// This is always true (the record always starts before the transition input)
					int bucket_idx = pt[1];
					bucket_is_active[bucket_idx] = 1;
					max_error = MAX(max_error, abs(required_bucket - bucket_idx));

					// The transition point is in the middle of record
					if (cursor[1] >= pt[0]) {
						int bucket_idx = pt[2];
						bucket_is_active[bucket_idx] = 1;
						max_error = MAX(max_error, abs(required_bucket - bucket_idx));
					}
				}
			}
		}

		// How many buckets should be active?
		scalar_t should_be_active_buckets = 0;
		for (uint32_t i=0; i<next_width; ++i) should_be_active_buckets+=bucket_should_be_active[i];

		// Calculate bucket coverage
		for (uint32_t i=0; i<next_width; ++i) coverage+=bucket_is_active[i];
		coverage /= should_be_active_buckets;

	} catch (const std::exception& e) {
		vector_list_free(matching_records);
		throw e;
	}

	// Lookup procedure requires error margin of 2
	max_error += 2;

	// Return coverage and maximum error
	vector_list_free(matching_records);
	return (scalar_pair_t){ (scalar_t)max_error, coverage };
}

/**
 * @brief Returns the width of a stage
 */
uint32_t rqrmi_tools_get_stage_width(rqrmi_probing_t* probe, uint32_t stage_idx) {
	return probe->stage_width[stage_idx];
}

/**
 * @brief Prints a responsibility to the screen (debug only)
 * @param probe An RQRMI probing data structure
 * @param stage_idx The required stage
 * @param bucket_idx The required bucket in stage
 */
void print_responsibility(rqrmi_probing_t* probe, uint32_t stage_idx, uint32_t bucket_index) {
#ifndef NDEBUG
	SimpleLogger::get() << SimpleLogger::lock();
	SimpleLogger::get().format("R<%u,%u>: ", stage_idx, bucket_index);
	vector_list_t* responsibility = probe->responsibilities[stage_idx][bucket_index];
	// In case the responsibility is empty
	if (vector_list_get_size(responsibility) == 0) {
		SimpleLogger::get().add("none");
	}
	// Print all responsibility intervals
	else {
		scalar_t* interval = (scalar_t*)vector_list_begin(responsibility);
		for(; interval; interval = (scalar_t*)vector_list_iterate(responsibility)) {
			SimpleLogger::get().format("[%f, %f) ", interval[0], SCALAR_NEXT(interval[1]));
		}
	}
	SimpleLogger::get() << SimpleLogger::endl() << SimpleLogger::release();
#endif
}

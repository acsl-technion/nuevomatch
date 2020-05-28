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

/** @file rqrmi_tools.h */

#pragma once

#include <basic_types.h>
#include <matrix_operations.h>
#include <vector_list.h>
#include <rqrmi_model.h>

#define RQRMI_PROBING_ERROR NULL

// Used to store probing datasets for fast operations
typedef struct rqrmi_probing rqrmi_probing_t;

#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief Generates an RQRMI probing data structure.
 * @param records The records the RQRMI index. Nx1 matrix. The key values / range start values. Should be sorted.
 * @param stage_num The number of stages in the RQRMI model
 * @param stages_width An array of integers with the width of each stage
 * @note On error returns RQRMI_PROBING_ERROR
 */
rqrmi_probing_t* rqrmi_tools_probe_new(matrix_t* records, uint32_t stage_num, uint32_t* stages_width);

/**
 * @brief Frees the allocated memory of rqrmi_probing_t structure
 */
void rqrmi_tools_probe_free(rqrmi_probing_t* probe);

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
vector_list_t* rqrmi_tools_calculate_transition_set(rqrmi_model_t* model, rqrmi_probing_t* probe, uint32_t stage_idx);

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
vector_list_t** rqrmi_tools_calculate_responsibility(rqrmi_model_t* model, rqrmi_probing_t* probe, uint32_t stage_idx);

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
matrix_t* rqrmi_tools_generate_dataset(rqrmi_probing_t* probe, uint32_t stage_idx, uint32_t bucket_idx,
		uint32_t num_of_samples, bool random, bool shuffle);

/**
 * @brief Calculates the maximum error and bucket coverage of a submodel
 * @param model An RQRMI model
 * @param probe An RQRMI probing data structure
 * @param stage_idx The required stage
 * @param bucket_idx The required submodel index
 * @returns A scalar pair of {maximum_error, coverage}. The error may be negative in case the submodels covers no records.
 * @note May throw exceptions
 */
scalar_pair_t rqrmi_tools_calculate_submodel_error(rqrmi_model_t* model, rqrmi_probing_t* probe, uint32_t stage_idx, uint32_t submodel_idx);

/**
 * @brief Returns the width of a stage
 */
uint32_t rqrmi_tools_get_stage_width(rqrmi_probing_t* probe, uint32_t stage_idx);

#ifdef __cplusplus
}
#endif

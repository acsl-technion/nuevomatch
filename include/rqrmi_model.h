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

/** @file rqrmi_model.h */

#pragma once

#include <matrix_operations.h>
#include <vector_list.h>

// Error return value
#define RQRMI_MODEL_ERROR NULL

// Used to identify read error
#define SUBMODEL_MAX_LAYERS 100
#define RQRMI_LAYER_MAX_WIDTH 1000

// Defines the number of input nodes of the RQRMI model
#define RQRMI_MODEL_INPUT_SIZE 1

typedef struct rqrmi_submodel rqrmi_submodel_t;
typedef struct rqrmi_stage rqrmi_stage_t;
typedef struct rqrmi_model rqrmi_model_t;

/**
 * @brief Holds all necessary information of an RQRMI submodel
 */
typedef struct {
	scalar_t w0;
	scalar_t b0;
	scalar_t w1[8];
	scalar_t b1[8];
	scalar_t w2[8];
	scalar_t b2;
	scalar_t output_factor;
	scalar_t output_min;
	scalar_t input_mean;
	scalar_t input_stddev;
	uint32_t error;
	uint32_t compiled;
} rqrmi_submodel_info_t;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Loads a RQRMI model
 * @param ptr A cursor for reading the data
 * @param size The size of the memory block pointed by cursor.
 * @returns A pointer to the RQRMI model, or MODEL_ERROR in case of an error
 * @note All numbers are stored as 32-bit little-endian encoding (integers & floats)
 */
rqrmi_model_t* rqrmi_load_model(void* ptr, int size);

/**
 * @brief Frees the memory allocated for the RQRMI model
 * @param rqrmi_model The RQRMI model
 */
void rqrmi_free_model(rqrmi_model_t* rqrmi_model);

/**
 * @brief Places an input to an RQRMI model (and all of its submodels)
 * @param rqrmi_model An RQRMI model
 * @param input The input to place
 * @note  This function is unsafe.
 */
void rqrmi_place_input(rqrmi_model_t* rqrmi_model, scalar_t input);

/**
 * @brief returns true iff the submodel is compiled
 * @param rqrmi_model The RQRMI model
 * @param stage_idx The stage index of the RQRMI model
 * @param submodel_idx The model index of the corresponding stage
 */
bool rqrmi_submodel_compiled(rqrmi_model_t* rqrmi_model, uint32_t stage_idx, uint32_t submodel_idx);

/**
 * @brief Returns the submodel corresponding to the specified stage and index
 * @param rqrmi_model The RQRMI model
 * @param stage The stage index of the RQRMI model
 * @param index The model index of the corresponding stage
 * @note Might throw exceptions
 */
rqrmi_submodel_t* rqrmi_get_submodel(rqrmi_model_t* rqrmi_model, uint32_t stage, uint32_t index);

/**
 * @brief Feed the specified model with the input
 * @param rqrmi_model The RQRMI model
 * @param model The model to feed
 * @returns The output of the model
 * @note Might throw exceptions
 */
matrix_t* rqrmi_evaluate_sub_model(rqrmi_model_t* rqrmi_model, rqrmi_submodel_t* model);

/**
 * @brief Feed the RQRMI model with an input and return the output
 * @param rqrmi_model The RQRMI model
 * @param input The input to feed to the model.
 * @note Might throw exceptions
 */
scalar_t rqrmi_evaluate_model(rqrmi_model_t* rqrmi_model, scalar_t input);

/**
 * @brief Calculate the trigger inputs of an RQRMI submodel
 * @param rqrmi_model the RQRMI model
 * @param stage_idx the required stage index
 * @param submodel_idx the required submodel index
 * @returns A matrix (Nx1) with sorted trigger inputs (scalar_t). On error, returns MATRIX_ERROR.
 * @note The matrix should be freed by the user
 */
matrix_t* rqrmi_calculate_trigger_inputs(rqrmi_model_t* rqrmi_model, uint32_t stage_idx, uint32_t submodel_idx);

/**
 * @brief Calculate the transition inputs of an RQRMI submodel
 * @param rqrmi_model the RQRMI model
 * @param stage_idx the required stage index
 * @param submodel_idx the required submodel index
 * @returns A vector list with sorted transition inputs. On error, returns VECTOR_LIST_ERROR.
 * 		    Format of rows: [ x, B(M(x-eps)), B(M(x+eps)) ]
 * @note The vector list should be freed by the user
 */
vector_list_t* rqrmi_calculate_transition_inputs(rqrmi_model_t* rqrmi_model,
		uint32_t stage_idx, uint32_t submodel_idx, uint32_t next_width);

// Getters & Setters

/**
 * @brief Returns number of stages in an RQRMI model
 */
uint32_t rqrmi_get_num_of_stages(rqrmi_model_t* rqrmi_model);

/**
 * @brief Returns number of submodel in an RQRMI stage
 */
uint32_t rqrmi_get_num_of_submodels(rqrmi_model_t* rqrmi_model, uint32_t stage);

/**
 * @brief Returns the RQRMI last submodel index after evaluation
 */
uint32_t rqrmi_get_last_sub_model_index(rqrmi_model_t* rqrmi_model);

/**
 * @brief Returns the RQRMI last error
 */
uint32_t rqrmi_get_last_error(rqrmi_model_t* rqrmi_model);

/**
 * @brief Returns the RQRMI error of a specific submodel
 * @param[out] list The error list
 * @param[out] size Number of elements in list
 */
void rqrmi_get_error_list(rqrmi_model_t* rqrmi_model, const uint32_t** list, uint32_t* size);

/**
 * @brief Get the requested submodel compilation state
 * @param rqrmi_model An RQRMI model
 * @param stage_idx The required stage's index
 * @param submodel_idx The required submodel's index
 * @returns 1 if the submodel is compiled, or 0 if the submodel is not available / not compiled
 */
uint32_t rqrmi_get_submodel_complie_state(rqrmi_model_t* rqrmi_model,
		uint32_t stage_idx, uint32_t submodel_idx);

/**
 * @brief Returns the input domain of an RQRMI model
 * @param rqrmi_model An RQRMI model
 * @return {min, max} scalar pair
 */
scalar_pair_t rqrmi_get_input_domain(rqrmi_model_t* rqrmi_model);

/**
 * @brief Sets number of stages in an RQRMI model
 */
void rqrmi_set_num_of_stages(rqrmi_model_t* rqrmi_model, uint32_t value);

/**
 * @brief Returns all the necessary information regarding an RQRMI submodel
 * @param rqrmi_model An RQRMI model
 * @param stage_idx The required stage's index
 * @param submodel_idx The required submodel's index
 * @param[out] output The output information
 * @returns 1 on success, 0 on error
 * @note RMI will not work here, only used for RQRMI models
 */
uint32_t rqrmi_get_submodel_info(rqrmi_model_t* rqrmi_model, uint32_t stage_idx,
		uint32_t submodel_idx, rqrmi_submodel_info_t* output);

#ifdef __cplusplus
}
#endif

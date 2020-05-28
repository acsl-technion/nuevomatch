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

#include <rqrmi_model.h>

#include <logging.h>
#include <matrix_operations.h>
#include <vector_list.h>

#define MAX_UINT 0xffffffff

#define MODEL_OP_ERROR 0
#define MODEL_OP_SUCCESS 1
// Scratchpads for fast evaluation
#define NUM_OF_SCRATCPADS 3

/**
 * @brief Macro for fast loading from memory and perform size checks
 * @param T Typename to read from buffer
 * @param buff The buffer to read from
 * @param size Size argument for buffer overflow check.
 */
#define safe_buffer_read(T, buf, size) \
	*(T*)buf; \
	if ((size-=sizeof(T))<0) { \
		throw error("cannot read " << size << " bytes from buffer"); \
	} \
	buf=(T*)buf+1;

// Model information is often not required for debugging.
// Set dedicated debugging flag for models
#ifdef DEBUG_MODEL
#  define model_info(...) info(__VA_ARGS__)
#else
#  define model_info(...)
#endif

struct rqrmi_submodel {
	uint32_t num_of_layers;
	matrix_t** biases;
	matrix_t** weights;
	unary_operation_t* activations;
	scalar_t input_stddev;
	scalar_t input_mean;
	scalar_t output_factor;
	scalar_t output_min;
};

struct rqrmi_stage {
	uint32_t num_of_models;
	rqrmi_submodel_t* models;
};

struct rqrmi_model {
	uint32_t num_of_stages;
	scalar_t input_domain_min;
	scalar_t input_domain_max;
	rqrmi_stage_t* stages;
	matrix_t* input_placeholder;

	uint32_t last_model;
	uint32_t* error_list;

	// Scratchpads for matrix evaluation
	void* scratchpads;
	uint32_t scratchpad_size;
};

// Local Methods
int load_submodel(rqrmi_model_t* rqrmi_model, uint32_t stage_index, uint32_t model_index, void** ptr, int* size);
void free_model(rqrmi_model_t* rqrmi_model, uint32_t stage_index, uint32_t model_index);
scalar_t rqrmi_submodel_preprocess_input(rqrmi_submodel_t* submodel, scalar_t input);
scalar_t rqrmi_submodel_reverse_process(rqrmi_submodel_t* submodel, scalar_t processed_input);


/**
 * @brief returns true iff the submodel is compiled
 * @param rqrmi_model The RQRMI model
 * @param stage_idx The stage index of the RQRMI model
 * @param submodel_idx The model index of the corresponding stage
 */
bool rqrmi_submodel_compiled(rqrmi_model_t* rqrmi_model, uint32_t stage_idx, uint32_t submodel_idx) {
	return rqrmi_model->stages[stage_idx].models[submodel_idx].biases != NULL;
}

/**
 * @brief Returns the NN model corresponding to the specified stage and index
 * @param rqrmi_model The RQRMI model
 * @param stage The stage index of the RQRMI model
 * @param index The model index of the corresponding stage
 * @note Might throw exceptions
 */
rqrmi_submodel_t* rqrmi_get_submodel(rqrmi_model_t* rqrmi_model, uint32_t stage, uint32_t index) {

// Perform checks on debug
#ifndef NDEBUG
	if (stage > rqrmi_model->num_of_stages) {
		throw error("requested to run submodel <" << stage << "," << index << "> with invalid index");
	}
#endif

	rqrmi_stage_t* stage_p = &rqrmi_model->stages[stage];

// Perform checks on debug
#ifndef NDEBUG
	if (index > stage_p->num_of_models) {
		throw error("requested to run submodel <" << stage << "," << index << "> with invalid index");
	}
#endif

	// Check whether the model is compiled
	if (!rqrmi_submodel_compiled(rqrmi_model, stage, index)) {
		throw error("requested to run submodel <" << stage << "," << index << "> which is not compiled. "
				"input: " << GET_SCALAR(rqrmi_model->input_placeholder,0,0) <<
				", last stage submodel: " << rqrmi_model->last_model);
	}
	return &stage_p->models[index];
}

/**
 * @brief Places an input to an RQRMI model (and all of its submodels)
 * @param rqrmi_model An RQRMI model
 * @param input The input to place
 * @note  This function is inline thus unsafe.
 */
void rqrmi_place_input(rqrmi_model_t* rqrmi_model, scalar_t input) {
	GET_SCALAR(rqrmi_model->input_placeholder, 0, 0) = input;
}

/**
 * @brief Feed the specified model with the input
 * @param rqrmi_model The RQRMI model
 * @param model The model to feed
 * @returns The output of the model
 * @note Might throw exceptions
 */
matrix_t* rqrmi_evaluate_sub_model(rqrmi_model_t* rqrmi_model, rqrmi_submodel_t* submodel) {

// Perform checks on debug
#ifndef NDEBUG
	if (submodel == RQRMI_MODEL_ERROR) {
		throw error("Cannot evaluate submodel: one of the inputs is invalid");
	}
#endif

	// References to scratchpads
	matrix_t* sp_0 = (matrix_t*)((char*)rqrmi_model->scratchpads);
	matrix_t* sp_1 = (matrix_t*)((char*)rqrmi_model->scratchpads + rqrmi_model->scratchpad_size);
	matrix_t* sp_2 = (matrix_t*)((char*)rqrmi_model->scratchpads + 2*rqrmi_model->scratchpad_size);

	matrix_t* current = rqrmi_model->input_placeholder;

	// Store old input value, as computation occur in-place
	scalar_t old_input_value = GET_SCALAR(current, 0, 0);

	// Preprocess the input to fit with the statistics
	GET_SCALAR(current, 0, 0) = rqrmi_submodel_preprocess_input(submodel, old_input_value);

	model_info("Input for submodel after preprocessing: " << GET_SCALAR(current, 0, 0));

	for (uint32_t i=0; i<submodel->num_of_layers; ++i) {
		// Multiply with weights
		mat_mul(current, submodel->weights[i], sp_0);

		// Add bias
		mat_op(sp_0, submodel->biases[i], op_add, sp_1);

		// In case of no activation
		if (submodel->activations[i] == OPERATION_BYPASS) {
			current = sp_1;
			continue;
		}

		// Enable activation
		mat_unary_op(sp_1, submodel->activations[i], sp_2);
		current = sp_2;
	}

	// Return the original input value
	GET_SCALAR(rqrmi_model->input_placeholder, 0, 0) = old_input_value;

	model_info("Submodel output before post-processing: " << GET_SCALAR(current, 0, 0));

	// Post-process output
	*(scalar_t*)get_element(current, 0, 0) =
			*(scalar_t*)get_element(current, 0, 0) * submodel->output_factor + submodel->output_min;
	return current;
}

/**
 * @brief Feed the RQRMI model with an input and return the output
 * @param rqrmi_model The RQRMI model
 * @param input The input to feed to the model.
 * @note Might throw exceptions
 */
scalar_t rqrmi_evaluate_model(rqrmi_model_t* rqrmi_model, scalar_t input) {
	matrix_t* current_output;
	scalar_t out_scalar = 0;

	// Set the input
	rqrmi_place_input(rqrmi_model, input);

	for(uint32_t i=0; i<rqrmi_model->num_of_stages; ++i) {

		// Set the model index of the current stage
		uint32_t next_submodel = (out_scalar < 0 ? 0 : out_scalar >= 1 ? 1-SCALAR_EPS : out_scalar) * rqrmi_model->stages[i].num_of_models;

		model_info("Evaluating stage " << i << "with model " << next_submodel);

		// Get the current model
		rqrmi_submodel_t* submodel = rqrmi_get_submodel(rqrmi_model, i, next_submodel);
		rqrmi_model->last_model = next_submodel;

		// Evaluate the model
		current_output = rqrmi_evaluate_sub_model(rqrmi_model, submodel);
		out_scalar = GET_SCALAR(current_output,0,0);

		model_info("Stage " << i << " output is " << GET_SCALAR(current_output, 0, 0) <<
				   ", corresponds to model index " << next_submodel);
	}

	out_scalar = out_scalar<0 ? 0 : out_scalar>=1 ? 1 - SCALAR_EPS : out_scalar;
	// No need to free current_output, it points to a scratchpad

	return out_scalar;
}

/**
 * @brief Loads a RQRMI model
 * @param ptr A cursor for reading the data
 * @param size The size of the memory block pointed by cursor.
 * @returns A pointer to the RQRMI model, or MODEL_ERROR in case of an error
 * @note  All numbers are stored as 32-bit little-endian encoding
 */
rqrmi_model_t* rqrmi_load_model(void* ptr, int size) {

#	ifndef NDEBUG
		void* orig_ptr = ptr;
#	endif

	// Allocate RQRMI memory
	rqrmi_model_t* rqrmi_model = (rqrmi_model_t*)malloc(sizeof(rqrmi_model_t));
	if (rqrmi_model == NULL) {
		warning("Cannot allocate memory for RQRMI model");
		return RQRMI_MODEL_ERROR;
	}

	info("Loading RQRMI model from address " << ptr << " with size " << size << " to " << rqrmi_model);

	try {

		// Initiate pointers and sizes
		rqrmi_model->input_placeholder = MATRIX_ERROR;
		rqrmi_model->scratchpad_size = 0;
		rqrmi_model->scratchpads = NULL;
		rqrmi_model->error_list = NULL;
		rqrmi_model->stages = MATRIX_ERROR;

		// Read the input domain
		rqrmi_model->input_domain_min = safe_buffer_read(scalar_t, ptr, size);
		rqrmi_model->input_domain_max = safe_buffer_read(scalar_t, ptr, size);

		uint32_t num_of_stages = safe_buffer_read(uint32_t, ptr, size);
		rqrmi_model->num_of_stages = num_of_stages;

		info("Input domain: [" << rqrmi_model->input_domain_min << "," <<
			 rqrmi_model->input_domain_max << "], num of stages: " << num_of_stages);

		// Allocate stages memory
		rqrmi_model->stages = (rqrmi_stage_t*)malloc(sizeof(rqrmi_stage_t) * num_of_stages);
		if (rqrmi_model->stages == NULL) {
			throw error("cannot allocate memory for RQRMI model stages");
		}

		// Set all stages to hold no data
		for (uint32_t i=0; i<num_of_stages; ++i) {
			rqrmi_model->stages[i].models = MATRIX_ERROR;
			rqrmi_model->stages[i].num_of_models=0;
		}

		// Allocate input placeholder
		rqrmi_model->input_placeholder = new_matrix(1,1);
		if (rqrmi_model->input_placeholder == MATRIX_ERROR) {
			throw error("cannot allocate memory for RQRMI model input placeholder");
		}

		// For each stage
		for (uint32_t i=0; i<num_of_stages; ++i) {

			info("Reading stage " << i << " out of " << num_of_stages-1);

			// Read stage's properties
			uint32_t num_of_models = safe_buffer_read(uint32_t, ptr, size);

			info("Number of models in stage: " << num_of_models);

			// Build the model array
			rqrmi_model->stages[i].num_of_models = num_of_models;
			rqrmi_model->stages[i].models = (rqrmi_submodel_t*)malloc(sizeof(rqrmi_submodel_t) * num_of_models);

			// Try to load each submodel
			for (uint32_t j=0; j<num_of_models; ++j) {
				if (load_submodel(rqrmi_model, i, j, &ptr, &size) == MODEL_OP_ERROR) {
					throw error("cannot load RQRMI model - error occurred on submodel <" << i << "," << j <<">");
				}
			}
		}

		// Allocate memory for error list
		uint32_t w = rqrmi_model->stages[num_of_stages-1].num_of_models;
		rqrmi_model->error_list = (uint32_t*)malloc(sizeof(uint32_t)*w);
		if (rqrmi_model->error_list == NULL) {
			throw error("cannot allocate memory for error list");
		}

		// Read the error list
		for (uint32_t i=0; i<w; ++i) {
			rqrmi_model->error_list[i] = safe_buffer_read(uint32_t, ptr, size);
		}

		// Print the error list
#ifndef NDEBUG
		SimpleLogger::get().add("Submodel error list: [");
		for (uint32_t i=0; i<w; ++i) {
			SimpleLogger::get() << SimpleLogger::format("%u", rqrmi_model->error_list[i]);
			if (i < w-1) SimpleLogger::get().add(", ");
		}
		SimpleLogger::get().add("]\n");
#endif

		// Allocate scratchpads for fast evaluation
		info("Allocating " << NUM_OF_SCRATCPADS << " scratchpads, each with size "
				<< std::hex << rqrmi_model->scratchpad_size << std::dec << " bytes.");
		rqrmi_model->scratchpads=malloc(NUM_OF_SCRATCPADS*rqrmi_model->scratchpad_size);
		if (rqrmi_model->scratchpads == NULL) {
			throw error("cannot allocate memory for scratchpads");
		}
		info("Scratchpads start at " << rqrmi_model->scratchpads << ", end at " <<
				(char*)rqrmi_model->scratchpads + NUM_OF_SCRATCPADS*rqrmi_model->scratchpad_size);

		// Set scratchpads elements reference
		for (int i=0; i<NUM_OF_SCRATCPADS;++i) {
			matrix_t* sp = (matrix_t*)( (char*)rqrmi_model->scratchpads + i*rqrmi_model->scratchpad_size);
			sp->elements = (void*)( (char*)sp + sizeof(matrix_t) );
		}

		// Log total bytes read
		info("Finished loading RQRMI model from memory. Total bytes read: " << ((uint64_t)ptr - (uint64_t)orig_ptr));

	} catch (const std::exception& e) {
		warning(e.what());
		rqrmi_free_model(rqrmi_model);
		rqrmi_model = RQRMI_MODEL_ERROR;
	}

	return rqrmi_model;
}


/**
 * @brief Loads a submodel
 * @param rqrmi_model The RQRMI model
 * @param stage_index The submodel's stage index
 * @param model_index The submodel's index within the stage
 * @param ptr A pointer to the cursor used for reading the data
 * @param size A pointer to the remaining buffer size.
 * 			   In case the parser excess that size, an error occurs
 * @returns MODEL_OP_SUCCESS on success, MODEL_OP_ERROR otherwise
 * @note  All numbers are stored as 32-bit little-endian encoding
 */
int load_submodel(rqrmi_model_t* rqrmi_model, uint32_t stage_index, uint32_t model_index, void** ptr, int* size) {

	info("Loading submodel <" << stage_index << "," << model_index << "> from address " << std::hex << *ptr << std::dec << "");

	// Check for size minimum requirement
	*size -= sizeof(uint8_t);
	if (*size < 0) {
		warning("Cannot read submodel <" << stage_index << "," << model_index << ">: buffer size smaller than required");
		return MODEL_OP_ERROR;
	}

	rqrmi_submodel_t* submodel = &rqrmi_model->stages[stage_index].models[model_index];
	if (submodel == MATRIX_ERROR) {
		warning("Cannot read submodel <" << stage_index << "," << model_index << ">: cannot acquire model");
		return MODEL_OP_ERROR;
	}

	// Read submodel's version
	uint8_t* char_ptr = *(uint8_t**)ptr;
	uint8_t model_version = *(char_ptr++);
	switch (model_version) {

		// The submodel is not compiled
		case 0:
			submodel->biases = MATRIX_ERROR;
			submodel->weights = MATRIX_ERROR;
			submodel->num_of_layers = 0;
			info("Submodel <" << stage_index << "," << model_index << "> is not compiled, skipping read");

			// Update the pointer location
			*ptr = char_ptr;
			return MODEL_OP_SUCCESS;

		// Default behavior
		case 1:
			*size -= 2*sizeof(scalar_t);
			if (*size < 0) {
				warning("Cannot read model <" << stage_index << "," << model_index <<
						">: buffer size smaller than required");
				return MODEL_OP_ERROR;
			}

			scalar_t* scalr_ptr = (scalar_t*)char_ptr;
			submodel->input_mean = *(scalr_ptr++);
			submodel->input_stddev = *(scalr_ptr++);
			submodel->output_factor = *(scalr_ptr++);
			submodel->output_min = *(scalr_ptr++);
			char_ptr = (uint8_t*)scalr_ptr;
			break;
	}

	info("Submodel <" << stage_index << "," << model_index << "> input mean: "
			<< submodel->input_mean << ", input stddev: " << submodel->input_stddev);

	// Check for size minimum requirement
	*size -= 2*sizeof(uint32_t);
	if (*size < 0) {
		error("Cannot read model <" << stage_index << "," << model_index <<
				">: buffer size smaller than required");
		return MODEL_OP_ERROR;
	}

	// Read the submodel's header
	uint32_t* header = (uint32_t*)char_ptr;

	// Build the layers
	uint32_t num_of_layers = *(header++);

	submodel->num_of_layers = num_of_layers;
	info("Loading model's data. Number of total layers: " << num_of_layers);

	// Validate no read errors
	if (num_of_layers > SUBMODEL_MAX_LAYERS) {
		warning("Read error for submodel <" << stage_index << "," << model_index <<
				">: number of stages is higher than maximum");
		return MODEL_OP_ERROR;
	}

	// Allocate memory for the matrices array
	submodel->biases = (matrix_t**)malloc(sizeof(matrix_t*) * (num_of_layers) );
	submodel->weights = (matrix_t**)malloc(sizeof(matrix_t*) * (num_of_layers) );
	submodel->activations = (unary_operation_t*)malloc(sizeof(unary_operation_t) * (num_of_layers) );
	if (submodel->biases == NULL || submodel->weights == NULL) {
		warning("Cannot allocate memory for model <" << stage_index << "," << model_index << ">");
		free(submodel->biases);
		free(submodel->weights);
		return MODEL_OP_ERROR;
	}

	// Generate the input layer
	uint32_t last_layer_width = RQRMI_MODEL_INPUT_SIZE;

	// Allocate memory for the matrices
	for (uint32_t k=0; k<num_of_layers; ++k) {
		uint32_t layer_width = *(header++);
		submodel->biases[k] = new_matrix(1, layer_width);
		submodel->weights[k] = new_matrix(last_layer_width, layer_width);
		last_layer_width = layer_width;

		// Set the required scratchpad size
		uint32_t requred_size = sizeof(matrix_t)+sizeof(scalar_t)*(last_layer_width*layer_width);
		rqrmi_model->scratchpad_size = requred_size > rqrmi_model->scratchpad_size ? requred_size : rqrmi_model->scratchpad_size;

		// Hidden layers have ReLU activation
		if (k>0 && k<num_of_layers-1) {
			submodel->activations[k] = op_relu;
		} else {
			submodel->activations[k] = OPERATION_BYPASS;
		}

		// Check for errors
		if (submodel->biases[k] == MATRIX_ERROR || submodel->weights[k] == MATRIX_ERROR) {
			warning("Cannot allocate memory for model <" << stage_index << "," << model_index << ">");
			return MODEL_OP_ERROR;
		}
	}

	// Read the model's data
	scalar_t* data = (scalar_t*)header;

	for (uint32_t k=0; k<num_of_layers; ++k) {
		matrix_t* mat_bias = submodel->biases[k];
		matrix_t* mat_weights = submodel->weights[k];
		info("Layer " << k << " - bias size: [" << mat_bias->rows << " x " << mat_bias->cols
				<< "]. weights size: [" << mat_weights->rows << " x " << mat_weights->cols << "]");

		// Check for size requirement
		*size -= (mat_bias->cols + mat_weights->rows*mat_weights->cols)*sizeof(scalar_t);
		if (*size < 0) {
			warning("Cannot read model <" << stage_index << "," << model_index <<
					">: buffer size smaller than required");
			return MODEL_OP_ERROR;
		}

		// Read biases
		for (uint32_t x=0; x<mat_bias->cols; ++x) {
			*(scalar_t*)get_element(mat_bias, 0, x) = *(data++);
		}
		// Read weights
		for (uint32_t y=0; y<mat_weights->rows; ++y) {
			for (uint32_t x=0; x<mat_weights->cols; ++x) {
				*(scalar_t*)get_element(mat_weights, y, x) = *(data++);
			}
		}
	}

	// Update the pointer location
	*ptr = data;
	return MODEL_OP_SUCCESS;
}

/**
 * @brief Calculate the trigger inputs of an RQRMI submodel
 * @param rqrmi_model the RQRMI model
 * @param stage_idx the required stage index
 * @param submodel_idx the required submodel index
 * @returns A matrix (Nx1) with sorted trigger inputs (scalar_t). On error, returns MATRIX_ERROR.
 * @note The matrix should be freed by the user
 */
matrix_t* rqrmi_calculate_trigger_inputs(rqrmi_model_t* rqrmi_model, uint32_t stage_idx, uint32_t submodel_idx) {

	scalar_t b0, w0;
	matrix_t* output = MATRIX_ERROR;
	vector_list_t* trigger_inputs = VECTOR_LIST_ERROR;
	rqrmi_submodel_t* submodel;
	scalar_t* cursor;

	try {

	    // Check inputs
		if ( rqrmi_model == RQRMI_MODEL_ERROR ||
			 stage_idx >= rqrmi_model->num_of_stages ||
			 submodel_idx >= rqrmi_model->stages[stage_idx].num_of_models)
		{
			throw error("invalid inputs for submodel <" << stage_idx << "," << submodel_idx << ">");
		}

	    // Get the relevant submodel
	    submodel = rqrmi_get_submodel(rqrmi_model, stage_idx, submodel_idx);
	    if (submodel == RQRMI_MODEL_ERROR) {
	    	throw error("cannot find submodel <" << stage_idx << "," << submodel_idx << ">");
	    }

	    // Allocate vecotr list
	    trigger_inputs = vector_list_create(1);
	    if (trigger_inputs == VECTOR_LIST_ERROR) {
	    	throw error("cannot allocate vector list");
	    }

	    // Analytically calculate trigger points
	    b0 = GET_SCALAR(submodel->biases[0], 0, 0);
	    w0 = GET_SCALAR(submodel->weights[0], 0, 0);
	    for (uint32_t i=0; i<submodel->biases[1]->cols; ++i) {
	    	scalar_t b1 = GET_SCALAR(submodel->biases[1], 0, i);
	    	scalar_t w1 = GET_SCALAR(submodel->weights[1], 0, i);
	    	// Calculate the NN input value that cause trigger
	    	scalar_t trg_value = -(b1/w0/w1) -(b0/w0);
	    	// Reverse submodel preprocessing
	    	trg_value = trg_value * submodel->input_stddev + submodel->input_mean;
	    	// Skip values outside input domain
	    	if (trg_value < rqrmi_model->input_domain_min || trg_value > rqrmi_model->input_domain_max) continue;
	    	cursor = (scalar_t*)vector_list_push_back_and_get(trigger_inputs);
	    	cursor[0] = trg_value;
	    }

	    // Add input domain boundaries
	    cursor = (scalar_t*)vector_list_push_back_and_get(trigger_inputs);
	    cursor[0] = rqrmi_model->input_domain_min;
	    cursor = (scalar_t*)vector_list_push_back_and_get(trigger_inputs);
		cursor[0] = rqrmi_model->input_domain_max;

	    // Allocate output matrix
	    output = vector_list_to_matrix(trigger_inputs);

	    // Sort the trigger inputs by value
	    qsort(output->elements, output->rows, sizeof(scalar_t), scalar_compare_asc);

	    // Debug Print
#ifndef NDEBUG
	    SimpleLogger::get().format("G<%u,%u>: {", stage_idx, submodel_idx);
			for (uint32_t i=0; i<output->rows; ++i) {

				// Calculate submodel output on point
				rqrmi_place_input(rqrmi_model, GET_SCALAR(output, i, 0));
				matrix_t* result = rqrmi_evaluate_sub_model(rqrmi_model, submodel);

				SimpleLogger::get() << SimpleLogger::format("%f (M=%f)", GET_SCALAR(output, i, 0), GET_SCALAR(result, 0, 0));
				if (i < output->rows-1) SimpleLogger::get().add(", ");
			}
			SimpleLogger::get().add("}\n");
#endif

	} catch (const std::exception& e) {
		warning(e.what());
		free_matrix(output);
		output=MATRIX_ERROR;
	}

	vector_list_free(trigger_inputs);
	return output;
}


/**
 * @brief Calculate the transition inputs of an RQRMI submodel
 * @param rqrmi_model the RQRMI model
 * @param stage_idx the required stage index
 * @param submodel_idx the required submodel index
 * @returns A vector list with sorted transition inputs. On error, returns VECTOR_LIST_ERROR.
 * 		    Format of rows: [ x, B(M(x-eps)), B(M(x+eps)) ]
 * @note The vector list should be freed by the user
 */
vector_list_t* rqrmi_calculate_transition_inputs(rqrmi_model_t* rqrmi_model, uint32_t stage_idx, uint32_t submodel_idx, uint32_t next_width) {

	rqrmi_submodel_t* submodel;
	matrix_t* submodel_output;
	matrix_t* trigger_inputs = MATRIX_ERROR;
	vector_list_t* transition_inputs = VECTOR_LIST_ERROR;

	// Check inputs
	if ( rqrmi_model == RQRMI_MODEL_ERROR ||
		 stage_idx >= rqrmi_model->num_of_stages ||
		 submodel_idx >= rqrmi_model->stages[stage_idx].num_of_models)
	{
		warning("invalid inputs for submodel <" << stage_idx << "," << submodel_idx << ">");
		return VECTOR_LIST_ERROR;
	}

	try {

		// Get the trigger inputs of the submodel m_{i,j}
		trigger_inputs = rqrmi_calculate_trigger_inputs(rqrmi_model, stage_idx, submodel_idx);
		if (trigger_inputs == MATRIX_ERROR) {
			throw error("Cannot calculate trigger inputs for submodel <" << stage_idx << "," << submodel_idx << ">");
		}

		// Get the relevant submodel
		submodel = rqrmi_get_submodel(rqrmi_model, stage_idx, submodel_idx);
		if (submodel == RQRMI_MODEL_ERROR) {
			throw error("cannot find submodel <" << stage_idx << "," << submodel_idx << ">");
		}

		// Hold output transition inputs
		transition_inputs = vector_list_create(3);
		if (transition_inputs == VECTOR_LIST_ERROR) {
			throw error("Cannot allocate memory for transition candidates list");
		}

		// Calculate the transition points of the submodel (lemma 3.3.2)
		for (uint32_t i=0; i<trigger_inputs->rows-1; ++i) {

			// Get two adjacent trigger inputs
			// (As the world sees them)
			scalar_t x_0 = GET_SCALAR(trigger_inputs, i, 0);
			scalar_t x_1 = GET_SCALAR(trigger_inputs, i+1, 0);

			// Evaluate submodel m_{i,j} on both x_0 and x_1
			rqrmi_place_input(rqrmi_model, x_0);
			submodel_output = rqrmi_evaluate_sub_model(rqrmi_model, submodel);
			scalar64_t M_0 = GET_SCALAR(submodel_output, 0, 0);

			rqrmi_place_input(rqrmi_model, x_1);
			submodel_output = rqrmi_evaluate_sub_model(rqrmi_model, submodel);
			scalar64_t M_1 = GET_SCALAR(submodel_output, 0, 0);

			// Both M_0 and M_1 are as the submodel sees its outputs

			// Calculate bucket functions of M_0 and M_1
			// Buckets are calculated according to the world (and not the submodel)
			uint32_t B_0 = (M_0 < 0 ? 0 : M_0 >= 1 ? 1 - SCALAR_EPS: M_0) * next_width;
			uint32_t B_1 = (M_1 < 0 ? 0 : M_1 >= 1 ? 1 - SCALAR_EPS: M_1) * next_width;

			// If no transition inputs in [x_0, x_1], continue
			if (B_0 == B_1) continue;

			// Get min & max bucket
			uint32_t B_min = MIN(B_0, B_1);
			uint32_t B_max = MAX(B_0, B_1);

			// Preprocess x_0 and x_1 (how the submodel sees them)
			scalar64_t xp_0 = rqrmi_submodel_preprocess_input(submodel, x_0);
			scalar64_t xp_1 = rqrmi_submodel_preprocess_input(submodel, x_1);

			// Add internal candidates between B_0 and B_1 (lemma 3.3.2 section 2)
			for (scalar64_t y=B_min+1; y<=B_max; ++y) {

				// Each candidate is the transition input between (y-1) to (y)

				// How the submodel sees the required output (for linear equation)
				scalar64_t raw_y = y / next_width;

				// Calculate candidate for M(x) from linear equation
				scalar64_t candidate = (raw_y-M_0) * (xp_1-xp_0) / (M_1-M_0) + xp_0;

				// The candidate must be between x_0 and x_1..
#				ifndef NDEBUG
					if (candidate < xp_0 || candidate > xp_1) {
						warning("candidate " << candidate << "for T<" << stage_idx << "," << submodel_idx <<
								"> is not between two trigger inputs: xp_0=" << xp_0 << ", xp_1=" << xp_1);
					}
#				endif

				// Convert candidate to the output world format.
				candidate = rqrmi_submodel_reverse_process(submodel, candidate);

				// The transition is set according to the linear behavior between x_0 and x_1
				uint32_t B_minus = M_0 < M_1 ? (y-1) : y;
				uint32_t B_plus = M_0 > M_1 ? (y-1) : y;

				// Add new transition input
				scalar_t* vector = (scalar_t*)vector_list_push_back_and_get(transition_inputs);
				vector[0] = candidate;
				vector[1] = B_minus;
				vector[2] = B_plus;
			}

			// Sort the transition inputs
			vector_list_sort(transition_inputs, 0, scalar_compare_asc);

		}
	} catch (const std::exception& e) {
		warning(e.what());
		vector_list_free(transition_inputs);
		transition_inputs = VECTOR_LIST_ERROR;
	}

	return transition_inputs;
}


/**
 * @brief Frees the memory allocated for the RQRMI model
 * @param rqrmi_model The RQRMI model
 */
void rqrmi_free_model(rqrmi_model_t* rqrmi_model) {
	if (rqrmi_model == NULL)
		return;

	// For each stage
	if (rqrmi_model->stages) {
		for (uint32_t i=0; i<rqrmi_model->num_of_stages; ++i) {
			for (uint32_t j=0; j<rqrmi_model->stages[i].num_of_models; ++j) {
				free_model(rqrmi_model, i, j);
			}
			free(rqrmi_model->stages[i].models);
		}
	}
	model_info("Free RQRMI fields");
	free_matrix(rqrmi_model->input_placeholder);
	free(rqrmi_model->stages);
	free(rqrmi_model->scratchpads);
	free(rqrmi_model->error_list);
	free(rqrmi_model);
}

/**
 * @brief Frees the memory allocated for a model
 * @param rqrmi_model The RQRMI model
 * @param stage_index The stage of the model
 * @param model_index The index of the model within the stage
 */
void free_model(rqrmi_model_t* rqrmi_model, uint32_t stage_index, uint32_t model_index) {

	rqrmi_submodel_t* model = &rqrmi_model->stages[stage_index].models[model_index];
	if (model == RQRMI_MODEL_ERROR) {
		warning("Cannot free model <" << stage_index << "," << model_index <<
				">: cannot acquire model");
		return;
	}
	model_info("Free submodel <" << stage_index << "," << model_index <<
			">: Number of layers: " << model->num_of_layers)

	// For each layer
	for (uint32_t k=0; k<model->num_of_layers; ++k) {
		free_matrix(model->biases[k]);
		free_matrix(model->weights[k]);
	}

	model_info("Free submodel <" << stage_index << "," << model_index << ">: Done");
}

/**
 * @brief Preprocesses the input of a submodel
 * @param submodel The submodel to be fed with the input
 * @param input The input scalar of the submodel
 * @returns Scalar after preprocessing
 */
scalar_t rqrmi_submodel_preprocess_input(rqrmi_submodel_t* submodel, scalar_t input) {
	return (input-submodel->input_mean)/submodel->input_stddev;
}

/**
 * @brief Reverses the preprocessing of a submodel input
 * @param submodel The submodel to be fed with the input
 * @param processed_input The processed input
 * @returns Scalar before preprocessing
 */
scalar_t rqrmi_submodel_reverse_process(rqrmi_submodel_t* submodel, scalar_t processed_input) {
	return processed_input*submodel->input_stddev+submodel->input_mean;
}


// Getters & Setters

/**
 * @brief Returns number of stages in an RQRMI model
 */
uint32_t rqrmi_get_num_of_stages(rqrmi_model_t* rqrmi_model) {
	if (rqrmi_model == RQRMI_MODEL_ERROR) return 0;
	return rqrmi_model->num_of_stages;
}

/**
 * @brief Returns number of submodel in an RQRMI stage
 */
uint32_t rqrmi_get_num_of_submodels(rqrmi_model_t* rqrmi_model, uint32_t stage) {
	return rqrmi_model->stages[stage].num_of_models;
}

/**
 * @brief Returns the RQRMI last submodel index after evaluation
 */
uint32_t rqrmi_get_last_sub_model_index(rqrmi_model_t* rqrmi_model) {
	return rqrmi_model->last_model;
}

/**
 * @brief Returns the RQRMI last error
 */
uint32_t rqrmi_get_last_error(rqrmi_model_t* rqrmi_model) {
	return rqrmi_model->error_list[rqrmi_model->last_model];
}

/**
 * @brief Returns the RQRMI error of a specific submodel
 * @param[out] list The error list
 * @param[out] size Number of elements in list
 */
void rqrmi_get_error_list(rqrmi_model_t* rqrmi_model, const uint32_t** list, uint32_t* size) {
	*list = rqrmi_model->error_list;
	*size = rqrmi_model->stages[rqrmi_model->num_of_stages-1].num_of_models;
}

/**
 * @brief Sets number of stages in an RQRMI model
 */
void rqrmi_set_num_of_stages(rqrmi_model_t* rqrmi_model, uint32_t value) {
	rqrmi_model->num_of_stages = value;
}

/**
 * @brief Returns the input domain of an RQRMI model
 * @param rqrmi_model An RQRMI model
 * @return {min, max} scalar pair
 */
scalar_pair_t rqrmi_get_input_domain(rqrmi_model_t* rqrmi_model) {
	return (scalar_pair_t){ .first = rqrmi_model->input_domain_min, .second = rqrmi_model->input_domain_max };
}

/**
 * @brief Get the requested submodel compilation state
 * @param rqrmi_model An RQRMI model
 * @param stage_idx The required stage's index
 * @param submodel_idx The required submodel's index
 * @returns 1 if the submodel is compiled, or 0 if the submodel is not available / not compiled
 */
uint32_t rqrmi_get_submodel_complie_state(rqrmi_model_t* rqrmi_model, uint32_t stage_idx, uint32_t submodel_idx) {

    // Check inputs
	if ( rqrmi_model == RQRMI_MODEL_ERROR ||
		 stage_idx >= rqrmi_model->num_of_stages ||
		 submodel_idx >= rqrmi_model->stages[stage_idx].num_of_models)
	{
		return 0;
	}

	return (rqrmi_model->stages[stage_idx].models[submodel_idx].biases == NULL) ? 0 : 1;
}

/**
 * @brief Returns all the necessary information regarding an RQRMI submodel
 * @param rqrmi_model An RQRMI model
 * @param stage_idx The required stage's index
 * @param submodel_idx The required submodel's index
 * @param[out] output The output information
 * @returns 1 on success, 0 on error
 * @note RMI will not work here, only used for RQRMI models
 */
uint32_t rqrmi_get_submodel_info(rqrmi_model_t* rqrmi_model, uint32_t stage_idx, uint32_t submodel_idx, rqrmi_submodel_info_t* output) {

	// Check inputs
	if ( rqrmi_model == RQRMI_MODEL_ERROR ||
		 stage_idx >= rqrmi_model->num_of_stages ||
		 submodel_idx >= rqrmi_model->stages[stage_idx].num_of_models)
	{
		info("invalid inputs");
		return 0;
	}

	// Get submodel
	rqrmi_submodel_t *submodel = &rqrmi_model->stages[stage_idx].models[submodel_idx];

	// Validate compilation
	if (submodel->biases == NULL) {
		// Dummy values, to avoid devision by zero or any other bugs
		output->output_factor = 1;
		output->output_min = 1;
		output->input_mean = 1;
		output->input_stddev = 1;
		output->error = 0;
		// Submodel is not compiled
		output->compiled = 0;
		return 1;
	}

	// Validate RQRMI size
	if ( submodel->num_of_layers != 3 ||
		 submodel->biases[0]->cols != 1 ||
		 submodel->biases[1]->cols != 8 ||
		 submodel->biases[2]->cols != 1 )
	{
		info("input is not valid RQRMI model (maybe RMI?)");
		return 0;
	}

	// Copy info
	output->compiled = 1;
	output->b0 = GET_SCALAR(submodel->biases[0], 0, 0);
	output->w0 = GET_SCALAR(submodel->weights[0], 0, 0);

	for (uint32_t k=0; k<8; ++k) {
		output->b1[k] = GET_SCALAR(submodel->biases[1], 0 ,k);
		output->w1[k] = GET_SCALAR(submodel->weights[1], 0 ,k);
		output->w2[k] = GET_SCALAR(submodel->weights[2], k ,0);
	}
	output->b2 = GET_SCALAR(submodel->biases[2], 0 ,0);
	output->output_factor = submodel->output_factor;
	output->output_min = submodel->output_min;
	output->input_mean = submodel->input_mean;
	output->input_stddev = submodel->input_stddev;

	// Update error
	if (stage_idx == rqrmi_model->num_of_stages -1) {
		output->error = rqrmi_model->error_list[submodel_idx];
	}
	return 1;
}

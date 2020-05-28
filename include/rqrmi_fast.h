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

#include <stdexcept>
#include <basic_types.h>
#include <rqrmi_model.h>
#include <simd_aux.h>

/**
 * @brief Used to pass information both float and integer vectors
 * @note Must be cached aligned, as AVX loads/stores must be aligned
 */
typedef union {
	scalar_t scalars[SIMD_WIDTH];
	uint32_t integers[SIMD_WIDTH];
} CACHE_ALIGNED wide_scalar_t;

class RQRMIFast {
private:

	/**
	 * @brief Used to store immediate vectors of RQRMI models
	 * @note RQRMI models have at most 8 scalars per vector
	 */
	typedef union {
		scalar_t scalars[8];
		uint32_t integers[8];
	} CACHE_ALIGNED rqrmi_vector_t;

	/**
	 * @brief Holds information of a single submodel
	 */
	typedef struct {
		uint8_t compiled;
		uint32_t error;
		scalar_t w0;
		scalar_t b0;
		rqrmi_vector_t w1;
		rqrmi_vector_t b1;
		rqrmi_vector_t w2;
		rqrmi_vector_t b2;
		scalar_t output_factor;
		scalar_t output_min;
		scalar_t input_mean;
		scalar_t input_stddev;
	} fast_submodel_t CACHE_ALIGNED;

	typedef fast_submodel_t* submodel_collection_t[SIMD_WIDTH];

	// Stage information
	uint32_t _num_of_stages;
	uint32_t *_stage_submodels;

	// Submodel information
	fast_submodel_t* _submodles;

public:

	/**
	 * @brief Returns the input width of the available SIMD engine for evaluation
	 */
	static constexpr uint32_t input_width() { return SIMD_WIDTH; }

	/**
	 * @brief Create new RQRMIFast instance from RQRMI model
	 * @throws std::runtime_error in case the model is not valid RQRMI model
	 */
	RQRMIFast(rqrmi_model_t *model);
	~RQRMIFast();

	/**
	 * @brief Evaluate fast RQRMI models using SIMD acceleration
	 * @param[in] inputs a vector of inputs
	 * @param[out] status a vector of output status (1 valid, 0 error)
	 * @param[out] output a vector of outputs
	 * @param[out] error a vector of error values
	 * @note The SIMD acceleration method should be set in compilation time (FMA, AVX512, AVX, SSE, NO_RQRMI_OPT)
	 */
	void evaluate(wide_scalar_t& inputs, wide_scalar_t& status, wide_scalar_t& output, wide_scalar_t& error) const;
};


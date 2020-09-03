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

#include <rqrmi_fast.h>
#include <logging.h>

// 0x11111111 float32 representation
#define BINARY_ONES 1.14437421e-28

// Model information is often not required for debugging.
// Set dedicated debugging flag for models
#ifdef DEBUG_MODEL
#  define model_info(...) info(__VA_ARGS__)
#else
#  define model_info(...)
#endif

/**
 * @brief Prints wide scalar type to screen
 */
void print_wide_scalar(wide_scalar_t& x) {
	for (uint32_t i=0; i<SIMD_WIDTH; ++i) {
		fprintf(stderr, "%f ", x.scalars[i]);
	}
	fprintf(stderr, "\n");
}

/**
 * @brief Create new RQRMIFast instance from RQRMI model
 * @throws std::runtime_error in case the model is not valid RQRMI model
 */
RQRMIFast::RQRMIFast(rqrmi_model_t *model) {
	_num_of_stages = rqrmi_get_num_of_stages(model);

	if (_num_of_stages > 3) {
		throw std::runtime_error("input model have more than three stages!");
	}

	model_info("Allocating data for %u stages", _num_of_stages);

	uint32_t total_submodels = 1;
	_stage_submodels = new uint32_t[_num_of_stages];

	// Initiate number of submodels
	for (uint32_t s=0; s<_num_of_stages; ++s) {
		_stage_submodels[s] = rqrmi_get_num_of_submodels(model, s);
		if (s > 0) {
			total_submodels += total_submodels*_stage_submodels[s];
		}
	}

	// Allocate memory
	model_info("Allocating data for %u submodels (in total)", total_submodels);
	_submodles = (fast_submodel_t*)aligned_alloc(64,
			sizeof(fast_submodel_t) * total_submodels);

	// Get submodel information
	uint32_t counter = 0;
	for (uint32_t s=0; s<_num_of_stages; ++s) {
		for (uint32_t m=0; m<_stage_submodels[s]; ++m) {
			rqrmi_submodel_info_t info;
			if (!rqrmi_get_submodel_info(model, s, m, &info)) {
				throw std::runtime_error("error while extracting information of a submodel");
			}
			// Copy info
			_submodles[counter].compiled = info.compiled;
			_submodles[counter].b0 = info.b0;
			_submodles[counter].w0 = info.w0;
			_submodles[counter].output_factor = info.output_factor;
			_submodles[counter].error = info.error;
			_submodles[counter].output_min = info.output_min;
			_submodles[counter].input_mean = info.input_mean;
			_submodles[counter].input_stddev = info.input_stddev;
			for (uint32_t i=0; i<8; ++i) {
				_submodles[counter].b1.scalars[i] = info.b1[i];
				_submodles[counter].w1.scalars[i] = info.w1[i];
				_submodles[counter].b2.scalars[i] = info.b2 / input_width(); // Menachem's tip for including this in reduce
				_submodles[counter].w2.scalars[i] = info.w2[i];
			}
			++counter;
		}
	}
}

RQRMIFast::~RQRMIFast() {
	delete[] _stage_submodels;
	free(_submodles);
}

/**
 * @brief Evaluate fast RQRMI models using SIMD acceleration
 * @param[in] inputs a vector of inputs
 * @param[out] status a vector of output status (1 valid, 0 error)
 * @param[out] output a vector of outputs
 * @param[out] error a vector of error values
 * @note The SIMD acceleration method should be set in compilation time (FMA, AVX512, AVX, SSE, NO_RQRMI_OPT)
 */
void RQRMIFast::evaluate(wide_scalar_t& inputs, wide_scalar_t& status, wide_scalar_t& output, wide_scalar_t& error) const {
	// Base index for submodel in array
	uint32_t base_idx = 0;
	// Next index of submodel in stage, per input
	wide_scalar_t next_idx;
	// A collection of submodels
	submodel_collection_t submodels {};

	// Holds a vector of submodels input normalization factors
	wide_scalar_t input_mean, input_stddev;

	// Holds a vector of submodels output post-processing factors
	wide_scalar_t output_factor, output_min;

	// Holds a vector of layer0 variables
	wide_scalar_t w0, b0;

#ifdef NO_RQRMI_OPT
	// Used for intermediate computations
	float reg0;

	// Used for ReLUs operations
	static float zeros = 0;
	static float ones = 1-SCALAR_EPS;

	// Initiate status
	status.scalars[0] = BINARY_ONES;
	next_idx.scalars[0] = 0;

#elif __AVX512F__
	// Used for intermediate computations
	__m512 reg0, reg1, reg2, reg3, reg4;

	// Used for ReLUs operations
	static const __m512 zeros = _mm512_setzero_ps();
	static const __m512 ones = _mm512_set1_ps(1-SCALAR_EPS);

	// Used for masked load operations
	static const __mmask16 mask = _mm512_int2mask(255);

	// Initiate status
	reg0 =_mm512_set1_ps(BINARY_ONES);
	_mm512_store_ps(status.scalars, reg0);
	_mm512_store_ps(next_idx.scalars, zeros);

#elif __AVX__
	// Used for intermediate computations
	__m256 reg0, reg1, reg2, reg3, reg4;

	// Used for ReLUs operations
	static const __m256 zeros = _mm256_setzero_ps();
	static const __m256 ones = _mm256_set1_ps(1-SCALAR_EPS);

	// Initiate status
	reg0 =_mm256_set1_ps(BINARY_ONES);
	_mm256_store_ps(status.scalars, reg0);
	_mm256_store_ps(next_idx.scalars, zeros);

#elif __SSE__
	// Used for intermediate computations
	__m128 reg0, reg1, reg2, reg3, reg4, reg5, reg6;

	// Used for ReLUs operations
	static const __m128 zeros = _mm_setzero_ps();
	static const __m128 ones = _mm_set1_ps(1-SCALAR_EPS);

	// Initiate status
	reg0 = _mm_set1_ps(BINARY_ONES); // 0x11111111 float32 representation
	_mm_store_ps(status.scalars, reg0);
	_mm_store_ps(next_idx.scalars, zeros);

#endif

	// For each stage
	for (uint32_t i=0; i<_num_of_stages; ++i) {

		// Get the address of each in the collection next submodels
		for (uint32_t j=0; j<SIMD_WIDTH; ++j) {
			submodels[j] = &_submodles[base_idx + (uint32_t)(_stage_submodels[i] * next_idx.scalars[j]) ];
			// Update status according to submodel compile status
			status.integers[j] &= submodels[j]->compiled;
			// Update input normalization factors
			input_mean.scalars[j] = submodels[j]->input_mean;
			input_stddev.scalars[j] = submodels[j]->input_stddev;
			// Update layer0
			w0.scalars[j] = submodels[j]->w0;
			b0.scalars[j] = submodels[j]->b0;
			// Update post-processing factors
			output_factor.scalars[j] = submodels[j]->output_factor;
			output_min.scalars[j] = submodels[j]->output_min;
		}

#ifdef NO_RQRMI_OPT
		// Preprocess input
		reg0 = (inputs.scalars[0] - input_mean.scalars[0]) / input_stddev.scalars[0];

		model_info("Input for submodel after preprocessing: %.12f", reg0);

		// Compute layer0
		FMA(reg0, w0.scalars[0], b0.scalars[0]);
		float result = 0;

#elif __AVX512F__
		// Load registers
		reg0 = _mm512_load_ps(inputs.scalars);
		reg1 = _mm512_load_ps(input_mean.scalars);
		reg2 = _mm512_load_ps(input_stddev.scalars);
		reg3 = _mm512_load_ps(w0.scalars);
		reg4 = _mm512_load_ps(b0.scalars);

		// Preprocess input collection
		reg0 = _mm512_sub_ps(reg0, reg1);
		reg0 = _mm512_div_ps(reg0, reg2);

		// Compute layer0
		FMA(reg0, reg3, reg4);

		// We wish to extract input per submodel
		wide_scalar_t result;
		_mm512_storeu_ps(result.scalars, reg0);

#elif __AVX__
		// Load registers
		reg0 = _mm256_load_ps(inputs.scalars);
		reg1 = _mm256_load_ps(input_mean.scalars);
		reg2 = _mm256_load_ps(input_stddev.scalars);
		reg3 = _mm256_load_ps(w0.scalars);
		reg4 = _mm256_load_ps(b0.scalars);

		// Preprocess input collection
		reg0 = _mm256_sub_ps(reg0, reg1);
		reg0 = _mm256_div_ps(reg0, reg2);

		// Compute layer0
		FMA(reg0, reg3, reg4);

		// We wish to extract input per submodel
		wide_scalar_t result;
		_mm256_storeu_ps(result.scalars, reg0);
#elif __SSE__
		// Load registers
		reg0 = _mm_load_ps(inputs.scalars);
		reg1 = _mm_load_ps(input_mean.scalars);
		reg2 = _mm_load_ps(input_stddev.scalars);
		reg3 = _mm_load_ps(w0.scalars);
		reg4 = _mm_load_ps(b0.scalars);

		// Preprocess input collection
		reg0 = _mm_sub_ps(reg0, reg1);
		reg0 = _mm_div_ps(reg0, reg2);

		// Compute layer0
		FMA(reg0, reg3, reg4);

		// We wish to extract input per submodel
		wide_scalar_t result;
		_mm_storeu_ps(result.scalars, reg0);
#endif

		// For each submodel
		for (uint32_t j=0; j<input_width(); ++j) {
#ifdef NO_RQRMI_OPT
			// Compute layer 1
			rqrmi_vector_t vector;
			for (int k=0; k<8; ++k) {
				vector.scalars[k] = (reg0 * submodels[j]->w1.scalars[k]) + submodels[j]->b1.scalars[k];
				if (vector.scalars[k] < 0) vector.scalars[k] = 0; // ReLU
			}

			model_info("Layer1 output: [%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f]",
					vector.scalars[0], vector.scalars[1], vector.scalars[2], vector.scalars[3],
					vector.scalars[4], vector.scalars[5], vector.scalars[6], vector.scalars[7]);

			// Compute layer2 & reduce output
			for (int k=0; k<8; ++k) {
				result += (vector.scalars[k] * submodels[j]->w2.scalars[k]);
			}
			result += submodels[j]->b2.scalars[0];

#elif __AVX512F__
			// Set current input
			reg1 = _mm512_set1_ps(result.scalars[j]);

			// Load registers
			reg2 = _mm512_mask_load_ps(reg2, mask, submodels[j]->w1.scalars);
			reg3 = _mm512_mask_load_ps(reg2, mask, submodels[j]->b1.scalars);

			// Compute layer1
			FMA(reg1, reg2, reg3);
			reg1 = _mm512_max_ps(reg1, zeros); // ReLU

			// Load registers
			reg2 = _mm512_mask_load_ps(reg2, mask, submodels[j]->w2.scalars);
			reg3 = _mm512_mask_load_ps(reg2, mask, submodels[j]->b2.scalars);

			// Compute layer2
			FMA(reg1, reg2, reg3);

			// Reduce output (https://stackoverflow.com/a/13222410/4103200)
			// reg1 = ( -, -, -, -, -, -, -, -, x7, x6, x5, x4, x3, x2, x1, x0 )
			const __m128 hiQuad = _mm512_extractf32x4_ps(reg1, 1); 	// ( x7, x6, x5, x4 )
			const __m128 loQuad = _mm512_castps512_ps128(reg1);	// ( x3, x2, x1, x0 )
			const __m128 sumQuad = _mm_add_ps(loQuad, hiQuad);		// ( x3 + x7, x2 + x6, x1 + x5, x0 + x4 )
			const __m128 loDual = sumQuad;							// ( -, -, x1 + x5, x0 + x4 )
			const __m128 hiDual = _mm_movehl_ps(sumQuad, sumQuad);  // ( -, -, x3 + x7, x2 + x6 )
			const __m128 sumDual = _mm_add_ps(loDual, hiDual);		// ( -, -, x1 + x3 + x5 + x7, x0 + x2 + x4 + x6 )
			const __m128 lo = sumDual;								// ( -, -, -, x0 + x2 + x4 + x6 )
			const __m128 hi = _mm_shuffle_ps(sumDual, sumDual, 0x1);// ( -, -, -, x1 + x3 + x5 + x7 )
			const __m128 sum = _mm_add_ss(lo, hi);					// ( -, -, -, x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 )
			result.scalars[j] = _mm_cvtss_f32(sum);
#elif __AVX__
			// Set current input
			reg1 = _mm256_set1_ps(result.scalars[j]);

			// Load registers
			reg2 = _mm256_load_ps(submodels[j]->w1.scalars);
			reg3 = _mm256_load_ps(submodels[j]->b1.scalars);

			// Compute layer1
			FMA(reg1, reg2, reg3);
			reg1 = _mm256_max_ps(reg1, zeros); // ReLU

			// Load registers
			reg2 = _mm256_load_ps(submodels[j]->w2.scalars);
			reg3 = _mm256_load_ps(submodels[j]->b2.scalars);

			// Compute layer2
			FMA(reg1, reg2, reg3);

			// Reduce output (https://stackoverflow.com/a/13222410/4103200)
			// reg1 = ( x7, x6, x5, x4, x3, x2, x1, x0 )
			const __m128 hiQuad = _mm256_extractf128_ps(reg1, 1); 	// ( x7, x6, x5, x4 )
			const __m128 loQuad = _mm256_castps256_ps128(reg1);	// ( x3, x2, x1, x0 )
			const __m128 sumQuad = _mm_add_ps(loQuad, hiQuad);		// ( x3 + x7, x2 + x6, x1 + x5, x0 + x4 )
			const __m128 loDual = sumQuad;							// ( -, -, x1 + x5, x0 + x4 )
			const __m128 hiDual = _mm_movehl_ps(sumQuad, sumQuad);  // ( -, -, x3 + x7, x2 + x6 )
			const __m128 sumDual = _mm_add_ps(loDual, hiDual);		// ( -, -, x1 + x3 + x5 + x7, x0 + x2 + x4 + x6 )
			const __m128 lo = sumDual;								// ( -, -, -, x0 + x2 + x4 + x6 )
			const __m128 hi = _mm_shuffle_ps(sumDual, sumDual, 0x1);// ( -, -, -, x1 + x3 + x5 + x7 )
			const __m128 sum = _mm_add_ss(lo, hi);					// ( -, -, -, x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 )
			result.scalars[j] = _mm_cvtss_f32(sum);

#elif __SSE__
			// Set current input
			reg1 = _mm_set1_ps(result.scalars[j]);
			reg2 = _mm_set1_ps(result.scalars[j]);

			// Load registers
			reg3 = _mm_load_ps(submodels[j]->w1.scalars);
			reg4 = _mm_load_ps(&submodels[j]->w1.scalars[4]);

			reg5 = _mm_load_ps(submodels[j]->b1.scalars);
			reg6 = _mm_load_ps(&submodels[j]->b1.scalars[4]);

			// Compute layer1
			FMA(reg1, reg3, reg5);
			FMA(reg2, reg4, reg6);

			reg1 = _mm_max_ps(reg1, zeros); // ReLU
			reg2 = _mm_max_ps(reg2, zeros); // ReLU

			// Load registers
			reg3 = _mm_load_ps(submodels[j]->w2.scalars);
			reg4 = _mm_load_ps(&submodels[j]->w2.scalars[4]);

			reg5 = _mm_load_ps(submodels[j]->b2.scalars);

			// Compute layer2
			FMA(reg1, reg3, reg5);
			FMA(reg2, reg4, zeros);

			// Reduce output (https://stackoverflow.com/a/13222410/4103200)
			// reg1 = ( x3, x2, x1, x0 )
			// reg2 = ( x7, x6, x5, x4 )
			const __m128 sumQuad = _mm_add_ps(reg1, reg2);		// ( x3 + x7, x2 + x6, x1 + x5, x0 + x4 )
			const __m128 loDual = sumQuad;							// ( -, -, x1 + x5, x0 + x4 )
			const __m128 hiDual = _mm_movehl_ps(sumQuad, sumQuad);  // ( -, -, x3 + x7, x2 + x6 )
			const __m128 sumDual = _mm_add_ps(loDual, hiDual);		// ( -, -, x1 + x3 + x5 + x7, x0 + x2 + x4 + x6 )
			const __m128 lo = sumDual;								// ( -, -, -, x0 + x2 + x4 + x6 )
			const __m128 hi = _mm_shuffle_ps(sumDual, sumDual, 0x1);// ( -, -, -, x1 + x3 + x5 + x7 )
			const __m128 sum = _mm_add_ss(lo, hi);					// ( -, -, -, x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 )
			result.scalars[j] = _mm_cvtss_f32(sum);
#endif
		}

#ifdef NO_RQRMI_OPT
		model_info("Submodel output before post-processing: %.12f", result);

		// Post-process outputs
		FMA(result, output_factor.scalars[0], output_min.scalars[0]);

		// Used for micro-debugging
		model_info("Final result for stage %u: %.12f", i, result);

		if (result < zeros) result = zeros;
		if (result > ones) result = ones;

		// Update next indices
		next_idx.scalars[0] = result;
#elif __AVX512F__
		// Re-load results to register0
		reg0 = _mm512_load_ps(result.scalars);

		// Load registers
		reg1 = _mm512_load_ps(output_factor.scalars);
		reg2 = _mm512_load_ps(output_min.scalars);

		// Post-process outputs
		FMA(reg0, reg1, reg2);

		reg0 = _mm512_max_ps(reg0, zeros);
		reg0 = _mm512_min_ps(reg0, ones);

		// Update next indices
		_mm512_storeu_ps(next_idx.scalars, reg0);

#elif __AVX__
		// Re-load results to register0
		reg0 = _mm256_load_ps(result.scalars);

		// Load registers
		reg1 = _mm256_load_ps(output_factor.scalars);
		reg2 = _mm256_load_ps(output_min.scalars);

		// Post-process outputs
		FMA(reg0, reg1, reg2);

		reg0 = _mm256_max_ps(reg0, zeros);
		reg0 = _mm256_min_ps(reg0, ones);

		// Update next indices
		_mm256_storeu_ps(next_idx.scalars, reg0);

#elif __SSE__
		// Re-load results to register0
		reg0 = _mm_load_ps(result.scalars);

		// Load registers
		reg1 = _mm_load_ps(output_factor.scalars);
		reg2 = _mm_load_ps(output_min.scalars);

		// Post-process outputs
		FMA(reg0, reg1, reg2);

		reg0 = _mm_max_ps(reg0, zeros);
		reg0 = _mm_min_ps(reg0, ones);

		// Update next indices
		_mm_storeu_ps(next_idx.scalars, reg0);
#endif
		// Update base index
		base_idx += _stage_submodels[i];
	}

	// Load the error vector
	for (uint32_t j=0; j<input_width(); ++j) {
		error.integers[j] = submodels[j]->error;
	}

	// Return result
	output = next_idx;

}

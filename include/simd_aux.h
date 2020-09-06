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

#include <sstream>

// Intrinsics
#include <x86intrin.h>

/**
 * @brief For clean code of fast multiple add
 */
#ifdef NO_RQRMI_OPT
	#define FMA(a,b,c) 				\
		a = (a*b)+c;
#elif __AVX512F__
	#define FMA(a,b,c) 				\
		a = _mm512_fmadd_ps(a, b, c);
#elif __FMA__
	#define FMA(a,b,c) 				\
		a = _mm256_fmadd_ps(a, b, c);
#elif __AVX__ // NO_SIMD
	#define FMA(a,b,c) 				\
		a = _mm256_mul_ps(a, b);	\
		a = _mm256_add_ps(a, c);
#elif __SSE__ // NO_SIMD
	#define FMA(a,b,c) 				\
		a = _mm_mul_ps(a, b);		\
		a = _mm_add_ps(a, c);
#endif // NO_SIMD

/**
 * @brief General defines that change according to available SIMD engine
 */
#ifdef NO_RQRMI_OPT
#	define SIMD_WIDTH 1
#	define PS_REG float
#	define EPU_REG unsigned int
#elif __AVX512F__ // NO_SIMD
#	define SIMD_WIDTH 16
#	define PS_REG __m512
#	define EPU_REG __m512i
#elif __AVX__ // NO_SIMD
#	define SIMD_WIDTH 8
#	define PS_REG __m256
#	define EPU_REG __m256i
#elif __SSE__ //NO_SIMD
#	define SIMD_WIDTH 4
#	define PS_REG __m128
#	define EPU_REG __m128i
#endif // NO_SIMD

/**
 * @brief Align to 64bytes
 */
#define CACHE_ALIGNED __attribute__ ((aligned (64)))

/**
 * @brief Returns a string representation of single-precision float vector
 */
std::string simd_ps_vector_logger(PS_REG v);

/**
 * @brief Returns a string representation of unsigned int vector
 */
std::string simd_epu_vector_logger(EPU_REG v);


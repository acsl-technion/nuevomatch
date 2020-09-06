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

#include <simd_aux.h>
#include <basic_types.h>

/**
 * @brief Returns a string representation of single-precision float vector
 */
std::string simd_ps_vector_logger(PS_REG v) {
#ifndef NO_RQRMI_OPT
	union {float f; int i;} fp;
#endif
	std::stringstream ss;
#ifdef NO_RQRMI_OPT
	ss<<"["<<v<<"]";
#elif __AVX512F__
	float f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,0),0); f0=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,0),1); f1=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,0),2); f2=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,0),3); f3=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,1),0); f4=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,1),1); f5=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,1),2); f6=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,1),3); f7=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,2),0); f8=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,2),1); f9=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,2),2); f10=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,2),3); f11=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,3),0); f12=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,3),1); f13=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,3),2); f14=fp.f;
	fp.i = _mm_extract_ps(_mm512_extractf32x4_ps(v,3),3); f15=fp.f;
	ss<<"["<<f0<<","<<f1<<","<<f2<<","<<f3<<","<<f4<<","
			<<f5<<","<<f6<<","<<f7<<","<<f8<<","<<f9<<","<<f10<<
			","<<f11<<","<<f12<<","<<f13<<","<<f14<<","<<f15<<"]";
#elif __AVX__
	float f0,f1,f2,f3,f4,f5,f6,f7;
	fp.i = _mm_extract_ps(_mm256_extractf128_ps(v,0),0); f0=fp.f;
	fp.i = _mm_extract_ps(_mm256_extractf128_ps(v,0),1); f1=fp.f;
	fp.i = _mm_extract_ps(_mm256_extractf128_ps(v,0),2); f2=fp.f;
	fp.i = _mm_extract_ps(_mm256_extractf128_ps(v,0),3); f3=fp.f;
	fp.i = _mm_extract_ps(_mm256_extractf128_ps(v,1),0); f4=fp.f;
	fp.i = _mm_extract_ps(_mm256_extractf128_ps(v,1),1); f5=fp.f;
	fp.i = _mm_extract_ps(_mm256_extractf128_ps(v,1),2); f6=fp.f;
	fp.i = _mm_extract_ps(_mm256_extractf128_ps(v,1),3); f7=fp.f;
	ss<<"["<<f0<<","<<f1<<","<<f2<<","<<f3<<","<<f4<<","
			<<f5<<","<<f6<<","<<f7<<"]";
#elif __SSE__
	float f0,f1,f2,f3;
	fp.i = _mm_extract_ps(v,0); f0=fp.f;
	fp.i = _mm_extract_ps(v,1); f1=fp.f;
	fp.i = _mm_extract_ps(v,2); f2=fp.f;
	fp.i = _mm_extract_ps(v,3); f3=fp.f;
	ss <<"["<<f0<<","<<f1<<","<<f2<<","<<f3<<"]";
#endif
	return ss.str();
}

/**
 * @brief Returns a string representation of unsigned int vector
 */
std::string simd_epu_vector_logger(EPU_REG v) {
#ifndef NO_RQRMI_OPT
	union {unsigned int f; int i;} fp;
#endif
	std::stringstream ss;
#ifdef NO_RQRMI_OPT
	ss<<"["<<v<<"]";
#elif __AVX512F__
	unsigned int f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,0),0); f0=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,0),1); f1=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,0),2); f2=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,0),3); f3=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,1),0); f4=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,1),1); f5=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,1),2); f6=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,1),3); f7=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,2),0); f8=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,2),1); f9=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,2),2); f10=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,2),3); f11=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,3),0); f12=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,3),1); f13=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,3),2); f14=fp.f;
	fp.i = _mm_extract_epi32(_mm512_extracti32x4_epi32(v,3),3); f15=fp.f;
	ss<<"["<<f0<<","<<f1<<","<<f2<<","<<f3<<","<<f4<<","
			<<f5<<","<<f6<<","<<f7<<","<<f8<<","<<f9<<","<<f10<<
			","<<f11<<","<<f12<<","<<f13<<","<<f14<<","<<f15<<"]";
#elif __AVX__
	unsigned int f0,f1,f2,f3,f4,f5,f6,f7;
	fp.i = _mm_extract_epi32(_mm256_extracti128_si256(v,0),0); f0=fp.f;
	fp.i = _mm_extract_epi32(_mm256_extracti128_si256(v,0),1); f1=fp.f;
	fp.i = _mm_extract_epi32(_mm256_extracti128_si256(v,0),2); f2=fp.f;
	fp.i = _mm_extract_epi32(_mm256_extracti128_si256(v,0),3); f3=fp.f;
	fp.i = _mm_extract_epi32(_mm256_extracti128_si256(v,1),0); f4=fp.f;
	fp.i = _mm_extract_epi32(_mm256_extracti128_si256(v,1),1); f5=fp.f;
	fp.i = _mm_extract_epi32(_mm256_extracti128_si256(v,1),2); f6=fp.f;
	fp.i = _mm_extract_epi32(_mm256_extracti128_si256(v,1),3); f7=fp.f;
	ss<<"["<<f0<<","<<f1<<","<<f2<<","<<f3<<","<<f4<<","
			<<f5<<","<<f6<<","<<f7<<"]";
#elif __SSE__
	unsigned int f0,f1,f2,f3;
	fp.i = _mm_extract_epi32(v,0); f0=fp.f;
	fp.i = _mm_extract_epi32(v,1); f1=fp.f;
	fp.i = _mm_extract_epi32(v,2); f2=fp.f;
	fp.i = _mm_extract_epi32(v,3); f3=fp.f;
	ss <<"["<<f0<<","<<f1<<","<<f2<<","<<f3<<"]";
#endif
	return ss.str();
}

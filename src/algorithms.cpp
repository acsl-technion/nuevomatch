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

#include <algorithms.h>

#include <iostream>
#include <random>
#include <cmath>

#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

using namespace std;

default_random_engine generator;
random_device rd;

/**
 * @brief Generates uniform random numbers between low and high
 * @param The type of the return type.
 * @param low Low bound (inclusive)
 * @param high High bound (inclusive)
 * @note Taken from here https://stackoverflow.com/questions/21096015/how-to-generate-64-bit-random-numbers
 */
uint32_t gen_uniform_random_uint32(uint32_t low, uint32_t high) {
    uniform_int_distribution<uint32_t> dist(low, high);
    return dist(rd);
}
scalar_t gen_uniform_random_scalar(scalar_t low, scalar_t high) {
	uniform_real_distribution<scalar_t> distribution(low, high);
	return distribution(generator);
}

/**
 * @brief Create random permutation of integers
 * @param dst Permutation destination location
 * @param size The size of the permutation
 * @note Unsafe: the user should allocate dst memory location.
 */
void random_permutation(uint32_t* dst, uint32_t size) {
	srand ( uint32_t ( time(0) ) );
	vector<uint32_t> vec;
	for (uint32_t i=0; i<size; ++i) {
		vec.push_back(i);
	}
	random_shuffle(vec.begin(), vec.end());
	for (uint32_t i=0; i<size; ++i) {
		dst[i] = vec[i];
	}
}



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

/** @file algorithms.h */

#pragma once

#include <set> // std::set
#include <vector> // std::vector
#include <basic_types.h>

/**
 * @brief Create random permutation of integers
 * @param dst Permutation destination location
 * @param size The size of the permutation
 * @note Unsafe: the user should allocate dst memory location.
 */
void random_permutation(uint32_t* dst, uint32_t size);

/**
 * @brief Generates uniform random unsigned integers between low and high
 * @param The type of the return type.
 * @param low Low bound (inclusive)
 * @param high High bound (inclusive)
 * @note Taken from here https://stackoverflow.com/questions/21096015/how-to-generate-64-bit-random-numbers
 */
uint32_t gen_uniform_random_uint32(uint32_t low, uint32_t high);

/**
 * @brief Generates uniform random scalars between low and high
 * @param The type of the return type.
 * @param low Low bound (inclusive)
 * @param high High bound (inclusive)
 * @note Taken from here https://stackoverflow.com/questions/21096015/how-to-generate-64-bit-random-numbers
 */
scalar_t gen_uniform_random_scalar(scalar_t low, scalar_t high);

/**
 * @brief Compute all possible combinations of all elements in an array of sets.
 * @param sets An array of sets, each with a list of unique elements
 * @returns A vector of all possible combinations of the sets' values
 * @note Complexity: O(n * (a*b*c*...)) for n number of sets, and {a,b,c,...} are set sizes
 * 		 This is also the total number of elements in the output array.
 */
template <typename T>
std::vector<std::vector<T>> calculate_all_combinations(const std::vector<std::set<T>> &sets) {

	// How many combinations are there?
	// Complexity: O(n) for n number of sets
	uint32_t total_combinations = 1;
	for (auto current_set : sets) {
		total_combinations *= current_set.size();
	}

	// Allocate memory for all combinations
	std::vector<std::vector<T>> output(total_combinations);

	// Store iterators for all sets
	// Complexity: O(n) for n number of sets
	std::vector<typename std::set<T>::const_iterator> iterator_array;
	uint32_t num_of_sets = sets.size();
	for (uint32_t j=0; j<num_of_sets; ++j) {
		iterator_array.push_back(sets[j].cbegin());
	}

	// Determine division factor per set
	// Complexity: O(n) for n number of sets
	std::vector<uint32_t> division_factor;
	uint32_t accumulator = 1;
	for (auto current_set : sets) {
		division_factor.push_back(accumulator);
		accumulator *= current_set.size();
	}

	// For each possible combinations
	// Complexity: O(n * (a*b*c*...)) for n number of sets, and {a,b,c,...} are set sizes
	// Note: this is also the total number of elements in the output array.
	for (uint32_t i=0; i<total_combinations; ++i) {
		// Initialize the current combination to hold S elements
		output[i].resize(num_of_sets);
		// For each set
		for (uint32_t j=0; j<num_of_sets; ++j) {
			// Populate the output with current iterator
			output[i][j] = *iterator_array[j];
			// Calculate the index divided by the current division factor
			uint32_t idx = i/division_factor[j];
			// Update iterator
			if ((idx + 1) % sets[j].size() == 0) {
				++iterator_array[j];
				if (iterator_array[j] == sets[j].end()) {
					iterator_array[j] = sets[j].cbegin();
				}
			}
		}
	}

	return output;
}



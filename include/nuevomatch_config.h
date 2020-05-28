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

#include <vector>

#include <generic_classifier.h>

/**
 * @brief Holds all configuration required for a NuevoMatch instance
 */
struct NuevoMatchConfig {

	/**
	 * @brief An instance of a valid classifier to perform the remainder lookup
	 */
	GenericClassifier* remainder_classifier = nullptr;

	/**
	 * @brief The queue size for concurrent workers
	 */
	uint32_t queue_size = 256;

	/**
	 * @brief Set the number of cores to run NuevoMatch
	 */
	uint32_t num_of_cores = 1;

	/**
	 * @brief If not-negative, limit the number of iSets available to classifier
	 */
	int max_subsets = -1;

	/**
	 * @brief Run iSet to start from index X
	 */
	uint32_t start_from_iset = 0;

	/**
	 * @brief Do not execute iSet classifiers
	 */
	bool disable_isets = false;

	/**
	 * @brief Do not execute remainder classifier
	 */
	bool disable_remainder = false;

	/**
	 * @brief Do not execute binary search on iSets
	 */
	bool disable_bin_search = false;

	/**
	 * @brief Do not execute the validation phase on iSets
	 */
	bool disable_validation_phase = false;

	/**
	 * @brief Used to test the batching and publishing time,
	 *        without performing any classification
	 */
	bool disable_all_classification = false;

	/**
	 * @brief Used to signal NuevoMatch that the remainder classifier was loaded
	 *        from any external source (e.g., file). If true, will make NuevoMatch
	 *        execute the remainder as it is, without any modifications!
	 */
	bool external_remainder = false;

	/**
	 * @brief Run NuevoMatch on any subset of fields of the original dataset
	 */
	std::vector<uint32_t> arbitrary_fields;
};

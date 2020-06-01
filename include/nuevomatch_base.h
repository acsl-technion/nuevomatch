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

#include <string>
#include <vector>
#include <sstream>

#include <basic_types.h>
#include <pipeline_thread.h>
#include <generic_classifier.h>

/**
 * @brief ActionBatch is a batch of N <action, priority> tuples.
 *        It is returned by the classify operation.
 */
template <uint32_t N>
using ActionBatch = WorkBatch<classifier_output_t, N>;

/**
 * @brief PacketBatch is a batch of N packets.
 *        A packet is a pointer to a continuous memory region
 *        of all header values. The number of headers is defined
 *        within the subsets of this, and should always be the same.
 *        PacketBatch is the input for the classify operation.
 *        Invalid packets can be set using invalid pointer (NULL).
 */
template <uint32_t N>
using PacketBatch = WorkBatch<const uint32_t*, N>;


/**
 * @brief An abstract class for a NuevoMatch subset.
 *        A subsets works on packet batches with N packets,
 *        and can return various information on the subset,
 *        main for performing core allocation
 * @tparam N Number of packets in batch
 */
template <uint32_t N>
class NuevoMatchSubset {
public:

	enum dynamic_type_t {ISET=0, REMAINDER};

	virtual ~NuevoMatchSubset() {};

	/**
	 * @brief Returns an integer by which subsets are allocated to cores
	 */
	virtual uint32_t get_size() const = 0;

	/**
	 * @brief Returns a string representation of this
	 */
	virtual std::string to_string() const = 0;

	/**
	 * @brief Returns the dynamic type of this
	 */
	virtual dynamic_type_t get_type() const = 0;
};

/**
 * @brief An adapter for any GenericClassifier to be used with NuevoMatchSubset
 * @tparam N The number of packets in a single batch
 */
template <uint32_t N>
class NuevoMatchRemainderClassifier : public NuevoMatchSubset<N> {
protected:
	GenericClassifier* _classifier;

public:

	NuevoMatchRemainderClassifier(GenericClassifier* classifier) : _classifier(classifier) {};
	virtual ~NuevoMatchRemainderClassifier() { delete _classifier; }

	/**
	 * @brief Invoke the classify method of the subset
	 * @param packets A batch of packet headers
	 * @returns The result of the subset
	 */
	ActionBatch<N> classify(PacketBatch<N>& packets, ActionBatch<N>& output) {
		for (uint32_t i=0; i<N; ++i) {
			if (packets[i] == nullptr) {
				continue;
			}
			int result = _classifier->classify_sync(packets[i], output[i].priority);
			output[i] = {result, result};
		}
		return output;
	}

	/**
	 * @brief Returns an integer by which subsets are allocated to cores
	 */
	virtual uint32_t get_size() const {
		return _classifier->get_size();
	}

	/**
	 * @brief Returns a string representation of this
	 */
	virtual std::string to_string() const {
		return _classifier->to_string();
	}

	/**
	 * @brief Returns the dynamic type of this
	 */
	virtual typename NuevoMatchSubset<N>::dynamic_type_t get_type() const {
		return NuevoMatchSubset<N>::dynamic_type_t::REMAINDER;
	}
};


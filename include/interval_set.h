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

// Intrinsics
#include <x86intrin.h>

#include <vector>
#include <sstream>

#include <basic_types.h>
#include <object_io.h>
#include <matrix_operations.h>
#include <rqrmi_model.h>
#include <rqrmi_fast.h>
#include <pipeline_thread.h>
#include <nuevomatch_base.h>
#include <simd_aux.h>

/**
 * @brief The basic information unit an iSet needs
 */
typedef struct {
	scalar_t rqrmi_input;
	scalar_t rqrmi_output;
	uint32_t rqrmi_error;
	const uint32_t* header;
	uint32_t valid;
} iset_info_t;

/**
 * @brief Information on batch of elements
 * @tparam N the number of packets in a single batch
 */
template<uint32_t N>
using IntervalSetInfoBatch = WorkBatch<iset_info_t, N>;

/**
 * @brief Abstract class. An interval set object, holds rule database and a Lookup object
 * @tparam N the number of packets in a single batch
 */
template<uint32_t N>
class IntervalSet : public NuevoMatchSubset<N> {
protected:

	// Number of scalars in vector
	constexpr static uint32_t scalars_in_vector = SIMD_WIDTH;

	// The lookup table of the secondary search
	// Divided to index (fast lookup) and values (validation)
	scalar_t* _index;

	// Validation database format:
	// Matrix of rows - validation phases, columns - header columns
	// Memory is arranged by rows (!!) and not by columns.
	// Reason: easier to perform validation using SIMD operations
	// and even enabling lazy evaluation!
	uint32_t* _validation_db;

	// The RQRMI model
	rqrmi_model_t* _model;
	RQRMIFast* _model_fast;

	// iSet index within classifier
	uint32_t _iset_index;

	// Additional information
	uint32_t _size;
	uint32_t _field_index;
	uint32_t _num_of_columns;
	uint32_t _num_of_validation_phases;
	uint32_t _size_kb;

	// Measure expected error
	double _search_counter, _error_counter;

	// Used for validation phase calculation
	uint32_t _F;
	uint32_t _vec_ops;
	uint32_t _remainder_size;
	uint32_t _rule_size;
	uint32_t _validation_low_size;
	uint32_t _remainder_mask[scalars_in_vector];

	/**
	 * @brief Updates all parameters required for the validation phase
	 */
	void update_vlidation_phase_params();

public:

	/**
	 * @brief Initialize new interval-set of specific index
	 * @param index The iSet index within classifier
	 */
	IntervalSet(uint32_t index);
	virtual ~IntervalSet();

	/**
	 * @brief Load an interval set from file
	 * @param object An object reader with binary information
	 * @throws In case the index is not ordered or in case of internal error
	 */
	void load(ObjectReader& object);

	/**
	 * @brief Extract the rules of this as a vector of OpenFlow rules
	 * @note Due to validation expansion, the number of output rules may be larger than expected
	 */
	std::vector<openflow_rule> extract_rules() const;

	/**
	 * @brief Search for packet within the interval set
	 * @param packets A batch of packets
	 * @returns A batch of tuples, each includes the matching rule's priority and action
	 */
	IntervalSetInfoBatch<N> rqrmi_search(PacketBatch<N>& packets) const;

	/**
	 * @brief Perform validation phase on packet header and a rule index
	 * @param packet A pointer to packet headers
	 * @param rule_idx The rule index to check
	 * @returns The output tuple of <priority, action>
	 */
	classifier_output_t do_validation(const uint32_t* packet, uint32_t rule_idx);

	/**
	 * @brief Rearranges the database of this to hold only a subset of the original indices
	 * @param indices A vector of field indices to keep
	 * @throws In case the field by which the iSet was created is not in the indices
	 */
	void rearrange_field_indices(const std::vector<uint32_t>& indices);

	/**
	 * @brief Returns the error list of the RQRMI model of this
	 */
	std::vector<uint32_t> get_error_list() const;

	/**
	 * @brief Returns the number of rules this holds
	 */
	uint32_t size() const { return _size; }

	/**
	 * @brief Returns the size of this (in KB)
	 */
	uint32_t get_size() const { return _size_kb; }

	/**
	 * @brief Returns the iSet index of this
	 */
	uint32_t get_iset_index() const { return _iset_index; }

	/**
	 * @brief Returns the field index of this
	 */
	uint32_t get_field_index() const { return _field_index; }

	/**
	 * @brief Returns a pointer to the database index of the iSet
	 */
	scalar_t get_index(uint32_t pos) const {
		return _index[pos];
	}

	/**
	 * @brief Returns the number of fields in the iSet
	 */
	uint32_t get_num_of_fields() const { return _num_of_columns / 2; }

	/**
	 * @brief Returns the number of columns in validation database
	 */
	uint32_t get_num_of_columns() const { return _num_of_columns; }

	/**
	 * @brief Returns the number of validation phases on this
	 */
	uint32_t get_num_of_validation_phases() const { return _num_of_validation_phases; }

	/**
	 * @brief Returns the expected error of this
	 */
	double get_expected_error() const { return _error_counter/ _search_counter; }

	/**
	 * @brief Returns a string representation of this
	 */
	virtual std::string to_string() const {
		std::stringstream ss;
		ss << "iSet-" << this->get_iset_index();
		return ss.str();
	}

	/**
	 * @brief Returns the dynamic type of this
	 */
	virtual typename NuevoMatchSubset<N>::dynamic_type_t get_type() const {
		return NuevoMatchSubset<N>::dynamic_type_t::ISET;
	}
};


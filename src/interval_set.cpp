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

#include <stdexcept>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <set>

#include <algorithms.h>
#include <rqrmi_model.h>
#include <rqrmi_fast.h>
#include <interval_set.h>

using namespace std;

/**
 * @brief A static register that holds 0xff in all of its bytes
 */
#ifndef NO_RQRMI_OPT
EPU_REG ff_vector;
#endif

/**
 * @brief Initialize new interval-set of specific index
 * @param index The iSet index within classifier
 * @param queue_size The lookup queue size
 */
template <uint32_t N>
IntervalSet<N>::IntervalSet(uint32_t index) :
		_index(nullptr), _validation_db(nullptr),
		_model(nullptr), _model_fast(nullptr),
		_iset_index(index), _size(0), _field_index(0),
		_num_of_columns(0), _num_of_validation_phases(1),
		_size_kb(0), _search_counter(0), _error_counter(0),
		_F(0), _vec_ops(0), _remainder_size(0), _rule_size(0),
		_validation_low_size(0)
{ }

template <uint32_t N>
IntervalSet<N>::~IntervalSet() {
	delete[] _index;
	delete[] _validation_db;
	delete _model_fast;
	rqrmi_free_model(_model);
}

/**
 * @brief Returns the error list of the RQRMI model of this
 */
template <uint32_t N>
std::vector<uint32_t> IntervalSet<N>::get_error_list() const {
	std::vector<uint32_t> output;
	const uint32_t *list;
	uint32_t size;
	rqrmi_get_error_list(_model, &list, &size);
	for (uint32_t i=0; i<size; ++i) {
		output.push_back(list[i]);
	}
	return output;
}

/**
 * @brief Load an interval set from file
 * @param object An object reader with binary information
 * @throws In case the index is not ordered or in case of internal error
 */
template <uint32_t N>
void IntervalSet<N>::load(ObjectReader& object) {

	// Get handlers for model and the lookup database
	ObjectReader model_handler = object.extract();
	ObjectReader index_db_handler = object.extract();
	ObjectReader validation_db_handler = object.extract();

	// Set size as the model size
	this->_size_kb = model_handler.size();

	// Load Objects
	this->_model = rqrmi_load_model(model_handler.buffer(), model_handler.size());
	matrix_t* database = load_matrix(index_db_handler.buffer(), index_db_handler.size());

	// Load RQRMIFast
	_model_fast = new RQRMIFast(this->_model);

	// Initiate the rule database
	this->_size=database->rows;
	this->_index = new scalar_t[this->_size];

	// Copy data. Note: database first column is always the one to index
	for (uint32_t i=0; i<this->_size; ++i) {
		this->_index[i] = GET_SCALAR(database, i ,0);

		// Check for order
		if (i>0 && this->_index[i] <= this->_index[i-1]) {
			free_matrix(database);
			throw error("index column in database is not ordered");
		}
	}
	free_matrix(database);

	// Read the validation database
	ObjectReader validation_reader(validation_db_handler.buffer(), validation_db_handler.size());

	uint32_t num_of_columns = validation_reader.read<uint32_t>();
	uint32_t num_of_validation_phases = 1;
	uint32_t version = 0;

	// In iSet versions later than 0 the first 4 bytes are always zero
	if (num_of_columns == 0) {
		version = validation_reader.read<uint32_t>();
	}

	// Act according to validation-db version
	switch (version) {
		case 0:
			// The number of packed column is 2*F+1
			// Remove the last column as it is the rule priority
			num_of_columns -= 1;
			break;
		case 1:
			num_of_validation_phases = validation_reader.read<uint32_t>();
			num_of_columns = validation_reader.read<uint32_t>();
			break;
		default:
			throw error("iSet validation database version not supported");
	}

	// Read validation-database header
	this->_num_of_columns = num_of_columns;
	this->_num_of_validation_phases = num_of_validation_phases;
	this->_field_index = validation_reader.read<uint32_t>();

	// Allocate the validation-database
	// Each rule has columns*num-phases values, +1 for rule priority
	uint32_t size_of_rule = this->_num_of_columns * this->_num_of_validation_phases + 1;
	uint32_t total_size = this->_size * size_of_rule;
	this->_validation_db = new uint32_t[total_size];

	loggerf("iSet size is %u, with %u columns, %u validation phases, and field index of %u. Total size: %u bytes",
			this->_size, this->_num_of_columns, this->_num_of_validation_phases,
			this->_field_index, total_size*sizeof(uint32_t));

	// Validate that the reader holds enough data
	if (validation_reader.size() != total_size*sizeof(uint32_t)) {
		throw errorf("Validation reader size (%u) does not contain enough data (%u)",
			validation_reader.size(), total_size*sizeof(uint32_t));
	}

	// Load the validation values
	for (uint32_t i=0; i<this->_size; ++i) {
		// Update rule validation values
		for (uint32_t p=0; p<this->_num_of_validation_phases; ++p) {
			for (uint32_t f=0; f<this->_num_of_columns; ++f) {
				uint32_t idx = i*size_of_rule+f*this->_num_of_validation_phases+p;
				this->_validation_db[idx] = validation_reader.read<uint32_t>();
			}
		}
		// Read rule priority
		uint32_t idx = (i+1)*size_of_rule-1;
		this->_validation_db[idx] = validation_reader.read<uint32_t>();
	}

	update_vlidation_phase_params();
}

/**
 * @brief Extract the rules of this as a vector of OpenFlow rules
 * @note Due to validation expansion, the number of output rules may be larger than expected
 */
template <uint32_t N>
std::vector<openflow_rule> IntervalSet<N>::extract_rules() const {

	list<openflow_rule> output_rules;
	uint32_t F = _num_of_columns / 2;

	uint32_t size_of_rule = this->_num_of_columns * this->_num_of_validation_phases + 1;

	// Go over all database
	for (uint32_t r=0; r<_size; ++r) {

		// Create a set of all valid values per field
		vector<set<range>> field_set(F);

		// Populate all sets with possible values from all phases
		for (uint32_t f=0; f<F; ++f) {
			for (uint32_t p=0; p<this->_num_of_validation_phases; ++p) {
				uint32_t lo_idx = r*size_of_rule+f*this->_num_of_validation_phases+p;
				uint32_t hi_idx = r*size_of_rule+(f+F)*this->_num_of_validation_phases+p;
				range current_range(
						this->_validation_db[lo_idx],
						this->_validation_db[hi_idx]);
				// Add range to set only if it's valid
				if (current_range.is_valid()) {
					field_set.at(f).insert(current_range);
				}
			}
		}

		// Read rule priority
		uint32_t idx = (r+1)*size_of_rule-1;
		uint32_t rule_prio = this->_validation_db[idx];

		// Get all possible combinations for all sets
		vector<vector<range>> combinations = calculate_all_combinations(field_set);

		// Update rules
		for (auto item : combinations) {
			openflow_rule new_rule;
			new_rule.fields.resize(F);
			for (uint32_t f=0; f<F; ++f) {
				new_rule.fields[f] = item[f];
			}
			new_rule.priority = rule_prio;
			output_rules.push_back(new_rule);
		}
	}
	return std::vector<openflow_rule>(output_rules.begin(), output_rules.end());
}

/**
 * @brief Rearranges the database of this to hold only a subset of the original indices
 * @param indices A vector of field indices to keep
 * @throws In case the field by which the iSet was created is not in the indices
 */
template <uint32_t N>
void IntervalSet<N>::rearrange_field_indices(const std::vector<uint32_t>& indices){

	// Find the new field index of this (as the fields may be reordered)
	uint32_t new_field_index = 0xffffffff;
	for (uint32_t i=0; i<indices.size(); ++i) {
		if (indices[i] == _field_index) {
			new_field_index = i;
			break;
		}
	}

	// Check that the field by which this was created is available
	if (new_field_index == 0xffffffff) {
		throw error("Cannot rearrange iSet by custom field subset: "
								 "The field by which the iSet was created is not available");
	}

	// Sizes of the current validation database
	uint32_t size_of_rule = this->_num_of_columns * this->_num_of_validation_phases + 1;
	uint32_t F = this->_num_of_columns / 2;

	// Sizes for the new validation database
	uint32_t new_num_of_columns = indices.size()*2;
	uint32_t new_F = (new_num_of_columns/2);
	uint32_t new_num_of_validation_phases = 1;

	// Holds sets of all values per field, per rule
	vector<vector<set<uint32_t>>> distinct_values_per_field_per_rule(_size);

	// Calculate the number of new validation phases
	for (uint32_t r=0; r<_size; ++r) {
		// Hold a set of all distinct values per field
		distinct_values_per_field_per_rule[r].resize(new_num_of_columns);
		// Go over all phases, populate the sets according to values of fields
		for (uint32_t p=0; p<this->_num_of_validation_phases; ++p) {
			for (uint32_t i=0; i<indices.size(); ++i) {
				// Get the range for the current field
				uint32_t lo_idx = r*size_of_rule+indices[i]*this->_num_of_validation_phases+p;
				uint32_t hi_idx = r*size_of_rule+(indices[i]+F)*this->_num_of_validation_phases+p;
				range current_range(
						this->_validation_db[lo_idx],
						this->_validation_db[hi_idx]);
				// Insert the range only if it's valid
				if (current_range.is_valid()) {
					distinct_values_per_field_per_rule[r][i].insert(current_range.low);
					distinct_values_per_field_per_rule[r][i+new_F].insert(current_range.high);
				}
			}
		}
		// Update the size of the new number of validation phases
		for (uint32_t c=0; c<new_F; ++c) {
			new_num_of_validation_phases = std::max(
					new_num_of_validation_phases,
					(uint32_t)distinct_values_per_field_per_rule[r][c].size());
		}
	}

	infof("Rearranging field indices in iSet %u. New number of columns: %u, new F: %u, new validation height: %u",
			this->_iset_index, new_num_of_columns, new_F, new_num_of_validation_phases);


	// Calculate the size of the new validation database, allocate it
	uint32_t new_size_of_rule = new_num_of_validation_phases * new_num_of_columns + 1;
	uint32_t new_size = _size * new_size_of_rule;
	uint32_t* new_validation_db = new uint32_t[new_size];

	// Populate new validation database
	for (uint32_t r=0; r<_size; ++r) {
		for (uint32_t c=0; c<new_F; ++c) {
			// Get the set iterators for the current column (low + high)
			auto low_set_iterator = distinct_values_per_field_per_rule[r][c].cbegin();
			auto low_iterator_end = distinct_values_per_field_per_rule[r][c].end();
			auto high_set_iterator = distinct_values_per_field_per_rule[r][c+new_F].cbegin();
			// Populate current column across all phases
			for (uint32_t p=0; p<new_num_of_validation_phases; ++p) {
				uint32_t lo_idx = r*new_size_of_rule+c*this->_num_of_validation_phases+p;
				uint32_t hi_idx = r*new_size_of_rule+(c+new_F)*this->_num_of_validation_phases+p;
				// In case there is no available value for column, set the column's range
				// of the current validation phase as invalid (0xffffffff to 0x00000000)
				if (low_set_iterator == low_iterator_end) {
					new_validation_db[lo_idx] = 0xffffffff;
					new_validation_db[hi_idx] = 0;
				}
				// Otherwise, get the next available value
				else {
					new_validation_db[lo_idx] = *(low_set_iterator++);
					new_validation_db[hi_idx] = *(high_set_iterator++);
				}
			}
		}
	}

	// Populate priority values for all rules
	for (uint32_t r=0; r<_size; ++r) {
		new_validation_db[(r+1)*new_size_of_rule-1] = this->_validation_db[(r+1)*size_of_rule-1];
	}

	// Delete old database, update this
	delete[] this->_validation_db;
	this->_validation_db = new_validation_db;
	this->_num_of_columns = new_num_of_columns;
	this->_num_of_validation_phases = new_num_of_validation_phases;
	this->_field_index = new_field_index;
	this->_F = new_F;

	update_vlidation_phase_params();
}

/**
 * @brief Updates all parameters required for the validation phase
 */
template <uint32_t N>
void IntervalSet<N>::update_vlidation_phase_params() {
	// Update remainder mask
	// Note: the validation phase is calculated column by column
	// the remainder is the extra rows that do not fit any vector
	_remainder_size = _num_of_validation_phases % scalars_in_vector;
	_rule_size = this->_num_of_validation_phases * this->_num_of_columns + 1;
	for (uint32_t i=0; i<scalars_in_vector; ++i) {
		_remainder_mask[i] = (_remainder_size > i) ? 0xffffffff : 0;
	}

#ifndef NDEBUG
	__m256i remainder_mask = _mm256_loadu_si256((__m256i const*)_remainder_mask);
	info("Remainder mask for iSet " << _iset_index << ": " << simd_epu_vector_logger(remainder_mask));
#endif

	// Initialize helper variables/registers for validation phase
#ifndef NO_RQRMI_OPT
	ff_vector = _mm256_set1_epi32(0xffffffff);
#endif

	// Calculate how many vector operations are required for validation
	_vec_ops = this->_num_of_validation_phases / scalars_in_vector;
	_F = _num_of_columns / 2;

	// Calculate the validation table total low bound elements
	_validation_low_size = _F * _num_of_validation_phases;
}

/**
 * @brief Search for packet within the interval set
 * @param packets A batch of packets
 * @returns A batch of tuples, each includes the matching rule's priority and action
 */
template <uint32_t N>
IntervalSetInfoBatch<N> IntervalSet<N>::rqrmi_search(PacketBatch<N>& packets) const {

	WorkBatch<iset_info_t, N> rqrmi_info;

	// RQRMIFast input/output parameters
	wide_scalar_t rqrmi_input, rqrmi_outputs, rqrmi_status, rqrmi_error;

	// Perform fast evaluation using vector of inputs

	// Main batch
	int i=0;
	for (i=0; i<(int)N-(int)RQRMIFast::input_width(); i+=(int)RQRMIFast::input_width()) {
		// Initiate SIMD inputs
		for (uint32_t k=0; k<RQRMIFast::input_width(); ++k) {
			if (packets[i+k] == nullptr) {
				rqrmi_input.scalars[k] = 0;
			} else {
				rqrmi_input.scalars[k] = packets[i+k][this->_field_index];
			}
		}

		// Perform SIMD inference
		this->_model_fast->evaluate(rqrmi_input, rqrmi_status, rqrmi_outputs, rqrmi_error);

		// Update RQRMI info
		for (uint32_t k=0; k<RQRMIFast::input_width(); ++k) {
			rqrmi_info[i+k].rqrmi_input = rqrmi_input.scalars[k];
			rqrmi_info[i+k].rqrmi_output = rqrmi_outputs.scalars[k];
			rqrmi_info[i+k].rqrmi_error = rqrmi_error.integers[k];
			rqrmi_info[i+k].valid = rqrmi_status.integers[k] & (packets[i+k] != nullptr);
			rqrmi_info[i+k].header = packets[i+k];
		}
	}

	// Remainder batch
	uint32_t counter = 0;
	int start = i;
	for (; i<(int)N; ++i) {
		if (packets[i] == nullptr) {
			rqrmi_input.scalars[counter++] = 0;
		} else {
			rqrmi_input.scalars[counter++] = packets[i][this->_field_index];
		}
	}

	// Perform SIMD inference
	this->_model_fast->evaluate(rqrmi_input, rqrmi_status, rqrmi_outputs, rqrmi_error);

	// Update RQRMI info
	for (uint32_t k=0; k<counter; ++k) {
		rqrmi_info[start+k].rqrmi_input = rqrmi_input.scalars[k];
		rqrmi_info[start+k].rqrmi_output = rqrmi_outputs.scalars[k];
		rqrmi_info[start+k].rqrmi_error = rqrmi_error.integers[k];
		rqrmi_info[start+k].valid = rqrmi_status.integers[k] & (packets[start+k] != nullptr);
		rqrmi_info[start+k].header = packets[start+k];
	}

	// Show debug messages
#ifndef NDEBUG
	infof("IntervalSet %u information for batch:", this->_field_index);
	for (uint32_t i=0; i<N; ++i) {
		infof("%u: input: %f, output: %.12f, error: %u, valid: %u, db_idx: %u",
				i, rqrmi_info[i].rqrmi_input, rqrmi_info[i].rqrmi_output,
				rqrmi_info[i].rqrmi_error, rqrmi_info[i].valid,
				(uint32_t)(rqrmi_info[i].rqrmi_output * this->_size));
	}
#endif

	return rqrmi_info;
}

/**
 * @brief Perform validation phase on packet header and a rule index
 * @param packet A pointer to packet headers
 * @param rule_idx The rule index to check
 * @returns The output tuple of <priority, action>
 */
template <uint32_t N>
classifier_output_t IntervalSet<N>::do_validation(const uint32_t* packet, uint32_t rule_idx) {

	uint32_t* cursor_lo = &this->_validation_db[rule_idx*_rule_size];
	uint32_t* cursor_hi = cursor_lo + _validation_low_size;

#if NO_RQRMI_OPT
	// For each column
	for (uint32_t f=0; f<_F; ++f) {
    if (packet[f] < *cursor_lo || packet[f] > *cursor_hi) return {-1, -1};
    cursor_lo++;
    cursor_hi++;
  }
#else
#ifndef __AVX2__
#error "NuevoMatch validation phase currently supports AVX2 extension only. Compile with -mavx2 with supported machines"
#endif
	uint32_t accumulated_result = 0;

	// TODO - write again for SSE and AVX512
	__m256i remainder_mask = _mm256_loadu_si256((__m256i const*)_remainder_mask);

	// For each column
	for (uint32_t f=0; f<_F; ++f) {
		// Set result to be false
		uint32_t column_result = 0xffffffff;
		// The current header
		__m256i header = _mm256_set1_epi32(packet[f]);

		// Divide phase to vector batches
		for (uint32_t i=0; i<_vec_ops; ++i) {
			__m256i vector_lo = _mm256_loadu_si256((const __m256i*)cursor_lo); // TODO - may cause segfault?
			__m256i vector_hi = _mm256_loadu_si256((const __m256i*)cursor_hi); // TODO - may cause segfault?
			// Which one is higher?
			__m256i result_lo = _mm256_max_epu32(header, vector_lo);
			__m256i result_hi = _mm256_max_epu32(header, vector_hi);
			// Act according to the column index:
			// Low bound we wish header >= vector_lo
			// High bound we wish vector_hi >= header
			result_lo = _mm256_cmpeq_epi32(result_lo, header);
			result_hi = _mm256_cmpeq_epi32(result_hi, vector_hi);
			// Result now hold scalars which are 0x00000000 (not good) or 0xffffffff (good)
			// We wish results that are valid in both low and high:
			// Element-wise min will return 0xffffffff
			__m256i result =  _mm256_min_epu32(result_lo, result_hi);
			// As we check across all validation phases, we wish an OR operation
			// So we want that at least one is not zero -> testz should return 0
			column_result &= _mm256_testz_si256(result, ff_vector);
			// Update cursor
			cursor_lo += scalars_in_vector;
			cursor_hi += scalars_in_vector;
		}
		// Perform calculation on the remaining validation-phases
		// Using the exact algorithm as above
		__m256i vector_lo = _mm256_loadu_si256((const __m256i*)cursor_lo); // TODO - may cause segfault?
		__m256i vector_hi = _mm256_loadu_si256((const __m256i*)cursor_hi); // TODO - may cause segfault?
		__m256i result_lo = _mm256_max_epu32(header, vector_lo);
		__m256i result_hi = _mm256_max_epu32(header, vector_hi);
		result_lo = _mm256_cmpeq_epi32(result_lo, header);
		result_hi = _mm256_cmpeq_epi32(result_hi, vector_hi);
		__m256i result =  _mm256_min_epu32(result_lo, result_hi);
		// Zero out all result elements what are not in mask
		result = _mm256_min_epu32(result, remainder_mask);
		column_result &= _mm256_testz_si256(result, ff_vector);
		// Update cursor
		cursor_lo += _remainder_size;
		cursor_hi += _remainder_size;
		accumulated_result |= column_result;
	}

	// In case the current column does not pass validation, stop
	if (accumulated_result) return {-1, -1};
#endif


	// At this point, all columns are valid
	int priority = *cursor_hi;
	return {priority, priority};
}



// Initiate templates
template class IntervalSet<1>;
template class IntervalSet<2>;
template class IntervalSet<4>;
template class IntervalSet<8>;
template class IntervalSet<16>;
template class IntervalSet<32>;
template class IntervalSet<64>;
template class IntervalSet<128>;
template class IntervalSet<256>;
template class IntervalSet<512>;

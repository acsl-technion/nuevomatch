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

#include <math.h>
#include <stdexcept>

#include <object_io.h>
#include <matrix_operations.h>
#include <lookup.h>
#include <cpu_core_tools.h>
#include <logging.h>



template <uint32_t N>
Lookup<N>::Lookup() : _index(nullptr), _values(nullptr),
	_size(0), _model(nullptr), _fast_model(nullptr), _own_model(false), _queue_size(0) {};

template <uint32_t N>
Lookup<N>::~Lookup() {
	delete[] _index;
	delete[] _values;
	delete _fast_model;
	if (_own_model) {
		rqrmi_free_model(_model);
	}
}

/**
 * @brief Loads the RQRMI model and database from objects
 * @param model The RQRMI object
 * @param database A database matrix (two columns)
 * @throws invalid_argument, domain_error
 */
template <uint32_t N>
void Lookup<N>::load(rqrmi_model_t* model, matrix_t* database) {
	// Set model
	_model = model;
	_fast_model = new RQRMIFast(_model);

	// Allocate members
	_size=database->rows;

	_index = new scalar_t[_size];
	_values = new scalar_t[_size];

	// Copy data. Note: database first column is always the one to index
	for (uint32_t i=0; i<_size; ++i) {
		_index[i] = GET_SCALAR(database, i ,0);
		_values[i] = GET_SCALAR(database, i ,1);

		// Check for order
		if (i>0 && _index[i] <= _index[i-1]) {
			free_matrix(database);
			throw domain_error("index column in database is not ordered");
		}
	}
}

/**
 * @brief Loads the RQRMI model and database from memory
 * @param object An object reader that contains the model
 * @throws invalid_argument, domain_error
 */
template <uint32_t N>
void Lookup<N>::load(ObjectReader object) {

	// Read both database and model
	ObjectReader model_handler = object.extract();
	ObjectReader database_handler = object.extract();

	// Load Objects
	_model = rqrmi_load_model(model_handler.buffer(), model_handler.size());
	matrix_t* database = load_matrix(database_handler.buffer(), database_handler.size());

	_fast_model = new RQRMIFast(_model);

	_fast_model = new RQRMIFast(_model);

	// Validate inputs
	if (database == MATRIX_ERROR) {
		throw invalid_argument("data matrix invalid");
	}
	if (_model == RQRMI_MODEL_ERROR){
		free_matrix(database);
		throw invalid_argument("invalid RQRMI model");
	}

	// Call load
	load(_model, database);

	// This owns the model
	_own_model = true;

	// Free resources
	free_matrix(database);
}

/**
 * @brief Adds listener object
 * @param listener A LookupListener object
 */
template <uint32_t N>
void Lookup<N>::add_listener(LookupListener& listener) {
	_listeners.push_back(&listener);
}

/**
 * @brief Perform binary search over the database
 * @param key The key to search
 * @param hint A hint to the possible location (RQRMI output)
 * @param error Maximum error bound for search
 * @throws out_of_range In case of lookup internal error
 */
template <uint32_t N>
void Lookup<N>::perform_binary_search(scalar_t key, scalar_t hint, uint32_t error) {

// Used for debugging
#ifdef NO_BIN_SEARCH
	for (auto it : _listeners) {
		it->on_new_result(key, 0, 1);
	}
	return;
#endif

// Used for debugging
#ifdef CUSTOM_ERROR
	error = CUSTOM_ERROR_VALUE;
#endif

	uint32_t pos = hint * _size;
	uint32_t u_bound = MIN((int)_size-1, (int)pos+(int)error);
	uint32_t l_bound = MAX(0, (int)pos-(int)error);

#ifndef LINEAR_SEARCH

	info("Performing binary search");

	// Perform binary search
	do {

		uint8_t current_value = _index[pos] <= key;
		uint8_t next_value = _index[pos+1] > key;

		if (current_value & next_value) {
			break;
		} else if (current_value) {
			l_bound = pos;
			pos=(l_bound+u_bound);
			pos=(pos>>1)+(pos&0x1); // Ceil
		} else {
			u_bound = pos;
			pos=(l_bound+u_bound)>>1; // Floor
		}

		error >>= 1;
	} while (error > 0);
#else
	info("Performing linear search");

	// Perform linear search
	for (pos=l_bound; pos<u_bound; ++pos) {
		uint8_t current_value = _index[pos] <= key;
		uint8_t next_value = _index[pos+1] > key;
		if (current_value & next_value) {
			break;
		}
	}
#endif

	info("Performing validation");

	// Perform validation
	int found = _values[pos] < key;
	for (auto it : _listeners) {
		it->on_new_result(key, pos, found);
	}
	return;
}

/**
 * @brief Returns a record from the database
 * @param index The requested index of the record
 * @returns A pair of scalars {interval_start_inclusive, interval_end_exclusive}
 * @throws overflow_error
 */
template <uint32_t N>
scalar_pair_t Lookup<N>::get_record(uint32_t index) {
	if (index >= _size) {
		throw overflow_error("record index our of range");
	}
	return (scalar_pair_t) { _index[index], _values[index] };
}

template <uint32_t N>
LookupCPU<N>::LookupCPU() : _worker_rqrmi(nullptr), _worker_lookup(nullptr) {}

template <uint32_t N>
LookupCPU<N>::~LookupCPU() {
	delete _worker_rqrmi;
	delete _worker_lookup;
}

/*
* @brief Loads the RQRMI model and database from objects
* @param model The RQRMI object
* @param database A database matrix (two columns)
* @throws invalid_argument, domain_error
*/
template <uint32_t N>
void LookupCPU<N>::load(rqrmi_model_t* model, matrix_t* database) {
	// Call super
	Lookup<N>::load(model, database);

	// Get current core index
	int current_core_idx = cpu_core_tools_get_index_of_current_thread();
	int rqrmi_core = cpu_core_tools_get_next_physical_core(current_core_idx);
	int lookup_core = cpu_core_tools_get_next_physical_core(rqrmi_core);

	// Check for performance hazards
	if (lookup_core == current_core_idx || rqrmi_core == lookup_core || current_core_idx == rqrmi_core) {
		warning("cannot access three different physical cores in system. performance may degenerate");
	}

	// Initiate worker threads
	_worker_rqrmi = new PipelineThread<batch_t>(this->_queue_size, rqrmi_core, consumer_rqrmi_thread, this);
	_worker_lookup = new PipelineThread<rqrmi_output_batch_t>(this->_queue_size, lookup_core, consumer_lookup_thread, this);
}

/**
 * @brief Perform lookup on the requested value
 * @param input The input to search
 * @returns true iff the search was performed
 */
template <uint32_t N>
bool LookupCPU<N>::search(batch_t input) {
	return _worker_rqrmi->produce(input);
}

/**
 * @brief Is invoked by the RQRMI thread when new job is available
 * @param input The new input
 * @param lookup_instance A pointer to lookup instance
 * @returns True iff the job was consumed
 */
template <uint32_t N>
bool LookupCPU<N>::consumer_rqrmi_thread(batch_t& input, void* lookup_instance) {
	LookupCPU* instance = (LookupCPU*)(lookup_instance);

	// Check new slot is available, fail fast (important for preventing deadlocks)
	if (!instance->_worker_lookup->available()) {
		return false;
	}

	// New slot is available
	rqrmi_output_batch_t rqrmi_produced_job;

	wide_scalar_t inputs, outputs, status, error;

	constexpr int first_batch = (int)N-(int)RQRMIFast::input_width();
	info("Performing RQRMIFast lookup on first batch with " << first_batch << " elements");

	// Perform lookup in batches for SIMD word
	int i;
	for (i=0; i<first_batch; i+=(int)RQRMIFast::input_width()) {
		for (uint32_t k=0; k<RQRMIFast::input_width(); ++k) {
			inputs.scalars[k] = input[i+k];
		}
		instance->_fast_model->evaluate(inputs, status, outputs, error);
		for (uint32_t k=0; k<RQRMIFast::input_width(); ++k) {
			rqrmi_produced_job[i+k].input = inputs.scalars[k];
			rqrmi_produced_job[i+k].result = outputs.scalars[k];
			rqrmi_produced_job[i+k].error = error.scalars[k];
			rqrmi_produced_job[i+k].valid = status.scalars[k];
		}
	}

	// Handle the remainder
	info("Handling the RQRMIFast remainder");

	uint32_t counter = 0;
	int start = i;
	for (; i<(int)N; ++i) {
		inputs.scalars[counter++] = input[i];
	}
	instance->_fast_model->evaluate(inputs, status, outputs, error);
	for (uint32_t k=0; k<counter; ++k) {
		rqrmi_produced_job[start+k].input = inputs.scalars[k];
		rqrmi_produced_job[start+k].result = outputs.scalars[k];
		rqrmi_produced_job[start+k].error = error.scalars[k];
		rqrmi_produced_job[start+k].valid = status.scalars[k];
	}

	info("Done with RQRMIFast");

	// Produce new lookup item
	if (!instance->_worker_lookup->produce(rqrmi_produced_job)) {
		throw error("cannot produce item in RQRMI consumer thread - should not get here!");
	}

	return true;
}

/**
 * @brief Is invoked by the lookup thread when new RQRMI result is available
 * @param input The RQRMI result
 * @param lookup_instance A pointer to lookup instance
 * @returns True iff the job was consumed
 */
template <uint32_t N>
bool LookupCPU<N>::consumer_lookup_thread(rqrmi_output_batch_t& input, void* lookup_instance) {
	LookupCPU* instance = (LookupCPU*)(lookup_instance);

	info("Lookup worker consume enter");

	// For each item in batch
	for (uint32_t i=0; i<N; ++i) {
		// Skip invalid inputs
		if (!input[i].valid) continue;
		// Perform the binary search
		// Invoke listeners
		instance->perform_binary_search(input[i].input, input[i].result, input[i].error);
	}

	// Item is always consumed
	return true;
}

/**
 * @brief Starts the performance measurement of this
 */
template <uint32_t N>
void LookupCPU<N>::start_performance_measurement() {
	_worker_lookup->start_performance_measurements();
	_worker_rqrmi->start_performance_measurements();
}

/**
 * @brief Stops the performance measurement of this
 */
template <uint32_t N>
void LookupCPU<N>::stop_performance_measurement() {
	_worker_lookup->stop_performance_measurements();
	_worker_rqrmi->stop_performance_measurements();
}

/**
 * @brief Prints statistics of this
 */
template <uint32_t N>
void LookupCPU<N>::print() const {
	SimpleLogger::get().format("RQRMI utilization: %.2f%% "
		   "(throughput: %.2f rpus, backpressure: %.2f rpus, average time per batch: %.2f us)\n"
		   "Lookup utilization: %.2f%% "
		   "(throughput: %.2f rpus, backpressure: %.2f rpus, average time per batch: %.2f us)\n",
			_worker_rqrmi->get_utilization(),   _worker_rqrmi->get_throughput(),
			_worker_rqrmi->get_backpressure(), _worker_rqrmi->get_average_work_time(),
			_worker_lookup->get_utilization(), _worker_lookup->get_throughput(),
			_worker_lookup->get_backpressure(), _worker_lookup->get_average_work_time());
}

// Explicit template Instantiation
template class Lookup<1>;
template class Lookup<4>;
template class Lookup<8>;
template class Lookup<16>;
template class Lookup<32>;
template class Lookup<64>;
template class Lookup<128>;
template class LookupCPU<1>;
template class LookupCPU<4>;
template class LookupCPU<8>;
template class LookupCPU<16>;
template class LookupCPU<32>;
template class LookupCPU<64>;
template class LookupCPU<128>;

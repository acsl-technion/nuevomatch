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
#include <exception>
#include <stdexcept>
#include <basic_types.h>
#include <object_io.h>
#include <matrix_operations.h>
#include <rqrmi_model.h>
#include <pipeline_thread.h>
#include <rqrmi_fast.h>

using namespace std;

/**
 * @brief Abstract class for lookup listeners
 */
class LookupListener {
public:
	virtual ~LookupListener() {};
protected:

	/**
	 * @brief Invoked by the lookup process when a new result is available
	 * @param input The search value
	 * @param index The lookup range-value pair index
	 * @param found 1 iff the value was found
	 */
	virtual void on_new_result(scalar_t input, uint32_t index, int found) = 0;

	// Lookup classes can call on_new_result
	template <uint32_t N>
	friend class Lookup;
};

/**
 * @brief Abstract class for performing Lookups using RQRMI models
 * @tparam Batch size per lookup
 */
template <uint32_t N>
class Lookup  {
public:
	/**
	 * @brief A batch job of input scalars
	 */
	typedef WorkBatch<scalar_t, N> batch_t;

protected:

	// The lookup table of the secondary search
	// Divided to index (fast lookup) and values (validation)
	scalar_t* _index;
	scalar_t* _values;
	uint32_t _size;

	// The RQRMI model
	rqrmi_model_t* _model;
	RQRMIFast* _fast_model;
	bool _own_model;

	// The queue size
	uint32_t _queue_size;

	// Holds listeners for search results
	vector<LookupListener*> _listeners;

	/**
	 * @brief Perform binary search over the database
	 * @param key The key to search
	 * @param hint A hint to the possible location (RQRMI output)
	 * @param error Maximum error bound for search
	 * @throws out_of_range In case of lookup internal error
	 */
	void perform_binary_search(scalar_t key, scalar_t hint, uint32_t error);

public:

	/**
	 * @brief Describe a lookup result
	 */
	typedef struct {
		scalar_t input;
		uint32_t index;
		uint32_t found;
	} result_t;

	Lookup();
	virtual ~Lookup();

	/**
	 * @brief Loads the RQRMI model and database from memory
	 * @param object An object reader that contains the model
	 * @throws invalid_argument, domain_error
	 */
	virtual void load(ObjectReader object);

	/**
	 * @brief Loads the RQRMI model and database from objects
	 * @param model The RQRMI object
	 * @param database A database matrix (two columns)
	 * @throws invalid_argument, domain_error
	 */
	virtual void load(rqrmi_model_t* model, matrix_t* database);

	/**
	 * @brief Adds listener object
	 * @param listener A LookupListener object
	 */
	void add_listener(LookupListener& listener);

	/**
	 * @brief Sets the size of the queue for the lookup procedure
	 * @param queue_size The queue size
	 */
	void set_queue_size(uint32_t queue_size) { _queue_size = queue_size; }

	/**
	 * @brief Perform lookup on the requested value
	 * @param input The input to search
	 * @returns true iff the search was performed
	 */
	virtual bool search(batch_t input) = 0;

	/**
	 * @brief Get the size of thid
	 */
	uint32_t get_size() const { return _size; }

	/**
	 * @brief Returns a record from the database
	 * @param index The requested index of the record
	 * @returns A pair of scalars {interval_start_inclusive, interval_end_exclusive}
	 * @throws overflow_error
	 */
	scalar_pair_t get_record(uint32_t index);

};


/**
 * @brief Performs lookup on CPU
 * @tparam Batch size per lookup
 */
template <uint32_t N>
class LookupCPU : public Lookup<N> {
public:

	/**
	 * @brief A batch job of input scalars
	 */
	typedef typename Lookup<N>::batch_t batch_t;

private:

	// RQRMI output pair of result and max-error
	typedef struct {
		scalar_t input;
		scalar_t result;
		uint32_t error;
		bool valid;
	} rqrmi_output_t;

	// RQRMI output as batch
	typedef WorkBatch<rqrmi_output_t, N> rqrmi_output_batch_t;

	// Worker threads
	PipelineThread<batch_t> *_worker_rqrmi;
	PipelineThread<rqrmi_output_batch_t> *_worker_lookup;

	// Worker methods

	/**
	 * @brief Is invoked by the RQRMI thread when new job is available
	 * @param input The new input
	 * @param lookup_instance A pointer to lookup instance
	 * @returns True iff the job was consumed
	 */
	static bool consumer_rqrmi_thread(batch_t& input, void* lookup_instance);

	/**
	 * @brief Is invoked by the lookup thread when new RQRMI result is available
	 * @param input The RQRMI result
	 * @param lookup_instance A pointer to lookup instance
	 * @returns True iff the job was consumed
	 */
	static bool consumer_lookup_thread(rqrmi_output_batch_t& input, void* lookup_instance);

public:

	LookupCPU();
	virtual ~LookupCPU();

	/*
	* @brief Loads the RQRMI model and database from objects
	* @param model The RQRMI object
	* @param database A database matrix (two columns)
	* @throws invalid_argument, domain_error
	*/
	virtual void load(rqrmi_model_t* model, matrix_t* database);

	// All overloads of load in super are visible here
	using Lookup<N>::load;

	/**
	 * @brief Perform lookup on the requested value
	 * @param input The input batch to search
	 * @returns true iff the search was performed
	 */
	virtual bool search(batch_t input);

	/**
	 * @brief Starts the performance measurement of this
	 */
	void start_performance_measurement();

	/**
	 * @brief Stops the performance measurement of this
	 */
	void stop_performance_measurement();

	/**
	 * @brief Prints statistics of this
	 */
	void print() const;
};


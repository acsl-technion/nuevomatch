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

// Used to skip in make
// Python's setup tools should compile this file
#ifndef NOPYTHON

#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <time.h>

// Include path set by makefile according to python-dev version
#include <Python.h>

#include <logging.h>
#include <basic_types.h>
#include <object_io.h>
#include <rqrmi_model.h>
#include <rqrmi_tools.h>
#include <lookup.h>

// RQRMI Capsules names
static const char* rqrmi_model_capsule_name = "RQRMI model";
static const char* rqrmi_probe_capsule_name = "RQRMI probe";
static const char* rqrmi_matrix_capsule_name = "RQRMI matrix";

// Used for lookup objects
static const char* lookup_cpu_capsule_name = "LookupCpu";
typedef struct {
	Lookup<1>* lookup;
	LookupListener* listener;
} lookup_pair_t;

// Listener for Lookup objects
class PythonLookupListener : public LookupListener {
protected:
	// Holds Python callback function pointer
	PyObject* _fun_ptr;
	void on_new_result(scalar_t input, uint32_t index, int found) {
		// This is a callback thread thus the GIL should be acquired
		// See: https://docs.python.org/3/c-api/init.html#non-python-created-threads
		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();

		// Create the argument list
		PyObject* arguments = PyTuple_New(3);
		PyTuple_SetItem(arguments, 0, PyFloat_FromDouble(input));
		PyTuple_SetItem(arguments, 1, PyLong_FromLong(index));
		PyTuple_SetItem(arguments, 2, PyLong_FromLong(found));

		// Call the python function with arguments
		PyObject_Call(_fun_ptr, arguments, NULL);

		// Release GIL
		PyGILState_Release(gstate);
	}

public:
	PythonLookupListener(PyObject* fun_ptr) : _fun_ptr(fun_ptr) {};
};

// Used to tack memory
static int print_memory_tracking = 0;
#define NEW_CAPSULE(NAME, P) if(print_memory_tracking) fprintf(stderr, "New %s object at %p\n", NAME, P);
#define REM_CAPSULE(NAME, P) if(print_memory_tracking) fprintf(stderr, "Delete %s object at %p\n", NAME, P);

#ifdef __cplusplus
extern "C" {
#endif

/*
 * @brief Is automatically called with RQRMI model capsules to deallocate their memory
 */
static void py_rqrmi_model_capsule_destructor(PyObject *args) {
	rqrmi_model_t* rqrmi_model = (rqrmi_model_t*)PyCapsule_GetPointer(args, rqrmi_model_capsule_name);
	// In case of an exception
	if (rqrmi_model==RQRMI_MODEL_ERROR) {
		return;
	}
	REM_CAPSULE(rqrmi_model_capsule_name, rqrmi_model);
	rqrmi_free_model(rqrmi_model);
}

/*
 * @brief Is automatically called with RQRMI probe capsules to deallocate their memory
 */
static void py_rqrmi_probe_capsule_destructor(PyObject *args) {
	rqrmi_probing_t* rqrmi_probe = (rqrmi_probing_t*)PyCapsule_GetPointer(args, rqrmi_probe_capsule_name);
	// In case of an exception
	if (rqrmi_probe==RQRMI_PROBING_ERROR) {
		return;
	}
	REM_CAPSULE(rqrmi_probe_capsule_name, rqrmi_probe);
	rqrmi_tools_probe_free(rqrmi_probe);
}

/*
 * @brief Is automatically called with RQRMI matrix capsules to deallocate their memory
 */
static void py_rqrmi_matrix_capsule_destructor(PyObject *args) {
	matrix_t* rqrmi_matrix = (matrix_t*)PyCapsule_GetPointer(args, rqrmi_matrix_capsule_name);
	// In case of an exception
	if (rqrmi_matrix==MATRIX_ERROR) {
		return;
	}
	REM_CAPSULE(rqrmi_matrix_capsule_name, rqrmi_matrix);
	free_matrix(rqrmi_matrix);
}

/*
 * @brief Is automatically called with LookupCpu capsules to deallocate their memory
 */
static void py_lookup_destructor(PyObject *args) {
	lookup_pair_t* object = (lookup_pair_t*)PyCapsule_GetPointer(args, lookup_cpu_capsule_name);
	// In case of an exception
	if (object==NULL) {
		return;
	}
	REM_CAPSULE(lookup_cpu_capsule_name, object);
	delete object->lookup;
	delete object->listener;
	delete object;
}

/**
 * @brief Creates an RQRMI lookup table object capsule for Python
 * @param An RQRMI matrix capsule
 * @param An RQRMI model capsule
 * @param Integer, queue size
 * @param Function pointer, listener to new results
 * @returns An RQRMI lookup table capsule
 * @throws RuntimeError In case the input is invalid or memory allocation error
 */
static PyObject* py_lookup_cpu_create(PyObject *self, PyObject *args) {

	// Parse arguments
	PyObject *rqrmi_matrix_capsule, *rqrmi_capsule, *listener;
	int queue_size;
	if (!PyArg_ParseTuple(args, "OOiO:lookup_cpu_create",
			&rqrmi_matrix_capsule, &rqrmi_capsule, &queue_size, &listener)) {
		return NULL;
	}

	// Get objects
	matrix_t* data = (matrix_t*)PyCapsule_GetPointer(rqrmi_matrix_capsule, rqrmi_matrix_capsule_name);
	rqrmi_model_t* rqrmi_model = (rqrmi_model_t*)PyCapsule_GetPointer(rqrmi_capsule, rqrmi_model_capsule_name);
	if (!data || !rqrmi_model) return NULL;

	// Validate last argument
	if (!PyCallable_Check(listener)) {
		PyErr_SetString(PyExc_RuntimeError, "listener function pointer is invalid");
		return NULL;
	}

	// Validate that the Python thread system is up and working
	// for the multithreaded environment (prior to python 3.7)
	// See: https://docs.python.org/3/c-api/init.html#c.PyEval_InitThreads
	Py_Initialize();
	PyEval_InitThreads();

	lookup_pair_t* output;
	try {
		// Create the new object
		output = new lookup_pair_t;
		output->lookup = new LookupCPU<1>();
		output->listener = new PythonLookupListener(listener);

		// Load and add listener
		output->lookup->set_queue_size(queue_size);
		output->lookup->load(rqrmi_model, data);
		output->lookup->add_listener(*output->listener);
	} catch (const char* ex) {
		PyErr_SetString(PyExc_RuntimeError, ex);
		return NULL;
	}

	// Wrap the pointer as Python capsule
	// See: https://docs.python.org/3.6/c-api/capsule.html
	PyObject *capsule = PyCapsule_New(output, lookup_cpu_capsule_name, py_lookup_destructor);
	if (capsule == NULL) {
		PyErr_SetString(PyExc_RuntimeError, "error creating capsule");
		return NULL;
	}

	NEW_CAPSULE(lookup_cpu_capsule_name, output);
	return capsule;
}

/**
 * @brief Adapter to Lookup search method
 * @param Lookup capsule
 * @param Float, what to search
 * @throws RuntimeError in case of invalid arguments or internal error
 */
static PyObject* py_lookup_search(PyObject *self, PyObject *args) {

	// Parse arguments
	PyObject *lookup_capsule;
	scalar_t what;
	if (!PyArg_ParseTuple(args, "Of:plookup_cpu_search", &lookup_capsule, &what)) {
		return NULL;
	}

	lookup_pair_t* object = (lookup_pair_t*)PyCapsule_GetPointer(lookup_capsule, lookup_cpu_capsule_name);
	if (object == NULL) {
		PyErr_SetString(PyExc_RuntimeError, "Invalid arguments");
		return NULL;
	}

	// Release GIL for listener thread will be able to work
	PyThreadState* state = PyEval_SaveThread();

	// Do until search is done
	while(!object->lookup->search( (Lookup<1>::batch_t){what} ));

	// Restore python thread
	PyEval_RestoreThread(state);

	// Return None
	Py_RETURN_NONE;
}

/**
 * @brief Creates an RQRMI matrix object capsule for Python
 * @param Packed Matrix <MxN> (scalar_t)
 * @returns An RQRMI matrix capsule
 * @throws RuntimeError In case the input is invalid or memory allocation error
 */
static PyObject* py_rqrmi_matrix_create(PyObject *self, PyObject *args) {

	// Parse the first argument as Python object
	PyObject *arg1;
	if (!PyArg_ParseTuple(args, "O:create_matrix", &arg1)) {
		return NULL;
	}

	// Request to read from the object as buffer
	Py_buffer buffer;
	if (PyObject_GetBuffer(arg1, &buffer, PyBUF_SIMPLE) < 0) {
		return NULL;
	}

	matrix_t* output = load_matrix(buffer.buf, buffer.len);

	// Allocate matrix
	if (output == MATRIX_ERROR) {
		PyErr_SetString(PyExc_RuntimeError, "Invalid packed matrix format");
		return NULL;
	}

	// Wrap the pointer as Python capsule
	// See: https://docs.python.org/3.6/c-api/capsule.html
	PyObject* capsule = PyCapsule_New(output, rqrmi_matrix_capsule_name, py_rqrmi_matrix_capsule_destructor);
	NEW_CAPSULE(rqrmi_matrix_capsule_name, output);
	return capsule;
}

/**
 * @brief Python adapter to call load_rqrmi_model method
 * @param Packed RQRMI Model
 * @returns RQRMI Model Capsule
 */
static PyObject* py_load_model(PyObject *self, PyObject *args) {

	// Example taken from here:
	// https://stackoverflow.com/a/15461457/4103200

	// Parse the first argument as Python object
	PyObject *arg1;
	if (!PyArg_ParseTuple(args, "O:load", &arg1)) {
		return NULL;
	}

	// Request to read from the object as buffer
	Py_buffer buffer;
	if (PyObject_GetBuffer(arg1, &buffer, PyBUF_SIMPLE) < 0) {
		return NULL;
	}

	// See Py_buffer object structure here:
	// https://docs.python.org/3/c-api/buffer.html#c.Py_buffer
	// Py_buffer.buf: For contiguous arrays points to the beginning of the memory block
	// Py_buffer.len: For contiguous arrays, this is the length of the underlying memory block.
	// Load the model from filename
	rqrmi_model_t* rqrmi_model = rqrmi_load_model(buffer.buf, buffer.len);
	PyObject* capsule;

	// Must release the buffer when finished, and see
	// https://docs.python.org/3/c-api/buffer.html#c.PyObject_GetBuffer
	PyBuffer_Release(&buffer);

	// Throw exception in case of an error
	// See: https://docs.python.org/3/c-api/exceptions.html
	if (rqrmi_model==RQRMI_MODEL_ERROR) {
		PyErr_SetString(PyExc_RuntimeError, SimpleLogger::get().get_buffer());
		return NULL;
	}

	// Wrap the pointer as Python capsule
	// See: https://docs.python.org/3.6/c-api/capsule.html
	capsule = PyCapsule_New(rqrmi_model, rqrmi_model_capsule_name, py_rqrmi_model_capsule_destructor);

	// Return a pointer to the RQRMI object
	NEW_CAPSULE(rqrmi_model_capsule_name, rqrmi_model);
	return capsule;
}

/**
 * @brief Python adapter for initiate_error_handling_module method
 * @param Integer, verbose
 */
static PyObject* py_initiate(PyObject *self, PyObject *args) {

	// Parse the first argument as integer
	int arg1;
	if (!PyArg_ParseTuple(args, "i:initiate", &arg1)) {
		return NULL;
	}

	// Initiate error handling
	SimpleLogger::get().set_sticky_force(arg1);
	print_memory_tracking = arg1;

	Py_RETURN_NONE;
}

/**
 * @brief Python adapter for evaluate_rqrmi_model method
 * @param RQRMI Model Capsule
 * @param RQRMI Matrix Capsule of scalar_t
 * @returns A tuple with two elements: (1) A tuple of (total_time, avg_time), and; (2) A list of outputs
 * @throws RuntimeError in case the RQRMI model is not built
 */
static PyObject* py_evaluate_model(PyObject *self, PyObject *args) {

	// Parse two inputs as capsules
	PyObject *rqrmi_model_capsule, *rqrmi_matrix_capsule;
	if (!PyArg_ParseTuple(args, "OO:evaluate_model", &rqrmi_model_capsule, &rqrmi_matrix_capsule)) {
		return NULL;
	}

	// Extract pointers from capsules
	rqrmi_model_t* rqrmi_model = (rqrmi_model_t*)PyCapsule_GetPointer(rqrmi_model_capsule, rqrmi_model_capsule_name);
	matrix_t* data = (matrix_t*)PyCapsule_GetPointer(rqrmi_matrix_capsule, rqrmi_matrix_capsule_name);

	// Throw exception in case RQRMI is not defined
	if (rqrmi_model==RQRMI_MODEL_ERROR || data==MATRIX_ERROR) {
		PyErr_SetString(PyExc_RuntimeError, "Cannot evaluate RQRMI model: invalid input");
		return NULL;
	}

	// Create output Python float32 list
	// See: https://docs.python.org/3/c-api/list.html
	uint32_t num_of_packets = data->rows;
	PyObject* result_list = PyList_New(num_of_packets);

	// Measure performance
	struct timespec start_time, end_time;

	// Catch model evaluation errors
	try {

		// Measure inference time
		clock_gettime(CLOCK_MONOTONIC, &start_time);

		// Perform the evaluation per input value
		for (uint32_t i=0; i<num_of_packets; ++i) {
			scalar_t input = GET_SCALAR(data, i, 0);
			// Evaluate current packet
			// Note: No AVX used since we don't know what is the submodels' size (RMI, for instance)
			// The Lookup object, which works with RQRMI only, uses AVX when possible.
			double out = (double)rqrmi_evaluate_model(rqrmi_model, input);
			// Set the output element
			// See: https://docs.python.org/3/c-api/float.html
			PyList_SetItem(result_list, i, PyFloat_FromDouble(out));
		}

		// Measure inference time
		clock_gettime(CLOCK_MONOTONIC, &end_time);

	} catch (const std::exception& e) {
		warning(e.what());
		PyErr_SetString(PyExc_RuntimeError, SimpleLogger::get().get_buffer());
		return NULL;
	}

	double last_eval_tot_ns = (end_time.tv_sec * 1e9 + end_time.tv_nsec) - (start_time.tv_sec * 1e9 + start_time.tv_nsec);
	double last_eval_avg_ns = last_eval_tot_ns/num_of_packets;

	// Create the measurements tuple (in micro seconds)
	// See: https://docs.python.org/3/c-api/tuple.html
	PyObject* measurement_list = PyTuple_New(2);
	PyTuple_SetItem(measurement_list, 0, PyFloat_FromDouble(last_eval_tot_ns / 1000));
	PyTuple_SetItem(measurement_list, 1, PyFloat_FromDouble(last_eval_avg_ns/ 1000));

	// Set the output tuple
	// See: https://docs.python.org/3/c-api/tuple.html
	PyObject* output_list = PyTuple_New(2);
	PyTuple_SetItem(output_list, 0, measurement_list);
	PyTuple_SetItem(output_list, 1, result_list);

	// Return the list
	return output_list;
}

/**
 * @brief Python adapter for rqrmi_tools_generate_dataset method
 * @param An RQRMI Probe Capsule
 * @param Boolean, randomize dataset
 * @param Boolean, shuffle dataset
 * @param Integer, stage index
 * @param Integer, bucket_index
 * @param Integer, The number of samples to sample from generated dataset.
 * @returns An RQRMI Matrix Capsule (scalar_t) generated dataset
 * @throws RuntimeError
 */
static PyObject* py_create_dataset(PyObject *self, PyObject *args) {
	PyObject *probe_capsule;
	int stage_idx, bucket_idx, samples;
	bool random, shuffle;

	// Parse input according to format
	if (!PyArg_ParseTuple(args, "Oppiii:create_dataset",
			&probe_capsule, &random, &shuffle, &stage_idx, &bucket_idx, &samples)) {
		return NULL;
	}

	// Extract pointers from capsules
	rqrmi_probing_t* rqrmi_probe = (rqrmi_probing_t*)PyCapsule_GetPointer(probe_capsule, rqrmi_probe_capsule_name);

	// Call function & check for errors
	matrix_t* output = rqrmi_tools_generate_dataset(rqrmi_probe, stage_idx, bucket_idx, samples, random, shuffle);
	if (output == MATRIX_ERROR) {
		PyErr_SetString(PyExc_RuntimeError, SimpleLogger::get().get_buffer());
		return NULL;
	}

	// Wrap the result as Python capsule
	// See: https://docs.python.org/3.6/c-api/capsule.html
	PyObject* capsule = PyCapsule_New(output, rqrmi_matrix_capsule_name, py_rqrmi_matrix_capsule_destructor);
	NEW_CAPSULE(rqrmi_matrix_capsule_name, output);
	return capsule;
}

/**
 * @brief Python adapter for rqrmi_tools_probe_new method
 * @param An RQRMI Matrix Capsule, records
 * @param A list of integers, each with the corresponding stage width
 * @returns An RQRMI Probe Capsule
 * @throws RuntimeError
 */
static PyObject* py_create_probe(PyObject *self, PyObject *args) {
	PyObject *records_capsule, *stage_list;

	// Parse input according to format
	if (!PyArg_ParseTuple(args, "OO:create_probe", &records_capsule, &stage_list)) {
		return NULL;
	}

	// Extract pointers from capsules
	matrix_t* records = (matrix_t*)PyCapsule_GetPointer(records_capsule, rqrmi_matrix_capsule_name);

	// Check that second argument is a list
	if (!PyList_Check(stage_list)) {
		PyErr_SetString(PyExc_ValueError, "second argument is not a valid list");
		return NULL;
	}

	// Generate integer array
	// List Documentation: https://docs.python.org/3.6/c-api/list.html
	uint32_t num_of_stages = PyList_Size(stage_list);
	uint32_t stage_width[num_of_stages];

	// Read list
	for (uint32_t i=0; i<num_of_stages; ++i) {
		// Long Documentation: https://docs.python.org/3.6/c-api/long.html
		PyObject* item = PyList_GetItem(stage_list, i);
		if (!PyLong_Check(item)) {
			PyErr_SetString(PyExc_ValueError, "list contains non integer element");
			return NULL;
		}
		stage_width[i] = PyLong_AsLong(item);
	}

	// Call function & check for errors
	rqrmi_probing_t* output = rqrmi_tools_probe_new(records, num_of_stages, stage_width);
	if (output == RQRMI_PROBING_ERROR) {
		PyErr_SetString(PyExc_RuntimeError, SimpleLogger::get().get_buffer());
		return NULL;
	}

	// Wrap the result as Python capsule
	// See: https://docs.python.org/3.6/c-api/capsule.html
	PyObject* capsule = PyCapsule_New(output, rqrmi_probe_capsule_name, py_rqrmi_probe_capsule_destructor);
	NEW_CAPSULE(rqrmi_probe_capsule_name, output);
	return capsule;
}

/**
 * @brief Python adapter for rqrmi_tools_calculate_transition_set method
 * @param An RQRMI Model Capsule
 * @param An RQRMI Probe Capsule
 * @param Integer, the required stage index
 * @returns A list of lists, each internal list is a transition input
 * @throws RuntimeError
 */
static PyObject* py_calculate_transition_set(PyObject *self, PyObject *args) {

	PyObject *rqrmi_capsule, *probe_capsule;
	int stage_idx;
	// Parse input according to format
	if (!PyArg_ParseTuple(args, "OOi:calculate_transition_set",
			&rqrmi_capsule, &probe_capsule, &stage_idx)) {
		return NULL;
	}

	// Extract pointers from capsules
	rqrmi_model_t* rqrmi_model = (rqrmi_model_t*)PyCapsule_GetPointer(rqrmi_capsule, rqrmi_model_capsule_name);
	rqrmi_probing_t* rqrmi_probe = (rqrmi_probing_t*)PyCapsule_GetPointer(probe_capsule, rqrmi_probe_capsule_name);

	// Calculate transition set
	vector_list_t* U_i = rqrmi_tools_calculate_transition_set(rqrmi_model, rqrmi_probe, stage_idx);
	if (U_i == VECTOR_LIST_ERROR) {
		PyErr_SetString(PyExc_RuntimeError, SimpleLogger::get().get_buffer());
		return NULL;
	}

	// Set the output  list
	uint32_t w = vector_list_get_size(U_i);
	uint32_t i=0;
	PyObject* output_list = PyList_New(w);
	scalar_t* cursor = (scalar_t*)vector_list_begin(U_i);
	for (; cursor; cursor = (scalar_t*)vector_list_iterate(U_i)) {
		PyObject* current = PyTuple_New(3);
		PyTuple_SetItem(current, 0, PyFloat_FromDouble(cursor[0]));
		PyTuple_SetItem(current, 1, PyFloat_FromDouble(cursor[1]));
		PyTuple_SetItem(current, 2, PyFloat_FromDouble(cursor[2]));
		PyList_SetItem(output_list, i++, current);
	}

	return output_list;
}

/**
 * @brief Python adapter for rqrmi_tools_calculate_responsibility method
 * @param An RQRMI Model Capsule
 * @param An RQRMI Probe Capsule
 * @param Integer, The required stage
 * @returns A list of responsibilities, each is a list of intervals, each is a tuple of two scalars (x,y)
 * @throws RuntimeError
 */
static PyObject* py_calculate_responsibility(PyObject *self, PyObject *args) {

	PyObject *rqrmi_capsule, *probe_capsule;
	int stage_idx;

	// Parse input according to format
	if (!PyArg_ParseTuple(args, "OOi:calculate_responsibility",
			&rqrmi_capsule, &probe_capsule, &stage_idx))
	{
		return NULL;
	}

	// Extract pointers from capsules
	rqrmi_model_t* rqrmi_model = (rqrmi_model_t*)PyCapsule_GetPointer(rqrmi_capsule, rqrmi_model_capsule_name);
	rqrmi_probing_t* rqrmi_probe = (rqrmi_probing_t*)PyCapsule_GetPointer(probe_capsule, rqrmi_probe_capsule_name);

	// Call function & check for errors
	vector_list_t** result = rqrmi_tools_calculate_responsibility(rqrmi_model, rqrmi_probe, stage_idx);
	if (result == VECTOR_LIST_ERROR) {
		PyErr_SetString(PyExc_RuntimeError, SimpleLogger::get().get_buffer());
		return NULL;
	}

	// Set the output responsibility list
	uint32_t w = rqrmi_tools_get_stage_width(rqrmi_probe, stage_idx);
	PyObject* output_list = PyList_New(w);
	for (uint32_t i=0; i<w; ++i) {
		// The current responsibility is a list
		PyObject* current_responsibility = PyList_New(vector_list_get_size(result[i]));
		uint32_t index = 0;
		scalar_t* cursor = (scalar_t*)vector_list_begin(result[i]);
		for (; cursor; cursor = (scalar_t*)vector_list_iterate(result[i])) {
			// Each item in the list is a tuple (start, end)
			PyObject* interval = PyTuple_New(2);
			PyTuple_SetItem(interval, 0, PyFloat_FromDouble(cursor[0]));
			PyTuple_SetItem(interval, 1, PyFloat_FromDouble(cursor[1]));
			PyList_SetItem(current_responsibility, index++, interval);
		}
		PyList_SetItem(output_list, i, current_responsibility);
	}

	// All is good, return the list
	return output_list;
}

/**
 * @brief Python adapter for rqrmi_tools_calculate_max_error_per_submodel method
 * @param An RQRMI model Capsule
 * @param An RQRMI Probe Capsule
 * @param Integer, stage index
 * @param Integer, submodel index
 * @returns A list with two integers: the submodel's maximum error and bucket coverage.
 * @throws RuntimeError
 */
static PyObject* py_calculate_submodel_error(PyObject *self, PyObject *args) {

	PyObject *rqrmi_capsule, *probe_capsule;
	int stage_idx, submodel_idx;
	// Parse input according to format
	if (!PyArg_ParseTuple(args, "OOii:calculate_submodel_error",
			&rqrmi_capsule, &probe_capsule, &stage_idx, &submodel_idx)) {
		return NULL;
	}

	// Extract pointers from capsules
	rqrmi_model_t* rqrmi_model = (rqrmi_model_t*)PyCapsule_GetPointer(rqrmi_capsule, rqrmi_model_capsule_name);
	rqrmi_probing_t* rqrmi_probe = (rqrmi_probing_t*)PyCapsule_GetPointer(probe_capsule, rqrmi_probe_capsule_name);

	PyObject* output_list;

	try{
		scalar_pair_t error_pair = rqrmi_tools_calculate_submodel_error(rqrmi_model, rqrmi_probe, stage_idx, submodel_idx);

		// Generate the output list
		output_list = PyList_New(2);
		PyList_SetItem(output_list, 0, PyFloat_FromDouble((scalar64_t)error_pair.first));
		PyList_SetItem(output_list, 1, PyFloat_FromDouble((scalar64_t)error_pair.second));

	} catch (const std::exception& e) {
		warning(e.what());
		PyErr_SetString(PyExc_RuntimeError, SimpleLogger::get().get_buffer());
		return NULL;
	}

	// Return the list
	return output_list;
}

/**
 * @brief Converts RQRMI Matrix Capsule to list of lists
 * @param An RQRMI Matrix Capsule
 * @returns A list of lists (scalars)
 * @throws RuntimeError
 */
static PyObject* py_matrix_to_list(PyObject *self, PyObject *args) {
	PyObject *matrixcapsule;

	// Parse input according to format
	if (!PyArg_ParseTuple(args, "O:matrix_to_list", &matrixcapsule)) {
		return NULL;
	}

	// Extract pointers from capsules
	matrix_t* mat = (matrix_t*)PyCapsule_GetPointer(matrixcapsule, rqrmi_matrix_capsule_name);

	// Create the output list
	PyObject* output_list = PyList_New(mat->rows);
	for (uint32_t r=0; r<mat->rows; ++r) {
		PyObject* current_row = PyList_New(mat->cols);
		for (uint32_t c=0; c<mat->cols; ++c) {
			scalar64_t val = GET_SCALAR(mat, r, c);
			PyList_SetItem(current_row, c, PyFloat_FromDouble(val));
		}
		PyList_SetItem(output_list, r, current_row);
	}

	return output_list;
}

/**
 * @brief Holds information regarding the exported methods
 */
static PyMethodDef module_methods[] = {
	{"initialize", py_initiate, METH_VARARGS,
			"Initializes the extenstion's logging module\n"
			"Args: \n"
			"\t verbose: Set to 1 in case the log should be printed to stderr\n"
			"This method will never fail.\n"
	},
    {"create_model", py_load_model, METH_VARARGS,
    		"Load an RQRMI model to memory for fast inference.\n"
    		"Args: \n"
    		"\t buffer: A byte-array object with packed RQRMI model.\n"
    		"Returns: \n"
    		"\t An RQRMI object"
    		"Throws: RuntimeError in case of any error in loading the model\n"
    },
	{"create_matrix", py_rqrmi_matrix_create, METH_VARARGS,
			"Returns a valid RQRMI matrix object\n"
			"Args: \n"
			"\t buff: Packed matrix as byte-array object (of scalars)\n"
			"Returns: \n"
			"\t An RQRMI matrix object \n"
			"Throws: RuntimeError in case of invalid format \n"
	},
	{"create_dataset", py_create_dataset, METH_VARARGS,
			"Generates dynamic dataset from records based on parameters \n"
			"Args: \n"
			"\t probe: An RQRMI Probe object \n"
			"\t random: True in case the dataset's values should be randomized \n"
			"\t shuffle: True in case the dataset should be shuffled \n"
			"\t stage_idx: The required stage \n"
			"\t bucket_idx: The reuquired bucket \n"
			"\t samples: The number of samples to sample from generated dataset. \n"
			"Returns: \n"
			"\t An RQRMI matrix object (scalars). Each row is a sample with format [input, expected_output]. \n"
			"Throws: RuntimeError with relevant message \n"
	},
	{"create_probe", py_create_probe, METH_VARARGS,
			"Generates an RQRMI probing data structure. \n"
			"Args: \n"
			"\t records: The records the RQRMI index RQRMI matrix object (scalars) \n"
			"\t stage_list: A list of integers, each is the corresponding stage width \n"
			"Returns: \n"
			"\t An RQRMI Probe object \n"
			"Throws: RuntimeError with relevant message \n"
	},
	{"create_lookup_cpu", py_lookup_cpu_create, METH_VARARGS,
			"Generates new LookupCpu instance. \n"
			"Args: \n"
			"\t data: A matrix with range-value pairs. Each row is a pair [start, end, value]. \n"
			"\t model: An RQRMI object \n"
			"\t queue_size: The size of lookup request buffer \n"
			"\t listener: Function pointer, listener to new results \n"
			"Returns: \n"
			"\t An RQRMI lookup table object \n"
			"Throws: RuntimeError with relevant message \n"
	},
	{"evaluate_model", py_evaluate_model, METH_VARARGS,
			"Evaluates previously loaded RQRMI model with an input stream \n"
			"Args: \n"
			"\t model: An RQRMI object \n"
			"\t input: An RQRMI matrix object with inputs of scalars\n"
			"Returns: \n"
			"\t A tuple with two elements: (1) A tuple of (total_time, avg_time), and; (2) A list of outputs \n"
			"Throws: RuntimeError in case the RQRMI model is not built \n"
	},
	{"lookup_search", py_lookup_search, METH_VARARGS,
			"Request to search a key \n"
			"Args: \n"
			"\t lookup: A Lookup object \n"
			"\t input: The input to search \n"
			"Note: \n"
			"\t The result will be broadcast to the listener method \n"
			"Throws: RuntimeError in case of invalid arguments or internal error \n"
	},
	{"calculate_transition_set", py_calculate_transition_set, METH_VARARGS,
			"Calculates an RQRMI last stage's properties \n"
			"Args: \n"
			"\t rqrmi: An RQRMI Model object \n"
			"\t probe: An RQRMI Probe object \n"
			"\t stage_idx: The required stage \n"
			"Returns: A list of lists, each internal list is a transition input [ x, B(S(x-eps)), B(S(x+eps)) ] \n"
			"Throws: RuntimeError with relevant message \n"
	},
	{"calculate_responsibility", py_calculate_responsibility, METH_VARARGS,
			"Calculate the responsibility of a stage within an RQRMI model \n"
			"Args: \n"
			 "\t model: An RQRMI model \n"
			 "\t probe: An RQRMI probing data structure (Modified) \n"
			 "\t required_stage: The stage index of which to calculate responsibility \n"
			 "\t Returns: A list of responsibilities, each is a list of intervals, each is a tuple of two scalars (x,y) \n"
			"Throws: RuntimeError with relevant message \n"
	},
	{"calculate_submodel_error", py_calculate_submodel_error, METH_VARARGS,
			"Calculate the maximum error and bucket coverage of a submodel\n"
			"Args: \n"
			"\t model: An RQRMI model \n"
			"\t probe: An RQRMI Probe object \n"
			"\t stage_idx: The required stage \n"
			"\t submodel_idx: The required submodel \n"
			"Returns: \n"
			"\t A list with two values: [maximum error, bucket coverage] \n"
			"Throws: RuntimeError with relevant message \n"
	},
	{"matrix_to_list", py_matrix_to_list, METH_VARARGS,
			"Converts an RQRMI Matrix object to list of lists \n"
			"Args: \n"
			"\t matrix: An RQRMI Matrix object \n"
			"Returns: \n"
			"\t List (rows) of lists (cols) representation of the matrix (scalars). \n"
			"Throws: RuntimeError with relevant message \n"
	},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

/**
 * @brief Holds information regarding the module
 */
static struct PyModuleDef module_info = {
    PyModuleDef_HEAD_INIT,
    "rqrmi",   	/* name of module */
    NULL, 		/* module documentation, may be NULL */
    -1,       	/* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
	module_methods
};

/**
 * @brief Initiate RQRMI module
 */
PyMODINIT_FUNC PyInit_rqrmi(void) {
    return PyModule_Create(&module_info);
}


#ifdef __cplusplus
}
#endif

#endif // NOPYTHON

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

#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/sysinfo.h>
#include <stdexcept>

#include <cpu_core_tools.h>
#include <logging.h>

/**
 * @brief Returns the number of processors currently available in the system
 */
int cpu_core_tools_get_core_count() {
	return get_nprocs();
}

/**
 * @brief Return the CPU index in which the current thread is running
 * @note  On error returns -1;
 */
int cpu_core_tools_get_index_of_current_thread() {

	pthread_t my_thread = pthread_self();
	int cpu_count = cpu_core_tools_get_core_count();

	int result;
	cpu_set_t cpu_set;

	// Get the current thread cpu_set
	result = pthread_getaffinity_np(my_thread, sizeof(cpu_set_t), &cpu_set);
	if (result) {
		warning("CPU set (thread " << pthread_self() << "): cannot acquire cpu_set affinity: " << strerror(result));
		return -1;
	}

	// Find the CPU id in which the current thread is running
	for (int i=0; i<cpu_count; ++i) {
		if(CPU_ISSET(i, &cpu_set)) {
			return i;
		}
	}

	// No CPU was found in set (should never get here)
	return -1;
}

/**
 * @brief Sets the thread affinity to run on a single CPU
 * @param thread The thread
 * @param cpu_index The desired CPU index
 * @throws In case of memory allocation error, or thread migration error
 */
void cpu_core_tools_set_thread_affinity(pthread_t thread, int cpu_index) {

	int cpu_count = cpu_core_tools_get_core_count();
	int dst_cpu_index = cpu_index % cpu_count;

	info("trying to migrate the thread (" << thread << ") to CPU " << dst_cpu_index << "...");
	// Allocate memory for cpu_set object
	cpu_set_t cpu_set;

	// Clear the set
	CPU_ZERO_S(sizeof(cpu_set), &cpu_set);

	// Set the CPU index
	CPU_SET_S(dst_cpu_index, sizeof(cpu_set), &cpu_set);

	// Set the thread affinity
	int result = pthread_setaffinity_np(thread, sizeof(cpu_set), &cpu_set);
	if (result) {
		throw error("cannot set affinity of thread (%" << thread << ") " << strerror(result));
	}

	info("thread (" << thread << ") was migrated to CPU " << dst_cpu_index << " out of " << cpu_count);
}

/**
 * @brief Returns the next available physical core index
 * @param core_idx The current core index
 */
int cpu_core_tools_get_next_physical_core(int core_idx) {
	// Important: make sure the hyperthreading is disabled!
	return (core_idx+1) % cpu_core_tools_get_core_count();
}

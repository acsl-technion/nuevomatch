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

/** @file cpu_core_tools.h */

#pragma once

#include <pthread.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Returns the number of processors currently available in the system
 */
int cpu_core_tools_get_core_count();

/**
 * @brief Return the CPU index in which the current thread is running
 * @note  On error returns -1;
 */
int cpu_core_tools_get_index_of_current_thread();

/**
 * @brief Returns the next available physical core index
 * @param core_idx The current core index
 */
int cpu_core_tools_get_next_physical_core(int core_idx);

/**
 * @brief Sets the thread affinity to run on a single CPU
 * @param thread The thread
 * @param cpu_index The desired CPU index
 * @throws In case of memory allocation error, or thread migration error
 */
void cpu_core_tools_set_thread_affinity(pthread_t thread, int cpu_index);

#ifdef __cplusplus
}
#endif

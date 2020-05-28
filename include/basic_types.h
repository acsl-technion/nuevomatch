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

/** @file basic_types.h */

#pragma once

#include <stdbool.h>
#include <float.h>

// Scalar limits
#define SCALAR_EPS FLT_EPSILON
#define SCALAR_MIN FLT_MIN
#define SCALAR_MAX FLT_MAX
#define SCALAR_NEXT(X)nextafterf(X, INFINITY)
#define SCALAR_PREV(X)nextafterf(X, -INFINITY)

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)>(b)?(b):(a))

// Types
typedef unsigned long uint64_t;
typedef unsigned short uint16_t;
typedef unsigned int uint32_t;
typedef unsigned char uint8_t;
typedef float  scalar32_t;
typedef double scalar64_t;

// RQRMI works with 32 bit precision
typedef scalar32_t scalar_t;

/**
 * @brief Holds pair of scalars
 */
typedef struct {
	scalar_t first;
	scalar_t second;
} scalar_pair_t;


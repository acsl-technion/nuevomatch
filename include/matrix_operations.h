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

/** @file matrix_operations.h */

#pragma once

#include <stdio.h>
#include <basic_types.h>

#define MATRIX_ERROR 		NULL
#define OPERATION_BYPASS 	NULL

/**
 * @brief A matrix data structure
 */
typedef struct matrix {
	uint32_t rows;
	uint32_t cols;
	void* elements;
} matrix_t;

typedef scalar_t(*binary_operation_t)(scalar_t *a, scalar_t *b);
typedef scalar_t(*unary_operation_t)(scalar_t *a);

#ifdef __cplusplus
extern "C" {
#endif

// Fast getting elements
#define GET_SCALAR(mat, row, col) *(scalar_t*)get_element(mat, row, col)
#define GET_UINT32(mat, row, col) *(uint32_t*)get_element(mat, row, col)


/**
 * @brief Creates new matrix
 * @param rows Number of rows
 * @param cols Number of columns
 * @note  Matrix memory should be handled explicitly.
 */
matrix_t* new_matrix(uint32_t rows, uint32_t cols);

/**
 * @brief Loads a matrix from buffer
 * @param buffer Byte array buffer
 * @param size The size of the buffer
 * @returns The matrix or MATRIX_ERROR on error
 */
matrix_t* load_matrix(void* buffer, uint32_t size);

/**
 * @brief Free previously allocated matrix
 * @param mat Previously allocated matrix
 */
void free_matrix(matrix_t *mat);

/**
 * @brief Performs matrix multiplication without checks and memory allocations
 * @param mat_a  Left hand side matrix
 * @param mat_b  Right hand side matrix
 * @param result Preallocated matrix with enough space to store results
 * @note  This method is unsafe, make sure to allocate enough space to result
 */
void mat_mul(matrix_t* mat_a, matrix_t* mat_b, matrix_t* result);

/**
 * @brief Performs an element-wise operation between two matrices without checks and memory allocations
 * @param mat_a  Left hand side matrix
 * @param mat_b  Right hand side matrix
 * @param op     Operation to perform element-wise
 * @param result Preallocated matrix with enough space to store results
 * @note  This method is unsafe, make sure to allocate enough space to result
 */
void mat_op(matrix_t* mat_a, matrix_t* mat_b, binary_operation_t op, matrix_t* result);

/**
 * @brief Performs an element-wise operation between two matrices without checks and memory allocations
 * @param mat_a  Left hand side matrix
 * @param scalar Right hand side scalar
 * @param op     Operation to perform element-wise
 * @param result Preallocated matrix with enough space to store results
 * @note  This method is unsafe, make sure to allocate enough space to result
 */
void mat_scalar_op(matrix_t* mat_a, scalar_t scalar, binary_operation_t op, matrix_t* result);

/**
 * @brief Perform an element-wise unary operation on a matrix without checks and memory allocations
 * @param mat    Left hand side matrix
 * @param op     Operation to perform element-wise
 * @param result Preallocated matrix with enough space to store results
 * @note  This method is unsafe, make sure to allocate enough space to result
 */
void mat_unary_op(matrix_t* mat, unary_operation_t op, matrix_t* result);

/**
 * @brief Prints the matrix to the destination file
 * @param mat The matrix to print
 * @param fd The file to print the matrix to
 */
void mat_print(matrix_t* mat, FILE* fd);

/**
 * @brief Returns a pointer to an element within matrix
 * @note  No overflows are checked, unsafe
 * @note  Each element is 4 bytes
 */
inline void* get_element(matrix_t* mat, uint32_t row, uint32_t col) {
	return (uint8_t*)mat->elements + (row * mat->cols + col)*sizeof(scalar_t);
}

/**
 * @brief Interval method for comparing two scalars for quicksort
 */
int scalar_compare_asc (const void * a, const void * b);

/**
 * @brief Interval method for comparing two scalars for quicksort in reverse
 */
int scalar_compare_dsc(const void * a, const void * b);

// Element-Wise Operations
scalar_t op_add(scalar_t *a, scalar_t *b);
scalar_t op_sub(scalar_t *a, scalar_t *b);
scalar_t op_mul(scalar_t *a, scalar_t *b);
scalar_t op_div(scalar_t *a, scalar_t *b);
scalar_t op_relu(scalar_t *a);

#ifdef __cplusplus
}
#endif

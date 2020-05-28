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
#include <stdio.h>

#include <matrix_operations.h>
#include <logging.h>

/**
 * @brief Macro for fast loading from memory and perform size checks
 * @param T Typename to read from buffer
 * @param buff The buffer to read from
 * @param size Size argument for buffer overflow check.
 */
#define safe_buffer_read(T, buf, size) \
	*(T*)buf; \
	if ((size-=sizeof(T))<0) { \
		throw error("cannot read " << size << " bytes from buffer"); \
	} \
	buf=(T*)buf+1;

// Matrix information is often not required for debugging.
// Set dedicated debugging flag for matrices
#ifdef DEBUG_MATRIX
#  define matrix_info(...) info(__VA_ARGS__)
#else
#  define matrix_info(...)
#endif

/**
 * @brief Creates new matrix
 * @param rows Number of rows
 * @param cols Number of columns
 * @note  Matrix memory should be handled explicitly.
 */
matrix_t* new_matrix(uint32_t rows, uint32_t cols) {

	matrix_t* out = (matrix_t*)malloc(sizeof(matrix_t));
	if (out == NULL) {
		warning("cannot allocate memory for matrix with size [" << rows << " x " << cols << "]");
		return MATRIX_ERROR;
	}

	out->rows = rows;
	out->cols = cols;
	out->elements = (scalar_t*)malloc(sizeof(scalar_t) * rows * cols);

	if (out->elements == NULL) {
		warning("cannot allocate memory for matrix with size [" << rows << " x " << cols << "]");
		return MATRIX_ERROR;
	}
	return out;
}

/**
 * @brief Free previously allocated matrix
 * @param mat Previously allocated matrix
 */
void free_matrix(matrix_t* mat) {
	if (mat == NULL) return;
	free(mat->elements);
	free(mat);
}

/**
 * @brief Loads a matrix from buffer
 * @param buffer Byte array buffer
 * @param size The size of the buffer
 * @returns The matrix or MATRIX_ERROR on error
 */
matrix_t* load_matrix(void* buffer, uint32_t size) {
	matrix_t* output = MATRIX_ERROR;
	try {
		// Read the matrix size
		uint32_t rows, cols;
		rows = safe_buffer_read(uint32_t, buffer, size);
		cols = safe_buffer_read(uint32_t, buffer, size);

		// Allocate matrix
		if ( (output=new_matrix(rows, cols)) == MATRIX_ERROR) {
			throw error("cannot allocate memory");
		}

		// Populate matrix
		for (uint32_t r=0; r<rows; ++r) {
			for (uint32_t c=0; c<cols; ++c) {
				GET_SCALAR(output, r, c) = safe_buffer_read(scalar_t, buffer, size);
			}
		}
	} catch (const std::exception& e) {
		warning(e.what());
		free_matrix(output);
		output = MATRIX_ERROR;
	}
	return output;
}

/**
 * @brief Performs matrix multiplication without checks and memory allocations
 * @param mat_a  Left hand side matrix
 * @param mat_b  Right hand side matrix
 * @param result Preallocated matrix with enough space to store results
 * @note  This method is unsafe, make sure to allocate enough space to result
 */
void mat_mul(matrix_t* mat_a, matrix_t* mat_b, matrix_t* result) {
	// For debug
	matrix_info(
			"mat_mul needs " << mat_a->rows*mat_b->cols*sizeof(scalar_t) << "bytes for result. "
			"Result points to " << result << ", elements point to " << result->elements);
	// Set the result size
	result->rows = mat_a->rows;
	result->cols = mat_b->cols;
	// Set the result elements
	for (uint32_t row=0; row<mat_a->rows; ++row) {
		for (uint32_t col=0; col<mat_b->cols; ++col) {
			scalar_t element = 0;
			for (uint32_t i=0; i<mat_a->cols; ++i) {
				element += (*(scalar_t*)get_element(mat_a, row, i)) * (*(scalar_t*)get_element(mat_b, i, col));
			}
			*(scalar_t*)get_element(result, row, col) = element;
		}
	}
}

/**
 * @brief Perform some operation element-wise between two matrices without checks and memory allocations
 * @param mat_a  Left hand side matrix
 * @param mat_b  Right hand side matrix
 * @param op     Operation to perform element-wise
 * @param result Preallocated matrix with enough space to store results
 * @note  This method is unsafe, make sure to allocate enough space to result
 */
void mat_op(matrix_t* mat_a, matrix_t* mat_b, binary_operation_t op, matrix_t* result) {
	// For debug
	matrix_info("fast_mat_op, needs " << mat_a->rows*mat_a->cols*sizeof(scalar_t) << " bytes for result. "
			"Result points to " << result << ", elements point to " << result->elements);
	// Set the result size
	result->rows = mat_a->rows;
	result->cols = mat_a->cols;
	// Set the result elements
	for (uint32_t row=0; row<mat_a->rows; ++row) {
		for (uint32_t col=0; col<mat_a->cols; ++col) {
			*(scalar_t*)get_element(result, row, col) =
					op( (scalar_t*)get_element(mat_a, row, col),
						(scalar_t*)get_element(mat_b, row, col));
		}
	}
}

/**
 * @brief Perform some operation element-wise between two matrices without checks and memory allocations
 * @param mat_a  Left hand side matrix
 * @param scalar Right hand side scalar
 * @param op     Operation to perform element-wise
 * @param result Preallocated matrix with enough space to store results
 * @note  This method is unsafe, make sure to allocate enough space to result
 */
void mat_scalar_op(matrix_t* mat_a, scalar_t scalar, binary_operation_t op, matrix_t* result) {
	// Set the result size
	result->rows = mat_a->rows;
	result->cols = mat_a->cols;
	// Set the result elements
	for (uint32_t row=0; row<mat_a->rows; ++row) {
		for (uint32_t col=0; col<mat_a->cols; ++col) {
			*(scalar_t*)get_element(result, row, col) =
					op( (scalar_t*)get_element(mat_a, row, col), &scalar );
		}
	}
}

/**
 * @brief Perform unary operation element-wise on a matrix without checks and memory allocations
 * @param mat    Left hand side matrix
 * @param op     Operation to perform element-wise
 * @param result Preallocated matrix with enough space to store results
 * @note  This method is unsafe, make sure to allocate enough space to result
 */
void mat_unary_op(matrix_t* mat, unary_operation_t op, matrix_t* result) {
	// For debug
	matrix_info("fast_mat_unary_op, need " << mat->rows*mat->cols*sizeof(scalar_t) << " bytes for result. "
				"Result points to " << result << ", elements point to " << result->elements);
	// Set the result size
	result->rows = mat->rows;
	result->cols = mat->cols;
	// Set the result elements
	for (uint32_t row=0; row<mat->rows; ++row) {
		for (uint32_t col=0; col<mat->cols; ++col) {
			*(scalar_t*)get_element(result, row, col) = op( (scalar_t*)get_element(mat, row, col) );
		}
	}
}

/**
 * @brief Prints the matrix to the destination file
 * @param mat The matrix to print
 * @param fd The file to print the matrix to
 */
void mat_print(matrix_t *mat, FILE* fd) {
	if (mat == MATRIX_ERROR) {
		fprintf(fd, "MATRIX_ERROR");
		return;
	}
	fprintf(fd, "Matrix Size: [%u x %u]\n", mat->rows, mat->cols);
	for (uint32_t y=0; y<mat->rows; ++y) {
		fprintf(fd, "| ");
		for (uint32_t x=0; x<mat->cols; ++x) {
			fprintf(fd, "%.12f ", *(scalar_t*)get_element(mat, y, x));
		}
		fprintf(fd, "|\n");
	}
}

/**
 * @brief Interval method for comparing two scalars for quicksort
 */
int scalar_compare_asc(const void * a, const void * b) {
  return ( *(scalar_t*)a < *(scalar_t*)b ) ? -1 : 1;
}

/**
 * @brief Interval method for comparing two scalars for quicksort in reverse
 */
int scalar_compare_dsc(const void * a, const void * b) {
  return ( *(scalar_t*)a > *(scalar_t*)b ) ? -1 : 1;
}

// Element-Wise Operations
scalar_t op_add(scalar_t *a, scalar_t *b) { return *a + *b; }
scalar_t op_sub(scalar_t *a, scalar_t *b) { return *a - *b; }
scalar_t op_mul(scalar_t *a, scalar_t *b) { return *a * *b; }
scalar_t op_div(scalar_t *a, scalar_t *b) { return *a / *b; }
scalar_t op_relu(scalar_t *a) { return *a > 0 ? *a : 0; }



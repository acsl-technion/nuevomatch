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

/** @file vector_list.h */

#pragma once

#include <basic_types.h>
#include <matrix_operations.h>

#define VECTOR_LIST_ERROR NULL

typedef struct vector_list vector_list_t;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Creates a new instance of vector list of width @a width (each element is 4 bytes)
 * @param width The width of each vector in the list
 * @returns The instance or VECTOR_LIST_ERROR in case of an error
 */
vector_list_t* vector_list_create(uint32_t width);

/**
 * @brief Frees a list
 * @param list a vector list
 */
void vector_list_free(vector_list_t* list);

/**
 * @brief Returns the size of the list
 * @param list A vector list instance
 */
uint32_t vector_list_get_size(vector_list_t* list);

/**
 * @brief Adds new empty vector to the list
 * @param list Pointer to the list
 * @returns On success returns 1, otherwise 0
 */
int vector_list_push_back(vector_list_t* list);

/**
 * @brief Removes node at position fro a vector list
 * @param pos The position of the vector (negative numbers allowed)
 * @note  In case an iterator was defined, it is moved to the previous element (might be NULL)
  * @note Might throw exceptions
 */
void vector_list_remove_at(vector_list_t* list, int pos);

/**
 * @brief Returns a pointer to a specific vector inside the list
 * @param pos The position of the vector (negative numbers allowed)
 * @returns A pointer to the vector or VECTOR_LIST_ERROR in case of an error
 */
void* vector_list_get(vector_list_t* list, int pos);

/**
 * @brief Returns a pointer to the last element of a list
 * @returns A pointer to the vector.
  * @note Might throw exceptions
 */
void* vector_list_get_last(vector_list_t* list);

/**
 * @brief Push back new element to the list and return a pointer to it
 * @param list The vector list
  * @note Might throw exceptions
 */
void* vector_list_push_back_and_get(vector_list_t* list);

/**
 * @brief Restart the list iterator to the beginning
 * @returns The first available element
  * @note Might throw exceptions
 */
void* vector_list_begin(vector_list_t* list);

/**
 * @brief Iterates over all list elements
 * @param list The vector list
 * @returns The next available element or VECTOR_LIST_ERROR
  * @note Might throw exceptions
 */
void* vector_list_iterate(vector_list_t* list);

/**
 * @brief Sort the vector list by a column value
 * @param list The vector list
 * @param col The column to sort by
 * @param compare Method to compare between elements
 * @note Might throw exceptions
 */
void vector_list_sort(vector_list_t* list, uint32_t col, int(*compare)(const void*,const void*));

/**
 * @brief Converts a matrix to a vector list
 * @param mat Matrix
 * @note Might throw exceptions
 */
vector_list_t* vector_list_from_matrix_rows(matrix_t* mat);

/**
 * @brief Converts a vector list to matrix
 * @param list A vector list
 * @note Might throw exceptions
 */
matrix_t* vector_list_to_matrix(vector_list_t* list);


#ifdef __cplusplus
}
#endif

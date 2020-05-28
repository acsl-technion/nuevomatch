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
#include <string.h>

#include <vector_list.h>

#include <logging.h>
#include <matrix_operations.h>

typedef struct vector_list_node {
	struct vector_list_node* prev;
	struct vector_list_node* next;
	void* elements;
} vector_list_node_t;

struct vector_list {
	uint32_t size;
	uint32_t width;
	vector_list_node_t* head;
	vector_list_node_t* tail;
	vector_list_node_t* cursor;
};

/**
 * @brief Creates a new instance of vector list of width @a width
 * @param width The width of each vector in the list
 * @returns The instance or VECTOR_LIST_ERROR in case of an error
 */
vector_list_t* vector_list_create(uint32_t width) {
	vector_list_t* list = (vector_list_t*)malloc(sizeof(vector_list_t));
	if (list == NULL) {
		return VECTOR_LIST_ERROR;
	}
	list->size=0;
	list->width=width;
	list->head=NULL;
	list->tail=NULL;
	list->cursor=NULL;
	return list;
}

/**
 * @brief Frees a list
 * @param list a vector list
 */
void vector_list_free(vector_list_t* list) {
	if (list == VECTOR_LIST_ERROR) return;
	vector_list_node_t* current = list->head;
	while(current != NULL) {
		free(current->elements);
		vector_list_node_t* prev = current;
		current=current->next;
		free(prev);
	}
}

/**
 * @brief Returns the size of the list
 * @param list A vector list instance
 */
uint32_t vector_list_get_size(vector_list_t* list) {
	if (list == VECTOR_LIST_ERROR) return 0;
	return list->size;
}

/**
 * @brief Adds new empty vector to the list
 * @param list Pointer to the list
 * @returns On success returns 1, otherwise 0
 */
int vector_list_push_back(vector_list_t* list) {
	if (list == VECTOR_LIST_ERROR) return 0;

	vector_list_node_t* node = (vector_list_node_t*)malloc(sizeof(vector_list_node_t));
	if (node == NULL) return 0;

	void* elements = (void*)malloc(sizeof(scalar_t)*list->width);
	if (elements == NULL) {
		free(node);
		return 0;
	}

	// Initialize all elements to be zero
	memset(elements, 0, list->width*sizeof(scalar_t));
	node->prev = list->tail;
	node->next = NULL;
	node->elements = elements;

	// Update list tail
	if (list->tail != NULL) {
		list->tail->next = node;
	}
	list->tail = node;

	// Update list head
	if (list->head == NULL) {
		list->head = node;
	}

	// Update list size
	++list->size;
	return 1;
}

/**
 * @brief Removes node at position fro a vector list
 * @param pos The position of the vector (negative numbers allowed)
 * @note  In case an iterator was defined, it is moved to the previous element (might be NULL)
 * @note Might throw exceptions
 */
void vector_list_remove_at(vector_list_t* list, int pos) {
	if (list == VECTOR_LIST_ERROR || list->size == 0) {
		throw error("invalid inputs");
	}

	// Get the relevant node
	vector_list_node_t* node = list->head;
	while(pos<0) pos+=list->size;
	pos %= list->size;
	for (int i=0; i<pos; ++i) {
		node=node->next;
	}

	if (node == NULL) return;

	// In case the node is the cursor, move cursor to previous
	if (list->cursor == node) list->cursor=list->cursor->prev;

	// Update connectivity to previous
	if (node == list->head) list->head = node->next;
	else node->prev->next = node->next;

	// Update connectivity to next
	if (node == list->tail) list->tail = node->prev;
	else node->next->prev = node->prev;

	// Update size
	--list->size;

	// Remove node
	free(node->elements);
	free(node);
}

/**
 * @brief Returns a pointer to a specific vector inside the list
 * @param pos The position of the vector (negative numbers allowed)
 * @returns A pointer to the vector.
 * @note Might throw exceptions
 */
void* vector_list_get(vector_list_t* list, int pos) {
	if (list == VECTOR_LIST_ERROR || list->size == 0) {
		throw error("invalid inputs");
	}

	// Get the relevant node
	vector_list_node_t* node = list->head;
	while(pos<0) pos+=list->size;
	pos %= list->size;
	for (int i=0; i<pos; ++i) {
		node=node->next;
	}

	// Get the element
	if (node == NULL) throw error("Cannot get vector with pos " << pos);
	return node->elements;
}

/**
 * @brief Returns a pointer to the last element of a list
 * @returns A pointer to the vector.
 * @note Might throw exceptions
 */
void* vector_list_get_last(vector_list_t* list) {
	if (list == VECTOR_LIST_ERROR || list->size == 0) {
		throw error("invalid inputs");
	}
	return list->tail->elements;
}

/**
 * @brief Push back new element to the list and return a pointer to it
 * @param list The vector list
 * @note Might throw exceptions
 */
void* vector_list_push_back_and_get(vector_list_t* list) {
	if (!vector_list_push_back(list)) {
		throw error("Error while pushing back new element");
	}
	return list->tail->elements;
}

/**
 * @brief Restart the list iterator to the beginning
 * @param list The vector list
 * @returns The first available element
 * @note Might throw exceptions
 */
void* vector_list_begin(vector_list_t* list) {
	if (list == VECTOR_LIST_ERROR) {
		throw error("Got VECTOR_LIST_ERROR: " << list);
	}
	list->cursor=list->head;
	if (list->cursor != NULL) {
		return list->cursor->elements;
	}
	return NULL;
}

/**
 * @brief Iterates over all list elements
 * @param list The vector list
 * @returns The next available element or VECTOR_LIST_ERROR
 * @note Might throw exceptions
 */
void* vector_list_iterate(vector_list_t* list) {
	if (list == VECTOR_LIST_ERROR) {
		throw error("Got VECTOR_LIST_ERROR: " << list);
	}
	// In case of invalid iterator
	if (list->cursor == NULL) return VECTOR_LIST_ERROR;
	// Get next element
	list->cursor = list->cursor->next;
	// In case of invalid iterator
	if (list->cursor == NULL) return VECTOR_LIST_ERROR;
	// Return the elements of the current node
	return list->cursor->elements;
}


/**
 * @brief Converts a matrix to a vector list
 * @param mat Matrix
 * @note Might throw exceptions
 */
vector_list_t* vector_list_from_matrix_rows(matrix_t* mat) {
	if (mat == MATRIX_ERROR) {
		throw error("invalid inputs");
	}
	vector_list_t* output = vector_list_create(mat->cols);
	if (output == VECTOR_LIST_ERROR) {
		throw error("cannot allocate memory for output list");
	}
	for (uint32_t r=0; r<mat->rows; ++r) {
		scalar_t* vec = (scalar_t*)vector_list_push_back_and_get(output);
		for (uint32_t c=0; c<mat->cols; ++c) {
			vec[c]=GET_SCALAR(mat,r,c);
		}
	}
	return output;
}


/**
 * @brief Converts a vector list to matrix
 * @param list A vector list
 * @note Might throw exceptions
 */
matrix_t* vector_list_to_matrix(vector_list_t* list) {
	if (list == VECTOR_LIST_ERROR) {
		throw error("invalid inputs");
	}
	matrix_t* mat = new_matrix(list->size, list->width);
	if (mat == MATRIX_ERROR) {
		throw error("cannot allocate memory for output matrix");
	}
	vector_list_node_t* cursor = list->head;
	for (uint32_t r=0; r<mat->rows; ++r) {
		for (uint32_t c=0; c<mat->cols; ++c) {
			GET_SCALAR(mat,r,c) = ((scalar_t*)cursor->elements)[c];
		}
		cursor=cursor->next;
	}
	return mat;
}


/**
 * @brief Sort the vector list by a column value
 * @param list The vector list
 * @param col The column to sort by
 * @param compare Method to compare between elements
 * @note Might throw exceptions
 */
void vector_list_sort(vector_list_t* list, uint32_t col, int(*compare)(const void*,const void*)) {

	// Validate inputs
	if (list == VECTOR_LIST_ERROR || list->width <= col) {
		throw error("vector_list_sort: Got invalid input");
	}

	// Create a temporary array with shuffled values
	scalar_t* arr = (scalar_t*)malloc(sizeof(scalar_t)*list->size*list->width);
	if (arr == NULL) throw error("Cannot allocate memory for array");
	int counter = 0;

	// Fill data to temporary array
	vector_list_node_t* node = list->head;
	for (uint32_t i=0; i<list->size; i++) {
		scalar_t* elements = (scalar_t*)node->elements;
		arr[counter++] = elements[col];
		for (uint32_t j=0; j<list->width; j++) {
			if (j == col) continue;
			arr[counter++] = elements[j];
		}
		node=node->next;
	}

	// Quicksort the temporary array
	qsort(arr, list->size, list->width*sizeof(scalar_t), compare);

	// Update list values
	counter = 0;
	node = list->head;
	for (uint32_t i=0; i<list->size; i++) {
		scalar_t* elements = (scalar_t*)node->elements;
		elements[col]=arr[counter++];
		for (uint32_t j=0; j<list->width; j++) {
			if (j == col) continue;
			elements[j]=arr[counter++];
		}
		node=node->next;
	}

	// Free temporary memory
	free(arr);
}


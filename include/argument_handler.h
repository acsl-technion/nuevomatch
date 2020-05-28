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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/**
 * Used to define and acquire command line arguments
 */
typedef struct {
	// Is set by the user
	const char* name;
	int required;
	int is_boolean;
	const char* value;
	const char* help;

	// Is set by parse_arguments method
	int available;

} argument_t;


/**
 * @brief Get argument by full-name
 */

#define ARG(x) get_argument_by_name(my_arguments, x)

/**
 * @brief Parse program arguments
 */
void parse_arguments(int argc, char** argv, argument_t* required_args);

/**
 * @brief Returns the relevant argument by its name
 */
inline argument_t* get_argument_by_name(argument_t* required_args, const char* name) {
	while (required_args->name != NULL && strcmp(name, required_args->name) != 0) ++required_args;
	if (required_args == NULL) {
		fprintf(stderr, "Fatal error: cannot get argument %s by name\n", name);
		exit(1);
	}
	return required_args;
}

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

#include <argument_handler.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

const char* help_0 = "-h";
const char* help_1 = "--help";

/**
 * @brief Parse program arguments
 */
void parse_arguments(int argc, char** argv, argument_t* required_args) {

	argument_t* current_arg = NULL;
	int return_value = 0;

	// Parse the arguments
	for(int idx=1; idx<argc; ++idx) {

		// Check which argument is it
		current_arg = required_args;
		while (current_arg->name != NULL) {

			// Check whether there is match
			if (strcmp(argv[idx], current_arg->name) == 0) {
				// Update the current argument
				current_arg->available=1;
				if (current_arg->is_boolean==0) {
					current_arg->value = argv[++idx];
				}
				break;
			}
			// Check whether help is requested
			else if ((strcmp(argv[idx], help_0) == 0) ||
					 (strcmp(argv[idx], help_1) == 0))
			{
				goto show_help;
			}

			// Go to the next argument
			++current_arg;
		}

		// In case the argument was found
		if (current_arg->available==1) {
			continue;
		}

		// No argument was found, show error
		printf("Argument %s is not defined\n\n", argv[idx]);
		return_value = 1;
		goto show_help;
	}

	// Check all the required arguments are available
	current_arg = required_args;
	while (current_arg->name != NULL) {
		if (current_arg->required && !current_arg->available) {
			printf("Argument %s is missing\n\n", current_arg->name);
			return_value = 1;
			goto show_help;
		}
		++current_arg;
	}

	return;

	// Show help
	show_help:

	// Try to print general description
	current_arg = required_args;
	while (current_arg->name != NULL) {
		++current_arg;
	};
	if (current_arg->help != NULL) {
		printf("%s\n", current_arg->help);
	}

	printf("Usage %s ", argv[0]);

	// Print mandatory arguments
	current_arg = required_args;
	while (current_arg->name != NULL) {
		if (current_arg->required) {
			printf("%s ", current_arg->name);
			if (!current_arg->is_boolean) {
				printf("<value> ");
			}
		}
		++current_arg;
	}

	// Print optional arguments
	current_arg = required_args;
	while (current_arg->name != NULL) {
		if (!current_arg->required) {
			printf("[%s", current_arg->name);
			if (!current_arg->is_boolean) {
				printf(" <value>] ");
			} else {
				printf("] ");
			}
		}
		++current_arg;
	}
	printf("\n");

	// Print show help
	printf("Use %s or %s to show this message\n", help_0, help_1);

	// Print help
	current_arg = required_args;
	while (current_arg->name != NULL) {
		// Build argument name
		char long_name[150] = {0};
		strcat(long_name, current_arg->name);
		if (!current_arg->is_boolean) strcat(long_name, " <value>");
		if (!current_arg->required) strcat(long_name, " (optional)");
		// Print spec
		printf("%s", long_name);
		if (strlen(long_name) > 25) {
			printf("\n%-25s", "");
		}
		else {
			for (int i=strlen(long_name); i<25; ++i) {
				printf(" ");
			}
		}
		printf(" %s", current_arg->help);
		if (current_arg->value) {
			printf(" (default: %s)", current_arg->value);
		}
		printf("\n");
		++current_arg;
	}
	printf("\n");

	exit(return_value);
}

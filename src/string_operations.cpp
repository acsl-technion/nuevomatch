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

#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>
#include <regex>
#include <list>

#include <string_operations.h>

namespace string_operations {

/**
 * @brief Splits a string base of delimiter char
 * @param str A string to split
 * @param delim A regex object
 * @return A vector of strings
 */
std::vector<std::string> split(const std::string& str, const std::regex& delim) {
	// Initiate regex iterators
    std::sregex_token_iterator it{str.begin(), str.end(), delim, -1}, end;

    // Read all string items from the original string, and parse them
    std::vector<std::string> items;
    while (it != end) {
    	std::string current = *it++;
    	// Skip empty strings
    	if (current.size() != 0) {
    		items.push_back(current);
    	}
    }

    return items;
}

/**
 * @brief Splits string to a vrctor of non-empty strings using a string delimiter
 * @param str The string to split
 * @param delim A string delimiter, each char is a delimiter
 * @returns A vector of strings
 */
std::vector<std::string> split(const std::string& str, const std::string& delim) {

	// Allocate new string, copy from the original
	char* new_string = strdup(str.c_str());
	int str_size = (int)str.size();
	int delim_size = (int)delim.size();
	int lst_idx = 0;

	std::vector<std::string> output;

	// Go over all chars
	for (int i=0; i<(int)str_size; ++i) {
		// Go over all delimiters
		// Convert all delimiter chars to \0
		// Add non-empty strings to output vector
		for (int j=0; j<(int)delim_size; ++j) {
			if (new_string[i] == delim[j]) {
				new_string[i] = 0;
				// In case the new string is not empty
				if (i-lst_idx > 0) {
					output.push_back(std::string(&new_string[lst_idx], &new_string[i]));
				}
				lst_idx=i+1;
				break;
			}
		}
	}

	// Add last field
	if (str_size-lst_idx > 0) {
		output.push_back(std::string(&new_string[lst_idx], &new_string[str_size]));
	}

	free(new_string);
	return output;
}

/**
 * @brief Converts an hex representation to 32bit integer
 * @param str A string with hex representation of a number
 * @note Taken from here: https://stackoverflow.com/a/39052987/4103200
 */
uint32_t hex2int(const std::string& str) {
    uint32_t val = 0;
    const char* hex = str.c_str();

    // Skip the chars 0x if exist
    if (*hex == '0' && *(hex+1) == 'x') hex+=2;

    while (*hex) {
        // get current character then increment
        uint8_t byte = *hex++;
        // transform hex character to the 4bit equivalent number, using the ascii table indexes
        if (byte >= '0' && byte <= '9') byte = byte - '0';
        else if (byte >= 'a' && byte <='f') byte = byte - 'a' + 10;
        else if (byte >= 'A' && byte <='F') byte = byte - 'A' + 10;
        // shift 4 to make space for new digit, and add the 4 bits of the new digit
        val = (val << 4) | (byte & 0xF);
    }
    return val;
}

/**
 * @brief Converts a string representation of integer to 32bit integer
 */
uint32_t str2int(const std::string& str) {
	return atoi(str.c_str());
}





};


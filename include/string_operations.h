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

#include <stdexcept>
#include <string>
#include <list>
#include <vector>
#include <regex>

namespace string_operations {

/**
 * @brief A function pointer type that gets a string and parses it to be of type T
 */
template<typename T>
using string_parser = T(*)(const std::string& item);

/**
 * @brief A function pointer type that gets a type T and returns its string representation
 */
template<typename T>
using convertor_to_string = std::string(*)(const T& item);

/**
 * @brief Splits a string base of delimiter char
 * @param str A string to split
 * @param delim A regex object
 * @param parser A function that parses the string elements to be of type T
 * @return A vector of T
 */
template <typename T>
std::vector<T> split(const std::string& str, const std::regex& delim, string_parser<T> parser) {
	// Initiate regex iterators
    std::sregex_token_iterator it{str.begin(), str.end(), delim, -1}, end;

    // Read all string items from the original string, and parse them
    std::vector<T> items;
    while (it != end) {
    	std::string current = *it++;
    	// Skip empty items
    	if (current.size() == 0) {
    		continue;
    	}
    	// Parse all items according to parser
    	items.push_back(parser(current));
    }

    return items;
}

/**
 * @brief Group a list of items to one string (using parser) with glue string between each
 * @param lst An list of items
 * @param glue String to glue elements with
 * @param parser A function pointer (or lambda expression) for converting each element to string
 */
template<typename T>
std::string join(const std::list<T>& lst, std::string glue, convertor_to_string<T> parser) {
	std::stringstream ss;
	uint32_t counter=0;
	for (auto it : lst) {
		ss << parser(it);
		if (counter++ != lst.size()-1) {
			ss << glue;
		}
	}
	return ss.str();
}

/**
 * @brief Splits a string base of delimiter char
 * @param str A string to split
 * @param delim A regex object
 * @return A vector of strings
 */
std::vector<std::string> split(const std::string& str, const std::regex& delim);

/**
 * @brief Splits string to a vrctor of non-empty strings using a string delimiter
 * @param str The string to split
 * @param delim A string delimiter, each char is a delimiter
 * @returns A vector of strings
 */
std::vector<std::string> split(const std::string& str, const std::string& delim);

/**
 * @brief Splits string to a vrctor of non-empty strings using a string delimiter
 * @param str The string to split
 * @param delim A string delimiter, each char is a delimiter
 * @param parser A function that parses the string elements to be of type T
 * @returns A vector of strings
 */
template <typename T>
std::vector<T> split(const std::string& str, const std::string& delim, string_parser<T> parser) {
	std::vector<std::string> vec = std::move(split(str, delim));
	std::vector<T> output;
	for (auto it : vec) {
		output.push_back(parser(it));
	}
	return output;
}

/**
 * @brief Converts an hex representation to 32bit integer
 * @param str A string with hex representation of a number
 * @note Taken from here: https://stackoverflow.com/a/39052987/4103200
 */
uint32_t hex2int(const std::string& str);

/**
 * @brief Converts a string representation of integer to 32bit integer
 */
uint32_t str2int(const std::string& str);

};

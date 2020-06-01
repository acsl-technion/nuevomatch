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

#include <rule_db.h>

/**
 * @brief Simple hash table for caching exact-match flows
 */
class ExactMatchTable {
public:

	ExactMatchTable(int size, bool bypass);

	~ExactMatchTable();

	/**
	 * @brief Checks whether the packet header is cached
	 * @returns The priority of the matched rule or -1 for miss.
	 */
	int lookup(trace_packet& packet);

	/**
	 * @brief Cache a trace packet and with its matching rule
	 */
	void add(trace_packet& packet, int priority);

	/**
	 * @brief Invalidates all entries
	 */
	void invalidate();

	/*
	 * @brief Returns the ratio of valid entries vs size
	 */
	double utilization() const;

private:


	static constexpr int num_of_ways = 2;

	// Cache way
	struct cache_way {
		trace_packet packet;
		int priority;
		bool valid;
	};

	// Table entry
	struct header_entry {
		cache_way way[num_of_ways];
		int lru;
	};

	// Hash table, header -> priority
	header_entry* table;

	// Size of table
	int size;
	bool bypass;

	int match(trace_packet& packet, header_entry* entry);
};

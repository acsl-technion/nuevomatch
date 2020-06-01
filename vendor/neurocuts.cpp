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

#include <assert.h>
#include <string>
#include <queue>
#include <sstream>
#include <iostream>
#include <queue>

#include <object_io.h>
#include <neurocuts.h>
#include <logging.h>

NeuroCuts::NeuroCuts() :
	_num_of_rules(0), _size(0), _build_time(0), _max_binth(0),
	_nodes(nullptr), _rules(nullptr), _root(nullptr),
	_is_clone(false) {}

NeuroCuts::~NeuroCuts() {
	if (!_is_clone) {
		delete[] _nodes;
		delete[] _rules;
	}
}

/**
 * @brief Creates this from a memory location
 * @param object An object-reader instance
 */
void NeuroCuts::load(ObjectReader& reader) {

    // Initiate packer
    uint32_t index;

	// Read build time
	reader >> _build_time;
	// Read all rules
	reader >> _num_of_rules;
	_rules = new matching_rule[_num_of_rules];
	for (uint32_t i=0; i<_num_of_rules; ++i) {
		reader >> _rules[i].priority;
		for (uint32_t f=0; f<FIVE_TUPLE_FIELDS; ++f) {
			reader >> _rules[i].field[f].low;
			reader >> _rules[i].field[f].high;
		}
	}

	// Reset statistics
	_max_binth = 0;

	// Read all nodes
	uint32_t num_of_nodes;
    reader >> num_of_nodes;
	_nodes = new node[num_of_nodes];
	for (uint32_t i=0; i<num_of_nodes; ++i) {

		// Read node's static attributes
		reader >> _nodes[i].id;
		reader >> _nodes[i].depth;

		_nodes[i].max_priority = -2;

		uint32_t item;
		reader >>  item;
		_nodes[i].is_partition = (item == 0) ? true : false;

		for (uint32_t f=0; f<FIVE_TUPLE_FIELDS; ++f) {
			reader >> _nodes[i].boundaries[f].low;
			reader >> _nodes[i].boundaries[f].high;
		}

		// Read all rules
		reader >> _nodes[i].num_of_rules;
		_nodes[i].rules = new matching_rule*[_nodes[i].num_of_rules];
		for (uint32_t j=0; j<_nodes[i].num_of_rules; ++j) {
			reader >> index;
			_nodes[i].rules[j] = &_rules[index];
		}

		// Read all children
		reader >> _nodes[i].num_of_children;
		_nodes[i].children = new node*[_nodes[i].num_of_children];
		for (uint32_t j=0; j<_nodes[i].num_of_children; ++j) {
			reader >> index;
			assert(index < num_of_nodes);
			_nodes[i].children[j] = &_nodes[index];
			_nodes[index].parent = &_nodes[i];
		}

		// Update max binth
		if (_nodes[i].num_of_children == 0) {
			_max_binth = std::max(_max_binth, _nodes[i].num_of_rules);
		}

	}

	// Set root node
	_root = &_nodes[0];
	_root->parent = nullptr;

	// Calculate max priority
	for (uint32_t i=0; i<num_of_nodes; ++i) {
		node* current = &_nodes[i];
		// Find max priority from the rules of this
		int max_priority = 0;
		for (uint32_t r=0; r<current->num_of_rules; ++r) {
			max_priority = std::max(max_priority, (int)current->rules[r]->priority);
		}
		current->max_priority = max_priority;
		// Update parents
		for (; current != nullptr; current=current->parent) {
			current->max_priority = std::max(current->max_priority, max_priority);
		}
	}

	// Calculate size
	for (uint32_t i=0; i<num_of_nodes; ++i) {
		node* current = &_nodes[i];
		_size += 1; // 1 Byte for action type
		// Leaf node
		if (current->num_of_children == 0) {
			// Rule pointer is 4 bytes
			_size += 4 * current->num_of_rules;
			// Number of rules is an integer
			_size += 4;
			continue;
		}
		// For fast calculation we store boundaries for both CUT and PARTITION nodes
		// each boundary is 8 bytes
		_size += FIVE_TUPLE_FIELDS * 8;
		// Each child pointer is 4 bytes
		_size += 4 * current->num_of_children;
		// Number of children is an integer
		_size += 4;
	}

	assert(reader.size() == 0);
}

/**
 * @brief Starts the performance measurement of this
 * @note Added by Alon Rashelbach
 */
void NeuroCuts::start_performance_measurement() {
	clock_gettime(CLOCK_MONOTONIC, &perf_start_time);
}

/**
 * @brief Stops the performance measurement of this
 * @note Added by Alon Rashelbach
 */
void NeuroCuts::stop_performance_measurement() {
	clock_gettime(CLOCK_MONOTONIC, &perf_end_time);
}

/**
 * @brief Helper recursive method for finding matching node
 * @param current The node to check
 * @param header The input packet header
 * @returns A pointer to a matching rule
 */
matching_rule* NeuroCuts::match(node* current, const uint32_t* header, int priority) {
	// In case the node is partition
	if (current->is_partition) {
		matching_rule* matches[current->num_of_children];

		// All children of this must be checked
		for(uint32_t i=0; i< current->num_of_children; ++i) {
			matches[i] = match(current->children[i], header, priority);
		}

		// Return the highest priority rule
		matching_rule* output = nullptr;
		for(uint32_t i=0; i< current->num_of_children; ++i) {
			if ( (matches[i]) && ((priority<0) || (priority > matches[i]->priority)) ) {
				if (output == nullptr) {
					output = matches[i];
				} else if (output->priority > matches[i]->priority) {
					output = matches[i];
				}
			}
		}
		return output;
	}
	// This node is CUT.
	// In case it has any children
	else if (current->num_of_children > 0) {
		// Return the first child that contains the packet
		for(uint32_t i=0; i< current->num_of_children; ++i) {
			bool child_match = true;
			node* child = current->children[i];
			for (uint32_t f=0; f<FIVE_TUPLE_FIELDS; ++f) {
				if (header[f] < child->boundaries[f].low || header[f] >= child->boundaries[f].high) {
					child_match = false;
					break;
				}
			}
			// In case the current child matches the packet
			if (child_match) {
				return match(child, header, priority);
			}
		}
		// No child matches the packet
		return nullptr;
	}
	// This node is CUT, but has no children
	else {
		for (uint32_t r=0; r<current->num_of_rules; ++r) {
			bool rule_match = true;
			matching_rule* rule = current->rules[r];
			for (uint32_t f=0; f<FIVE_TUPLE_FIELDS; ++f) {
				if (header[f] < rule->field[f].low || header[f] >= rule->field[f].high) {
					rule_match = false;
					break;
				}
			}
			// The current rule matches the packet
			if ( (rule_match) && ((priority<0) || (priority > rule->priority)) ) {
				return rule;
			}
		}
		// No rule matches the packet
		return nullptr;
	}
}

/**
 * @brief Start a synchronous process of classification an input packet.
 * @param header An array of 32bit integers according to the number of supported fields.
 * @returns The matching rule action/priority (or 0xffffffff if not found)
 */
uint32_t NeuroCuts::classify_sync(const uint32_t* header, int priority) {
	matching_rule* rule = match(_root, header, priority);
	return (rule == nullptr) ? priority : rule->priority;
}

/**
 * @brief Start an asynchronous process of classification for an input packet.
 * @param header An array of 32bit integers according to the number of supported fields.
 * @returns A unique id for the packet
 */
uint32_t NeuroCuts::classify_async(const uint32_t* header, int priority) {
	uint32_t packet_id = 0xffffffff;
	matching_rule* rule = nullptr;

	// Perform lookup
	// Skip invalid packets
	if (header) {
		rule = match(_root, header, priority);
		packet_id = _packet_counter++;
	}

	// In case the rule was found
	if (rule != nullptr) {
		priority = rule->priority;
	}

	// Broadcast result
	for (auto it : _listeners) {
		it->on_new_result(packet_id, priority, priority, _additional_args);
	}

	return packet_id;
}

/**
 * @brief Computes the memory access of this
 * @param n The node from which start calculating. At normal cases, should be root.
 */
uint32_t NeuroCuts::compute_mem_access(node* n) const {
	// Leaf
	if (n->num_of_children == 0) {
		return 1;
	}
	// Partition - sum of all children
	else if (n->is_partition) {
		uint32_t sum = 0;
		for (uint32_t i=0; i< n->num_of_children; ++i) {
			sum += compute_mem_access(n->children[i]);
		}
		return sum;
	}
	// Cut - max of all children
	else {
		uint32_t max = 0;
		for (uint32_t i=0; i< n->num_of_children; ++i) {
			max =  std::max(max, compute_mem_access(n->children[i]));
		}
		return 1+max;
	}
}


/**
 * @brief Prints debug information
 * @param verbose Set the verbosity level of printing
 */
void NeuroCuts::print(uint32_t verbose) const {
	// Measure performance
	uint32_t total_usec = (perf_end_time.tv_sec * 1e6 + perf_end_time.tv_nsec / 1e3) -
							  (perf_start_time.tv_sec * 1e6 + perf_start_time.tv_nsec / 1e3);

	messagef("Performance: total time %u usec. Average time: %.3f usec per packet.", total_usec, (double)total_usec / _packet_counter);

	// Medium Verbosity
	if (verbose > 1) {
		// Compute memory access
		messagef("Required memory accesses (worst-case): %u, binth: %u", compute_mem_access(_root), _max_binth);
	}
}

/**
 * @brief Used for debugging. Prints a NeuroCuts node
 */
void NeuroCuts::print_node(node* n) {
	std::stringstream s;
	s << "ID: " << n->id << ", "
	  << "Action: " << (n->is_partition ? "partition" : " cut") << " "
	  << "Depth: " << n->depth
	  << "Range: [";

	for (uint32_t f=0; f<FIVE_TUPLE_FIELDS; ++f) {
		s << n->boundaries[f].low << "," << n->boundaries[f].high;
		if (f < FIVE_TUPLE_FIELDS-1) s << "; ";
	}

	s << "]" << std::endl
	  << "Children: ";
	for(uint32_t i=0; i< n->num_of_children; ++i) {
		s << n->children[i]->id << " ";
	}
	s << std::endl;

	s << "Rules: " << std::endl;
	for (uint32_t r=0; r<n->num_of_rules; ++r) {
		for (uint32_t f=0; f<FIVE_TUPLE_FIELDS; ++f) {
			s << n->rules[r]->field[f].low << ", " << n->rules[r]->field[f].high;
			if (f < FIVE_TUPLE_FIELDS-1) s << "; ";
		}
		s << std::endl;
	}
	std::cout << s.str();
}

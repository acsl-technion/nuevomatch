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

#include <iostream>
#include <stdexcept> // runtime error
#include <regex> // regex
#include <fstream> // file stream
#include <algorithm> //copy, find
#include <sstream> // string stream
#include <set>

#include <object_io.h>
#include <rule_db.h>
#include <array_operations.h>
#include <string_operations.h>


using namespace std;

/**
 * @brief Returns a string representation of this
 */
string trace_packet::to_string() const {
	stringstream ss;
	for (auto f : this->header) {
		ss << f << "\t";
	}
	ss << this->match_priority;
	return ss.str();
}

/**
 * @brief Returns a string representation of this
 */
std::string range::to_string() const {
	stringstream ss;
	ss << low << "-" << high;
	return ss.str();
}

/**
 * @brief Create new OpenFlow rule from a classic 5-tuple rule
 */
openflow_rule::openflow_rule(const matching_rule& rule) {
	priority = rule.priority;
	for (int i=0; i<5; ++i) {
		fields.push_back(rule.field[i]);
	}
}

/**
 * @brief Returns a string representation of this
 */
std::string openflow_rule::to_string() const {
	stringstream ss;
	for (auto f : fields) {
		ss << f.to_string() << ",";
	}
	ss << "prio: " << this->priority;
	return ss.str();
}

/**
 * @brief Converts this to a 5-tuple rule
 * @throws error In case the rule cannot be converted
 */
matching_rule openflow_rule::convert_to_five_tuple() const {
	if (fields.size() != 5) {
		throw error("rule does not have exactly 5 fields");
	}
	matching_rule rule;
	for (int i=0; i<5; ++i) {
		rule.field[i] = fields[i];
	}
	rule.priority=priority;
	return rule;
}

/**
 * @brief Used for sorting rules in STL containers
 */
bool openflow_rule::operator<(const openflow_rule& rhs) const {
	return this->priority < rhs.priority;
}

/**
 * @brief Loads rules database from memory
 * @param reader An object-reader with binary data
 * @return A list of open-flow rules
 */
list<openflow_rule> load_rule_database(ObjectReader& reader) {

	list<openflow_rule> output;

	uint32_t num_of_rules, num_of_fields, version_number;

	// Read the number of rules
	reader >> num_of_rules;
	num_of_fields = FIVE_TUPLE_FIELDS;
	version_number = 0;

	// This is an extension for newer versions - ugly patch
	if (num_of_rules == 0) {
		// Read version number
		reader >> version_number;
	}

	// Act according to version number
	switch(version_number) {
	// Nothing new, original version of rule-db
	case 0: break;
	// Add support for arbitrary number of fields
	case 1:
		reader >> num_of_rules >> num_of_fields;
		break;
	}

	// Read the database from file
	for (uint32_t i=0; i<num_of_rules; ++i) {
		openflow_rule rule;
		reader >> rule.priority;
		rule.fields.resize(num_of_fields);
		// Set the rule's fields
		for (uint32_t j=0; j<num_of_fields; ++j) {
			reader >> rule.fields[j].low >> rule.fields[j].high;
		}
		output.push_back(rule);
	}

	// Return the output
	return output;
}

/**
 * @brief Reads a ruleset indices file and return a set with all indices available in file
 * @param reader An object-reader with binary data
 */
std::set<uint32_t> read_indices_file(ObjectReader& reader) {
	std::set<uint32_t> output;
	uint32_t size = reader.read<uint32_t>();
	for (uint32_t i=0; i<size; ++i) {
		uint32_t current = reader.read<uint32_t>();
		output.insert(current);
	}
	return output;
}

/**
 * @brief Truncate the ruleset according to a set of indices
 * @param ruleset Original ruleset
 * @param indices A set of indices
 * @param reverse If true, select the rules with indices not in the indices set
 */
list<openflow_rule> apply_indices_on_ruleset(std::list<openflow_rule>& ruleset, std::set<uint32_t>& indices, bool reverse) {
	list<openflow_rule> output;
	auto it = ruleset.cbegin();
	for (uint32_t i=0; i<ruleset.size(); ++i) {
		bool in_set = indices.find(i) != indices.end();
		if ((in_set && (!reverse)) || ((!in_set) && reverse)) {
			output.push_back(*it);
		}
		++it;
	}
	return output;
}

/**
 * @brief Reads a textual trace file into memory
 * @param[in] trace_filename The textual trace filename
 * @param[in] indices A vector of custom fields to look at
 * @param[out] num_of_packets The number of packet in trace
 * @returns An array of trace packet headers, or NULL in case of an error
 */
trace_packet* read_trace_file(const char* trace_filename, const std::vector<uint32_t>& indices, uint32_t* num_of_packets) {
	// Open file
	fstream fs;
	fs.open(trace_filename, fstream::in);

	// Check whether file exists
	if (!fs.good()) {
		throw error("cannot open file for reading trace: file does not exist");
	}

	// Output list
	list<trace_packet> output;

	trace_packet packet;
	char buffer[25], c;
	int current = 0;

	// Flags
	bool field_end = false, check_packet = false;

	do {
		// Get char
		c = fs.get();

		// In case of number
		if ((c >= '0') && (c <= '9')) {
			buffer[current++]=c;
		}
		// In case of field delimiter
		else if ( (c == ' ') || (c == '\t')) {
			field_end = true;
		}
		// In case of new line or EOF
		else if ((c == '\n') || (fs.eof())) {
			field_end = true;
			check_packet = true;
		}

		// The current field was ended
		if (field_end) {
			if (current > 0) {
				buffer[current++] = '\0';
				packet.header.push_back(atoi(buffer));
			}
			field_end = false;
			current = 0;
		}

		// The whole packet as ended
		if (check_packet) {
			// In case the current packet is not empty
			if (packet.header.size() > 0) {

				// Remove last header, set as match index
				packet.match_priority = packet.header.back();
				packet.header.pop_back();

				// Shuffle header according to indices
				if (indices.size() > 0) {
					std::vector<uint32_t> new_header(indices.size());
					for (uint32_t i=0; i<indices.size(); ++i) {
						if (indices[i] > packet.header.size()) {
							throw errorf("Cannot extract field %d from trace packet, as it has only %d fields",
									indices[i], packet.header.size());
						}
						new_header[i] = packet.header[indices[i]];
					}
					packet.header = new_header;
				}

				// Add packet to output list
				output.push_back(std::move(packet));
			}
			check_packet = false;
		}
	} while (!fs.eof());

	*num_of_packets = output.size();
	return list_to_array(output);
}

/**
 * @brief Parse an IPv4-mask string (xxx.xxx.xxx.xxx/xx)
 * @param ip_address The IP-mask string
 * @return The range (in 32bit space) as {start, end}
 */
vector<uint32_t> parse_ip_mask_address(const string& ip_address) {

	// Split string to numeric components
	static const regex delim("\\.|\\/");

	vector<uint32_t> parts = string_operations::split(
			ip_address, delim, string_operations::str2int);

	if (parts.size() != 5) {
		throw error("IP/mask string is invalid");
	}

	// Mask
	if (parts[4]>0) parts[4] = 0xffffffff << (32-parts[4]) & 0xffffffff;
	else parts[4]=0;

	uint32_t ip_start = (parts[0] << 24 | parts[1] << 16 | parts[2] << 8 | parts[3]) & parts[4];
	uint32_t ip_end   = ip_start | ~parts[4];
	return {ip_start, ip_end};
}

/**
 * @brief Returns an IP string prefix length
 */
uint32_t get_ip_prefix_length(const string& ip_address) {
	// Split string to numeric components
	static const regex delim("\\.|\\/");
	vector<uint32_t> parts = string_operations::split(
			ip_address, delim, string_operations::str2int);
	return parts[4];
}

/**
 * @brief Parses a MAC address and returns its least significant 32 bits
 * @param mac_address A string representation of a MAC address
 */
uint32_t parse_mac_address(const string& mac_address) {
	auto parts = string_operations::split(mac_address, ":");
	uint32_t out = 0;
	for (auto p : parts) {
		out = (out << 8) | string_operations::hex2int(p);
	}
	return out;
}

/**
 * @brief Parses protocol range (0xXXXX/0xXXXX)
 */
vector<uint32_t> parse_hex_range(const string& str) {
	auto vals = string_operations::split(str, "/", string_operations::hex2int);
	if (vals[1] == 0xff) {
		return {vals[0], vals[0]};
	} else{
		return {0, 255};
	}
}

/**
 * @brief Reads a textual indices file
 * @param filename Path to Classbench file
 * @return A list of integers
 */
std::list<uint32_t> read_indices_file(const char* filename) {
	std::list<uint32_t> output;
	try{
		// Open file
		fstream fs;
		fs.open(filename, fstream::in);



	} catch (exception& e) {
		warningf("Cannot load indices file: %s", e.what());
	}
	return output;
}

/**
 * @brief Reads Classbench file into rule_table_t data structure
 * @param filename Path to Classbench file
 * @return A list of open-flow rules
 */
std::list<openflow_rule> read_classbench_file(const char* filename) {

	list<openflow_rule> output;
	openflow_rule rule;

	try{
		// Open file
		fstream fs;
		fs.open(filename, fstream::in);

		while (1) {

			// Read next line
			char line_buffer[2048];
			fs.getline(line_buffer, 2048);
			string line(line_buffer);

			// Stop reading file
			if (fs.eof()) break;

			// Skip empty lines
			if (line.size() == 0) {
				continue;
			}

			// Split line according to delimiters
			auto fields = string_operations::split(line, "@ \t");

			// Validate there are 6 fields:
			if (fields.size() != 10) {
				throw error("Classbench line has illegal number of fields: %lu", fields.size());
			}

			// Validate the fields 3 and 6 are ":" according to the Classbench format
			if (fields[3].compare(":") || fields[6].compare(":")) {
				throw error("Classbench line: field 3 is '%s', field 6 is '%s'; both should be ':'",
						fields[3].c_str(), fields[6].c_str());
			}

			// Create new rule
			rule.fields.resize(5);
			rule.fields[0] = parse_ip_mask_address(fields[0]); // src-ip
			rule.fields[1] = parse_ip_mask_address(fields[1]); // dst-ip
			rule.fields[2].low  = string_operations::str2int(fields[2]); // src-port from
			rule.fields[2].high = string_operations::str2int(fields[4]); // src-port to
			rule.fields[3].low  = string_operations::str2int(fields[5]); // dst-port from
			rule.fields[3].high = string_operations::str2int(fields[7]); // dst-port to
			rule.fields[4] = parse_hex_range(fields[8]);	   // protocol
			rule.priority = output.size();

			// Get the rule prefixes, for hash-based algoes such as TM
			rule.prefixes.resize(5);
			rule.prefixes[0] = get_ip_prefix_length(fields[0]);
			rule.prefixes[1] = get_ip_prefix_length(fields[1]);
			rule.prefixes[2] = (rule.fields[2].low == rule.fields[2].high) ? 32 : 16;
			rule.prefixes[3] = (rule.fields[3].low == rule.fields[3].high) ? 32 : 16;
			rule.prefixes[4] = (rule.fields[4].low == rule.fields[4].high) ? 32 : 24;

			output.push_back(std::move(rule));
		}
	} catch (exception& e) {
		warningf("Cannot load Classbench file: %s", e.what());
	}
	return output;
}

/**
 * @brief Reads Classbench-ng file with OpenFlow rules
 * @param filename Path to Classbench-ng file
 * @return A list of open-flow rules
 */
std::list<openflow_rule> read_classbench_ng_file(const char* filename) {
	try {
		// Open file
		fstream fs;
		fs.open(filename, fstream::in);

		// Create the regex for delimiters
		regex field_delim(",\\s+"), kv_delim("\\=");

		// Read line by line
		list<openflow_rule> output;
		while (1) {

			// Read next line
			char line_buffer[2048];
			fs.getline(line_buffer, 2048);
			string line(line_buffer);

			if (fs.eof()) break;

			// Initiate new rule with zeros
			openflow_rule rule;
			rule.fields.resize(9);
			for (auto i=0u; i<rule.fields.size(); ++i) {
				rule.fields[i].low = rule.fields[i].high = 0;
			}
			rule.priority = output.size();

			// Get all fields from the current line
			auto elements = string_operations::split(line, field_delim);

			// For all fields
			for (auto item : elements) {
				// Split to key and value
				auto parts = string_operations::split(item, kv_delim);

				// Act according to the key
				if (parts[0]=="dl_src") { // Layer 2
					rule.fields[0].low = rule.fields[0].high = parse_mac_address(parts[1]);
				} else if (parts[0]=="dl_dst") { // Layer 2
					rule.fields[1].low = rule.fields[1].high = parse_mac_address(parts[1]);
				} else if (parts[0]=="eth_type") { // Layer 2
					rule.fields[2].low = rule.fields[2].high = string_operations::hex2int(parts[1]);
				} else if (parts[0]=="in_port") { // Layer 1
					rule.fields[3].low = rule.fields[3].high = string_operations::str2int(parts[1]);
				} else if (parts[0]=="nw_src") { // Layer 3
					vector<uint32_t> vec = parse_ip_mask_address(parts[1]);
					rule.fields[4].low = vec[0];
					rule.fields[4].high = vec[1];
				} else if (parts[0]=="nw_dst") { // Layer 3
					vector<uint32_t> vec = parse_ip_mask_address(parts[1]);
					rule.fields[5].low = vec[0];
					rule.fields[5].high = vec[1];
				} else if (parts[0]=="nw_proto") { // Layer 3
					rule.fields[6].low = rule.fields[6].high = string_operations::str2int(parts[1]);
				} else if (parts[0]=="tp_dst") { // Layer 4
					rule.fields[7].low = rule.fields[7].high = string_operations::str2int(parts[1]);
				} else if (parts[0]=="tp_src") { // Layer 4
					rule.fields[8].low = rule.fields[8].high = string_operations::str2int(parts[1]);
				}
			}

			// Add rule to output
			output.push_back(rule);
		}

		return output;
	} catch (...) {
		throw error("cannot open Calssbench-ng file");
	}
}

/**
 * @brief Read a textual rule-set file and identify its origins.
 * @param filename Path to the rule-set file
 * @return The type of the rule-set
 */
ruleset_type_t classify_ruleset_file(const char* filename) {
	try {
		// Open file
		fstream fs;
		fs.open(filename, fstream::in);

		// Read first first line
		char line[2048];
		fs.getline(line, 2048);

		// Check if binary (version 1 or later)
		if (line[0] == 0) return BINARY;

		// The file contains text
		string first_line(line);

		// Check whether it is Classbench

		if (first_line[0] == '@') {
			return CLASSBENCH;
		}

		// Check whether it is Classbench-ng
		regex re("(\\s*\\w+=([^,]+),?)+");
		smatch match;
		if (regex_match(first_line, match, re)) {
			return CLASSBENCHNG;
		}

	} catch (exception& e) {
		warningf("Cannot open rule-set in path %s: %s", filename, e.what());
	}

	return UNKNOWN;
}

/**
 * @brief Returns the upper bound for a field of a rule-set
 * @param field_idx The index of the field within the rule-db
 * @param type The rule-set type
 */
uint64_t get_field_bound(int field_idx, ruleset_type_t type) {
	// The classic 5-tuple bounds
	if (type==CLASSBENCH) {
		static uint32_t bounds[5] = {0xffffffff, 0xffffffff, 0xffff, 0xffff, 0xff};
		if (field_idx > 4) {
			throw error("Classbench has at most five fields");
		}
		return bounds[field_idx];
	}
	// Classbench-NG with 9 OpenFlow fields
	else if (type==CLASSBENCHNG) {
		// dl_src, dl_dst, eth_type, in_port, nw_src, nw_dst, nw_proto, tp_dst, tp_src
		static uint32_t bounds[9] = {0xffffffff, 0xffffffff, 0xffff, 0xffffffff, 0xffffffff, 0xffffffff, 0xff, 0xffff, 0xffff};
		if (field_idx > 9) {
			throw error("Classbench-ng has at most nine fields");
		}
		return bounds[field_idx];
	}
	else if (type==BINARY) {
		// All fields are between 0 and 0xffffffff
		return 0xffffffff;
	}
	else {
		throw error("Unsupported rule-set type");
	}
}

/**
 * @brief Prints a rule-db to stdout
 * @param rule_db The rule database
 */
void print_rule_db(std::list<openflow_rule> rules) {
	uint32_t i = 0;
	for(auto rule : rules) {
		cout << i++ << ": ";
		for (auto field : rule.fields) {
			cout << field.low << "-" << field.high << "\t";
		}
		cout << "<p: " << rule.priority << ">" << endl;
	}
}


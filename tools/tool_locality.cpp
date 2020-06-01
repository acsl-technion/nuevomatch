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

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>

#include <argument_handler.h>
#include <logging.h>
#include <zipf.h>

using namespace std;

// Holds arguments information
static argument_t my_arguments[] = {
		// Name,			         Required,	IsBoolean,	Default,	Help
		{"--N",		      		0,			0,			"450000",	"Number of valid locality entries."},
		{"--percent",           0,			0,			"0.03",		"Measure the generated locality for X percent of the packets. in [0,1]"},

		/* Input mode */
		{"--zipf",	      		0,			1,			NULL,		"Use Zipf locality"},
		{"--from-file",    		0,			0,			NULL,		"Extract the locality according to the uniqueness of text lines."},
		{"--from-locality",		0,			0,			NULL,		"Read the locality from a text file with integers."},

		/* Zipf mode */
		{"--n",		      		0,			0,			"700000",	"Zipf, number of packets to generate"},
		{"--alpha",	      		0,			0,			"1",		"Zipf, parameter alpha"},
		
		{"--trace-file",	      0,			0,			NULL,		"Input trace-file to shuffle according to the locality"},

		/* Output mode */
		{"--output-trace",		0,			0,			NULL,		"(Output mode) Output to trace file. Filename."},
		{"--output-locality",	0,			0,			NULL,		"(Output mode) Output to locality text file. Filename."},

		{NULL,                  0,			0,			NULL,		"Applies locality to existing packet traces. The locality can be "
																"generated from Zipf distribution or extracted from text files "
																"(such as parsed PCAP files)." } /* Sentinel */
};


/**
 * @brief Prints progres to the screen
 */
void print_progress(int counter, const char* message, size_t size) {
	if ( (size ==0) || (counter < 0) ) {
		fprintf(stderr, "\r%s... Done   \n", message);
	} else {
		int checkpoint = size < 100 ? 1 : size/100;
		if (counter%checkpoint==0) {
			fprintf(stderr, "\r%s... (%u%%)", message, counter/checkpoint);
		}
	}
}

/**
 * @brief Generates a vector of 'n' integers with Zipf(N,alpha) distribution.
 */
vector<int> generate_zipf_locality() {

	int n = atoi(ARG("--n")->value);
	int N = atoi(ARG("--N")->value);
	double alpha = atof(ARG("--alpha")->value);

	messagef("Generating %d packets with Zipf locality with N=%d and alpha=%lf", N, alpha);
	vector<int> output;

	for (int i=0; i<n; ++i) {
		output.push_back(zipf(alpha, N));
		print_progress(i, "Generating zipf distribution", n);
	}

	print_progress(0, "Generating zipf distribution", 0);
	return output;
}

/**
 * @brief Reads an external file and extract its locality according
 * to the uniqueness of its lines.
 */
vector<int> extract_locality_from_file() {

	const char* filename = ARG("--from-file")->value;
	ifstream file_in (filename);

	if (!file_in.is_open()) {
		throw errorf("Cannot read file \"%s\"", filename);
	}

	messagef("Extracting the locality from \"%s\"...", filename);
	vector<int> output;

	map<string, int> flows;
	string current;
	int counter = 0;
	int N = atoi(ARG("--N")->value);
	
	// For generating uniform random
	random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> distrib(0, N-1);

	while (getline(file_in, current)) {
		if (flows.find(current) == flows.end()) {
			// In case the counter is in range
			if (counter < N) {
				flows[current] = counter++;
			}
			// Counter is not in range. Randomize value for flow.
			// (This models the case when various flows match a single
			// matching-rule).
			else {
				flows[current] = distrib(gen);
			}
		}
		output.push_back(flows[current]);
	}

	file_in.close();
	return output;
}

/**
 * @brief Reads the locality from a text file with 
 * integers, s.t each line is the locality.
 */
vector<int> read_locality_from_file() {

	const char* filename = ARG("--from-locality")->value;
	ifstream file_in (filename);

	if (!file_in.is_open()) {
		throw errorf("Cannot read file \"%s\"", filename);
	}

	messagef("Reading locality from \"%s\"...", filename);
	vector<int> output;

	string current;

	while (getline(file_in, current)) {
		output.push_back(atoi(current.c_str()));
	}

	file_in.close();
	return output;

}

/**
 * @brief Reads a trace file and returns a vector of lines
 */
vector<string> read_trace_file() {
	const char* filename = ARG("--trace-file")->value;
	messagef("Reading trace file \"%s\"...", filename);
	ifstream file_in (filename);

	if (!file_in.is_open()) {
		throw errorf("Cannot read file \"%s\"", filename);
	}

	vector<string> trace_lines;
	string current;
	while (getline(file_in, current)) {
		trace_lines.push_back(current);
	}

	// Shuffle trace lines. This is important for removing possible bias
	auto rng = std::default_random_engine {};
	shuffle(begin(trace_lines), end(trace_lines), rng);
	file_in.close();

	return trace_lines;
}

/**
 * @brief Analyzes the locality percent for X percent of the traffic
 */
void analyze_locality_percent(vector<int>& locality) {
	double traffic_percent = atof(ARG("--percent")->value);
	int N = atoi(ARG("--N")->value);
	int max_bound = N * traffic_percent;
	double counter = 0;
	for (auto x : locality) {
		if (x<=max_bound) {
			counter++;
		}
	}
	messagef("Locality: %.0lf%% most frequent flows hold %.0lf%% of the traffic (%d available flows, %ld traffic size)", 
		traffic_percent*100, counter/locality.size()*100, N, locality.size());
}

/**
 * @brief Writes the output trace file
 */
void write_output_trace_file(vector<int>& locality, vector<string>& trace_lines) {
	const char* filename = ARG("--output-trace")->value;
	messagef("Writing output trace to file \"%s\"...", filename);
	ofstream file_out (filename);
	for (size_t pos : locality) {
		file_out << trace_lines[pos] << endl;
	}
	file_out.close();
}

/**
 * @brief Writes the output locality file
 */
void write_output_locality_file(vector<int>& locality) {
	const char* filename = ARG("--output-locality")->value;
	messagef("Writing output locality to file \"%s\"...", filename);
	ofstream file_out (filename);
	for (size_t pos : locality) {
		file_out << pos << endl;
	}
	file_out.close();
}

int main(int argc, char** argv) {
	try {
		// Parse arguments
		parse_arguments(argc, argv, my_arguments);

		// Generate locality according to arguments
		vector<int> locality;
		if (ARG("--zipf")->available) {
			locality = generate_zipf_locality();
		} else if (ARG("--from-file")->available) {
			locality = extract_locality_from_file();
		} else if (ARG("--from-locality")->available) {
			locality = read_locality_from_file();
		} else {
			throw error("Input mode must be set by arguments.");
		}
		
		// Analyze the locality
		analyze_locality_percent(locality);
		
		// Do we write a locality output?
		if (ARG("--output-locality")->available) {
			// Write output locality file
			write_output_locality_file(locality);
		}
		// Do we write a trace output?
		else if (ARG("--output-trace")->available) {
			// Reads a trace file from argumens
			vector<string> trace_lines = read_trace_file();
	
			// Write the output file
			write_output_trace_file(locality, trace_lines);
		}
		// All other options are invalid
		else {
			throw error ("Output mode must be set by arguments.");
		}

	} catch (exception& e) {
		messagef("%s", e.what());
		return 1;
	}
}

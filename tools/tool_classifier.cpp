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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <regex.h>
#include <list>
#include <set>
#include <vector>
#include <sys/mman.h> // mmap

#include <logging.h>
#include <argument_handler.h>
#include <object_io.h>
#include <cut_split.h>
#include <efficuts.h>
#include <neurocuts.h>
#include <tuple_merge.h>
#include <nuevomatch.h>
#include <rule_db.h>
#include <nuevomatch_config.h>
#include <parallel_classifier.h>
#include <string_operations.h>
#include <em_table.h>

// Internal methods
list<openflow_rule> read_rule_db();

// Holds arguments information
static argument_t my_arguments[] = {
		// Name,						Required,	IsBoolean,	Default,	Help
		{"-m",							1,			0,			NULL,		"Set classifier type. Valid options are: [neurocuts, cutsplit, efficuts, tuplemerge, nuevomatch]."},
		{"-in",							0,			0,			NULL,		"Input file. Meaning changes across modes."},

		/* Create Mode */
		{"-c",							0,			1,			NULL,		"(Create Mode) Generate new classifier."},
		{"-out",						0,			0,			NULL,		"(Create Mode) Name of output packed classifier."},
		{"--indices",					0,			0,			NULL,		"(Create Mode) Override rule indices when creating classifier."},

		/* Load Mode */
		{"-l",							0,			1,			NULL,		"(Load Mode) Load packed classifier from file"},

		/* Print Mode */
		{"-v",							0,			0,			"0",		"Print classifier properties. Set verbosity=[0,1,2,3,4,5,6]."},

		/* CutSplit / EfiiCuts Modes */
		{"--binth",						0,			0,			"8",		"(EffiCuts / CutSplit Mode) Binth value."},
		{"--threshold",					0,			0,			"24",		"(CutSplit Mode) Threshold value."},

		/* TupleMerge Mode */
		{"--collision-limit",			0,			0,			"10",		"(TupleMerge) Collision-limit value"},

		/* NuevoMatch Mode */
		{"--disable-isets",				0,			1,			NULL,		"(NuevoMatch Mode) Disable iSets when performing classification."},
		{"--disable-remainder",			0,			1,			NULL,		"(NuevoMatch Mode) Disable remainder when performing classification."},
		{"--disable-bin-search",		0,			1,			NULL,		"(NuevoMatch Mode) Disable binary search in all iSets."},
		{"--disable-validation",		0,			1,			NULL,		"(NuevoMatch Mode) Disable validation phase in all iSets."},
		{"--disable-classification",	0,			1,			NULL,		"(NuevoMatch Mode) Disable all classification. "
																			"Used to measure performance of batching and publishing."},
		{"--max-subsets",				0,			0,			"-1",		"(NuevoMatch Mode) Maximum iSets for NuevoMatch. Note: iSet limitation "
																			"causes remainder classifier generation on the fly."},
		{"--start-from-iset",			0,			0,			"0",		"(NuevoMatch Mode) The index of iSet to start NuevoMatch with" },
		{"--arbitrary-fields",			0,			0,			NULL,		"(NuevoMatch Mode) Run NuevoMatch on any subset of fields."
																			"Usage: --arbitrary-fields \"1,2,6\""},
		{"--arbitrary-core-allocation",	0,			0,			"",			"(NuevoMatch Mode) Manually set the subset-core allocation."
																					"Example: use \"0=1,3;2=4,5\" to allocate subsets {1,3}"
																					"to core 0, and subset {4,5} to core 2. Remainder subset is always 0."},

		{"--force-remainder-build",		0,			1,			NULL,		"(NuevoMatch Mode) Force rebuilding the remainder classifier on the fly."},
		{"--remainder-type",			0,			0,			"cutsplit",	"(NuevoMatch Mode) Set NuevoMatch remainder classifier."
																			"Valid options are: [cutsplit, efficuts, nuerocuts, tuplemerge]."
																			"Note: flags for CutSplit / EffiCuts are valid for the remainder classifier as well."},
		{"--external-remainder",		0,			0,			NULL,		"(NuevoMatch Mode) Load the remainder classifier from a file."
																			"Note: requires the '--remainder-type' flag."
																			"Usage: --external-remainder FILENAME"},

		/* Caching */
		{"--cache",						0,			1,			NULL,		"Add small & simple exact-match cache before classifier."},

		/* Parallel Mode */
		{"--parallel",					0,			0,			"1",		"(Parallel Mode) Start any classifier with X parallel threads. "},
		{"--queue-size",				0,			0,			"256",		"(Parallel Mode) Inter-core messages queue size."},
		{"--batch-size",				0,			0,			"128",		"(Parallel Mode) Packet batch size."},

		/* Trace benchmark */
		{"--trace",						0,			0,			NULL,		"(Trace Mode) Activate trace mode. Set trace filename."},
		{"--trace-from",				0,			0,			"0",		"(Trace Mode) Limit the number of packets in trace. Start packet"},
		{"--trace-to",					0,			0,			"-1",		"(Trace Mode) Limit the number of packets in trace. End packet"},
		{"--trace-no-warm",				0,			1,			NULL,		"(Trace Mode) Disable cache warm for trace."},
		{"--trace-silent",				0,			1,			NULL,		"(Trace Mode) Do not print classification errors to screen"},
		{"--trace-repeat",				0,			0,			"1",		"(Trace Mode) Repeat the experiment multiple times"},
		{"--trace-fail-fast",			0,			1,			NULL,		"(Trace Mode) Fail on classification error"},

		{NULL,							0,			0,			NULL,		"Classifier generation and benchmark tool."} /* Sentinel */
};

/**
 * @brief Holds the trace packets to test
 */
trace_packet* trace_packets;
uint32_t start_packet, end_packet;
volatile bool fail_fast;

/**
 * @brief Exact match table to help with skewed traffic
 */
ExactMatchTable* em_table;

class BenchmarkListener : public GenericClassifierListener {
public:

	volatile uint32_t num_of_results;
	bool silent;
	BenchmarkListener(bool silent) : num_of_results(0), silent(silent) {};

	/**
	 * @brief Is invoked by the classifier when new result is available
	 * @param id A unique packet id
	 * @param priority The priority of the matching rule
	 * @param action The action to take on the packet
	 * @param args Any additional information
	 */
	virtual void on_new_result(unsigned int id, int priority, int action, void* args) {
		// Skip invalid results
		if (id != 0xffffffff) {
			// Skip invalid packets
			if (num_of_results < (end_packet-start_packet)) {
				// Cache packet
				em_table->add(trace_packets[start_packet+id], priority);
				// Check result match trace
				if (!silent && (uint32_t)action != trace_packets[start_packet+id].match_priority) {
					warningf("packet %u does not match!. Got: %u, expected: %u",
							id, action, trace_packets[start_packet+id].match_priority);
					// In case of fail fast argument was set, throw an exception
					if (fail_fast) {
						throw error("Classification error");
					}
				}
			}
		}
		++num_of_results;
	}
};

/**
 * @brief Work in CutSplit mode
 */
CutSplit* mode_cutsplit() {

	uint32_t binth = 	 atoi( get_argument_by_name(my_arguments,"--binth")->value );
	uint32_t threshold = atoi( get_argument_by_name(my_arguments,"--threshold")->value );

	// Modes
	bool mod_generate = ARG("-c")->available;
	bool mod_load = ARG("-l")->available;

	// Inputs and outputs
	argument_t* input_arg = ARG("-in");
	argument_t* output_arg = ARG("-out");

	CutSplit* output = nullptr;

	// In case of mode generate
	if (mod_generate) {
		messagef("Generating new CutSplit classifier. Input: Classbench rule database. Output: classifier file");

		if (!input_arg->available || !output_arg->available) {
			throw error("Cannot generate classifier: -in and -out arguments are required");
		}

		// Load the rule-db
		list<openflow_rule> rule_db = read_rule_db();

		// Generate CutSplit classifier
		messagef("Generating CutSplit classifier...");
		output = new CutSplit(threshold, binth);
		output->build(rule_db);

		// Pack the classifier to bytes
		messagef("Packing classifier to %s...", output_arg->value);
		uint8_t *buffer;
		uint32_t size;
		output->pack().pack(&buffer, &size);

		// Write the output file
		FILE* file = fopen(output_arg->value, "w");
		if (!file) {
			throw error("cannot open output filename");
		}
		fwrite(buffer, sizeof(uint8_t), size, file);
		fclose(file);

		// Load the classifier from buffer, should be faster
		free(output);
		output = new CutSplit(threshold, binth);
		ObjectReader classifier_handler(buffer, size);
		output->load(classifier_handler);

		messagef("Done.");
	}
	// Mode load
	else if (mod_load) {
		messagef("Loading CutSplit classifier from file. Input: classifier filename");

		if (!input_arg->available) {
			throw error("-in Argument is required with classifier filename");
		}

		// Read classifier file to memory
		ObjectReader classifier_handler(input_arg->value);

		// Load the classifier
		output = new CutSplit(threshold, binth);
		output->load(classifier_handler);
	}

	return output;
}

/**
 * @brief Work in TupleMerge mode
 */
TupleMerge* mode_tuplemerge() {

	// Modes
	bool mod_generate = ARG("-c")->available;
	bool mod_load = ARG("-l")->available;

	// Inputs and outputs
	argument_t* input_arg = ARG("-in");
	argument_t* output_arg = ARG("-out");

	TupleMerge* output = nullptr;

	// In case of mode generate
	if (mod_generate) {
		messagef("Generating new TupleMerge classifier. Input: Classbench rule database. Output: classifier file");

		if (!input_arg->available || !output_arg->available) {
			throw error("Cannot generate classifier: -in and -out arguments are required");
		}

		// Load the rule-db
		list<openflow_rule> rule_db = read_rule_db();

		// Generate CutSplit classifier
		messagef("Generating TupleMerge classifier...");
		uint32_t limit = atoi(ARG("--collision-limit")->value);
		output = new TupleMerge(limit);
		output->build(rule_db);

		// Pack the classifier to bytes
		messagef("Packing classifier to %s...", output_arg->value);
		uint8_t *buffer;
		uint32_t size;
		output->pack().pack(&buffer, &size);

		// Write the output file
		FILE* file = fopen(output_arg->value, "w");
		if (!file) {
			throw error("cannot open output filename");
		}
		fwrite(buffer, sizeof(uint8_t), size, file);
		fclose(file);

		// Load the classifier from buffer, should be faster
		free(output);
		output = new TupleMerge();
		ObjectReader classifier_handler(buffer, size);
		output->load(classifier_handler);

		messagef("Done.");
	}
	// Mode load
	else if (mod_load) {
		messagef("Loading TupleMerge classifier from file. Input: classifier filename");

		if (!input_arg->available) {
			throw error("-in Argument is required with classifier filename");
		}

		// Read classifier file to memory
		ObjectReader classifier_handler(input_arg->value);

		// Load the classifier
		output = new TupleMerge();
		output->load(classifier_handler);
	}

	return output;
}

/**
 * @brief Work in EffiCuts mode
 */
EffiCuts* mode_efficuts() {

	uint32_t binth = atoi( get_argument_by_name(my_arguments,"--binth")->value );

	// Modes
	bool mod_generate = ARG("-c")->available;
	bool mod_load = ARG("-l")->available;

	// Inputs and outputs
	argument_t* input_arg = ARG("-in");
	argument_t* output_arg = ARG("-out");

	EffiCuts* output = nullptr;

	if (mod_generate) {
			messagef("Generating new EffiCuts classifier. Input: Classbench rule database. Output: classifier file");

			if (!input_arg->available || !output_arg->available) {
				throw error("Cannot generate classifier: -in and -out arguments are required");
			}

			// Load the rule-db
			list<openflow_rule> rule_db = read_rule_db();

			// Generate CutSplit classifier
			messagef("Generating EffiCuts classifier...");
			output = new EffiCuts(binth);
			output->build(rule_db);

			// Pack the classifier to bytes
			messagef("Packing classifier to %s...", output_arg->value);
			uint8_t *buffer;
			uint32_t size;
			output->pack().pack(&buffer, &size);

			// Write the output file
			FILE* file = fopen(output_arg->value, "w");
			if (!file) {
				throw error("cannot open output filename");
			}
			fwrite(buffer, sizeof(uint8_t), size, file);
			fclose(file);

			messagef("Done.");
	}
	// Mode load
	else if (mod_load) {
		messagef("Loading EffiCuts classifier from file. Input: classifier filename");

		if (!input_arg->available) {
			throw error("-in Argument is required with classifier filename");
		}

		// Read classifier file to memory
		ObjectReader classifier_handler(input_arg->value);

		// Load the classifier
		output = new EffiCuts(binth);
		output->load(classifier_handler);
	}

	return output;
}


/**
 * @brief Work in NeuroCuts mode
 */
NeuroCuts* mode_neurocuts() {

	// Modes
	bool mod_generate = ARG("-c")->available;
	bool mod_load = ARG("-l")->available;

	// Inputs and outputs
	argument_t* input_arg = ARG("-in");

	NeuroCuts* output = nullptr;

	if (mod_generate) {
			throw error("NeuoCuts mode is not supported with Create mode. Use Python script to convert from NeuroCuts output to loadable classifier");
	}
	// Mode load
	else if (mod_load) {
		messagef("Loading NeuroCuts classifier from file. Input: classifier filename");

		if (!input_arg->available) {
			throw error("-in Argument is required with classifier filename");
		}

		// Read classifier file to memory
		ObjectReader classifier_handler(input_arg->value);

		// Load the classifier
		output = new NeuroCuts();
		output->load(classifier_handler);
	}

	return output;
}

/**
 * @brief Work in nuevomatch mode
 */
GenericClassifier* mode_nuevomatch() {

	// Modes
	bool mod_generate = ARG("-c")->available;

	// Inputs and outputs
	argument_t* input_arg = ARG("-in");
	argument_t* output_arg = ARG("-out");

	// Create new NuevoMatch object
	GenericClassifier* output;

	// Set configuration for NuevoMatch
	NuevoMatchConfig config;
	config.queue_size = atoi( get_argument_by_name(my_arguments,"--queue-size")->value );
	config.num_of_cores = MAX(1, atoi( ARG("--parallel")->value ));
	config.max_subsets = atoi( ARG("--max-subsets")->value );
	config.start_from_iset = atoi( ARG("--start-from-iset")->value );
	config.disable_isets = ARG("--disable-isets")->available;
	config.disable_remainder = ARG("--disable-remainder")->available;
	config.disable_bin_search = ARG("--disable-bin-search")->available;
	config.disable_validation_phase = ARG("--disable-validation")->available;
	config.disable_all_classification = ARG("--disable-classification")->available;
	config.arbitrary_subset_clore_allocation = ARG("--arbitrary-core-allocation")->value;
	config.force_rebuilding_remainder = ARG("--force-remainder-build")->available;

	// Arbitrary field argument
	if (ARG("--arbitrary-fields")->available) {
		static regex re(",");
		config.arbitrary_fields = string_operations::split(
				ARG("--arbitrary-fields")->value, re, string_operations::str2int);
	}

	// Read configuration for the remainder classifier
	uint32_t binth = 	 atoi( get_argument_by_name(my_arguments,"--binth")->value );
	uint32_t threshold = atoi( get_argument_by_name(my_arguments,"--threshold")->value );

	// Get the remainder type. Default is CutSplit
	const char* remainder_type = ARG("--remainder-type")->value;
	config.remainder_type = remainder_type;

	// Build new remainder classifier according to type
	if (!strcmp(remainder_type, "cutsplit")) {
		config.remainder_classifier = new CutSplit(binth, threshold);
	} else if (!strcmp(remainder_type, "efficuts")) {
		config.remainder_classifier = new EffiCuts(binth);
	} else if (!strcmp(remainder_type, "neurocuts")) {
		config.remainder_classifier = new NeuroCuts();
	} else if (!strcmp(remainder_type, "tuplemerge")) {
		config.remainder_classifier = new TupleMerge();
		 // Note: this is due to licb critical error when trying to build base on previous data and failing!
		config.force_rebuilding_remainder = true;
	} else {
		throw errorf("Remainder classifier type is not valid. Got '%s'.", remainder_type);
	}

	// In case the remainder classifier should be loaded from an external file
	if (ARG("--external-remainder")->available) {
		const char* external_remainder = ARG("--external-remainder")->value;
		loggerf("Setting remainder classifier as external, loading from file %s...", external_remainder);
		ObjectReader reader(external_remainder);
		// Load the remainder from file
		config.remainder_classifier->load(reader);
		config.external_remainder = true;
	}

	// Set the batch size
	switch (atoi( ARG("--batch-size")->value )) {
	case 512:
		output = new NuevoMatch<512>(config);
		break;
	case 256:
		output = new NuevoMatch<256>(config);
		break;
	case 128:
		output = new NuevoMatch<128>(config);
		break;
	case 64:
		output = new NuevoMatch<64>(config);
		break;
	case 32:
		output = new NuevoMatch<32>(config);
		break;
	case 16:
		output = new NuevoMatch<16>(config);
		break;
	case 8:
		output = new NuevoMatch<8>(config);
		break;
	case 4:
		output = new NuevoMatch<4>(config);
		break;
	case 2:
		output = new NuevoMatch<2>(config);
		break;
	default:
		throw error("--batch-size is valid only with values [2, 4, 8, 16, 32, 64, 128, 256, 512].");
	}

	if (!input_arg->available) {
		throw error("-in Argument is required with classifier filename");
	}

	// Read classifier file to memory
	ObjectReader classifier_handler(input_arg->value);


	messagef("Loading nuevomatch with batch size of %s...", ARG("--batch-size")->value);

	// Load nuevomatch
	// This will work for both classifiers without remainder classifier set
	// and classifiers with remainder classifiers set
	output->load(classifier_handler);

	// In case of generate object, write the output file
	if (mod_generate) {
		if (!output_arg->available) {
			throw error("-out Argument is required with output filename");
		}

		uint8_t* buffer;
		uint32_t size;
		output->pack().pack(&buffer, &size);;

		// Write the output file
		FILE* file = fopen(output_arg->value, "w");
		if (!file) {
			throw error("cannot open output filename");
		}
		fwrite(buffer, sizeof(uint8_t), size, file);
		fclose(file);
	}

	return output;
}


/**
 * @brief Main entry point
 */
int main(int argc, char** argv) {
 try{
 	// Parse arguments
 	parse_arguments(argc, argv, my_arguments);
 
 	// Get working mode
 	const char* mode = ARG("-m")->value;
 
 	// Get classifier according to mode
 	GenericClassifier* classifier = nullptr;
 	bool nuevomatch_enabled = false;
 
 	if (strcmp(mode, "cutsplit") == 0) {
 		classifier = mode_cutsplit();
 	} else if (strcmp(mode, "nuevomatch") == 0) {
 		classifier = mode_nuevomatch();
 		nuevomatch_enabled = true;
 	} else if (strcmp(mode, "efficuts") == 0) {
 		classifier = mode_efficuts();
 	} else if (strcmp(mode, "neurocuts") == 0) {
 		classifier = mode_neurocuts();
 	} else if (strcmp(mode, "tuplemerge") == 0) {
 		classifier = mode_tuplemerge();
 	} else {
 		throw error("mode is invalid");
 	}
 
 	// Print classifier attributes
 	if (classifier == nullptr) {
 		throw error("classifier was not initialized");
 	}
 
 	messagef("Classifier attributes:");
 	messagef("Total rules: %u", classifier->get_num_of_rules());
 	messagef("Total size (bytes): %u", classifier->get_size());
 	messagef("Build time (ms): %u", classifier->get_build_time());
 
 	// In case of parallel classifier
 	if (!nuevomatch_enabled && ARG("--parallel")->available) {
 
 		// One core for pushing, other cores for processing
 		uint32_t num_of_classifiers = atoi(ARG("--parallel")->value);
 		GenericClassifier** classifiers = new GenericClassifier*[num_of_classifiers];
 
 		// Clone classifiers
 		classifiers[0] = classifier;
 		for (uint32_t i=1; i<num_of_classifiers; ++i) {
 			classifiers[i] = classifier->clone();
 		}
 
 		// Get arguments
 		uint32_t queue_size = atoi(ARG("--queue-size")->value);
 
 		// Set the batch size
 		switch (atoi( ARG("--batch-size")->value )) {
 		case 512:
 			classifier = new ParallelClassifier<512>(queue_size, num_of_classifiers, classifiers);
 			break;
 		case 256:
 			classifier = new ParallelClassifier<256>(queue_size, num_of_classifiers, classifiers);
 			break;
 		case 128:
 			classifier = new ParallelClassifier<128>(queue_size, num_of_classifiers, classifiers);
 			break;
 		case 64:
 			classifier = new ParallelClassifier<64>(queue_size, num_of_classifiers, classifiers);
 			break;
 		case 32:
 			classifier = new ParallelClassifier<32>(queue_size, num_of_classifiers, classifiers);
 			break;
 		case 16:
 			classifier = new ParallelClassifier<16>(queue_size, num_of_classifiers, classifiers);
 			break;
 		case 8:
 			classifier = new ParallelClassifier<8>(queue_size, num_of_classifiers, classifiers);
 			break;
 		case 4:
 			classifier = new ParallelClassifier<4>(queue_size, num_of_classifiers, classifiers);
 			break;
 		case 2:
 			classifier = new ParallelClassifier<2>(queue_size, num_of_classifiers, classifiers);
 			break;
 		default:
 			throw error("--batch-size is valid only with values [2, 4, 8, 16, 32, 64, 128, 256, 512].");
 		}
 	}
 
 	// Register a listener for classifier
 	BenchmarkListener listener( ARG("--trace-silent")->available );
 	classifier->add_listener(listener);
 
 	// In case of mode trace - initialize trace engine
 	bool mod_trace = ARG("--trace")->available;
 	if (mod_trace) {
  	
	 	// Create exact-match table 
	 	em_table = new ExactMatchTable(8192, !ARG("--cache")->available);

 		// Fail fast?
 		fail_fast = ARG("--trace-fail-fast")->available;
 
 		// Arbitrary field argument
 		vector<uint32_t> arbitrary_fields;
 		if (ARG("--arbitrary-fields")->available) {
 			static regex re(",");
 			arbitrary_fields = string_operations::split(
 					ARG("--arbitrary-fields")->value, re, string_operations::str2int);
 		}
 
 		// Read the textual trace file
 		messagef("Reading trace file...");
 		const char* trace_filename = ARG("--trace")->value;
 		uint32_t num_of_packets;
 		trace_packets = read_trace_file(trace_filename, arbitrary_fields, &num_of_packets);
 		if (!trace_packets) {
 			throw error("error while reading trace file");
 		}
 		messagef("Total %u packets in trace", num_of_packets);
 
 		// Limit the number of packets
 		start_packet = atoi( get_argument_by_name(my_arguments,"--trace-from")->value );
 		end_packet   = atoi( get_argument_by_name(my_arguments,"--trace-to")->value );
 		if (end_packet > num_of_packets) end_packet = num_of_packets;
 
 		// Warm cache
 		uint32_t warm_repetitions = 5;
 		if (ARG("--trace-no-warm")->available) {
 			messagef("Skipping cache warm");
 			warm_repetitions = 0;
 		} else {
 			messagef("Warming cache...");
 		}
 		for (uint32_t r=0; r<warm_repetitions; ++r) {
 			messagef("Iteration %u...", r);
 			for (uint32_t i=start_packet; i<end_packet; ++i) {
 				classifier->classify_async(trace_packets[i].get(), -1);
 			}
 			// Request to process remaining packets
 			classifier->classify_async(nullptr, -1);
 			while(listener.num_of_results < (end_packet-start_packet));
 			// Reset counters
 			listener.num_of_results=0;
 			classifier->reset_counters();
 			em_table->invalidate();
 		}
 	}

	// Chace hit rate
	double cache_hit=0;

 	// Perform the experiment, repeat X times
 	uint32_t time_to_repeat = atoi( ARG("--trace-repeat")->value );
 	messagef("Repeating experiment %u times", time_to_repeat);
 
 	for (uint32_t i=0; i<time_to_repeat; ++i) {
 		if (mod_trace) {
 			messagef("Starting trace test for classifier with %u packets...", (end_packet-start_packet));
 
 			// Reset counters
 			classifier->reset_counters();
 			listener.num_of_results=0;
 			em_table->invalidate();

 			classifier->start_performance_measurement();
 			// Run the lookup
 			for (uint32_t i=start_packet; i<end_packet; ++i) {
 				// Check cache for hit
 				int hit_priority = em_table->lookup(trace_packets[i]);
 				if (hit_priority != -1) {
					listener.on_new_result(i-start_packet, hit_priority, hit_priority, nullptr);
					classifier->advance_counter();
					++cache_hit;
					continue;
 				}
 				// On cache miss
 				classifier->classify_async(trace_packets[i].get(), -1);
 			}
 
 			// Request to process remaining packets
 			classifier->classify_async(nullptr, -1);
 
 			// Wait for results
 			while(listener.num_of_results < (end_packet-start_packet));
 			classifier->stop_performance_measurement();
 		}
 
 		bool mod_print = ARG("-v")->available;
 		if (mod_print) {
			messagef("Cache hit rate: %.2f, utilization: %.2f",
				cache_hit/(end_packet-start_packet), em_table->utilization());
 			messagef("Classifier Information:");
 			classifier->print( atoi(ARG("-v")->value) );
 		}
 	}
 
 	delete classifier;
 	messagef("done.");
  return 0;
 } catch (std::exception& e) {
  messagef("%s", e.what());
  return 1;
 }
}

/**
 * @brief Reads the input rule-set file and return a list of OpenFlow rules.
 */
list<openflow_rule> read_rule_db() {
	list<openflow_rule> rule_db;
	const char* filename = ARG("-in")->value;

	// Which type is the rule-db?
	ruleset_type_t ruleset_type = classify_ruleset_file(filename);
	if (ruleset_type == UNKNOWN) {
		throw error("Cannot parse rule-set file: file format not recognized");
	} else if (ruleset_type == CLASSBENCH) {
		messagef("Recognized rule-set as Classbench text file");
		rule_db = read_classbench_file(filename);
	} else if (ruleset_type == CLASSBENCHNG) {
		messagef("Recognized rule-set as Classbench-ng text file");
		rule_db = read_classbench_ng_file(filename);
	}
	// Check size of rule-set
	if (rule_db.size() == 0) {
		throw error("Rule-set has zero rules");
	}

	// Should the priorities be overridden?
	filename = ARG("--indices")->value;
	if (filename != NULL) {
		ObjectReader reader(filename);
		set<uint32_t> indices = read_indices_file(reader);
		if (indices.size() != rule_db.size()) {
			throw error("Number of indices in indices-file mismatch number of rules");
		}
		list<openflow_rule>::iterator rule_it = rule_db.begin();
		set<uint32_t>::iterator indices_it = indices.begin();
		for (uint32_t i=0; i<rule_db.size(); ++i) {
			rule_it->priority = *indices_it;
			rule_it++;
			indices_it++;
		}
	}

	return rule_db;
}

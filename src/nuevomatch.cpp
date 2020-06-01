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

#include <set>
#include <string>
#include <vector>
#include <list>
#include <queue>
#include <algorithm>
#include <string.h>

#include <object_io.h>
#include <cpu_core_tools.h>

#include <array_operations.h>
#include <string_operations.h>
#include <nuevomatch.h>

/**
 * @brief Initiate a new NuevoMatch instance.
 * @param config The configuration for NuevoMatch
 */
template <uint32_t N>
NuevoMatch<N>::NuevoMatch(NuevoMatchConfig config) :
	_configuration(config), _isets(nullptr),
	_worker_serial(nullptr), _workers_parallel(nullptr),
	_last_iset_idx(0),
	_num_of_isets(0),_num_of_rules(0),
	_size(0), _build_time(0),
	_pack_buffer(nullptr), _pack_size(0),
	_next_batch_items(0), _batch_counter(0), _reducer(nullptr) { };

template <uint32_t N>
NuevoMatch<N>::~NuevoMatch() {
	for (uint32_t i=0; i<_configuration.num_of_cores-1; ++i) {
		delete _workers_parallel[i];
	}
	delete[] _workers_parallel;
	delete _worker_serial;
}

/**
 * @brief Creates this from a memory location
 * @param object An object-reader instance
 */
template <uint32_t N>
void NuevoMatch<N>::load(ObjectReader& reader) {

	// Used for packing
	_pack_buffer = new uint8_t[reader.size()];
	_pack_size = reader.size();
	memcpy(_pack_buffer, reader.buffer(), reader.size());

	// Read static information
	reader 	>> _num_of_isets >> _num_of_rules
			>> _size >> _build_time;

	// The size is measured by the iSets, and not by what was packed
	// Reason: support dynamic size for dynamic iSets
	_size = 0;

	// Show general information
	if (_configuration.disable_bin_search) {
		loggerf("Disabling binary search in all iSets");
	}
	if (_configuration.disable_remainder) {
		loggerf("Disabling remainder classifier");
	}
	if (_configuration.disable_validation_phase) {
		loggerf("Disabling validation phase in all iSets");
	}
	if (_configuration.disable_all_classification) {
		loggerf("Disabling classification");
	}

	// Check configuration error
	if (!_configuration.disable_remainder && !_configuration.remainder_classifier) {
		throw error("Remainder classifier is enabled but is not set");
	}

	// Load all subsets from file
	load_subsets(reader);

	// Load the remainder classifier
	load_remainder(reader);

	// Group subsets to groups, initialize workers
	group_subsets_to_cores();

	// Initialize the reducer
	_reducer = new reducer_job_t[_configuration.queue_size];
	for (uint32_t i=0; i<_configuration.queue_size; ++i) {
		_reducer[i].lock = 0;
	}
}

/**
 * @brief Packs this to byte array
 * @returns An object-packer with the binary data
 */
template <uint32_t N>
ObjectPacker NuevoMatch<N>::pack() const {

	// Pack the remainder classifier
	ObjectPacker remainder_packer = _configuration.remainder_classifier->pack();

	// Pack this
	ObjectPacker output;
	output.push(_pack_buffer, _pack_size);
	output << remainder_packer;

	return output;
}

/**
 * @brief Resets the all classifier counters
 */
template <uint32_t N>
void NuevoMatch<N>::reset_counters() {
	GenericClassifier::reset_counters();
	_batch_counter = 0;
	_next_batch_items = 0;
}

/**
 * @brief Advance the packet counter. Should be used when skipping 
 * classification of packets, such as with caches.
 */
template <uint32_t N>
void NuevoMatch<N>::advance_counter() { 
	_packet_counter++;
}

/**
 * @brief Start an asynchronous process of classification for an input packet.
 * @param header An array of 32bit integers according to the number of supported fields.
 * @returns A unique id for the packet
 */
template <uint32_t N>
uint32_t NuevoMatch<N>::classify_async(const uint32_t* header, int priority) {

	// Build next batch
	if (header != NULL) {
		uint32_t batch_modulo = _batch_counter & (_configuration.queue_size - 1);
		_next_batch[_next_batch_items] = header;
		_reducer[batch_modulo].packet_id[_next_batch_items] = _packet_counter;
		_next_batch_items++;
	}

	if (header == NULL || _next_batch_items == N) {
		process_batch();
	}

	return _packet_counter++;
}

/**
 * @brief Process a new batch of packets.
 */
template <uint32_t N>
void NuevoMatch<N>::process_batch() {

	// Reset reducer of batch
	uint32_t batch_modulo = _batch_counter & (_configuration.queue_size - 1);

	// Invalidate remaining packets in batch
	for (uint32_t i=_next_batch_items; i<N; ++i) {
		_next_batch[i] = nullptr;
		_reducer[batch_modulo].packet_id[i] = -1;
	}
		
	_reducer[batch_modulo].counter = 0;
	_reducer[batch_modulo].valid_items = _next_batch_items;
	
	// Produce next batch in all parallel workers
	for (uint32_t i=1; i<_configuration.num_of_cores; ++i) {
	 while(!_workers_parallel[i-1]->classify(batch_modulo, _next_batch));
	}
	
	// Do serial work
	while(!_worker_serial->classify(batch_modulo, _next_batch));
	
	// Produce the reducer
	infof("Produced batch with %u packets starting from id %u" ,_next_batch_items, packet_id);
	
	// Update counters
	++_batch_counter;
	_next_batch_items = 0;
}

/**
 * @brief Callback. Invoked by the iSet on result
 * @param info The batch information generated by the iSet
 * @param iset_index The iSet index
 * @param batch_id A unique id for the batch
 */
template <uint32_t N>
void NuevoMatch<N>::on_new_result(WorkBatch<classifier_output_t, N> info, uint32_t iset_index, uint32_t batch_id) {

	infof("subset %u trying to acquire lock for batch %u", iset_index, batch_id);

	// Get the reduce instance of current batch
	reducer_job_t* reduce = &this->_reducer[batch_id];

	// Lock on reduce
	while (__sync_val_compare_and_swap(&reduce->lock, 0, 1));

	infof("subset %u acquired lock for batch %u", iset_index, batch_id);

	// Update priority per output
	for (uint32_t i=0; i<N; ++i) {
		// Update reducer result in case the current result is better than the last
		// Flags that relevant for overriding result
		bool first_result  = (reduce->counter == 0);
		bool result_better = (uint32_t)info[i].priority < (uint32_t)reduce->results[i].priority;
		// Override result
		if (first_result || result_better) {
			reduce->results[i] = info[i];
		}
	}

	// Update counter
	++reduce->counter;

	// Publish in case of last output
	if (reduce->counter == _configuration.num_of_cores) {
		infof("subset %u publish batch %u", iset_index, batch_id);

		for (uint32_t i=0; i<N; ++i) {
			// Do not publish invalid, virtual packets
			if (i < reduce->valid_items) {
				// Publish current result
				for (auto it : _listeners) {
					it->on_new_result(
							reduce->packet_id[i],
							reduce->results[i].priority,
							reduce->results[i].action,
							_additional_args);
				}
			}	
		}
	}

	// Release lock
	reduce->lock = 0;

	infof("subset %u released lock for batch %u", iset_index, batch_id);
}

/**
 * @brief Starts the performance measurement of this
 */
template <uint32_t N>
void NuevoMatch<N>::start_performance_measurement() {
	for (uint32_t i=1; i<_configuration.num_of_cores; ++i) {
		_workers_parallel[i-1]->start_performance_measurements();
	}
	_worker_serial->start_performance_measurements();
	clock_gettime(CLOCK_MONOTONIC, &start_time);
}

/**
 * @brief Stops the performance measurement of this
 */
template <uint32_t N>
void NuevoMatch<N>::stop_performance_measurement() {
	clock_gettime(CLOCK_MONOTONIC, &end_time);
	for (uint32_t i=1; i<_configuration.num_of_cores; ++i) {
		_workers_parallel[i-1]->stop_performance_measurements();
	}
	_worker_serial->stop_performance_measurements();
}


/**
 * @brief Prints statistical information
 * @param verbose Set the verbosity level of printing
 */
template <uint32_t N>
void NuevoMatch<N>::print(uint32_t verbose) const {

	// High verbosity
	if (verbose > 2) {

		// Print the errors of all RQRMI
		for (uint32_t i=_configuration.start_from_iset; i<_last_iset_idx; ++i) {
			SimpleLogger::get() << "Error list for iSet " << i << ": [";
			auto&& error_list = _isets[i]->get_error_list();
			bool first = true;
			for (auto it : error_list) {
				if (!first) SimpleLogger::get() << ", ";
				SimpleLogger::get() << it;
				first = false;
			}
			SimpleLogger::get() << "]" << SimpleLogger::endl();
		}
		// Print expected errors
		for (uint32_t i=_configuration.start_from_iset; i<_last_iset_idx; ++i) {
			message_s("Expected error for iSet " << i << ": "
				 << _isets[i]->get_expected_error());
		}
	}

	// Measure performance
	double total_usec = (double)((end_time.tv_sec * 1e9 + end_time.tv_nsec) -
						  (start_time.tv_sec * 1e9 + start_time.tv_nsec)) / 1e3;

	messagef("Performance: total time %.3lf usec. Average time: %.3lf usec per packet.",
			total_usec, total_usec / _packet_counter);

	// Medium verbosity
	if (verbose > 1) {
		messagef("Serial worker 0 total time: %.3lf used, avg time per batch: %.3lf usec, publish time: %.3f us",
				_worker_serial->get_work_time(),
				_worker_serial->get_work_time() / _batch_counter,
				_worker_serial->get_publish_time());

		for (uint32_t i=1; i<_configuration.num_of_cores; ++i) {
			messagef("Parallel worker %u statistics: utilization: %.2lf%%, throughput: %.2lf rpus, "
					"backpressure: %.2lf rpus, avg time per batch: %.3lf us, publish time: %.3f us", i,
					_workers_parallel[i-1]->get_utilization(), _workers_parallel[i-1]->get_throughput(),
					_workers_parallel[i-1]->get_backpressure(), _workers_parallel[i-1]->get_average_work_time(),
					_workers_parallel[i-1]->get_publish_time());
		}

		if (!_configuration.disable_remainder) {
			messagef("Remainder classifier total size: %u bytes", this->_configuration.remainder_classifier->get_size());
		}
	}

	// Max verbosity
	if (verbose > 3 && !_configuration.disable_remainder) {
		messagef("Remainder classifier information");
		this->_configuration.remainder_classifier->print(verbose-1);
	}
}

/**
 * @brief Loads all subsets (iSets/Remainder) from file
 * @param reader An object-reader with binary data
 */
template <uint32_t N>
void NuevoMatch<N>::load_subsets(ObjectReader& reader) {

	// Lists to populate available iSets and any remainder rules
	_remainder_rules.clear();

	// Statistics
	uint32_t iset_rule_count=0;

	_isets = new IntervalSet<N>*[_num_of_isets];

	// Populate lists based on configuration
	for (uint32_t i=0; i<_num_of_isets; ++i) {

		// Get the handler of the next stored iSet
		ObjectReader sub_reader;
		reader >> sub_reader;

		// Read the current iSet
		IntervalSet<N>* iset = new IntervalSet<N>(i);
		iset->load(sub_reader);

		bool skip_current_iset =
				// Skip the current iSet in case the maximum number of iSets is limited
				((_configuration.max_subsets >= 0) && ((uint32_t)_configuration.max_subsets <= i)) ||
				// Skip the current iSet in case the minimal iSet number if limited
				(_configuration.start_from_iset > i) ||
				// The iSet field index should be skipped
				((_configuration.arbitrary_fields.size() > 0) &&
						(std::find(	_configuration.arbitrary_fields.begin(),
									_configuration.arbitrary_fields.end(),
									iset->get_field_index())
						 == _configuration.arbitrary_fields.end())
				);

		// In case the current iSet is valid but should not run
		if (!skip_current_iset && _configuration.disable_isets) {
			auto rules = iset->extract_rules();
			loggerf("Created a disabled iSet (%u) with %lu rules.", i, rules.size());
			_isets[i] = nullptr;
			delete iset;
		}
		// In case the current iSet should be skipped
		else if (skip_current_iset) {

			auto rules = iset->extract_rules();
			_remainder_rules.insert(_remainder_rules.end(), rules.begin(), rules.end());

			_isets[i] = nullptr;
			delete iset;

			loggerf("Skipping iSet %u. Extracted %lu rules.", i, rules.size());
		}
		// In case the current iSet is valid, and the disable-isets option is disabled
		else {
			_isets[i] = iset;

			// In case the field list is reconfigured
			if (_configuration.arbitrary_fields.size() > 0) {
				iset->rearrange_field_indices(_configuration.arbitrary_fields);
			}

			// Update statistics
			iset_rule_count += iset->size();
			_size += iset->get_size();
		}
	}

	// Read the predefined remainder rule-set, add to remainder
	ObjectReader db_reader(reader.buffer(), reader.size()); // TODO This is ugly, change packing to be using sub-reader
	std::list<openflow_rule> predefined_remainder = load_rule_database(db_reader);
	_remainder_rules.insert(_remainder_rules.end(), predefined_remainder.begin(), predefined_remainder.end());

	// Sort remainder rules by priority
	_remainder_rules.sort();
	uint32_t net_total_rules = (iset_rule_count+_remainder_rules.size());
	loggerf("Total rules after removing validation phase duplicates: %u", net_total_rules);

	// Print iSet coverage status
	for (uint32_t i=0; i<_num_of_isets; ++i) {
		if (_isets[i] == nullptr) continue;
		loggerf("iSet %u holds %u rules (coverage: %.2f) for field %u with RQRMI size of %u bytes",
				i, _isets[i]->size(), (scalar_t)_isets[i]->size() / net_total_rules * 100,
				_isets[i]->get_field_index(), _isets[i]->get_size());
	}

	// Print coverage status
	loggerf("NuevoMatch total coverage: %.2f%%", (double)iset_rule_count/net_total_rules*100);
}

/**
 * @brief Loads the remainder classifier based on subset configuration and input file
 * @param buffer The input buffer
 * @param size The buffer size in bytes
 */
template <uint32_t N>
void NuevoMatch<N>::load_remainder(ObjectReader& reader) {

	ObjectReader sub_reader;

	// In case the remainder classifier should be avoided
	if (_configuration.disable_remainder) {
		delete _configuration.remainder_classifier;
		_configuration.remainder_classifier = nullptr;
		return;
	}

	// In case the remainder classifier is external, do not change it
	if (_configuration.external_remainder) {
		if (_configuration.remainder_classifier == nullptr) {
			throw error("Remainder classifier was set as external, but is not available");
		}
		return;
	}

	// In case at least on iSet is missing, the classifier should be built
	bool rebuild_remainder = _configuration.force_rebuilding_remainder;
	for (uint32_t i=0; i<_num_of_isets;++i) {
		if (_isets[i] == nullptr) {
			rebuild_remainder = true;
			break;
		}
	}

        // Build remainder classifier from temporary rule-set
        if (rebuild_remainder) {
                sub_reader = build_remainder();
        }
        // Load the sub-reader from reader
        else {
                try {
                        reader >> sub_reader;
                } catch (const exception& e) {
                        throw error("Error while extracting remainder classifier: " << e.what());
                }
        }

        // Load classifier from sub-reader
        try {
                _configuration.remainder_classifier->load(sub_reader);
                return;
        } catch (const exception& e) {
                warning("Error while loading remainder classifier: " << e.what());
        }

        // Try to recover
        loggerf("Recovering by rebuilding remainder classifier");
        sub_reader = build_remainder();
        try {   
                _configuration.remainder_classifier->load(sub_reader);
                return;
        } catch (const exception& e) {
                error("Error while loading remainder classifier: " << e.what());
        }
}


/**
 * @brief Manually build remainder classifier
 */
template <uint32_t N>
ObjectReader NuevoMatch<N>::build_remainder() {
        loggerf("Manually building remainder classifier (remainder holds %lu rules)", _remainder_rules.size());
        // Building new classifier might thrash cash.
        // Therefore, the building is done using a temporary object
        GenericClassifier* gc;
        if (_configuration.remainder_type == "cutsplit") {
                gc = new CutSplit(24, 8);
        } else if (_configuration.remainder_type == "tuplemerge") {
                gc = new TupleMerge();
        } else {
                throw errorf("NuevoMatch cannot rebuild a remainder classifier of type %s", _configuration.remainder_type);
        }

        gc->build(_remainder_rules);
        // Pack classifier into reader
        ObjectReader output = ObjectReader(gc->pack());
        delete gc;
        return output;
}


/**
 * @brief Group the subsets based on their size (load-balance), and assign them to cores
 * @throws In case no valid subsets are available
 */
template <uint32_t N>
void NuevoMatch<N>::group_subsets_to_cores() {

	// Create a list of all subset classifiers based on availability
	std::vector<NuevoMatchSubset<N>*> subsets;

	// Add iSets
	for (uint32_t i=0; i<_num_of_isets;++i) {
		if (_isets[i] != nullptr) {
			subsets.push_back(_isets[i]);
		}
	}

	// Add remainder classifier
	if (_configuration.remainder_classifier != nullptr) {
		subsets.push_back(new NuevoMatchRemainderClassifier<N>(_configuration.remainder_classifier));
	}

	if (subsets.size() == 0) {
		throw error("Classifier has no valid subsets");
	}

	// Sort subsets based on the number of rules they hold (high to low)
	std::sort(subsets.begin(), subsets.end(),
			[](const NuevoMatchSubset<N>* a, const NuevoMatchSubset<N>* b) {
				return a->get_size() > b->get_size();
			});

	// Load balance between all classifiers and workers
	std::list<NuevoMatchSubset<N>*> classifier_list[_configuration.num_of_cores];

	// Is an arbitrary core allocation was set in the configuration?
	if (_configuration.arbitrary_subset_clore_allocation != "") {

		// Make sure the remainder classifier is first
		std::sort(subsets.begin(), subsets.end(),
			[](const NuevoMatchSubset<N>* a, const NuevoMatchSubset<N>* b) {
				return (a->get_type() != NuevoMatchSubset<N>::dynamic_type_t::ISET);
			});

		bool state_core = true;
		uint32_t current = 0;
		for (auto it : _configuration.arbitrary_subset_clore_allocation) {
			if (state_core && it == '=') state_core = false;
			else if (!state_core && it == ';') state_core = true;
			else if (!state_core && it == ',') continue;
			else if ((it < '0') || (it > '9')) throw errorf("Cannot parse char %c when allocating subsets to cores", it);
			else if (state_core) {
				current = it - '0';
				current = current >= _configuration.num_of_cores ? _configuration.num_of_cores - 1 : current;
			}
			else if (!state_core) {
				uint32_t idx = it - '0';
				if (idx >= subsets.size()) continue;
				classifier_list[current].push_back(subsets[idx]);
			}

		}

	}
	// There is not an arbitrary allocation, allocate based on size
	else {
		// Store the size in bytes used in each core
		uint32_t core_size[_configuration.num_of_cores];
		for (uint32_t i=0; i<_configuration.num_of_cores; ++i) {
			core_size[i] = 0;
		}

		for (auto it : subsets) {
			uint32_t current = 0, size_min=core_size[0];
			// Choose the core with minimum size
			for (uint32_t i=0; i<_configuration.num_of_cores; ++i) {
				if (core_size[i] < size_min) {
					current = i;
					size_min = core_size[i];
				}
			}
			// Add the current subset to the core
			classifier_list[current].push_back(it);
			core_size[current] += it->get_size();
		}
	}

	// Hold information of available CPUs
	vector<int> cpu_vec;
	int core_idx = cpu_core_tools_get_index_of_current_thread();
	cpu_vec.push_back(core_idx);

	// The current thread will run a serial worker
	_worker_serial = new NuevoMatchWorkerSerial<N>(0, this->_configuration);
	_worker_serial->add_listener(*this);
	for (auto it : classifier_list[0]) {
		_worker_serial->add_subset(*it);
	}

	// All other threads will run a parallel worker
	_workers_parallel = new NuevoMatchWorkerParallel<N>*[_configuration.num_of_cores-1];
	for (uint32_t i=1; i<_configuration.num_of_cores; ++i) {
		// Get next core index
		core_idx = cpu_core_tools_get_next_physical_core(core_idx);
		if (std::find(cpu_vec.begin(), cpu_vec.end(), core_idx) != cpu_vec.end()) {
			warningf("No available free CPU core for worker %u. Performance will degenerate", i);
		}
		cpu_vec.push_back(core_idx);

		// Build parallel worker
		_workers_parallel[i-1] = new NuevoMatchWorkerParallel<N>(i, this->_configuration, core_idx);
		_workers_parallel[i-1]->add_listener(*this);
		for (auto it : classifier_list[i]) {
			_workers_parallel[i-1]->add_subset(*it);
		}
	}

	// Print status of all workers
	for (uint32_t i=0; i<_configuration.num_of_cores; ++i) {

		// Calculate KB for the current worker
		uint32_t size = 0;
		for (auto it : classifier_list[i]) {
			size += it->get_size();
		}

		string_operations::convertor_to_string<NuevoMatchSubset<N>*> classifier_to_string =
				[](NuevoMatchSubset<N>* const& item) -> string { return item->to_string(); };
		string subset_string = string_operations::join(classifier_list[i], " ", classifier_to_string);

		// Print status of serial worker
		logger("NuevoMatch worker 0 on CPU " << cpu_vec[i] << " holds: {" << subset_string << "} of total " << size << " KB.");
	}

}

// Initiate template with custom batch sizes
template class NuevoMatch<1>;
template class NuevoMatch<2>;
template class NuevoMatch<4>;
template class NuevoMatch<8>;
template class NuevoMatch<16>;
template class NuevoMatch<32>;
template class NuevoMatch<64>;
template class NuevoMatch<128>;
template class NuevoMatch<256>;
template class NuevoMatch<512>;

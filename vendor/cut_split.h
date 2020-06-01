/**
 * Original Author:  Wenjun Li (Peking University, email: wenjunli@pku.edu.cn)
 * Updates By: 		 Alon Rashelbach (The Technion - Israel Institute of Technology, email: alonrs@campus.technion.ac.il)
 */

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <queue>
#include <vector>
#include <time.h>

#include <array_operations.h>
#include <generic_classifier.h>
#include <hyper_split.h>
#include <rule_db.h>

#define MAXCUTS  16
#define PTR_SIZE 4
#define LEAF_NODE_SIZE 4
#define TREE_NODE_SIZE 8

using namespace std;

struct field_length {
	uint32_t length[5];
	uint32_t size[5];
	uint32_t flag_smallest[4];
};

class CutSplitTrie {
private:

	typedef enum {CUT = 0, SPLIT} node_type_t;

	typedef struct {
		bool is_leaf;
		uint32_t num_of_rules;
		uint32_t num_of_cuts;
		uint32_t layer_number;
		node_type_t flag;
		// Arrays
		uint32_t *rule_indices;
		uint32_t *child_indices;
		int max_priority;
		range field[DIM];
		// HyperSplit SubTree
		hs_node_t* rootnode;
	} node_t;

	// Used for statistics
	uint32_t num_of_rules;
	uint32_t binth;
	uint32_t threshold;
	uint32_t dimension;

	// Pointer to rule database
	matching_rule	*rule_db;
	uint32_t *root_rules;

	// Used for statistics
	uint32_t total_rules;
	uint32_t total_leaves;
	uint32_t total_non_leaves;
	uint32_t max_depth;
	uint32_t max_layer;
	uint32_t total_nodes;
	float total_hs_memory_in_KB;

	// The set of nodes
	node_t* node_set;

	// Used for lookup
	int field_width[DIM];

	// Performance
	double work_time_ns, split_work_time_ns, linear_rule_time;
	bool measure_performance;
	uint32_t num_of_packets, split_lookups, max_linear_rules;

	// Private methods
	int  count_np_ficut(node_t*);
	void createtrie();
	uint32_t get_nbits(unsigned int n);
	uint32_t get_pow(unsigned int n);

public:

	typedef enum { SA = 0, DA } trie_type_t;

	CutSplitTrie(uint32_t num_of_rules, uint32_t binth, matching_rule* rule_db, uint32_t* root_rules, uint32_t threshold, trie_type_t tree_type);
	~CutSplitTrie();

	/**
	 * @brief Perform packet lookup
	 * @param header The packet header
	 */
	int lookup(const uint32_t* header, int priority);

	/**
	 * @brief Packs the trie to byte-array
	 */
	ObjectPacker pack();

	/**
	 * @brief Unpack trie from buffer
	 * @return 1 on success, otherwise 0
	 */
	int unpack(ObjectReader& packer);

	/**
	 * @brief Prints statistics of the current tree
	 */
	void print();

	/**
	 * @brief Returns the memory size of this in bytes
	 */
	uint32_t get_size() const;

	/**
	 * @brief Returns the maximum depth of this
	 */
	uint32_t get_max_gepth() const {
		return max_depth;
	}

	/**
	 * @brief Starts the performance measurement of this
	 */
	void start_performance_measurement();

	/**
	 * @brief Stops the performance measurement of this
	 */
	void stop_performance_measurement();
};


/**
 * @brief A CutSplit object matching the vector interface
 */
class CutSplit : public GenericClassifier {
private:

	// Rule DB
	uint32_t num_of_rules;
	matching_rule* rule_db;
	field_length* rule_field_length;

	// Hyper Parameters
	uint32_t threshold;
	uint32_t binth;

	// Performance
	struct timespec start_time, end_time;
	double big_tree_work_time_ns;
	bool measure_performance;

	// Statistics
	uint32_t size;
	uint32_t build_time;
	hs_result_t hs_result;

	// Hold the tries
	CutSplitTrie *tree_sa, *tree_da;
	bool has_big_tree;
	hs_node_t tree_big;

	/**
	 * @brief record length of field and corresponding size
	 */
	void count_length();

	/**
	 * @brief partition rule-set into subsets based on address field(2 dim.)
	 */
	void partition_v1(uint32_t* subset[3], uint32_t num_subset[3], uint32_t threshold_value[2]);

public:

	CutSplit(uint32_t threshold, uint32_t binth);

	/**
	 * @brief Build the classifier data structure
	 * @returns 1 On success, 0 on fail
	 */
	int build(const std::list<openflow_rule>& rule_db);

	/**
	 * @brief Packs this to byte array
	 * @returns An object-packer with the binary data
	 */
	ObjectPacker pack() const;

	/**
	 * @brief Creates this from a memory location
	 * @param object An object-reader instance
	 */
	void load(ObjectReader& object);

	/**
	 * @brief Returns the number of rules
	 */
	uint32_t get_num_of_rules() const { return num_of_rules; }

	/**
	 * @brief Returns the memory size of this in bytes
	 */
	uint32_t get_size() const { return size; }

	/**
	 * @brief Returns the building time of this in milliseconds
	 */
	uint32_t get_build_time() const { return build_time; }

	/**
	 * @brief Starts the performance measurement of this
	 */
	void start_performance_measurement();

	/**
	 * @brief Stops the performance measurement of this
	 */
	void stop_performance_measurement();

	/**
	 * @brief clones this to another instance
	 */
	virtual GenericClassifier* clone() {
		return new CutSplit(*this);
	}

	/**
	 * @brief Start an asynchronous process of classification for an input packet.
	 * @param header An array of 32bit integers according to the number of supported fields.
	 * @param priority The priority of a previous matching rule.
	 * Stops classifiying when there is no potential better priority
	 * @returns A unique id for the packet
	 */
	virtual unsigned int classify_async(const unsigned int* header, int priority);

	/**
	 * @brief Start a synchronous process of classification an input packet.
	 * @param header An array of 32bit integers according to the number of supported fields.
	 * @param priority The priority of a previous matching rule.
	 * Stops classifiying when there is no potential better priority
	 * @returns The matching rule action/priority (or 0xffffffff if not found)
	 */
	virtual unsigned int classify_sync(const unsigned int* header, int priority);

	/**
	 * @brief Prints statistical information
	 * @param verbose Set the verbosity level of printing
	 */
	virtual void print(uint32_t verbose=1) const;

	/**
	 * @brief Returns a string representation of this
	 */
	virtual const std::string to_string() const { return "CutSplit"; }

	/**
	 * @brief Returns the maximum supported number of fields this can classify
	 */
	virtual const unsigned int get_supported_number_of_fields() const { return FIVE_TUPLE_FIELDS; }
};

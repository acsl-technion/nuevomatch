/**
 *  Name:			trie.cpp for CutSplit
 *  Description:	trie construction for CutSplit: Pre-Cutting + Post-Splitting
 *  Version:		2.0 (release)
 *  Author:  		Wenjun Li (Peking University, email: wenjunli@pku.edu.cn)
 *  Updates: 		Alon Rashelbach (The Technion - Israel Institute of Technology, email: alonrs@campus.technion.ac.il)
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <queue>
#include <stack>
#include <list>
#include <sys/time.h>

#include <logging.h>
#include <object_io.h>
#include <cut_split.h>

#define MAX_UINT 0xffffffff

#ifndef MAX
#    define MAX(a,b) (a) > (b) ? (a) : (b);
#endif

/**
 * @brief This feature enables comprehensive performance tests,
 * but affects the performance itself. Use only for debugging.
 */
// #define MICRO_PERFORMANCE_TEST

using namespace std;

// Local methods
rule_set_t* get_hyper_split_rules(uint32_t num_of_rules, matching_rule* rule_db, uint32_t* indices);


CutSplitTrie::CutSplitTrie(uint32_t num_of_rules, uint32_t binth, matching_rule* rule_db, uint32_t* root_rules,
		uint32_t threshold, trie_type_t tree_type)
	: num_of_rules(num_of_rules), binth(binth), threshold(threshold),
	  dimension(tree_type), rule_db(rule_db), root_rules(root_rules),
	  total_rules(0), total_leaves(0), total_non_leaves(0), max_depth(0), max_layer(0), total_nodes(0),
	  total_hs_memory_in_KB(0), node_set(nullptr),
	  work_time_ns(0), split_work_time_ns(0), linear_rule_time(0),
	  num_of_packets(0), split_lookups(0), max_linear_rules(0)
{
	// The cut threshold is in power of 2
	this->threshold = pow(2, threshold);

	// Used for 5-tuple rules
	//sip:32 dip:32 sport:16 dport:16 protocol:8
    for(int i = 0; i< DIM; i++){
        if(i == 4) field_width[i] = 8;
        else if(i >= 2) field_width[i] = 16;
        else field_width[i] = 32;
    }

    createtrie();
}

CutSplitTrie::~CutSplitTrie() {
	// TODO memory leak. fix.
  delete [] this->node_set;
  delete [] this->rule_db;
}

/**
 * @brief count np_ficut in pre-cutting stage (using FiCut[from HybridCuts] for one small field)
 */
int CutSplitTrie::count_np_ficut(node_t *node) {
   int nump=0;
   // Check the value of the field in the tree's dimension (SA / DA)
   if(node->field[dimension].high == node->field[dimension].low) {
	   nump=1;
   } else {
	   nump=2;
   }

   // As long an additional cut is available according to threshold
   while( (nump < MAXCUTS) && (node->field[dimension].calc() > threshold) ) {
	   nump*=2;
   }
   return nump;
}

void CutSplitTrie::createtrie() {
   // Holds values of the current constructed rule
   node_t* new_node = nullptr;
   uint32_t node_idx = 0;

   // Holds temporary vectors for building the nodes
   vector<uint32_t> child_nodes;
   vector<uint32_t> rule_id;
   vector<node_t*> nodes;

   // The number of possible partitions
   uint32_t num_partitions = 0;

   // A queue for node relationships
   queue<uint32_t> node_queue;

   // Always create root node
   new_node = new node_t;
   new_node->is_leaf = false;
   new_node->num_of_rules = num_of_rules;
   new_node->num_of_cuts = 0;
   new_node->layer_number = 1;
   new_node->flag = CUT;
   new_node->rootnode = nullptr;
   new_node->max_priority = -1;

   // Reset root fields to all rule-space
   // TODO - this is a bit awkward. Fix.
   for (int i=0; i<DIM; i++){
	   new_node->field[i].low = 0;
       if (i<2) new_node->field[i].high = 0xffffffff;
       else if (i==4) new_node->field[i].high = 255;
       else new_node->field[i].high = 65535;
   }

   // Initiate the root with the subset
   rule_id.clear();
   for (uint32_t i=0; i<num_of_rules; ++i) {
	   rule_id.push_back(root_rules[i]);
   }
   new_node->rule_indices = vector_to_array(rule_id);

   // Push the current node
   nodes.push_back(new_node);
   node_queue.push(node_idx++);
   while(!node_queue.empty()){
	   // Pop first node
	   node_t* current_node = nodes[node_queue.front()];
	   node_queue.pop();

	   // In case the state is cut (FiCuts)
	   if (current_node->flag == CUT) {
		   // Count the number of possible partitions
		   num_partitions = count_np_ficut(current_node);
		   // Check whether to change to split
		   if(num_partitions < MAXCUTS) {
			   current_node->flag = SPLIT;
		   }
	   }

	   // In case of a leaf node
	   if (current_node->num_of_rules <= binth || num_partitions == 1) {
		   current_node->is_leaf = true;
		   total_rules += current_node->num_of_rules;
		   ++total_leaves;
		   // Update max depth
		   uint32_t current_depth = (current_node->num_of_rules+current_node->layer_number);
		   max_depth=MAX(max_depth, current_depth);
		   continue;
	   }

	   // Current node is not leaf
	   ++total_non_leaves;

	   // In case the state is cut (FiCuts)
	   if (current_node->flag == CUT) {

		   current_node->num_of_cuts = num_partitions;

		   // Initiate children
		   child_nodes.clear();

		   uint32_t stride  = current_node->field[dimension].calc() / num_partitions;
		   uint32_t l_bound, h_bound = current_node->field[dimension].low - 1;

		   for (uint32_t i=0; i<num_partitions; ++i) {
			   // Update bounds
			   l_bound = h_bound+1;
			   h_bound = l_bound + stride;

			   // Check how many rules exist inside current partition
			   rule_id.clear();
			   for (uint32_t r=0; r<current_node->num_of_rules; ++r) {
				   // Get the field of the current rule
				   uint32_t index = current_node->rule_indices[r];
				   range* rule_field = &rule_db[index].field[dimension];
				   // Check boundaries
				   if (	(rule_field->low  >= l_bound && rule_field->low  <= h_bound) ||
						(rule_field->high >= l_bound && rule_field->high <= h_bound) ||
						(rule_field->low  <= l_bound && rule_field->high >= h_bound) )
				   {
					   rule_id.push_back(index);
				   }
			   }

			   // In case the current partition is empty of rules
			   if (rule_id.size() == 0) {
				   child_nodes.push_back(MAX_UINT);
				   continue;
			   }

			   // The partition is not empty

			   // Create new node
			   new_node = new node_t;
			   new_node->num_of_rules = rule_id.size();
			   new_node->layer_number = current_node->layer_number+1;
			   new_node->is_leaf = rule_id.size() < binth;
			   new_node->flag = (num_partitions < MAXCUTS) ? SPLIT : CUT;
			   new_node->rule_indices = vector_to_array(rule_id);
			   new_node->max_priority = -1;
			   new_node->rootnode = nullptr;

			   // Update max layers
			   if (!new_node->is_leaf) {
				   max_layer = MAX(max_layer, new_node->layer_number);
			   }

			   // Initiate new node's fields
			   for (uint32_t f=0; f<DIM; f++){
				   if (f != dimension) {
					   new_node->field[f] = current_node->field[f];
				   } else {
					   new_node->field[f] = (range){ l_bound, h_bound };
				   }
			   }

			   // Set the new node as child of current node
			   child_nodes.push_back(node_idx);

			   // Push the new node
			   nodes.push_back(new_node);
			   node_queue.push(node_idx++);
		   }

		   // Set the children of the current node
		   current_node->child_indices = vector_to_array(child_nodes);
		   assert(child_nodes.size() == num_partitions);
	   }
	   // In case the current node is SPLIT (HyperSplit)
	   else {

		   // Convert the relevant rules to HyperSplit representation
		   rule_set_t* rule_set = get_hyper_split_rules(current_node->num_of_rules, rule_db, current_node->rule_indices);

		   // Create HyperSplit trie
		   current_node->rootnode = (hs_node_t*) malloc(sizeof(hs_node_t));
		   HyperSplitTrie hyper_split_tree(rule_set, binth, current_node->rootnode);

		   // Update statistics
		   total_hs_memory_in_KB += hyper_split_tree.result.total_mem_kb;
		   current_node->layer_number += hyper_split_tree.result.wst_depth;
		   max_depth=MAX(max_depth, current_node->layer_number);
       }
   }

   // Calculate the max priorities of all nodes (DFS)
   stack<int> node_stack;
   node_stack.push(0);
   while (!node_stack.empty()) {
	   node_t* current = nodes[node_stack.top()];
	   // In case the current node has already priority
	   if (current->max_priority > 0) {
		   node_stack.pop();
		   continue;
	   }
	   // In case this node is leaf
	   else if (current->is_leaf) {
		   int max_priority = 0x7fffffff;
		   for (int i=0; i<current->num_of_rules; ++i) {
			   matching_rule* rule = &rule_db[current->rule_indices[i]];
			   max_priority = std::min(max_priority, (int)rule->priority);
		   }
		   current->max_priority = max_priority;
		   node_stack.pop();
	   }
	   // In case the child nodes of this should be processed
	   else if (current->max_priority == -1) {
		   // In case the current node is CUT
		   if (current->flag == CUT) {
				for (uint32_t j=0; j<current->num_of_cuts; ++j) {
					uint32_t child_idx = current->child_indices[j];
					if (child_idx != MAX_UINT) node_stack.push(child_idx);
				}
				current->max_priority = -2;
		   }
		   // In case of split
		   else {
				if (current->rootnode != nullptr) {					
					current->max_priority = current->rootnode->max_priority;
				} else {
					current->max_priority = 0;
				}
			   node_stack.pop();
		   }
	   }
	   // The childs of this were already processed
	   else if (current->max_priority == -2) {
		   int max_priority = 0x7fffffff;
		   for (uint32_t j=0; j<current->num_of_cuts; ++j) {
			   uint32_t child_idx = current->child_indices[j];
			   if (child_idx != MAX_UINT) {
				   node_t* child = nodes[current->child_indices[j]];
				   max_priority = std::min(max_priority, child->max_priority);
			   }
		   }
		   current->max_priority = max_priority;
		   node_stack.pop();
	   }
   }

   // Convert the node vector to array
   node_set = vector_to_array(nodes);
   total_nodes = node_idx;
}

/**
 * @brief Prints statistics of the current tree
 */
void CutSplitTrie::print() {
	messagef("***%s Subset Tree (using FiCuts + HyperSlit):***", dimension == 0 ? "SA" : "DA:");
	messagef(" number of rules:%d", num_of_rules);
	messagef(" worst case tree level: %d", max_layer);
	messagef(" worst case tree depth: %d", max_depth);
	messagef(" total memory (all) (Pre-Cutting + Post_Splitting): %u(bytes)", get_size());
#ifdef MICRO_PERFORMANCE_TEST
	messagef(" Performance: average work time per packet: %.3lf us. "
			"average split work time per packet: %.3lf., "
			"average linear work time per packet: %.3lf., "
			"max linear rules: %u ,"
			"Total split lookups: %u",
			work_time_ns / 1e3 / num_of_packets,
			split_work_time_ns / 1e3 / split_lookups,
			linear_rule_time / 1e3 / (num_of_packets-split_lookups),
			max_linear_rules,
			split_lookups);
#endif
}

/**
 * @brief Returns the memory size of this in bytes
 */
uint32_t CutSplitTrie::get_size() const {

	// Total FiCuts nodes memory
	float total_ficuts_memory_bytes=
			total_rules*PTR_SIZE+
			total_nodes*PTR_SIZE+
			total_leaves*LEAF_NODE_SIZE+
			total_non_leaves*TREE_NODE_SIZE;

	return total_ficuts_memory_bytes + total_hs_memory_in_KB * 1024;
}

uint32_t CutSplitTrie::get_nbits(unsigned int n) {
    int k = 0;
    while (n >>= 1) k++;
    return k;
}

uint32_t CutSplitTrie::get_pow(unsigned int n) {
    int k = 1;
    while (n--) k*=2;
    return	k;
}

/**
 * @brief Starts the performance measurement of this
 */
void CutSplitTrie::start_performance_measurement() {
	num_of_packets = 0;
	work_time_ns = 0;
	split_work_time_ns = 0;
	split_lookups = 0;
	linear_rule_time = 0;
	max_linear_rules = 0;
	measure_performance = true;
}

/**
 * @brief Stops the performance measurement of this
 */
void CutSplitTrie::stop_performance_measurement() {
	measure_performance = false;
}

/**
 * @brief Converts a set of rules to HyperSplit rules
 * @param num_of_rules Total number of rules
 * @param indices The indices of the rules in the rule_db (or NULL, in case rule_db should be taken as is)
 */
rule_set_t* get_hyper_split_rules(uint32_t num_of_rules, matching_rule* rule_db, uint32_t* indices) {

	// Copy from local rule representation to HyperSplit rule representation
	rule_t* hyper_rule_db = (rule_t *)malloc(num_of_rules*sizeof(rule_t));
	for(uint32_t r=0; r<num_of_rules; ++r) {
		// Get current rule
		uint32_t index = indices ? indices[r] : r;
		matching_rule* current_rule = &rule_db[index];

		hyper_rule_db[r].pri = current_rule->priority;
		for(uint32_t f=0; f < DIM; f++){
		   hyper_rule_db[r].range[f][0] = current_rule->field[f].low;
		   hyper_rule_db[r].range[f][1] = current_rule->field[f].high;
		}
	}

	// Set the output
	rule_set_t* rule_set = new rule_set_t;
	rule_set->num = num_of_rules;
	rule_set->ruleList = hyper_rule_db;
	return rule_set;
}

/**
 * @brief Perform packet lookup
 * @param header The packet header
 */
int CutSplitTrie::lookup(const uint32_t* header, int priority) {
	// Not found is the default
	int result = priority;

	// Initiate the bits in field
	int current_bit = field_width[dimension];

    // Get three root
    node_t* current_node = &node_set[0];
    // Traverse until reaching a leaf
    while(!current_node->is_leaf){

    	// Stop in case the priority is higher
    	// than the maximum of the current node
    	if ( (priority >= 0) && (priority < current_node->max_priority) ) return priority;

    	uint32_t num_of_bits = get_nbits(current_node->num_of_cuts);

    	// Get the relevant child index
        uint32_t child_index = 0;
        for(uint32_t i = current_bit; i > current_bit-num_of_bits; i--){
            if((header[dimension] & 1<<(i-1)) != 0) {
                child_index += (int)get_pow(i-current_bit+num_of_bits-1);
            }
        }

        // Set the child index
        child_index = current_node->child_indices[child_index];

		// In case no relevant child, return not-found
        if (child_index == MAX_UINT) {
        	return priority;
        }

        // Update the current node
        current_node = &node_set[child_index];

        // In case the current node is SPLIT (HyperSplit) and not a leaf
        if (current_node->flag == SPLIT && !current_node->is_leaf) {
        	result = LookupHSTree(current_node->rootnode, header, priority);
            // Return the result
            return result;
        }

        // Update the current bit
        current_bit-= num_of_bits;
    }
    // Go over all rules
    for(uint32_t i=0; i<current_node->num_of_rules; ++i){
    	matching_rule* current_rule = &rule_db[current_node->rule_indices[i]];
    	// Get first rule that covers header
    	int cover = 1;
    	for(uint32_t j=0; j < DIM; j++){
			if(current_rule->field[j].low > header[j] || current_rule->field[j].high < header[j]){
				cover = 0;
				break;
			}
		}
    	// Return the rule index
    	if (cover) {
    		result = current_rule->priority;
    		break;
    	}
    }
    // Return the result
    return result;
}


/**
 * @brief Packs the trie to byte-array
 */
ObjectPacker CutSplitTrie::pack() {
	ObjectPacker output;

	// Write internal information
	output << num_of_rules;
	output << binth;
	output << threshold;
	output << dimension;
	output << total_rules;
	output << total_leaves;
	output << total_non_leaves;
	output << max_depth;
	output << max_layer;
	output << total_nodes;
	output << total_hs_memory_in_KB;

	// Pack the field width
	for (uint32_t f=0; f<DIM; ++f) {
		output << field_width[f];
	}

	// Pack the nodes
	for (uint32_t i=0; i<total_nodes; ++i) {
		output << node_set[i].num_of_rules;
		output << node_set[i].num_of_cuts;
		output << node_set[i].layer_number;
		output << (node_set[i].flag == CUT ? 0 : 1);
		output << node_set[i].is_leaf;

		for (uint32_t j=0; j<node_set[i].num_of_rules; ++j) {
			output << node_set[i].rule_indices[j];
		}

		if (node_set[i].flag == CUT && !node_set[i].is_leaf) {
			for (uint32_t j=0; j<node_set[i].num_of_cuts; ++j) {
				output << node_set[i].child_indices[j];
			}
		} else if (node_set[i].flag == SPLIT) {
			output << hstrie_pack(node_set[i].rootnode);
		}
	}

	return output;
}


/**
 * @brief Unpack trie from buffer
 * @return 1 on success, otherwise 0
 */
int CutSplitTrie::unpack(ObjectReader& reader) {
	uint32_t item;

	// Read static internal information
	reader >> this->num_of_rules;
	reader >> this->binth;
	reader >> this->threshold;
	reader >> this->dimension;
	reader >> this->total_rules;
	reader >> this->total_leaves;
	reader >> this->total_non_leaves;
	reader >> this->max_depth;
	reader >> this->max_layer;
	reader >> this->total_nodes;
	reader >> this->total_hs_memory_in_KB;

	// Read field width
	for (uint32_t f=0; f<DIM; ++f) {
		reader >> field_width[f];
	}

	// Allocate and read nodes
	node_set = new node_t[total_nodes];
	for (uint32_t i=0; i<total_nodes; ++i) {

		// Read node information
		reader >> node_set[i].num_of_rules;
		reader >> node_set[i].num_of_cuts;
		reader >> node_set[i].layer_number;
		reader >> item;
		node_set[i].flag = (item == 0 ? CUT : SPLIT);
		reader >> node_set[i].is_leaf;
		node_set[i].rootnode = nullptr;
		node_set[i].child_indices = nullptr;
		node_set[i].max_priority = -1;

		// Read the rules of the current node
		node_set[i].rule_indices = new uint32_t[node_set[i].num_of_rules];
		for (uint32_t j=0; j<node_set[i].num_of_rules; ++j) {
			reader >> node_set[i].rule_indices[j];
		}

		// In case the node is CUT and has children
		if (node_set[i].flag == CUT && !node_set[i].is_leaf) {
			// Read children of current node
			node_set[i].child_indices = new uint32_t[node_set[i].num_of_cuts];
			for (uint32_t j=0; j<node_set[i].num_of_cuts; ++j) {
				reader >> node_set[i].child_indices[j];
			}
		}
		// In case the node ise SPLIT
		else if (node_set[i].flag == SPLIT) {
			// Unpack internal object
			ObjectReader sub_reader;
			reader >> sub_reader;
			node_set[i].rootnode = HyperSplitTrie_unpack(sub_reader);
		}
	}

	// Calculate the max priorities of all nodes (DFS)
	stack<int> node_stack;
	node_stack.push(0);
	while (!node_stack.empty()) {
	   node_t* current = &node_set[node_stack.top()];
	   // In case the current node has already priority
	   if (current->max_priority > 0) {
		   node_stack.pop();
		   continue;
	   }
	   // In case this node is leaf
	   else if (current->is_leaf) {
		   int max_priority = 0x7fffffff;
		   for (int i=0; i<current->num_of_rules; ++i) {
			   matching_rule* rule = &rule_db[current->rule_indices[i]];
			   max_priority = std::min(max_priority, (int)rule->priority);
		   }
		   current->max_priority = max_priority;
		   node_stack.pop();
	   }
	   // In case the child nodes of this should be processed
	   else if (current->max_priority == -1) {
		   // In case the current node is CUT
		   if (current->flag == CUT) {
				for (uint32_t j=0; j<current->num_of_cuts; ++j) {
					uint32_t child_idx = current->child_indices[j];
					if (child_idx != MAX_UINT) node_stack.push(child_idx);
				}
				current->max_priority = -2;
		   }
		   // In case of split
		   else {
			  	if (current->rootnode != nullptr) {					
					current->max_priority = current->rootnode->max_priority;
				} else {
					current->max_priority = 0;
				}
			   node_stack.pop();
		   }
	   }
	   // The childs of this were already processed
	   else if (current->max_priority == -2) {
		   int max_priority = 0x7fffffff;
		   for (uint32_t j=0; j<current->num_of_cuts; ++j) {
			   uint32_t child_idx = current->child_indices[j];
			   if (child_idx != MAX_UINT) {
				   node_t* child = &node_set[current->child_indices[j]];
				   max_priority = std::min(max_priority, child->max_priority);
			   }
		   }
		   current->max_priority = max_priority;
		   node_stack.pop();
	   }
	}

	return 1;
}

/**
 * @brief record length of field and corresponding size
 */
void CutSplit::count_length() {
   unsigned temp_size=0;
   unsigned temp_value=0;
   for(uint32_t i=0;i<num_of_rules;i++) {
       for(uint32_t j=0;j<DIM;j++) {  //record field length in rule_field_length[i]
          rule_field_length[i].length[j]=rule_db[i].field[j].calc();
          if(rule_field_length[i].length[j]==0xffffffff) {
             rule_field_length[i].size[j]=32; //for address *
          } else {
             temp_size=0;
             temp_value=rule_field_length[i].length[j]+1;
             while((temp_value=temp_value/2)!=0) {
                temp_size++;
             }
             //for port number
             if((rule_field_length[i].length[j]+1 - pow(2,temp_size))!=0) {
               temp_size++;
             }
             rule_field_length[i].size[j]=temp_size;
          }
       }
   }
}

/**
 * @brief partition rule-set into subsets based on address field(2 dim.)
 */
void CutSplit::partition_v1(uint32_t* subset[3], uint32_t num_subset[3], uint32_t threshold_value[2]) {
  int num_small_tmp[num_of_rules];
  for(uint32_t i=0;i<num_of_rules;i++){
      num_small_tmp[i]=0;
      for(uint32_t k=0;k<2;k++)
         if(rule_field_length[i].size[k] <= threshold_value[k])
            num_small_tmp[i]++;
  }

  // Big tree - all rules that their ip-address (src/dst)
  // is larger than threshold
  int count_big=0;
  for(uint32_t i=0;i<num_of_rules;i++)
     if(num_small_tmp[i]==0)
        subset[0][count_big++]=i;
  num_subset[0]=count_big;

  int count_sa=0;
  int count_da=0;
  for(uint32_t i=0;i<num_of_rules;i++){

	  // All rules that are large in only one field are partitioned between
	  // tree SA (source address) and tree DA (destination address)
      if((num_small_tmp[i]==1)&&(rule_field_length[i].size[0]<=threshold_value[0]))
         subset[1][count_sa++]=i;
      if((num_small_tmp[i]==1)&&(rule_field_length[i].size[1]<=threshold_value[1]))
         subset[2][count_da++]=i;

      // Small in both fields (src+dst ip-address)
      if(num_small_tmp[i]==2) {
    	 // Partition based on the smaller field
         if(rule_field_length[i].size[0]<rule_field_length[i].size[1])
            subset[1][count_sa++]=i;
         else if(rule_field_length[i].size[0]>rule_field_length[i].size[1])
            subset[2][count_da++]=i;
         // The fields sizes are equal, add to the smaller tree
         else if(count_sa<=count_da)
            subset[1][count_sa++]=i;
         else
            subset[2][count_da++]=i;
      }
  }
   num_subset[1]=count_sa;
   num_subset[2]=count_da;
}

CutSplit::CutSplit(uint32_t threshold, uint32_t binth) :
		num_of_rules(0), rule_db(nullptr), rule_field_length(nullptr),
		threshold(threshold), binth(binth),
		big_tree_work_time_ns(0), measure_performance(false),
		size(0), build_time(0), tree_sa(nullptr), tree_da(nullptr),
		has_big_tree(false) {}

/**
 * @brief Build the classifier data structure
 * @returns 1 On success, 0 on fail
 */
int CutSplit::build(const std::list<openflow_rule>& rules) {

	// Convert the rule database to classic 5-tuple rule array
	uint32_t num_of_rules = rules.size(), counter = 0;
	matching_rule* rule_db = new matching_rule[num_of_rules];
	for (auto rule : rules) {
		rule_db[counter++] = rule.convert_to_five_tuple();
	}

	// Set the properties of this
	this->num_of_rules = num_of_rules;
	this->rule_db = rule_db;

	// Build this
	this->rule_field_length= new field_length[num_of_rules];
	count_length();

	// Partition rule-set into subsets based on address field(2 dim.)
	uint32_t* subset_3[3];
	for(int n=0;n<3;n++) {
	  subset_3[n]=(uint32_t*)malloc(num_of_rules*sizeof(uint32_t));
	}

	uint32_t num_subset_3[3]={0,0,0};
	uint32_t threshold_value_3[2]={threshold, threshold};
	partition_v1(subset_3, num_subset_3, threshold_value_3);

	// Process trees
	struct timespec	gStartTime,gEndTime;
	clock_gettime(CLOCK_MONOTONIC, &gStartTime);
	tree_sa = new CutSplitTrie(num_subset_3[1], binth, rule_db, subset_3[1], threshold, (CutSplitTrie::trie_type_t)0);
	tree_da = new CutSplitTrie(num_subset_3[2], binth, rule_db, subset_3[2], threshold, (CutSplitTrie::trie_type_t)1);

	// Process big rules using HyperSplit
	if(num_subset_3[0] > 0){
		rule_set_t* rule_set = get_hyper_split_rules(num_subset_3[0], rule_db, subset_3[0]);
		HyperSplitTrie trie(rule_set, binth, &tree_big);
		has_big_tree = true;
		// Update statistics
		hs_result = trie.result;
		size += trie.result.total_mem_kb * 1024;
	}
	clock_gettime(CLOCK_MONOTONIC, &gEndTime);

	// Set statistics
	build_time = (gEndTime.tv_sec  - gStartTime.tv_sec) * 1e3 +
				 (gEndTime.tv_nsec - gStartTime.tv_nsec) / 1e6;

	size = tree_sa->get_size() + tree_da->get_size();

	// All is well
	return 1;
}

/**
 * @brief Start a synchronous process of classification an input packet.
 * @param header An array of 32bit integers according to the number of supported fields.
 * @returns The matching rule action/priority (or 0xffffffff if not found)
 */
uint32_t CutSplit::classify_sync(const uint32_t* header, int priority) {

	int match_id = -1;
	int match_sa = tree_sa->lookup(header, priority);
	int match_da = tree_da->lookup(header, priority);
	int match_big = -1;


	if (has_big_tree) {
#ifdef MICRO_PERFORMANCE_TEST
		struct timespec big_tree_start_time, big_tree_end_time;
		bool measure = measure_performance;
		// Measure time for big tree classification
		if (measure) {
			clock_gettime(CLOCK_MONOTONIC, &big_tree_start_time);
		}
#endif

		// Perform classification
		match_big = LookupHSTree(&tree_big, header, priority);

#ifdef MICRO_PERFORMANCE_TEST
		// Measure time for big tree classification
		if (measure) {
			clock_gettime(CLOCK_MONOTONIC, &big_tree_end_time);
			big_tree_work_time_ns += (big_tree_end_time.tv_sec - big_tree_start_time.tv_sec) * 1e9 +
					(big_tree_end_time.tv_nsec - big_tree_start_time.tv_nsec);
		}
#endif
	}

	infof("match_sa = %d   match_da = %d   match_big = %d", match_sa, match_da, match_big);

	// Get highest prio rule
	if(match_sa != -1) match_id = match_sa;
	if((match_id == -1) || (match_da != -1 && match_da < match_id)) match_id = match_da;
	if((match_id == -1) || (match_big != -1 && match_big < match_id)) match_id = match_big;
	return match_id;
}

/**
 * @brief Start an asynchronous process of classification for an input packet.
 * @param header An array of 32bit integers according to the number of supported fields.
 * @returns A unique id for the packet
 */
uint32_t CutSplit::classify_async(const unsigned int* header, int priority) {
	uint32_t match_id = -1;
	uint32_t packet_id = 0xffffffff;

	// Lookup only valid packets
	if (header != nullptr) {
		// Perform lookup
		match_id = classify_sync(header, priority);
		packet_id = _packet_counter++;
	}

	// Notify all listeners
	for (auto it : _listeners) {
		it->on_new_result(packet_id, match_id, match_id, _additional_args);
	}
	return packet_id;
}

/**
 * @brief Prints debug information
 * @param verbose Set the verbosity level of printing
 */
void CutSplit::print(uint32_t verbose) const {

	// Maximum verbosity
	if (verbose > 2) {
		messagef("CutSplit internal tries:");
		tree_sa->print();
		tree_da->print();
	}

	// Calculate required memory accesses
	uint32_t total_mem_access = tree_sa->get_max_gepth() + tree_da->get_max_gepth();

	// Medium verbosity
	if (verbose > 1) {
		if (has_big_tree) {
			messagef("***Big rules(using HyperSlit):***");
			messagef(" number of rules: %d", hs_result.num_rules);
			messagef(" number of children: %d", hs_result.num_childnode);
			messagef(" worst case tree depth: %d", hs_result.wst_depth);
			messagef(" average tree depth: %f", hs_result.avg_depth);
			messagef(" number of tree nodes:%d", hs_result.num_tree_node);
			messagef(" number of leaf nodes:%d", hs_result.num_leaf_node);
			messagef(" total memory: %f(bytes)", hs_result.total_mem_kb * 1024);
			messagef(" Performance: average work time per packet: %.3lf us", big_tree_work_time_ns / 1e3 / _packet_counter);
			total_mem_access += hs_result.wst_depth;
		}
		messagef("Required memory accesses (worst-case): %u, binth: %u", total_mem_access, binth);
	}

	// Measure performance
	uint32_t total_usec = (end_time.tv_sec * 1e6 + end_time.tv_nsec / 1e3) -
						  (start_time.tv_sec * 1e6 + start_time.tv_nsec /1e3);
	messagef("Performance: total time %u usec. Average time: %.3f usec per packet.", total_usec, (double)total_usec / _packet_counter);

}

/**
 * @brief Packs this to byte array
 * @returns An object-packer with the binary data
 */
ObjectPacker CutSplit::pack() const {
	ObjectPacker output;

	// Write my own properties
	output << this->num_of_rules;
	output << this->threshold;
	output << this->binth;
	output << this->_packet_counter;
	output << this->size;
	output << this->build_time;
	output << this->has_big_tree;

	// Write the rule field length for this
	for (uint32_t i=0; i<num_of_rules; ++i) {
		for (int j=0; j<4; j++){
			output << this->rule_field_length[i].flag_smallest[j];
		}
		for (int j=0; j<DIM; j++){
			output << this->rule_field_length[i].length[j];
			output << this->rule_field_length[i].size[j];
		}
		output << rule_db[i].priority;
		for (uint32_t f=0; f<DIM; ++f) {
			output << rule_db[i].field[f].low;
			output << rule_db[i].field[f].high;
		}
	}

	// Write big trie results
	if (has_big_tree) {
		output << this->hs_result.avg_depth;
		output << this->hs_result.num_childnode;
		output << this->hs_result.num_leaf_node;
		output << this->hs_result.num_rules;
		output << this->hs_result.num_tree_node;
		output << this->hs_result.total_mem_kb;
		output << this->hs_result.wst_depth;
	}

	// Pack tries
	output << tree_da->pack();
	output << tree_sa->pack();

	// Pack HS trie
	if (has_big_tree) {
		output << hstrie_pack(&tree_big);
	}

	return output;
}

/**
 * @brief Creates this from a memory location
 * @param object An object-reader instance
 */
void CutSplit::load(ObjectReader& reader) {

	// Read my own properties
	reader >> this->num_of_rules;
	reader >> this->threshold;
	reader >> this->binth;
	reader >> this->_packet_counter;
	reader >> this->size;
	reader >> this->build_time;
	reader >> this->has_big_tree;

	// Write the rule field length for this
	this->rule_field_length= new field_length[num_of_rules];
	this->rule_db = new matching_rule[num_of_rules];
	for (uint32_t i=0; i<num_of_rules; ++i) {
		for (int j=0; j<4; j++){
			reader >> rule_field_length[i].flag_smallest[j];
		}
		for (int j=0; j<5; j++){
			reader >> rule_field_length[i].length[j];
			reader >> rule_field_length[i].size[j];
		}
		reader >> rule_db[i].priority;
		for (uint32_t f=0; f<DIM; ++f) {
			reader >> rule_db[i].field[f].low;
			reader >> rule_db[i].field[f].high;
		}
	}

	// Write big trie results
	if (has_big_tree) {
		reader >> this->hs_result.avg_depth;
		reader >> this->hs_result.num_childnode;
		reader >> this->hs_result.num_leaf_node;
		reader >> this->hs_result.num_rules;
		reader >> this->hs_result.num_tree_node;
		reader >> this->hs_result.total_mem_kb;
		reader >> this->hs_result.wst_depth;
	}

	// Read trees
	ObjectReader sub_reader;

	tree_da = new CutSplitTrie(0, binth, rule_db, nullptr, threshold, (CutSplitTrie::trie_type_t)1);
	reader >> sub_reader;
	tree_da->unpack(sub_reader);

	tree_sa = new CutSplitTrie(0, binth, rule_db, nullptr, threshold, (CutSplitTrie::trie_type_t)0);
	reader >> sub_reader;
	tree_sa->unpack(sub_reader);

	if (has_big_tree) {
		reader >> sub_reader;
		tree_big = *HyperSplitTrie_unpack(sub_reader);
	}

	// Validate all bytes were read
	assert(reader.size() == 0);
}

/**
 * @brief Starts the performance measurement of this
 */
void CutSplit::start_performance_measurement() {
	big_tree_work_time_ns = 0;
#ifdef MICRO_PERFORMANCE_TEST
	measure_performance = true;
	tree_sa->start_performance_measurement();
	tree_da->start_performance_measurement();
#endif
	clock_gettime(CLOCK_MONOTONIC, &start_time);
}

/**
 * @brief Stops the performance measurement of this
 */
void CutSplit::stop_performance_measurement() {
	clock_gettime(CLOCK_MONOTONIC, &end_time);
#ifdef MICRO_PERFORMANCE_TEST
	tree_sa->stop_performance_measurement();
	tree_da->stop_performance_measurement();
	measure_performance = false;
#endif
}

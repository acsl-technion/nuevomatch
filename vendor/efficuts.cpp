/**
 * The EffiCuts program
 * Original program written by Balajee Vamanan, Gwendolyn Voskuilen and T. N. Vijaykumar
 * of Purdue University
 * Changes made by James Daly as marked
 * Several tests was added by kun
 * The entire program was merged into a single file by Alon Rashelbach
 * Many changes were made by Alon to support building large classifiers
 */

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <list>
#include <set>
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <assert.h>

#include <time.h>

#include <efficuts.h> // Added by Alon
#include <logging.h> // Added by Alon

#define DIMS_PER_REP 5
#define NUM_JUNK 5

#define KUN_TEST

#define MAX_ALLOWED_LEVELS 200

#define PTR_SIZE 4
#define HEADER_SIZE 4
#define BOUNDARY_SIZE 16
#define INTERVAL_SIZE 2

#define INTERNAL_NODE_SIZE (HEADER_SIZE + BOUNDARY_SIZE + PTR_SIZE)
#define LEAF_NODE_SIZE HEADER_SIZE

#define INTERNAL_COST_64    3
#define INTERNAL_COST_128   2
#define RULE_COST_64    3
#define RULE_COST_128   2

#define TOO_MUCH 16

using namespace std;

uint32_t Percents[NUM_JUNK] =
{
  83,
  10,
  4,
  2,
  1
};

vector<uint32_t> UpperBounds;

uint32_t Cutoffs[NUM_JUNK];

int numReps = 1;

uint32_t bucketSize = 16;
double spfac = 8.0;
FILE *fpr;
FILE *fpt; //Added by kun for trace file
int hypercuts = 1;
int compressionON = 1; //kun: compress contiguous sibling nodes or not
int binningON = 1; //kun: separable trees
int mergingON = 1; //kun: selective tree merging
int fineOn = 1; //kun: use equi-dense cut or not
uint32_t num_intervals = 7; //kun: upperbound of a node's children number
//kun: identifying separable rules
double bin = 0.5; //for non-ip dimension
double IPbin = 0.05; //for ip dimension
int thirtyone = 0; //thirtyone == 1 --> make a difference between rules with one and no wildcards

int trace_rule_num; //Added by kun

int numTrees = 0;

// tree related
list <matching_rule> classifier;
int numrules=0;
node *root;

int rulelists[31];
list<matching_rule*> bigrules[5];
list<matching_rule*> kindabigrules[10];
list<matching_rule*> mediumrules[10];
list<matching_rule*> littlerules[5];
list<matching_rule*> smallrules;

int Num_Partitions;
int Avg_Degree;
int Max_Degree;
uint32_t Max_WorklistSize;
// Statistics
// live records
int Max_Depth;
int Max_Levels;
int Max_Cuts;
int Max_Access64Bit;
int Max_Access128Bit;
int Rules_at_the_Leaf;
int Rules_along_path;
uint32_t Total_Rule_Size;
uint32_t Total_Rules_Moved_Up;
uint32_t Total_Array_Size;
uint32_t Node_Count;
uint32_t Problematic_Node_Count;
uint32_t NonLeaf_Node_Count;
uint32_t Compressed_NonLeaf_Node_Count;
uint32_t Uncompressed_NonLeaf_Node_Count;
map <unsigned,uint32_t> interval_per_node;
map <unsigned,uint32_t> cuts_per_node;

int updateReads = 0;
int updateWrites = 0;

struct TreeStat
{
  // independent vars
  int Id;
  int No_Rules;
  int Max_Depth;
  int Max_Levels;
  int Max_Cuts;
  int Max_Access64Bit;
  int Max_Access128Bit;
  int Rules_at_the_Leaf;
  int Rules_along_path;
  uint32_t Total_Rule_Size;
  uint32_t Total_Rules_Moved_Up;
  uint32_t Total_Array_Size;
  uint32_t Node_Count;
  uint32_t Problematic_Node_Count;
  uint32_t NonLeaf_Node_Count;
  uint32_t Compressed_NonLeaf_Node_Count;
  uint32_t Uncompressed_NonLeaf_Node_Count;
  map <unsigned,uint32_t> interval_per_node;
  map <unsigned,uint32_t> cuts_per_node;
  // dependent vars
  uint32_t ruleptr_memory;
  uint32_t array_memory;
  uint32_t leaf_node_memory;
  uint32_t compressed_int_node_memory;
  uint32_t uncompressed_int_node_memory;
  uint32_t total_memory;
  uint32_t total_memory_in_KB;
};

struct MemBin
{
  int Max_Depth;
  int Max_Levels;
  int Max_Access64Bit;
  int Max_Access128Bit;
  uint32_t total_memory;
  list<TreeStat*> Trees;
};

// accumulated records
int treecount = 0;
TreeStat* p_record = nullptr;
list <TreeStat*> Statistics;
list<TreeDetails> efficuts_trees;

/*
 * Function declarations added by James Daly
 * Mostly for my own sanity
 */
//Comparers
bool mycomparison(matching_rule* first,matching_rule* second);
bool myequal(matching_rule* first,matching_rule* second);
bool mystatsort(TreeStat* first,TreeStat* second);
bool mymemsort(MemBin* first,MemBin* second);
//Checkers
int CheckIPBounds(range fld);
int CheckPortBounds(range fld);
int CheckProtoBounds(range fld);

bool IsPowerOfTwo(int x);

matching_rule get_bound(node *curr_node, uint32_t diemnsion, uint32_t index);
bool is_present(matching_rule boundary,matching_rule *rule);
bool is_equal(matching_rule rule1,matching_rule rule2, matching_rule boundary);

void calc_dimensions_to_cut(node *curr_node,int *select_dim);
list<node*> calc_num_cuts_1D(node *root,int dim);
list<node*> calc_num_cuts_2D(node *root,int *dim);
list<node*> calc_cuts(node *curr_node);
void createBoundary(node *a,node *b,node *c);
bool NodeCompress(list <node*> &node_list);

void InitStats(int No_Rules);
void NodeStats(node *curr_node);

void InterValHist(map <unsigned,uint32_t> interval_per_node);
void CutsHist(map <unsigned,uint32_t> cuts_per_node);

void PrintStatRecord(TreeStat *p_record);
void PrintStats();
void RecordTreeStats();

int samerules(node * r1, node * r2);
list<node*> merge_children(node * curr_node);
void regionCompaction(node * curr_node);
void create_tree(list <matching_rule*> & p_classifier);
//void node* CreateTreeFromRuleList(list<matching_rule*> p_classifier);

void binRules(list<matching_rule> & ruleList);
void MergeTrees();
void ComputeCutoffs();

// Functions written by J. Daly
bool DoRulesIntersect(matching_rule* r1, matching_rule* r2);

uint32_t BestSplitPoint(node * current, int dim);
list<node*> CalcMultiEquiCuts1D(node *root, int dim);
list<node*> CalcEquiCuts1D(node *root, int dim);
list<node*> CalcEquiCuts2D(node *root, int* dim);
node* SpawnChild(node* parent);

list<node*> CutNode(node *curr_node);
void StatNode(node* currNode);

void PrintRule(matching_rule* rule);
void PrintNode(node* node);
void PrintRuleList(list<matching_rule*> &rules);

int ColorOfList(list<matching_rule*> rules, const uint32_t *pt);
int ColorOfTree(node* tree, const uint32_t * pt);
int ColorOfTrees(list<TreeDetails>& trees, const uint32_t *pt);

void setNumReps(int reps);

void setNumReps(int reps) {
    numReps = reps;
    UpperBounds.resize(CLASSIFIER_FIELDS);
    for (int i = 0; i < reps; i++) {
        int index = i * DIMS_PER_REP;
        UpperBounds[index] = 4294967295;
        UpperBounds[index + 1] = 4294967295;
        UpperBounds[index + 2] = 65535;
        UpperBounds[index + 3] = 65535;
        UpperBounds[index + 4] = 255;
    }
}

// makes a = a + b
void createBoundary(node *a,node *b,node *c) {
	// This method was rewritten by AR
    for (int i = 0;i < CLASSIFIER_FIELDS;i++) {
        c->boundary.field[i].low = std::min(a->boundary.field[i].low, b->boundary.field[i].low);
        c->boundary.field[i].high = std::max(a->boundary.field[i].high, b->boundary.field[i].high);
    }

}

/**
 * @brief Compress adjacent nodes
 * @param node_list A list of nodes to compress
 * @note Changes made by Alon. Original version was both computation & memory inefficient
 */
bool NodeCompress(list <node*> &node_list) {
    uint32_t merge_possible = compressionON;
    uint32_t original_size = node_list.size();

    while (merge_possible) {
        merge_possible = 0;
        // Find the node with the largest number of rules
        uint32_t max=0;
        for (auto node : node_list) {
        	max = node->classifier.size() > max ? node->classifier.size() : max;
        }

        // Try to merge adjacent nodes
        auto first_it = node_list.begin();
        while (first_it != node_list.end()) {

        	// Get first node
        	node* current_item = *first_it;

        	// Get next node
        	auto next_it = first_it;
        	++next_it;
        	node* next_item = *next_it;
			if (next_it == node_list.end()) break;

        	// Try to merge both nodes
        	list<matching_rule*> merged_rules = current_item->classifier;
        	merged_rules.insert(merged_rules.begin(), next_item->classifier.begin(), next_item->classifier.end());
        	merged_rules.sort(mycomparison);
        	merged_rules.unique(myequal);

        	// What is the larger rule list?
        	uint32_t maxsize = std::max(current_item->classifier.size(), next_item->classifier.size());

        	// In case the merge is valid
        	if ( (merged_rules.size() <= bucketSize) ||
        		((merged_rules.size() <= maxsize) && (merged_rules.size() < max)) )
        	{
        		node tmp_node;
        		createBoundary(current_item, next_item, &tmp_node);
        		current_item->classifier=merged_rules;
        		current_item->boundary=tmp_node.boundary;

        		// Next node is no longer required

        		delete next_item;
        		node_list.erase(next_it);
        		merge_possible = 1;
        	}

        	++first_it;
        }
    }

    return node_list.size() < original_size;
}

/**
 * @brief Merge redundant children of a node
 * @note This method was rewritten by alon as the original version was messy
 */
list<node*> merge_children(node * curr_node) {
    list<node*> newlist = curr_node->children;
    auto first_it = newlist.begin();

    while (first_it != newlist.end()) {

    	auto second_it = first_it;
    	++second_it;

        while(second_it != newlist.end()) {
        	// Merge nodes that have the same rules
            if (samerules(*first_it, *second_it)) {
            	// The merged node should have the tighter boundaries between the two
                for (int i = 0; i < CLASSIFIER_FIELDS; i++) {
                    if ((*first_it)->boundary.field[i].low > (*second_it)->boundary.field[i].low) {
                        (*first_it)->boundary.field[i].low = (*second_it)->boundary.field[i].low;
                    }
                    if ((*first_it)->boundary.field[i].high > (*second_it)->boundary.field[i].high) {
                        (*first_it)->boundary.field[i].high = (*second_it)->boundary.field[i].high;
                    }
                }
                // The second node was merged into the first, so delete it
                delete(*second_it);
                second_it = newlist.erase(second_it);
            } else {
                second_it++;
            }
        }
        ++first_it;
    }
    return newlist;
}

void regionCompaction(node * curr_node) {
	// This method was rewritten by alon
	for (auto rule : curr_node->classifier) {
		for (uint32_t f=0; f<(uint32_t)CLASSIFIER_FIELDS; ++f) {
			curr_node->boundary.field[f].low = std::min(rule->field[f].low, curr_node->boundary.field[f].low);
			curr_node->boundary.field[f].high = std::max(rule->field[f].high, curr_node->boundary.field[f].high);
		}
	}
}


void binRules(list<matching_rule> &ruleList) {
    int min, wild;
    double field[5];

    int count = 0;
    for (list<matching_rule>::iterator itr = ruleList.begin(); itr != ruleList.end(); itr++) {
        count++;
        wild = 0;
        field[0] = ((double) ((*itr).field[0].high - (*itr).field[0].low))/0xFFFFFFFF;
        field[1] = ((double) ((*itr).field[1].high - (*itr).field[1].low))/0xFFFFFFFF;
        field[2] = ((double) ((*itr).field[2].high - (*itr).field[2].low))/65535;
        field[3] = ((double) ((*itr).field[3].high - (*itr).field[3].low))/65535;
        if (((*itr).field[4].low == 0) && ((*itr).field[4].high == 0xFF)) {
            field[4] = 1;
            wild++;
        } else {
            field[4] = 0;
        }
        min = 0;
        if (field[0] >= IPbin) { wild++; }
        if (field[1] >= IPbin) { wild++; }
        if (field[2] >= bin) { wild++; }
        if (field[3] >= bin) { wild++; }
        for (int i = 0; i < 4; i++) {
            if (field[i] < field[min]) {
                min = i;
            }
        }
        if (wild >= 4) {
            if ((field[0] > IPbin) && (field[1] > IPbin) && (field[2] > bin) && (field[3] > bin) && (field[4] != 1))
                bigrules[4].push_back(&(*itr));
            else {
                bigrules[min].push_back(&(*itr));
            }
        } else if (wild == 3) {
            if ((field[0] < IPbin) && (field[1] < IPbin)) {
                kindabigrules[9].push_back(&(*itr));    // wc except 0 and 1
            } else if ((field[0] < IPbin) && (field[2] < bin)){
                kindabigrules[8].push_back(&(*itr));    // wc except 0 and 2
            } else if ((field[0] < IPbin) && (field[3] < bin)) {
                kindabigrules[7].push_back(&(*itr));    // wc except 0 and 3
            } else if ((field[0] < IPbin) && (field[4] < bin)) {
                kindabigrules[6].push_back(&(*itr));    // wc except 0 and 4
            } else if ((field[1] < IPbin) && (field[2] < bin)) {
                kindabigrules[5].push_back(&(*itr));    // wc except 1 and 2
            } else if ((field[1] < IPbin) && (field[3] < bin)) {
                kindabigrules[4].push_back(&(*itr));    // wc except 1 and 3
            } else if ((field[1] < IPbin) && (field[4] < bin)) {
                kindabigrules[3].push_back(&(*itr));    // wc except 1 and 4
            } else if ((field[2] < bin) && (field[3] < bin)) {
                kindabigrules[2].push_back(&(*itr));    // wc except 2 and 3
            } else if ((field[2] < bin) && (field[4] < bin)) {
                kindabigrules[1].push_back(&(*itr));    // wc except 2 and 4
            } else if ((field[3] < bin) && (field[4] < bin)) {
                kindabigrules[0].push_back(&(*itr));    // wc except 3 and 4
            } else {
                printf("ERROR: Rule had 3 wc but did not match any of the bins!\n");
            }
        } else if (wild == 2) {
            if ((field[0] < IPbin) && (field[1] < IPbin) && (field[2] < bin)) {
                mediumrules[9].push_back(&(*itr));  // wc except 0, 1 and 2
            } else if ((field[0] < IPbin) && (field[1] < IPbin) && (field[3] < bin)){
                mediumrules[8].push_back(&(*itr));  // wc except 0, 1 and 3
            } else if ((field[0] < IPbin) && (field[1] < IPbin) && (field[4] < bin)) {
                mediumrules[7].push_back(&(*itr));  // wc except 0, 1 and 4
            } else if ((field[0] < IPbin) && (field[2] < bin) && (field[3] < bin)) {
                mediumrules[6].push_back(&(*itr));  // wc except 0, 2 and 3
            } else if ((field[0] < IPbin) && (field[2] < bin) && (field[4] < bin)) {
                mediumrules[5].push_back(&(*itr));  // wc except 0, 2 and 4
            } else if ((field[0] < IPbin) && (field[3] < bin) && (field[4] < bin)) {
                mediumrules[4].push_back(&(*itr));  // wc except 0, 3 and 4
            } else if ((field[1] < IPbin) && (field[2] < bin) && (field[3] < bin)) {
                mediumrules[3].push_back(&(*itr));  // wc except 1, 2 and 3
            } else if ((field[1] < IPbin) && (field[2] < bin) && (field[4] < bin)) {
                mediumrules[2].push_back(&(*itr));  // wc except 1, 2 and 4
            } else if ((field[1] < IPbin) && (field[3] < bin) && (field[4] < bin)) {
                mediumrules[1].push_back(&(*itr));  // wc except 1, 3 and 4
            } else if ((field[2] < bin) && (field[3] < bin) && (field[4] < bin)) {
                mediumrules[0].push_back(&(*itr));  // wc except 2, 3 and 4
            } else {
                printf("ERROR: Rule had 2 wc but did not match any of the bins!: %lf, %lf, %lf, %lf, %lf\n",field[0],field[1],field[2],field[3],field[4]);
            }
        } else {
            if (thirtyone) {
                if (wild == 1) {
                    if (field[0] >= IPbin) {
                        littlerules[0].push_back(&(*itr));
                        //printf("littlerules[0]\n");
                    } else if (field[1] >= IPbin){
                        //printf("littlerules[1]\n");
                        littlerules[1].push_back(&(*itr));
                    } else if (field[2] >= IPbin) {
                        littlerules[2].push_back(&(*itr));
                        //printf("littlerules[1]\n");
                    } else if (field[3] >= IPbin) {
                        //printf("littlerules[3]\n");
                        littlerules[3].push_back(&(*itr));
                    } else if (field[4] >= IPbin) {
                        //printf("littlerules[4]\n");
                        littlerules[4].push_back(&(*itr));
                    } else {
                        printf("ERROR: Rule had 1 wc but did not match any of the bins!\n");
                    }
                } else {
                    smallrules.push_back(&(*itr));
                }
            } else {
                smallrules.push_back(&(*itr));
            }
        }
    }
    numTrees = 0;

    for (int i = 0; i < 5; i++) {
        if (bigrules[i].size() > 0) {
            numTrees++;
            rulelists[i] = 1;
        } else {
            rulelists[i] = 0;
        }
    }
    for (int j = 0; j < 10; j++) {
        if (kindabigrules[j].size() > 0) {
            numTrees++;
            rulelists[(j+5)] = 1;
        } else {
            rulelists[(j+5)] = 0;
        }
    }
    for (int k = 0; k < 10; k++) {
        if (mediumrules[k].size() > 0) {
            numTrees++;
            rulelists[k+15] = 1;
        } else {
            rulelists[k+15] = 0;
        }
    }
    for (int l = 0; l < 5; l++) {
        if (littlerules[l].size() > 0) {
            numTrees++;
        }
    }
    if (smallrules.size() > 0) {
        numTrees++;
        rulelists[25] = 1;
    } else {
        rulelists[25] = 0;
    }
}

/*
 *  Method to merge trees together
 *  Will try to merge trees that have no more than one field that is not overlapping (i.e. where one tree is WC and one tree is not)
 */
void MergeTrees() {
#ifndef KUN_TEST
    printf("Number of trees before merge: %d\n",numTrees);
#endif
    int merged[26]; // array - if the value is 0 than that try is not merged, if it is 1 it has been and is NOT a candidate for merging anymore!
    for (int i = 0; i < 26; i++) { merged[i] = 0; } // make sure array is initialized to 0

    // try to merge any of the 1* into a 2* if it exists
    if (rulelists[0] == 1) {
        if (rulelists[11] == 1 && !(merged[11])) {
            bigrules[0].merge(kindabigrules[6],mycomparison);
            rulelists[11] = 0;
            merged[0] = 1;
            numTrees--;
        } else if (rulelists[12] == 1 && !(merged[12])) {
            bigrules[0].merge(kindabigrules[7],mycomparison);
            rulelists[12] = 0;
            merged[0] = 1;
             numTrees--;
        } else if (rulelists[13] == 1 && !(merged[13])) {
            bigrules[0].merge(kindabigrules[8],mycomparison);
            rulelists[13] = 0;
            merged[0] = 1;
            numTrees--;
        } else if (rulelists[14] == 1 && !(merged[14])) {
            bigrules[0].merge(kindabigrules[9],mycomparison);
            rulelists[14] = 0;
            merged[0] = 1;
            numTrees--;
        }
    }

    if (rulelists[1] == 1) {
         if (rulelists[8] == 1 && !(merged[8])) {
            bigrules[1].merge(kindabigrules[3],mycomparison);
            rulelists[8] = 0;
            merged[1] = 1;
            numTrees--;
        } else if (rulelists[9] == 1 && !(merged[9])) {
            bigrules[1].merge(kindabigrules[4],mycomparison);
            rulelists[9] = 0;
            merged[1] = 1;
            numTrees--;
        } else if (rulelists[10] == 1 && !(merged[10])) {
            bigrules[1].merge(kindabigrules[5],mycomparison);
            rulelists[10] = 0;
            merged[1] = 1;
            numTrees--;
        } else if (rulelists[14] == 1 && !(merged[14])) {
            bigrules[1].merge(kindabigrules[9],mycomparison);
            rulelists[14] = 0;
            merged[1] = 1;
            numTrees--;
        }
    }
     if (rulelists[2] == 1) {
         if (rulelists[6] == 1 && !(merged[6])) {
            bigrules[2].merge(kindabigrules[1],mycomparison);
            rulelists[6] = 0;
            merged[2] = 1;
            numTrees--;
        } else if (rulelists[7] == 1 && !(merged[7])) {
            bigrules[2].merge(kindabigrules[2],mycomparison);
            rulelists[7] = 0;
            merged[2] = 1;
            numTrees--;
        } else if (rulelists[10] == 1 && !(merged[10])) {
            bigrules[2].merge(kindabigrules[5],mycomparison);
            rulelists[10] = 0;
            merged[2] = 1;
            numTrees--;
        } else if (rulelists[13] == 1 && !(merged[13])) {
            bigrules[2].merge(kindabigrules[8],mycomparison);
            rulelists[13] = 0;
            merged[2] = 1;
            numTrees--;
        }
    }
    if (rulelists[3] == 1) {
         if (rulelists[5] == 1 && !(merged[5])) {
            bigrules[3].merge(kindabigrules[0],mycomparison);
            rulelists[5] = 0;
            merged[3] = 1;
            numTrees--;
        } else if (rulelists[7] == 1 && !(merged[7])) {
            bigrules[3].merge(kindabigrules[2],mycomparison);
            rulelists[7] = 0;
            merged[3] = 1;
            numTrees--;
        } else if (rulelists[9] == 1 && !(merged[9])) {
            bigrules[3].merge(kindabigrules[4],mycomparison);
            rulelists[9] = 0;
            merged[3] = 1;
            numTrees--;
        } else if (rulelists[12] == 1 && !(merged[12])) {
            bigrules[3].merge(kindabigrules[7],mycomparison);
            rulelists[12] = 0;
            merged[3] = 1;
            numTrees--;
        }
    }
    if (rulelists[4] == 1) {
        if (rulelists[5] == 1 && !(merged[5])) {
            bigrules[4].merge(kindabigrules[0],mycomparison);
            rulelists[5] = 0;
            merged[4] = 1;
            numTrees--;
        } else if (rulelists[6] == 1 && !(merged[6])) {
            bigrules[4].merge(kindabigrules[1],mycomparison);
            rulelists[6] = 0;
            merged[4] = 1;
            numTrees--;
        } else if (rulelists[8] == 1 && !(merged[8])) {
            bigrules[4].merge(kindabigrules[3],mycomparison);
            rulelists[8] = 0;
            merged[4] = 1;
            numTrees--;
        } else if (rulelists[11] == 1 && !(merged[11])) {
            bigrules[4].merge(kindabigrules[6],mycomparison);
            rulelists[11] = 0;
            merged[4] = 1;
            numTrees--;
        }
    }
    if (rulelists[5] == 1) {
        if (rulelists[15] == 1 && !(merged[15])) {
            kindabigrules[0].merge(mediumrules[0],mycomparison);
            rulelists[15] = 0;
            merged[5] = 1;
            numTrees--;
        } else if (rulelists[16] == 1 && !(merged[16])) {
            kindabigrules[0].merge(mediumrules[1],mycomparison);
            rulelists[16] = 0;
            merged[5] = 1;
            numTrees--;
        } else if (rulelists[19] == 1 && !(merged[19])) {
            kindabigrules[0].merge(mediumrules[4],mycomparison);
            rulelists[19] = 0;
            merged[5] = 1;
            numTrees--;
        }
    }
    if (rulelists[6] == 1) {
        if (rulelists[15] == 1 && !(merged[15])) {
            kindabigrules[1].merge(mediumrules[0],mycomparison);
            rulelists[15] = 0;
            merged[6] = 1;
            numTrees--;
        } else if (rulelists[17] == 1 && !(merged[17])) {
            kindabigrules[1].merge(mediumrules[2],mycomparison);
            rulelists[17] = 0;
            merged[6] = 1;
            numTrees--;
        } else if (rulelists[20] == 1 && !(merged[20])) {
            kindabigrules[1].merge(mediumrules[5],mycomparison);
            rulelists[20] = 0;
            merged[6] = 1;
            numTrees--;
        }
    }
    if (rulelists[7] == 1) {
        if (rulelists[15] == 1 && !(merged[15])) {
            kindabigrules[2].merge(mediumrules[0],mycomparison);
            rulelists[15] = 0;
            merged[7] = 1;
            numTrees--;
        } else if (rulelists[18] == 1 && !(merged[18])) {
            kindabigrules[2].merge(mediumrules[3],mycomparison);
            rulelists[18] = 0;
            merged[7] = 1;
            numTrees--;
        } else if (rulelists[21] == 1 && !(merged[21])) {
            kindabigrules[2].merge(mediumrules[6],mycomparison);
            rulelists[21] = 0;
            merged[7] = 1;
            numTrees--;
        }
    }
    if (rulelists[8] == 1) {
        if (rulelists[16] == 1 && !(merged[16])) {
            kindabigrules[3].merge(mediumrules[1],mycomparison);
            rulelists[16] = 0;
            merged[8] = 1;
            numTrees--;
        } else if (rulelists[17] == 1 && !(merged[17])) {
            kindabigrules[3].merge(mediumrules[2],mycomparison);
            rulelists[17] = 0;
            merged[8] = 1;
            numTrees--;
        } else if (rulelists[22] == 1 && !(merged[22])) {
            kindabigrules[3].merge(mediumrules[7],mycomparison);
            rulelists[22] = 0;
            merged[8] = 1;
            numTrees--;
        }
    }
    if (rulelists[9] == 1) {
        if (rulelists[16] == 1 && !(merged[16])) {
            kindabigrules[4].merge(mediumrules[1],mycomparison);
            rulelists[16] = 0;
            merged[9] = 1;
            numTrees--;
        } else if (rulelists[18] == 1 && !(merged[18])) {
            kindabigrules[4].merge(mediumrules[3],mycomparison);
            rulelists[18] = 0;
            merged[9] = 1;
            numTrees--;
        } else if (rulelists[23] == 1 && !(merged[23])) {
            kindabigrules[4].merge(mediumrules[8],mycomparison);
            rulelists[23] = 0;
            merged[9] = 1;
            numTrees--;
        }
    }
    if (rulelists[10] == 1) {
        if (rulelists[17] == 1 && !(merged[17])) {
            kindabigrules[5].merge(mediumrules[2],mycomparison);
            rulelists[17] = 0;
            merged[10] = 1;
            numTrees--;
        } else if (rulelists[18] == 1 && !(merged[18])) {
            kindabigrules[5].merge(mediumrules[3],mycomparison);
            rulelists[18] = 0;
            merged[10] = 1;
            numTrees--;
        } else if (rulelists[24] == 1 && !(merged[24])) {
            kindabigrules[5].merge(mediumrules[9],mycomparison);
            rulelists[24] = 0;
            merged[10] = 1;
            numTrees--;
        }
    }
    if (rulelists[11] == 1) {
        if (rulelists[19] == 1 && !(merged[19])) {
            kindabigrules[6].merge(mediumrules[4],mycomparison);
            rulelists[19] = 0;
            merged[11] = 1;
            numTrees--;
        } else if (rulelists[20] == 1 && !(merged[20])) {
            kindabigrules[6].merge(mediumrules[5],mycomparison);
            rulelists[20] = 0;
            merged[11] = 1;
            numTrees--;
        } else if (rulelists[22] == 1 && !(merged[22])) {
            kindabigrules[6].merge(mediumrules[7],mycomparison);
            rulelists[22] = 0;
            merged[11] = 1;
            numTrees--;
        }
    }
    if (rulelists[12] == 1) {
        if (rulelists[19] == 1 && !(merged[19])) {
            kindabigrules[7].merge(mediumrules[4],mycomparison);
            rulelists[19] = 0;
            merged[12] = 1;
            numTrees--;
        } else if (rulelists[21] == 1 && !(merged[21])) {
            kindabigrules[7].merge(mediumrules[6],mycomparison);
            rulelists[21] = 0;
            merged[12] = 1;
            numTrees--;
        } else if (rulelists[23] == 1 && !(merged[23])) {
            kindabigrules[7].merge(mediumrules[8],mycomparison);
            rulelists[23] = 0;
            merged[12] = 1;
            numTrees--;
        }
    }
    if (rulelists[13] == 1) {
        if (rulelists[20] == 1 && !(merged[20])) {
            kindabigrules[8].merge(mediumrules[5],mycomparison);
            rulelists[20] = 0;
            merged[13] = 1;
            numTrees--;
        } else if (rulelists[21] == 1 && !(merged[21])) {
            kindabigrules[8].merge(mediumrules[6],mycomparison);
            rulelists[21] = 0;
            merged[13] = 1;
            numTrees--;
        } else if (rulelists[24] == 1 && !(merged[24])) {
            kindabigrules[8].merge(mediumrules[9],mycomparison);
            rulelists[24] = 0;
            merged[13] = 1;
            numTrees--;
        }
    }
    if (rulelists[14] == 1) {
        if (rulelists[22] == 1 && !(merged[22])) {
            kindabigrules[9].merge(mediumrules[7],mycomparison);
            rulelists[22] = 0;
            merged[14] = 1;
            numTrees--;
        } else if (rulelists[23] == 1 && !(merged[23])) {
            kindabigrules[9].merge(mediumrules[8],mycomparison);
            rulelists[23] = 0;
            merged[14] = 1;
            numTrees--;
        } else if (rulelists[24] == 1 && !(merged[24])) {
            kindabigrules[9].merge(mediumrules[9],mycomparison);
            rulelists[24] = 0;
            merged[14] = 1;
            numTrees--;
        }
    }
#if 0
    for (int i = 0; i < 9; i++) {
        if (rulelists[i+15] == 1 && rulelists[25] == 1) {
            mediumrules[i].merge(smallrules,mycomparison);
            merged[i+15] = 1;
            rulelists[25] = 0;
            numTrees--;
            break;
        }
    }
#endif
#ifndef KUN_TEST
    printf("Number of trees after merge: %d\n",numTrees);
#endif
}

/**
 * Method for creating arbirary 2D cuts
 * Written by James Daly
 */
list<node*> CalcEquiCuts2D(node *root, int* dims) {
    const int numSplits = 2;

    node current = *root;
    list<node*> output;
    uint32_t splits[numSplits];

    for (int i = 0; i < numSplits; i++) {
        splits[i] = BestSplitPoint(&current, dims[i]);
        root->cuts[i] = 2;
    }

    for (int i = 0; i < numSplits; i++)
    {
        for (int j = 0; j < numSplits; j++)
        {

            node* child = SpawnChild(&current);

            if (i == 0)
            {
                child->boundary.field[dims[0]].high = splits[0];
            }
            else
            {
                child->boundary.field[dims[0]].low = splits[0] + 1;
            }

            if (j == 0)
            {
                child->boundary.field[dims[1]].high = splits[1];
            }
            else
            {
                child->boundary.field[dims[1]].low = splits[1] + 1;
            }

            // TODO : set bounds
            for (auto rule : current.classifier)
            {
                if (is_present(child->boundary, rule)) {
                    child->classifier.push_back(rule);
                }
            }
            output.push_back(child);
        }
    }

	return output;
}

/**
 * Written for creating wider arbitrary splits
 */
list<node*> CalcMultiEquiCuts1D(node *root, int dim) {
	list<node*> output;
    list<node*> clone_list = CalcEquiCuts1D(root, dim);

    for (auto nodex : clone_list) {
        if (nodex->classifier.size() > bucketSize) {
            list<node*> children = CalcEquiCuts1D(nodex, dim);
            while (!children.empty()) {
                node* n = children.front();
                children.pop_front();
                n->depth--;
                output.push_back(n);
            }
        } else {
        	output.push_back(nodex);
        }
    }
    root->cuts[dim] = output.size();
    return output;
}

/**
 * Method for creating arbitrary splits
 * Written by James Daly
 */
list<node*> CalcEquiCuts1D(node *root, int dim)
{
    list<node*> results;

    node current = *root;

    uint32_t split = BestSplitPoint(&current, dim);

    // Split

    node* lowerNode = SpawnChild(&current);
    node* upperNode = SpawnChild(&current);
    lowerNode->boundary.field[dim].high = split;
    upperNode->boundary.field[dim].low = split + 1;
    // TODO : initialize nodes
    for (auto rule : current.classifier)
    {
        if (rule->field[dim].low <= split)
        {
            lowerNode->classifier.push_back(rule);
        }
        if (rule->field[dim].high > split)
        {
            upperNode->classifier.push_back(rule);
        }
    }

    if (lowerNode->classifier.size() == current.classifier.size())
    {
        cout << "size problems:" << split << endl;
        delete lowerNode;
        delete upperNode;
        results.push_back(root);
    } else {
        results.push_back(lowerNode);
        results.push_back(upperNode);
    }

    return results;
}

/**
 * Helper method for creating child nodes
 * Memory is allocated by this method that the client is responsible for
 * Written by James Daly
 */
node* SpawnChild(node* parent) {
    node* child = new node;
    child->depth = parent->depth + 1;
    child->boundary = parent->boundary;
    return child;
}

/**
 * Helper method for calculating the best split point
 * Written by James Daly
 */
uint32_t BestSplitPoint(node * current, int dim)
{
    vector<uint32_t> ubounds;
    vector<int> endsByCount;
    vector<int> crossesCount;

    // Find all of the unique upper bounds
    for (list <matching_rule*>::iterator rule = current->classifier.begin();
            rule != current->classifier.end(); ++rule)
    {
        uint32_t upper = min((*rule)->field[dim].high, current->boundary.field[dim].high);
        bool contains = false;
        for (vector<uint32_t>::iterator item = ubounds.begin();
                item != ubounds.end(); ++item)
        {
            if (*item == upper)
            {
                contains = true;
                break;
            }
        }
        if (!contains)
        {
            ubounds.push_back(upper);
            endsByCount.push_back(0);
            crossesCount.push_back(0);
        }
    }
    sort(ubounds.begin(), ubounds.end());

    // tally the number of rules that end at certain spots
    // or that cross each spot
    for (list <matching_rule*>::iterator rule = current->classifier.begin();
            rule != current->classifier.end(); ++rule)
    {
        for (int i = 0; i < (int)ubounds.size(); i++)
        {
            if ((*rule)->field[dim].high <= ubounds[i])
            {
                endsByCount[i] = endsByCount[i] + 1;
                break;
            }
            else
            {
                if ((*rule)->field[dim].low < ubounds[i])
                {
                    crossesCount[i] = crossesCount[i] + 1;
                }
            }
        }
    }

    // find the best split point
    // Try to get half the rules on either side of the split
    // and no crossings
    // Should probably weight crossings higher
    // so rules don't get replicated
    uint32_t split = 0;
    int worstCost = numeric_limits<int>::max();
    int tally = 0;//-current->classifier.size() / 2;
    for (int i = 0; i < (int)ubounds.size() - 1; i++)
    {
        tally += endsByCount[i];
        //int cost = (tally < 0) ? -tally : tally;
        //cost += crossesCount[i];// * 2;
        int cost = 2 * tally - current->classifier.size() + crossesCount[i];
        cost *= (cost < 0) ? -1 : 1;
        cost += crossesCount[i];
        //cout << cost << " " << tally << " " << crossesCount[i] << endl;
        if (cost < worstCost)
        {

            worstCost = cost;
            split = ubounds[i];
        }
    }

    //cout << worstCost << " " << split << endl;

    if (split >= ubounds.back())
    {
        cout << "Bad split: " << split << " " << worstCost << endl;
        for (int i = 0; i < (int)ubounds.size(); i++)
        {
            cout << endsByCount[i] << " " << crossesCount[i] << endl;
        }
        exit(1);
    }
    return split;
}


void PrintRule(matching_rule* rule)
{
    for (int i = 0; i < CLASSIFIER_FIELDS; i++)
    {
        printf("[%010u,%010u]\t", (unsigned int)rule->field[i].low, (unsigned int)rule->field[i].high);
    }
    printf("%u\n", rule->priority);
}

void PrintNode(node *curr_node)
{
    printf("Node at depth %u\n", curr_node->depth);
    PrintRule(&curr_node->boundary);
    printf("%u rules\n", (unsigned int)curr_node->classifier.size());
    PrintRuleList(curr_node->classifier);
}

void PrintRuleList(list<matching_rule*> &rules)
{
    for (list<matching_rule*>::iterator iter = rules.begin();
            iter != rules.end(); iter++)
    {
        PrintRule(*iter);
    }
}

matching_rule get_bound(node *curr_node, uint32_t diemnsion, uint32_t index) {

    matching_rule boundary;
    uint32_t interval;

    // This loop was changed by alon
    for (uint32_t i = 0;i < (uint32_t)CLASSIFIER_FIELDS;i++) {
    	// Get the interval of the boundary for each dimension
        interval = curr_node->boundary.field[i].high - curr_node->boundary.field[i].low + 1;
        interval = interval / curr_node->cuts[i]; // cuts can be either 1 or 2

        // Set the boundary on current field
        uint32_t offset = (i==diemnsion) ? index : 0;
        boundary.field[i].low = curr_node->boundary.field[i].low + offset * interval;
        boundary.field[i].high = std::min(curr_node->boundary.field[i].high, boundary.field[i].low + interval - 1);
    }

    // check the node's bounds
    if (CheckIPBounds(boundary.field[0]))
    {
        printf("Error: get_bound bounds check for 0 failed\n");
        printf("[%u - %u] => [%u - %u] @ %d\n",
            curr_node->boundary.field[0].low,curr_node->boundary.field[0].high,
            boundary.field[0].low,boundary.field[0].high,curr_node->cuts[0]);
        exit(1);
    }
    if (CheckIPBounds(boundary.field[1]))
    {
        printf("Error: get_bound bounds check for 1 failed\n");
        printf("[%u - %u] => [%u - %u] @ %d\n",
            curr_node->boundary.field[1].low,curr_node->boundary.field[1].high,
            boundary.field[1].low,boundary.field[1].high,curr_node->cuts[1]);
        exit(1);
    }
    if (CheckPortBounds(boundary.field[2]))
    {
        printf("Error: get_bound bounds check for 2 failed\n");
        printf("[%u - %u] => [%u - %u] @ %d\n",
            curr_node->boundary.field[2].low,curr_node->boundary.field[2].high,
            boundary.field[2].low,boundary.field[2].high,curr_node->cuts[2]);
        exit(1);
    }
    if (CheckPortBounds(boundary.field[3]))
    {
        printf("Error: get_bound bounds check for 3 failed\n");
        printf("[%u - %u] => [%u - %u] @ %d\n",
            curr_node->boundary.field[3].low,curr_node->boundary.field[3].high,
            boundary.field[3].low,boundary.field[3].high,curr_node->cuts[3]);
        exit(1);
    }
    if (CheckProtoBounds(boundary.field[4]))
    {
        printf("Error: get_bound bounds check for 4 failed\n");
        printf("[%u - %u] => [%u - %u] @ %d\n",
            curr_node->boundary.field[4].low,curr_node->boundary.field[4].high,
            boundary.field[4].low,boundary.field[4].high,curr_node->cuts[4]);
        exit(1);
    }
    return boundary;
}

void ComputeCutoffs()
{
    if (binningON == 0)
    {
        Cutoffs[0] = numrules;
    }
    for (int i = 0;i < NUM_JUNK;i++)
    {
        Cutoffs[i] = numrules * Percents[i] / 100;
#ifndef KUN_TEST
        printf("Cutoffs[%d] = %u\n",i,Cutoffs[i]); //Kun: may not include all rules!!!
#endif
    }
}


void InitStats(int No_Rules) {
    p_record = new TreeStat;
    p_record->Id = treecount++;
    p_record->No_Rules = No_Rules;

    Max_Depth = 0;
    Max_Levels = 0;
    Max_Cuts = 0;
    Max_Access64Bit = 0;
    Max_Access128Bit = 0;
    Rules_at_the_Leaf = 0;
    Rules_along_path = 0;
    Max_WorklistSize = 0;
    Node_Count = 0;
    Problematic_Node_Count = 0;
    NonLeaf_Node_Count = 0;
    Compressed_NonLeaf_Node_Count = 0;
    Uncompressed_NonLeaf_Node_Count = 0;
    Total_Array_Size = 0;
    Total_Rule_Size = 0;
    Total_Rules_Moved_Up = 0;
    Max_Degree = 0;
    Avg_Degree = 0;
    Num_Partitions = 0;

    interval_per_node.clear();
    cuts_per_node.clear();
}

void NodeStats(node *curr_node)
{

    if (curr_node->problematic == 1)
        Problematic_Node_Count++;


    Num_Partitions = curr_node->cuts[0] * curr_node->cuts[1] * curr_node->cuts[2] * curr_node->cuts[3] * curr_node->cuts[4];

    int mcuts = curr_node->cuts[0];
    for (int di = 1;di < CLASSIFIER_FIELDS;di++)
        if (curr_node->cuts[di] > mcuts)
            mcuts = curr_node->cuts[di];

    if (mcuts > Max_Cuts)
        Max_Cuts = mcuts;

    // checks
    if ( curr_node->classifier.size() > bucketSize &&
            curr_node->children.size() == 0 && curr_node->problematic == 0) {
        printf("Error: This node is not cut further!\n");
        printf("\tIt has %u rules!\n", (unsigned int)curr_node->classifier.size());
        printf("\tactual-children: %u\n", (unsigned int)curr_node->actual_children.size());
        PrintNode(curr_node);
        exit(1);
    }

    if (curr_node->problematic == 1 && curr_node->classifier.size() > TOO_MUCH)
    {
        printf("Error: This problematic node has %d rules!\n",(int)curr_node->classifier.size());
        // Edit: JED 13/10/2014
        // Allow nodes to stop with too many children if they can't split further
        ////exit(1);
    }

    if (Num_Partitions != (int)curr_node->children.size()
            && curr_node->children.size() != 0 && compressionON == 0)
    {
        printf("Error: num children != partitions!(%d != %d)\n",(int)curr_node->children.size(),(int)Num_Partitions);
        exit(1);
    }

    for (int i = 0;i < CLASSIFIER_FIELDS;++i)
        if (IsPowerOfTwo(curr_node->cuts[i]) == false && !fineOn)
        {
            printf("Error: ncuts[%d] = %d is not a power of 2!\n",i,curr_node->cuts[i]);
            exit(1);
        }

    // check the node's bounds
    if (CheckIPBounds(curr_node->boundary.field[0]))
    {
        printf("Error: NodeStat bounds check for 0 failed\n");
        exit(1);
    }
    if (CheckIPBounds(curr_node->boundary.field[1]))
    {
        printf("Error: NodeStat bounds check for 1 failed\n");
        exit(1);
    }
    if (CheckPortBounds(curr_node->boundary.field[2]))
    {
        printf("Error: NodeStat bounds check for 2 failed\n");
        exit(1);
    }
    if (CheckPortBounds(curr_node->boundary.field[3]))
    {
        printf("Error: NodeStat bounds check for 3 failed\n");
        exit(1);
    }
    if (CheckProtoBounds(curr_node->boundary.field[4]))
    {
        printf("Error: NodeStat bounds check for 4 failed\n");
        exit(1);
    }

    // stats
    Node_Count++;

    if (curr_node->children.size() != 0)
    {
        NonLeaf_Node_Count++;
        if (curr_node->is_compressed == true)
            Compressed_NonLeaf_Node_Count++;
        else
            Uncompressed_NonLeaf_Node_Count++;
    }

    if (curr_node->is_compressed == true && curr_node->children.size() == 0)
    {
        printf("Error: How the heck is leaf node compressed, exiting..\n");
        exit(1);
    }


    int Actual_Curr_Level = curr_node->depth;


    if (curr_node->is_compressed != true && !fineOn)
        Actual_Curr_Level++;

    if (Actual_Curr_Level > Max_Levels)
    {

        Max_Levels = Actual_Curr_Level;
        //printf("[Tree %d] currently at level %d ...with %d children\n",p_record->Id,Max_Levels,curr_node->children.size());
        if (Max_Levels > MAX_ALLOWED_LEVELS)
        {
            printf("Error: [Tree %d] more that %d levels!\n",
                    p_record->Id,MAX_ALLOWED_LEVELS);
            exit(1);
        }

        if (curr_node->children.empty())
            Rules_at_the_Leaf = curr_node->classifier.size();

        if (curr_node->node_has_rule == 1)
            Rules_along_path++;
    }

    int depth = curr_node->depth + (curr_node->children.empty() ? curr_node->classifier.size() : 0);
    int cost64 = curr_node->depth * INTERNAL_COST_64 +
        (curr_node->children.empty() ? curr_node->classifier.size() * RULE_COST_64 : 0);
    int cost128 = curr_node->depth * INTERNAL_COST_128 +
        (curr_node->children.empty() ? curr_node->classifier.size() * RULE_COST_128 : 0);


    if (depth > Max_Depth)
        Max_Depth = depth;
    if (cost64 > Max_Access64Bit) {
        Max_Access64Bit = cost64;
    }
    if (cost128 > Max_Access128Bit)
        Max_Access128Bit = cost128;

    if (curr_node->children.size() != 0)
        Total_Array_Size += curr_node->children.size();
    else
    {
        Total_Rule_Size += curr_node->classifier.size();
    }

    if ((int)curr_node->children.size() > Max_Degree)
        Max_Degree = curr_node->children.size();

    Avg_Degree += curr_node->children.size();

    // intervals per node
    if (curr_node->children.size() != 0 && curr_node->is_compressed == true)
    {
        map<unsigned,uint32_t>::iterator iter = interval_per_node.find(curr_node->children.size());
        if (iter != interval_per_node.end())
        {
            uint32_t count = iter->second;
            count++;
            interval_per_node[curr_node->children.size()] = count;
        }
        else
        {
            interval_per_node[curr_node->children.size()] = 1;
        }
    }

    // cuts per node
    if (curr_node->children.size() != 0)
    {
        map<unsigned,uint32_t>::iterator iter = cuts_per_node.find(Num_Partitions);
        if (iter != cuts_per_node.end())
        {
            uint32_t count = iter->second;
            count++;
            cuts_per_node[Num_Partitions] = count;
        }
        else
        {
            cuts_per_node[Num_Partitions] = 1;
        }
    }

}

void InterValHist(map <unsigned,uint32_t> interval_per_node)
{
    for (map<unsigned,uint32_t>::iterator iter = interval_per_node.begin();
            iter != interval_per_node.end();++iter)
    {
        printf("I %u,%u\n",(*iter).first,(*iter).second);
    }
}
void CutsHist(map <unsigned,uint32_t> cuts_per_node)
{
    for (map<unsigned,uint32_t>::iterator iter = cuts_per_node.begin();
            iter != cuts_per_node.end();++iter)
    {
        printf("C %u,%u\n",(*iter).first,(*iter).second);
    }
}

void PrintStatRecord(TreeStat *p_record)
{
    printf("******************************************\n");
    printf("Tree: %d\n",p_record->Id);
    printf("******************************************\n");
    printf("Rules: %d\n",p_record->No_Rules);
    printf("Cost64: %d\n",p_record->Max_Access64Bit);
    printf("Cost128: %d\n",p_record->Max_Access128Bit);
    printf("Depth: %d\n",p_record->Max_Depth);
    printf("Levels: %d\n",p_record->Max_Levels);
    printf("Cuts: %d\n",p_record->Max_Cuts);
    printf("Rules_at_the_Leaf: %d\n",p_record->Rules_at_the_Leaf);
    printf("Rules_along_path: %d\n",p_record->Rules_along_path);
    printf("Total_Rule_Size: %u\n",p_record->Total_Rule_Size);
    printf("Total_Rules_Moved_Up: %u\n",p_record->Total_Rules_Moved_Up);
    printf("Total_Array_Size: %u\n",p_record->Total_Array_Size);
    printf("Node_Count: %u\n",p_record->Node_Count);
    printf("Problematic_Node_Count: %u\n",p_record->Problematic_Node_Count);
    printf("NonLeaf_Node_Count: %u\n",p_record->NonLeaf_Node_Count);
    printf("Compressed_NonLeaf_Node_Count: %u\n",p_record->Compressed_NonLeaf_Node_Count);
    printf("Uncompressed_NonLeaf_Node_Count: %u\n",p_record->Uncompressed_NonLeaf_Node_Count);
    printf("------------------------------------------\n");
    printf("ruleptr_memory: %u\n",p_record->ruleptr_memory);;
    printf("array_memory: %u\n",p_record->array_memory);;
    printf("leaf_node_memory: %u\n",p_record->leaf_node_memory);;
    printf("compressed_int_node_memory: %u\n",p_record->compressed_int_node_memory);;
    printf("uncompressed_int_node_memory: %u\n",p_record->uncompressed_int_node_memory);;
    printf("total_memory: %u\n",p_record->total_memory);;
    printf("total_memory_in_KB: %u\n",p_record->total_memory_in_KB);;
    printf("------------------------------------------\n");
    InterValHist(p_record->interval_per_node);
    printf("------------------------------------------\n");
    CutsHist(p_record->cuts_per_node);
}

void PrintStats()
{
    uint32_t OVERALL_MEMORY = 0;
    int OVERALL_DEPTH = 0;
    int OVERALL_LEVELS = 0;
    int No_Rules = 0;
    int sumCost64 = 0;
    int sumCost128 = 0;

    for (list<TreeStat*>::iterator iter = Statistics.begin();
                iter != Statistics.end();iter++)
    {
        PrintStatRecord(*iter);
        OVERALL_MEMORY += (*iter)->total_memory_in_KB;
        No_Rules += (*iter)->No_Rules;
        if ((*iter)->Max_Depth > OVERALL_DEPTH)
            OVERALL_DEPTH = (*iter)->Max_Depth;
        if ((*iter)->Max_Levels > OVERALL_LEVELS)
            OVERALL_LEVELS = (*iter)->Max_Levels;
        sumCost64 += (*iter)->Max_Access64Bit;
        sumCost128 += (*iter)->Max_Access128Bit;
    }
    printf("******************************************\n");
    printf("OVERALL_MEMORY: %u\n",OVERALL_MEMORY);
    printf("OVERALL_DEPTH: %d\n",OVERALL_DEPTH);
    printf("OVERALL_LEVELS: %d\n",OVERALL_LEVELS);

    printf("SUM_COST64: %d\n", sumCost64);
    printf("SUM_COST128: %d\n", sumCost128);

    printf("Update Reads: %d\n", updateReads);
    printf("Update Writes: %d\n", updateWrites);

    // some final checks
    if (No_Rules != (int)classifier.size())
    {
        printf("Error: Some rules got dropped while binning!\n");
        ////exit(1);
    }

}


void PrintTree(node* currNode)
{
    printf("Node: depth %u, %u children\n", (unsigned int)currNode->depth, (unsigned int)currNode->children.size());
    for (list<node*>::iterator iter = currNode->children.begin();
        iter != currNode->children.end(); iter++)
    {
        PrintTree(*iter);
    }
}

void StatNode(node* currNode)
{
    //printf("node size: %u rules %u children %u actual\n", currNode->classifier.size(), currNode->children.size(), currNode->actual_children.size());
    NodeStats(currNode);

    if (currNode->count > 0)
    {
        printf("Node has been visited %u times before!\n", currNode->count);
    }

    currNode->count++;

    for (list<node*>::iterator iter = currNode->actual_children.begin();
            iter != currNode->actual_children.end(); iter++)
    {
        StatNode(*iter);
    }
}

void RecordTreeStats()
{
    p_record->Max_Depth = Max_Depth;

    p_record->Max_Access64Bit = Max_Access64Bit;
    p_record->Max_Access128Bit = Max_Access128Bit;

    p_record->Max_Levels = Max_Levels;

    p_record->Max_Cuts = Max_Cuts;

    p_record->Rules_at_the_Leaf = Rules_at_the_Leaf;

    p_record->Rules_along_path = Rules_along_path;

    p_record->Total_Rule_Size = Total_Rule_Size;

    p_record->Total_Rules_Moved_Up = Total_Rules_Moved_Up;

    p_record->Total_Array_Size = Total_Array_Size;

    p_record->Node_Count = Node_Count;

    p_record->Problematic_Node_Count = Problematic_Node_Count;

    p_record->NonLeaf_Node_Count = NonLeaf_Node_Count;

    p_record->Compressed_NonLeaf_Node_Count = Compressed_NonLeaf_Node_Count;

    p_record->Uncompressed_NonLeaf_Node_Count = Uncompressed_NonLeaf_Node_Count;

    p_record->interval_per_node = interval_per_node;

    p_record->cuts_per_node = cuts_per_node;

    p_record->ruleptr_memory =  PTR_SIZE * Total_Rule_Size;

    p_record->array_memory = PTR_SIZE * Total_Array_Size;

    p_record->leaf_node_memory = LEAF_NODE_SIZE * (Node_Count - NonLeaf_Node_Count);

    p_record->compressed_int_node_memory = (INTERNAL_NODE_SIZE + INTERVAL_SIZE * num_intervals) *
                                                                                Compressed_NonLeaf_Node_Count;

    p_record->uncompressed_int_node_memory = INTERNAL_NODE_SIZE * Uncompressed_NonLeaf_Node_Count;

    p_record->total_memory = p_record->ruleptr_memory + p_record->array_memory + p_record->leaf_node_memory
                                                    + p_record->compressed_int_node_memory + p_record->uncompressed_int_node_memory;

    p_record->total_memory_in_KB = p_record->total_memory / 1024;

    Statistics.push_back(p_record);

}

bool DoesRuleContainPoint(matching_rule* rule, const uint32_t* pt)
{
    for (int i = 0; i < CLASSIFIER_FIELDS; i++)
    {
        if (rule->field[i].low > pt[i] || rule->field[i].high < pt[i])
        {
            return false;
        }
    }
    return true;
}

bool DoRulesIntersect(matching_rule* r1, matching_rule* r2)
{
    for (int i = 0; i < CLASSIFIER_FIELDS; i++)
    {
        if (r1->field[i].high < r2->field[i].low || r1->field[i].low > r2->field[i].high)
        {
            return false;
        }
    }
    return true;
}

bool mycomparison(matching_rule* first,matching_rule* second)
{
  return (first->priority < second->priority);
}

bool myequal(matching_rule* first,matching_rule* second)
{
  return (first->priority == second->priority);
}

bool mystatsort(TreeStat* first,TreeStat* second)
{
  if (first->Max_Depth > second->Max_Depth)
  {
    return true;
  }
  else
  {
    if (first->Max_Depth == second->Max_Depth)
    {
      if (first->total_memory > second->total_memory)
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }
}

bool mymemsort(MemBin* first,MemBin* second)
{
  if (first->Max_Depth < second->Max_Depth)
  {
    return true;
  }
  else
  {
    if (first->Max_Depth == second->Max_Depth)
    {
      if (first->total_memory < second->total_memory)
      {
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      return false;
    }
  }
}

int CheckIPBounds(range fld)
{
  if (fld.low > 0xFFFFFFFF)
  {
    printf("Error: IPRange is buggy!(%u)\n",fld.low);
    return 1;
  }
  if (fld.high > 0xFFFFFFFF)
  {
    printf("Error: IPRange is buggy!(%u)\n",fld.high);
    return 1;
  }
  if (fld.low > fld.high)
  {
    printf("Error: IPRange is buggy!(%u - %u)\n",fld.low,fld.high);
    return 1;
  }
  return 0;
}

int CheckPortBounds(range fld)
{
  if (fld.low > 0xFFFF)
  {
    printf("Error: PortRange is buggy!(%u)\n",fld.low);
    return 1;
  }
  if (fld.high > 0xFFFF)
  {
    printf("Error: PortRange is buggy!(%u)\n",fld.high);
    return 1;
  }
  if (fld.low > fld.high)
  {
    printf("Error: PortRange is buggy!(%u - %u)\n",fld.low,fld.high);
    return 1;
  }
  return 0;
}

int CheckProtoBounds(range fld)
{
  if (fld.low > 0xFF)
  {
    printf("Error: ProtoRange is buggy!(%u)\n",fld.low);
    return 1;
  }
  if (fld.high > 0xFF)
  {
    printf("Error: ProtoRange is buggy!(%u)\n",fld.high);
    return 1;
  }
  if (fld.low > fld.high)
  {
    printf("Error: ProtoRange is buggy!(%u - %u)\n",fld.low,fld.high);
    return 1;
  }
  return 0;
}

bool IsPowerOfTwo(int x)
{
  return (x & (x - 1)) == 0;
}

bool is_present(matching_rule boundary,matching_rule *rule)
{
  if ( ((rule->field[0].low  <= boundary.field[0].low  && rule->field[0].high >= boundary.field[0].low)  ||  // cuts to the left of range
        (rule->field[0].high >= boundary.field[0].high && rule->field[0].low  <= boundary.field[0].high) ||  // cuts to the right of range
        (rule->field[0].low  >= boundary.field[0].low  && rule->field[0].high <= boundary.field[0].high)) && // completely inside the range
      ((rule->field[1].low  <= boundary.field[1].low  && rule->field[1].high >= boundary.field[1].low)  ||  // cuts to the left of range
       (rule->field[1].high >= boundary.field[1].high && rule->field[1].low  <= boundary.field[1].high) ||  // cuts to the right of range
       (rule->field[1].low  >= boundary.field[1].low  && rule->field[1].high <= boundary.field[1].high)) && // completely inside the range
      ((rule->field[2].low  <= boundary.field[2].low  && rule->field[2].high >= boundary.field[2].low)  ||  // cuts to the left of range
       (rule->field[2].high >= boundary.field[2].high && rule->field[2].low  <= boundary.field[2].high) ||  // cuts to the right of range
       (rule->field[2].low  >= boundary.field[2].low  && rule->field[2].high <= boundary.field[2].high)) && // completely inside the range
      ((rule->field[3].low  <= boundary.field[3].low  && rule->field[3].high >= boundary.field[3].low)  ||  // cuts to the left of range
       (rule->field[3].high >= boundary.field[3].high && rule->field[3].low  <= boundary.field[3].high) ||  // cuts to the right of range
       (rule->field[3].low  >= boundary.field[3].low  && rule->field[3].high <= boundary.field[3].high)) && // completely inside the range
      ((rule->field[4].low  <= boundary.field[4].low  && rule->field[4].high >= boundary.field[4].low)  ||  // cuts to the left of range
       (rule->field[4].high >= boundary.field[4].high && rule->field[4].low  <= boundary.field[4].high) ||  // cuts to the right of range
       (rule->field[4].low  >= boundary.field[4].low  && rule->field[4].high <= boundary.field[4].high)) )  // completely inside the range
  {
    return true;
  }
  else
  {
    return false;
  }

}

bool is_equal(matching_rule rule1,matching_rule rule2, matching_rule boundary)
{
  int count = 0;
  range r1, r2;
  for (int i = 0;i < CLASSIFIER_FIELDS;i++)
  {
    if (rule1.field[i].low > boundary.field[i].low) {
      r1.low = rule1.field[i].low;
    } else {
      r1.low = boundary.field[i].low;
    }
    if (rule1.field[i].high < boundary.field[i].high) {
      r1.high = rule1.field[i].high;
    } else {
      r1.high = boundary.field[i].high;
    }
    if (rule2.field[i].low > boundary.field[i].low) {
      r2.low = rule2.field[i].low;
    } else {
      r2.low = boundary.field[i].low;
    }
    if (rule2.field[i].high < boundary.field[i].high) {
      r2.high = rule2.field[i].high;
    } else {
      r2.high = boundary.field[i].high;
    }
    if (r1.low <= r2.low && r1.high >= r2.high)
    {
      count++;
    }
  }

  if (count == CLASSIFIER_FIELDS)
    return true;
  else
    return false;
}

/**
 * @brief Checks whether two nodes have the same rules
 * @note Was rewritten by Alon as original version was messy
 */
int samerules(node * r1, node * r2) {

	if (r1->classifier.empty() || r2->classifier.empty()) {
		return 0;
	}
	if (r1->classifier.size() != r2->classifier.size()) {
		return 0;
	}

	// Check that each rule in r1 is also found in r2
	for (auto rule_a : r1->classifier) {
		bool found = false;
		for (auto rule_b : r2->classifier) {
			if (rule_a->priority == rule_b->priority) {
				found = true;
				break;
			}
		}
		if (!found) {
			return 0;
		}
	}
	return 1;
}

int ColorOfList(list<matching_rule*> rules, const uint32_t *pt)
{
    for (list<matching_rule*>::iterator iter = rules.begin(); iter != rules.end(); iter++)
    {
        bool isMatch = true;
        for (int d = 0; d < CLASSIFIER_FIELDS; d++)
        {
            if (pt[d] < (*iter)->field[d].low || pt[d] > (*iter)->field[d].high)
            {
                isMatch = false;
                break;
            }
        }
        if (isMatch)
        {
            return (*iter)->priority;
        }
    }

    // No match
    return -1;
}

int ColorOfTree(node* tree, const uint32_t * pt)
{
    for (int dim = 0; dim < CLASSIFIER_FIELDS; dim++)
    {
        if (pt[dim] < tree->boundary.field[dim].low || pt[dim] > tree->boundary.field[dim].high)
            return -1; // Out of bounds
    }

    if (tree->actual_children.size() == 0)
    {
        // Check list
        return ColorOfList(tree->classifier, pt);
    }
    else
    {
        // Check children

        int color = -1;

        for (list<node*>::iterator iter = tree->actual_children.begin();
                iter != tree->actual_children.end(); iter++)
        {
            if (DoesRuleContainPoint(&(*iter)->boundary, pt))
            {
                int c = ColorOfTree(*iter, pt);
                if (color < 0 || (c < color && c >= 0))
                    color = c;
            }
        }

        return color;
    }

    // No match
    return -1;
}

int ColorOfTrees(list<TreeDetails>& trees, const uint32_t *pt)
{
    int color = -1;

    for (list<TreeDetails>::iterator iter = trees.begin(); iter != trees.end(); iter++)
    {
        int c = ColorOfTree(iter->root, pt);
        //cout << "Found color: " << c << endl;
        if (color < 0 || (c < color && c >= 0))
            color = c;
    }

    return color;
}


void calc_dimensions_to_cut(node *curr_node,int *select_dim)
{
  int unique_elements[CLASSIFIER_FIELDS];
  double average = 0;
  //int average = 0;
  range check;
  for (int i = 0;i < CLASSIFIER_FIELDS;++i)
  {
    list <range> rangelist;
    rangelist.clear();
    for (list<matching_rule*>::iterator rule = curr_node->classifier.begin();
        rule != curr_node->classifier.end();++rule)
    {
      int found = 0;
      if ((*rule)->field[i].low > curr_node->boundary.field[i].low) {
        check.low = (*rule)->field[i].low;
      } else {
        check.low = curr_node->boundary.field[i].low;
      }
      if ((*rule)->field[i].high < curr_node->boundary.field[i].high) {
        check.high = (*rule)->field[i].high;
      } else {
        check.high = curr_node->boundary.field[i].high;
      }
      for (list <range>::iterator range = rangelist.begin();
          range != rangelist.end();++range)
      {
        if (check.low == (*range).low && check.high == (*range).high)
        {
          found = 1;
          break;
        }
      }
      if (!found)
        rangelist.push_back(check);
    }
    unique_elements[i] = rangelist.size();
      //printf("unique_elements[%d] = %d\n",i,unique_elements[i]);

  }

  int dims_cnt = 0;
  for (int i = 0;i < CLASSIFIER_FIELDS;++i)
  {
    if (curr_node->boundary.field[i].high > curr_node->boundary.field[i].low)
    {
      average += unique_elements[i];
      dims_cnt++;
    }
  }
  average = average / dims_cnt;

  int max = -1;
  for (int i = 0;i < CLASSIFIER_FIELDS;++i)
  {
    if (curr_node->boundary.field[i].high > curr_node->boundary.field[i].low)
      if (unique_elements[i] > max)
        max = unique_elements[i];
  }

  // Daly: made it so only dimensions strictly greater than the average get picked
  // Such a dimension must exist unless all are equal
  // So we detect such a case
  // This encourages slightly higher draws
  bool areEqual = true;
  for (int i = 0;i < CLASSIFIER_FIELDS;++i)
  {
    select_dim[i] = 0;
    if (unique_elements[i] != average && curr_node->boundary.field[i].high > curr_node->boundary.field[i].low)
        areEqual = false;
  }

  if (areEqual)
  {
      cout << "all are equal" << endl;
      for (int i = 0; i < CLASSIFIER_FIELDS; i++)
      {
          if (curr_node->boundary.field[i].high > curr_node->boundary.field[i].low)
          {
              cout << i << " " << unique_elements[i] << endl;
              select_dim[i] = 1;
              return;
          }
      }
  }

  int dim_count = 0;
  for (int i = 0;i < CLASSIFIER_FIELDS;++i)
  {
    if (curr_node->boundary.field[i].high > curr_node->boundary.field[i].low)
    {
      if (hypercuts)
      {
        if (unique_elements[i] > average) // Daly: changed to strictly greater
        {
          select_dim[i] = 1;
          dim_count++;
          // don't cut on more than 2 dimensions
          if (dim_count == 2)
            break;
        }
      }
      else
      {
        if (unique_elements[i] == max)
        {
          select_dim[i] = 1;
          break;
        }
      }
    }
  }

}

list<node*> calc_num_cuts_1D(node *root,int dim) {
  node *curr_node;

  int nump = 0;

  int spmf = int(floor(root->classifier.size() * spfac));
  int sm = 0;

  int prev_depth = -1;

  int index = 0;

  list<node*> output;
  node* top = new node(*root);

  output.push_back(top);

  while (!output.empty()) {

    curr_node = output.front();

    if (prev_depth != curr_node->depth)
    {
      if (sm < spmf)
      {
        nump++;
        sm = 1 << nump;
        prev_depth = curr_node->depth;
        index = 0;
      }
      else
        break;
    }

    for (int k = 0;k < 2;k++)
    {
      curr_node->cuts[dim] = 2;

      node* child = new node;
      child->depth = curr_node->depth + 1;

      child->boundary = get_bound(curr_node, dim, k);
      child->children.clear();

      for (int i=0;i < CLASSIFIER_FIELDS;i++)
        child->cuts[i] = 1;

      for (list <matching_rule*>::iterator rule = curr_node->classifier.begin();rule != curr_node->classifier.end();
          ++rule)
      {
        if (is_present(child->boundary,(*rule)) == true)
        {
          child->classifier.push_back(*rule);
        }
      }

      child->index = index++;

      child->is_compressed = false;

      sm += child->classifier.size();
      if (child->boundary.field[0].low == child->boundary.field[0].high &&
          child->boundary.field[1].low == child->boundary.field[1].high &&
          child->boundary.field[2].low == child->boundary.field[2].high &&
          child->boundary.field[3].low == child->boundary.field[3].high &&
          child->boundary.field[4].low == child->boundary.field[4].high )
        if (child->classifier.size() > 1)
        {
          printf("Error: Box 1X1X1X1X1 cannot contain more than 1 rule!\n");
          exit(1);
        }

      output.push_back(child);

    }

    output.pop_front();

    delete curr_node;
  }

  root->cuts[dim] = 1 << nump;
  return output;
}

list<node*> calc_num_cuts_2D(node *root,int *dim) {
  root->row = 0;
  root->column = 0;

  node *curr_node;

  int nump[2];
  for (int i=0;i<2;i++)
    nump[i] = 0;

  int spmf = int(floor(root->classifier.size() * spfac));
  int sm = 0;

  int prev_depth = -1;

  unsigned short chosen = 1;

  list<node*> child_list;
  node* top = new node(*root);

  child_list.push_back(top);

  while (!child_list.empty())
  {
    curr_node = child_list.front();

    if (prev_depth != curr_node->depth)
    {
      if (sm < spmf)
      {
        chosen = chosen ^ 1;
        nump[chosen]++;
        sm = 1 << (nump[0] + nump[1]);
        prev_depth = curr_node->depth;
      }
      else
        break;
    }

    for (int k = 0;k < 2;k++)
    {
      curr_node->cuts[dim[chosen]] = 2;
      curr_node->cuts[dim[chosen ^ 1]] = 1;

      node* child = new node;
      child->depth = curr_node->depth + 1;
      child->boundary = get_bound(curr_node, dim[chosen], k);
      child->children.clear();

      for (int i=0;i < CLASSIFIER_FIELDS;i++)
        child->cuts[i] = 1;

      for (list <matching_rule*>::iterator rule = curr_node->classifier.begin();rule != curr_node->classifier.end();
          ++rule)
      {
        if (is_present(child->boundary,*rule) == true)
        {
          child->classifier.push_back(*rule);
        }
      }

      if (chosen == 0)
      {
        child->row = (2 * curr_node->row) + k;
        child->column = curr_node->column;
      }
      else
      {
        child->column = (2 * curr_node->column) + k;
        child->row = curr_node->row;
      }

      child->is_compressed = false;

      //    printf("[%d,%d] ",child->row,child->column);

      sm += child->classifier.size();
      if (child->boundary.field[0].low == child->boundary.field[0].high &&
          child->boundary.field[1].low == child->boundary.field[1].high &&
          child->boundary.field[2].low == child->boundary.field[2].high &&
          child->boundary.field[3].low == child->boundary.field[3].high &&
          child->boundary.field[4].low == child->boundary.field[4].high )
        if (child->classifier.size() > 1)
        {
          printf("Box 1X1X1X1X1 cannot contain more than 1 rule!\n");
          exit(1);
        }
      child_list.push_back(child);

    }

    child_list.pop_front();

	delete curr_node;
  }

  root->cuts[dim[0]] = ( 1 << nump[0]);
  root->cuts[dim[1]] = ( 1 << nump[1]);

  if (compressionON) {
	// Linearize children
	for (auto item : child_list) {
	item->index = item->column * root->cuts[dim[0]] + item->row;
	}

	// Sort children
    auto comapre_method = [] (node* a, node* b) {
    	return a->index < b->index;
    };
    child_list.sort(comapre_method);
  }
  return child_list;
}

list<node*> calc_cuts(node *curr_node) {
	int select_dim[CLASSIFIER_FIELDS];
	int chosen_dim[2];
	int chosen_cnt = 0;

	calc_dimensions_to_cut(curr_node,select_dim);

	for (int i = 0;i < CLASSIFIER_FIELDS;++i)
	if (select_dim[i])
	  chosen_dim[chosen_cnt++] = i;

	if (chosen_cnt > 2) {
	printf("Error: More than 2 dimensions are cut!\n");
	}

	if (chosen_cnt > 1 && hypercuts == 0) {
	printf("Error: Hicut: More than 1 dimensions are cut!\n");
	}

	if (chosen_cnt == 0) {
	printf("Error: Atleast 1 dimension needs to be cut!\n");
	}

	list<node*> output;
	if (chosen_cnt == 2 && !fineOn) {
		output = calc_num_cuts_2D(curr_node,chosen_dim);
	} else if (chosen_cnt == 2 && fineOn) {
		output = CalcEquiCuts2D(curr_node, chosen_dim);
	} else if (chosen_cnt == 1 && !fineOn) {
		output = calc_num_cuts_1D(curr_node,chosen_dim[0]);
	} else if (chosen_cnt == 1 && fineOn) {
		output = CalcMultiEquiCuts1D(curr_node, chosen_dim[0]);
	}

	return output;
}

node* CreateRootNode(list<matching_rule*> p_classifier)
{
    // create a currNode node, put all rules in it.
    node* currNode = new node;
    currNode->depth = 1;
    for (int i = 0;i < CLASSIFIER_FIELDS;++i) {
        currNode->boundary.field[i].low = 0;
        if (i < 2)
            currNode->boundary.field[i].high = 0xffffffff;
        else if (i < 4)
            currNode->boundary.field[i].high = 0xffff;
        else
            currNode->boundary.field[i].high = 0xff;
    }
    currNode->children.clear();
    for (int i=0;i < CLASSIFIER_FIELDS;i++)
        currNode->cuts[i] = 1;

    for (list <matching_rule*>::iterator i = p_classifier.begin();i != p_classifier.end();++i)
    {
        currNode->classifier.push_back((*i));
    }

    int count = (int)currNode->classifier.size();
    if (count != (int)currNode->classifier.size()) {
        cout << "Redundant rules removed!" << endl;
    }

    return currNode;
}

void create_tree(list <matching_rule*> & p_classifier) {
	list <node*> worklist;

	root = CreateRootNode(p_classifier);
	if (root->classifier.size() > bucketSize) {
		worklist.push_back(root);
	} else {
		root->problematic = 0;
		NodeStats(root);
	}

	messagef("starting create_tree");

	while (!worklist.empty()) {

		node *curr_node = worklist.back();
		worklist.pop_back();

		list<node*> topush = CutNode(curr_node);

		for (auto item : topush) {
			if (item->classifier.size() > bucketSize)
			{
				if (item->boundary.field[0].low == curr_node->boundary.field[0].low
						&& item->boundary.field[0].high == curr_node->boundary.field[0].high
						&& item->boundary.field[1].low == curr_node->boundary.field[1].low
						&& item->boundary.field[1].high == curr_node->boundary.field[1].high
						&& item->boundary.field[2].low == curr_node->boundary.field[2].low
						&& item->boundary.field[2].high == curr_node->boundary.field[2].high
						&& item->boundary.field[3].low == curr_node->boundary.field[3].low
						&& item->boundary.field[3].high == curr_node->boundary.field[3].high
						&& item->boundary.field[4].low == curr_node->boundary.field[4].low
						&& item->boundary.field[4].high == curr_node->boundary.field[4].high
						&& item->classifier.size() == curr_node->classifier.size())
				{
					printf("Warning: parent and child are identical with %d rules!\n",(int)curr_node->classifier.size());
					item->problematic = 1;
					NodeStats(item);
				} else {

					worklist.push_back(item);
					if (worklist.size() > Max_WorklistSize) {
						Max_WorklistSize = (int)worklist.size();
						if (Max_WorklistSize % 100 == 0)
							printf("Worklist.size() = %u\n",Max_WorklistSize);
					}
				}

			} else {
				if (!item->classifier.empty()) {
					item->problematic = 0;
					NodeStats(item);
				}
			}
		}

		curr_node->problematic = 0;
		NodeStats(curr_node);
	}
}

list<node*> CutNode(node* curr_node) {

    if (hypercuts) {
        regionCompaction(curr_node);
    }

    list<node*> new_nodes = calc_cuts(curr_node);
	curr_node->children.insert(curr_node->children.end(), new_nodes.begin(), new_nodes.end());

    for (auto item : curr_node->children) {
    	item->depth = curr_node->depth + 1;
        for (int i=0;i<CLASSIFIER_FIELDS;i++) {
        	item->cuts[i] = 1;
        }
    }

    if (compressionON) {

        // backup the number of children, in case compression can't fit in!
        list <node*> original_children;
        for (auto item : curr_node->children) {
            // create a new node and make a copy
            node *child_copy = new node(*item);
            original_children.push_back(child_copy);
        }

        curr_node->is_compressed = NodeCompress(curr_node->children);

        // The compression was not successful, roll-back
        if (curr_node->children.size() > num_intervals) {

			// Remove compressed nodes
			for (auto node : curr_node->children) {
				delete node;
			}

        	curr_node->children = original_children;
            curr_node->is_compressed = false;
        }
        // The compression was successful, no backup is required
        else {
            for (auto node : original_children) {
                delete node;
            }
        }
    }

    // HEURISTIC 1 - create a list of nodes that should actually exist - both leaf and non-leaf
	new_nodes = merge_children(curr_node);
    curr_node->actual_children.insert(curr_node->actual_children.end(), new_nodes.begin(), new_nodes.end());

    return curr_node->actual_children;
}


EffiCuts::EffiCuts(uint32_t binth) :
		_num_of_rules(0), _binth(binth), _size(0), _build_time(0), rule_db(nullptr) {}

EffiCuts::~EffiCuts() {
	delete[] rule_db;
}

/**
 * @brief Build the classifier data structure
 * @returns 1 On success, 0 on fail
 * @note  Merged by Alon
 */
int EffiCuts::build(const std::list<openflow_rule>& rules) {

	// Convert the rule database to classic 5-tuple rule array
	uint32_t num_of_rules = rules.size(), counter = 0;
	matching_rule* rule_db = new matching_rule[num_of_rules];
	for (auto rule : rules) {
		rule_db[counter++] = rule.convert_to_five_tuple();
	}

	// Hard-coded arguments
	spfac=8;
	hypercuts=1;
	compressionON=1;
	binningON=1;
	mergingON=1;
	fineOn=1;

	// Tuned arguments
	bucketSize = _binth;
	_num_of_rules = num_of_rules;

	// Perform checks
	bool ok = 1;
	if(bucketSize <= 0){
		warning("bucketSize should be > 0\n");
		ok = 0;
	}
	if(spfac < 0){
		warning("space factor should be >= 0\n");
		ok = 0;
	}
	if (compressionON > 1) {
		warning("c can be only 0 - no compress, 1 - linear\n");
		ok = 0;
	}
	if (binningON > 2) {
		warning("g can be only 0 - no binning, 1 - binning 2 - static\n");
		ok = 0;
	}
	if (hypercuts > 1) {
		warning("m can be only 0 - hicut, 1 - hypercut\n");
		ok = 0;
	}
	if(num_intervals < 0){
		warning("num_intervals should be >= 0\n");
		ok = 0;
	}
	if (!ok) {
		return 0;
	}

	this->rule_db = rule_db;

	// What is this, I don't know
	setNumReps(1);

	// Convert from matching rule to rule_pc (EffiCuts representation)
	for (uint32_t i=0; i<num_of_rules; ++i) {
		matching_rule rule;
		rule.priority = rule_db[i].priority;
		for (int f=0; f<CLASSIFIER_FIELDS; ++f) {
			rule.field[f].low = rule_db[i].field[f].low;
			rule.field[f].high = rule_db[i].field[f].high;
		}
		classifier.push_back(rule);
	}

	// Record building time
	struct timespec start_time, end_time;
	clock_gettime(CLOCK_MONOTONIC, &start_time);

	// Start computing
	ComputeCutoffs();

	if (binningON == 1)
	{
		messagef("bin rules for classifier");
		binRules(classifier);

		messagef("merge rules");
		if (mergingON == 1)
			MergeTrees();

		messagef("creating tree for big rules");
        for (int i = 0; i < 5; i++) {
			if (!(bigrules[i].empty())) {
				InitStats(bigrules[i].size());
				create_tree(bigrules[i]);
				//Store tree, added by kun
				TreeDetails details;
				details.root = root;
				efficuts_trees.push_back(details);
				RecordTreeStats();
				bigrules[i].clear();
			}
		}

		messagef("creating tree for kindabig rules");
		for (int j = 0; j < 10; j++) {
			if (!(kindabigrules[j].empty())) {
				InitStats(kindabigrules[j].size());
				create_tree(kindabigrules[j]);
				//Store tree, added by kun
				TreeDetails details;
				details.root = root;
				efficuts_trees.push_back(details);
				RecordTreeStats();
				kindabigrules[j].clear();
			}
		}

		messagef("creating tree for medium rules");
		for (int k = 0; k < 10; k++) {
			if (!(mediumrules[k].empty())) {
				InitStats(mediumrules[k].size());
				create_tree(mediumrules[k]);
				//Store tree, added by kun
				TreeDetails details;
				details.root = root;
				efficuts_trees.push_back(details);
				RecordTreeStats();
				mediumrules[k].clear();
			}
		}

		messagef("creating tree for little rules");
		for (int l = 0; l < 5; l++) {
			if (!(littlerules[l].empty())) {
				InitStats(littlerules[l].size());
				create_tree(littlerules[l]);
				//Store tree, added by kun
				TreeDetails details;
				details.root = root;
				efficuts_trees.push_back(details);
				RecordTreeStats();
				mediumrules[l].clear();
			}
		}

		messagef("creating tree for small rules");
		if (!(smallrules.empty())) {
			InitStats(smallrules.size());
			create_tree(smallrules);
			//Store tree, added by kun
			TreeDetails details;
			details.root = root;
			efficuts_trees.push_back(details);
			RecordTreeStats();
			smallrules.clear();
		}

	}

	clock_gettime(CLOCK_MONOTONIC, &end_time);

	// Count size
    _size = 0;
    for (auto iter : Statistics) {
    	_size += iter->total_memory;
    }

    // Copy trees to local copy
    _trees = efficuts_trees;

    // Calculate build time
    _build_time = (end_time.tv_sec - start_time.tv_sec) * 1e3 + (end_time.tv_nsec - start_time.tv_nsec) / 1e6;
    return 1;
}

/**
 * @brief Start a synchronous process of classification an input packet.
 * @param header An array of 32bit integers according to the number of supported fields.
 * @returns The matching rule action/priority (or 0xffffffff if not found)
 * @note Added by Alon Rashelbach
 */
uint32_t EffiCuts::classify_sync(const uint32_t* header, int priority) {
	return ColorOfTrees(_trees, header);
}

/**
 * @brief Start an asynchronous process of classification for an input packet.
 * @param header An array of 32bit integers according to the number of supported fields.
 * @returns A unique id for the packet
 * @note Added by Alon Rashelbach
 */
uint32_t EffiCuts::classify_async(const uint32_t* header, int priority) {
	int result = -1;
	uint32_t packet_id = 0xffffffff;

	// Lookup only valid packets
	if (header != nullptr) {
		result = classify_sync(header, priority);
		packet_id = _packet_counter++;
	}

	// Broadcast
	for(auto it : _listeners) {
		it->on_new_result(packet_id, result ,result, _additional_args);
	}
	return packet_id;
}

/**
 * @brief Starts the performance measurement of this
 * @note Added by Alon Rashelbach
 */
void EffiCuts::start_performance_measurement() {
	clock_gettime(CLOCK_MONOTONIC, &perf_start_time);
}

/**
 * @brief Stops the performance measurement of this
 * @note Added by Alon Rashelbach
 */
void EffiCuts::stop_performance_measurement() {
	clock_gettime(CLOCK_MONOTONIC, &perf_end_time);
}

/**
 * @brief Prints debug information
  * @param verbose Set the verbosity level of printing
 * @note Added by Alon Rashelbach
 */
void EffiCuts::print(uint32_t verbose) const {
	// Measure performance
	uint32_t total_usec = (perf_end_time.tv_sec * 1e6 + perf_end_time.tv_nsec / 1e3) -
							  (perf_start_time.tv_sec * 1e6 + perf_start_time.tv_nsec / 1e3);

	messagef("Performance: total time %u usec. Average time: %.3f usec per packet.", total_usec, (double)total_usec / _packet_counter);
}

/**
 * @brief Packs this to byte array
 * @returns An object-packer with the binary data
 * @note Added by Alon Rashelbach
 */
ObjectPacker EffiCuts::pack() const {
	ObjectPacker output;
	uint32_t item;

	// Write global parameters
	output << CLASSIFIER_FIELDS;
	output << this->_num_of_rules;
	output << this->_binth;
	output << this->_size;
	output << this->_build_time;

	// We wish to convert node/rule pointers to indices
	// for fast packing / unpacking
	map<node*, uint32_t> node_ptr_to_index;
	queue<node*> nodes_to_write;
	map<matching_rule*, uint32_t> rule_ptr_to_index;
	queue<matching_rule*> rules_to_write;

	// For each element in _trees
	for (auto it : _trees) {

		// Go over all nodes, index node & rule pointers
		queue<node*> node_queue;

		node_queue.push(it.root);

		while (!node_queue.empty()) {
			// Get current node
			node* current_node = node_queue.front();
			node_queue.pop();

			// In case the current node is not indexed
			if (node_ptr_to_index.find(current_node) == node_ptr_to_index.end()) {
				// index current node
				node_ptr_to_index[current_node] = nodes_to_write.size();
				nodes_to_write.push(current_node);
			}

			// Add the children of this to queue
		    for (auto node_it : current_node->actual_children) {
				node_queue.push(node_it);
			}

			// index the rules of this
			for (auto rule_it : current_node->classifier) {
				// In case the current rule is not indexed
				if (rule_ptr_to_index.find(rule_it) == rule_ptr_to_index.end()) {
					// index current node
					rule_ptr_to_index[rule_it] = rules_to_write.size();
					rules_to_write.push(rule_it);
				}
			}
		}

	}

	// Write all rules
	output << (uint32_t)rules_to_write.size();
	while (!rules_to_write.empty()) {
		// Get current rule
		matching_rule* current_rule = rules_to_write.front();
		rules_to_write.pop();

		output << current_rule->priority;
		for (uint32_t f=0; f<CLASSIFIER_FIELDS; ++f) {
			output << current_rule->field[f].low;
			output << current_rule->field[f].high;
		}
	}

	// Write all nodes
	output << (uint32_t)nodes_to_write.size();
	while (!nodes_to_write.empty()) {
		// Get current node
		node* current_node = nodes_to_write.front();
		nodes_to_write.pop();
        infof("Nodes left to write: %lu. Buffer size: %u bytes",
                        nodes_to_write.size(), output.size());

		// Write node attributes
		output << current_node->depth;
		output << current_node->problematic;
		output << current_node->node_has_rule;
		output << current_node->row;
		output << current_node->column;
		output << current_node->index;
		output << current_node->is_compressed;
		output << current_node->count;

		// Write node's rule
		output << current_node->boundary.priority;
		for (uint32_t i=0; i<CLASSIFIER_FIELDS; ++i) {
			output << current_node->boundary.field[i].low;
			output << current_node->boundary.field[i].high;
		}

		// Write the node's cuts
		output << (uint32_t)current_node->cuts.size();
		for (auto cut_it : current_node->cuts) {
			output << cut_it;
		}

		// Write node's all other rules indices
		output << (uint32_t)current_node->classifier.size();
		for (auto rule_it : current_node->classifier) {
			output << rule_ptr_to_index[rule_it];
		}

		// Write node's children
		output << (uint32_t)current_node->actual_children.size();
		for (auto node_it : current_node->actual_children) {
			output << node_ptr_to_index[node_it];
		}
	}

	// Write all EffiCut Trees
	item = (uint32_t)_trees.size();
	output << item;
	for (auto it : _trees) {
		// Write the wide-fields
		item = it.wideFields.size();
		output << item;
		for (auto wide_it : it.wideFields) {
			item = wide_it;
			output << item;
		}

		// Write the current tree root node
		item = node_ptr_to_index[it.root];
		output << item;
	}

	return output;
}


/**
 * @brief Creates this from a memory location
 * @param object An object-reader instance
 */
void EffiCuts::load(ObjectReader& reader) {

	// Read my own properties
	reader >> _num_of_rules; // Dummy
	reader >> _num_of_rules;
	reader >> _binth;
	reader >> _size;
	reader >> _build_time;

	// Read all rules
	uint32_t num_of_rules;
	reader >> num_of_rules;

	matching_rule* rules = new matching_rule[num_of_rules];
	for (uint32_t i=0; i<num_of_rules; ++i) {
		reader >> rules[i].priority;
		for (uint32_t f=0; f<CLASSIFIER_FIELDS; ++f) {
			reader >> rules[i].field[f].low;
			reader >> rules[i].field[f].high;
		}
	}

	// Read all nodes
	uint32_t num_of_nodes;
	reader >> num_of_nodes;
	node* nodes = new node[num_of_nodes];
	for (uint32_t i=0; i<num_of_nodes; ++i) {

		// Read static info
		reader >> nodes[i].depth;
		reader >> nodes[i].problematic;
		reader >> nodes[i].node_has_rule;
		reader >> nodes[i].row;
		reader >> nodes[i].column;
		reader >> nodes[i].index;
		reader >> nodes[i].is_compressed;
		reader >> nodes[i].count;

		// Read node's boundary rule
		reader >> nodes[i].boundary.priority;
		for (uint32_t f=0; f<CLASSIFIER_FIELDS; ++f) {
			reader >> nodes[i].boundary.field[f].low;
			reader >> nodes[i].boundary.field[f].high;
		}

		// Read node's cuts
		uint32_t num_of_cuts;
		reader >> num_of_cuts;
		for (uint32_t j=0; j<num_of_cuts; ++j) {
			int item;
			reader >> item;
			nodes[i].cuts.push_back(item);
		}

		// Read node's rules
		reader >> num_of_rules;
		for (uint32_t j=0; j<num_of_rules; ++j) {
			uint32_t index;
			reader >> index;
			nodes[i].classifier.push_back(&rules[index]);
		}

		// Read node's children
		uint32_t num_of_children;
		reader >> num_of_children;
		for (uint32_t j=0; j<num_of_children; ++j) {
			uint32_t index;
			reader >> index;
			nodes[i].actual_children.push_back(&nodes[index]);
		}
	}

	// Read all EffiCuts trees
	uint32_t num_of_trees;
	reader >> num_of_trees;
	for (uint32_t i=0; i<num_of_trees; ++i) {
		TreeDetails tree;

		// Read wide fields
		uint32_t num_of_wide_fields;
		reader >> num_of_wide_fields;
		for (uint32_t j=0; j<num_of_wide_fields; ++j) {
			uint32_t item ;
			reader >> item;
			tree.wideFields.push_back((bool)item);
		}

		// Read root
		uint32_t index;
		reader >> index;
		tree.root = &nodes[index];

		// Push current tree
		_trees.push_back(tree);
	}

	// Validate all bytes were read
	assert(reader.size() == 0);
}

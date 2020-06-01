/*-----------------------------------------------------------------------------
 *
 *  Name:			hs.h
 *  Description:	hypersplit packet classification algorithm
 *  Version:		1.0 (release)
 *  Author:			Wenjun Li (wenjunli@pku.edu.cn)
 *  Updates: 		Alon Rashelbach (The Technion - Israel Institute of Technology, email: alonrs@campus.technion.ac.il)
 *  Date:			5/3/2019
 *
 *  *** Wenjun Li Comments: ***
 *  1:		Modified work based on the code of hypersplit from Yaxuan Qi
 *  2:		We simply use hypersplit for big rules in this version. However, its performance may be very bad for big rule sets(e.g., 100K).
 *               A simple way to resolve this problem: increase threshold value to reduce the number of big rules(e.g. from 2^16 to 2^24,2^28).
 *               Another more general way is to introduce a few cuttings before splitting to reduce set size, which will be given in future versions.
 *
 *  *** Alon Rashelbach Comments: ***
 *  1: My updates are commented as <AR>
 *
 *-----------------------------------------------------------------------------*/

#pragma once

#include <object_io.h> // <AR>

/* for 5-tuple classification */
#define DIM			5

//#define	DEBUG

/* for function return value */
#define SUCCESS		1
#define FAILURE		0
#define TRUE		1
#define FALSE		0

/* for bitmap */
#define MAXFILTERS	115536 /* support 64K rules */
#define WORDLENGTH	32	/* for 32-bit system */
#define BITMAPSIZE	256 /* MAXFILTERS/WORDLENGTH */

/*-----------------------------------------------------------------------------
 *  structure
 *-----------------------------------------------------------------------------*/
struct FILTER
{
	unsigned int cost;
	unsigned int dim[DIM][2];
	unsigned char act;
};

struct FILTSET
{
	unsigned int	numFilters;
	struct FILTER	filtArr[MAXFILTERS];
};


struct TPOINT
{
	unsigned int value;
	unsigned char flag;
};

struct FRAGNODE
{
	unsigned int start;
	unsigned int end;
	struct FRAGNODE *next;
};

struct FRAGLINKLIST
{
	unsigned int fragNum;
	struct FRAGNODE *head;
};

struct TFRAG
{
	unsigned int value;						// end point value
	unsigned int cbm[BITMAPSIZE];					// LENGTH * SIZE bits, CBM
};
						// released after tMT[2] is generated

struct FRAG
{
	unsigned int value;
};


struct CES
{
	unsigned short eqID;					// 2 byte, eqID;
	unsigned int  cbm[BITMAPSIZE];
	struct	CES *next;								// next CES
};

struct LISTEqS
{
	unsigned short	nCES;					// number of CES
	struct			CES *head;								// head pointer of LISTEqS
	struct			CES *rear;								// pointer to end node of LISTEqS
};


struct PNODE
{
	unsigned short	cell[65536];			// each cell stores an eqID
	struct			LISTEqS listEqs;					// list of Eqs
};

typedef struct {
	unsigned int	pri;
	unsigned int	range[DIM][2];
} rule_t;

typedef struct rule_set_s
{
	unsigned int	num; /* number of rules in the rule set */
	rule_t*	ruleList; /* rules in the set */
} rule_set_t;

typedef	struct seg_point_s
{
	unsigned int	num;	/* number of segment points */
	unsigned int*	pointList;	/* points stores here */
} seg_point_t;

typedef struct segments_s
{
	unsigned int	num;		/* number of segment */
	unsigned int	range[2];	/* segment */
} segments_t;

typedef	struct search_space_s
{
	unsigned int	range[DIM][2];
} search_space_t;

typedef struct hs_node_s
{
    unsigned int		d2s;		/* dimension to split, 2bit is enough */
    unsigned int		depth;		/* tree depth of the node, x bits supports 2^(2^x) segments */
	unsigned int		thresh;		/* thresh value to split the current segments */
	int max_priority; // <AR>
	rule_set_t* ruleset;
	struct hs_node_s*	child[2];	/* pointer to child-node, 2 for binary split */
} hs_node_t;

// The Results of HyperSplit trie
typedef struct {
	int num_rules;
	int num_childnode;
	int wst_depth;
	float avg_depth;
	int num_tree_node;
	int num_leaf_node;
	float total_mem_kb;
} hs_result_t;


class HyperSplitTrie {
public:
	unsigned int	binth;
	unsigned int	gChildCount;
	unsigned int	gNumTreeNode;
	unsigned int	gNumLeafNode;
	unsigned int	gWstDepth;
	unsigned int	gAvgDepth;
	unsigned int	gNumNonOverlappings[DIM];
	unsigned long long	gNumTotalNonOverlappings;

	hs_result_t result;

	HyperSplitTrie(rule_set_t* rule_set, uint32_t binth, hs_node_t* node);

	/* build hyper-split-tree */
	int BuildHSTree(rule_set_t* ruleset, hs_node_t* node, unsigned int depth); /* main */
};

/* lookup hyper-split-tree */
int LookupHSTree(hs_node_t* rootnode, const uint32_t* header, int priority);

/**
 * @brief Export hs_node_t array to byte array
 * @param node_array The nodes to pack
 * @returns The packed object
 * @note  <AR> Modified by Alon Rashelbach
 */
ObjectPacker hstrie_pack(const hs_node_t* node_array);

/**
 * @brief Unpacks byte array to hs node array
* @note  <AR> Created by Alon Rashelbach
 */
hs_node_t* HyperSplitTrie_unpack(ObjectReader& packer);

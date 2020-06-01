/*-----------------------------------------------------------------------------
 *
 *  Name:			hs.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <stdexcept>
#include <queue>
#include <stack>

#include <object_io.h> // <AR>
#include <hyper_split.h>

using namespace std;

HyperSplitTrie::HyperSplitTrie(rule_set_t* rule_set, uint32_t binth1, hs_node_t* node)
{
  gChildCount = 0;
  gNumTreeNode = 0;
  gNumLeafNode = 0;
  gWstDepth = 0;
  gAvgDepth = 0;
  gNumTotalNonOverlappings = 0;
  binth = binth1;

  // build hyper-split tree
  //hs_node_t rootnode;
//  printf("\n\n>>Building HyperSplit tree (%u rules, 5-tuple)", ruleset.num);

  BuildHSTree(rule_set, node, 0);
/*#ifdef	LOOKUP
  LookupHSTree(&ruleset, node);
#endif*/

  //result of hypersplit tree
  result.num_rules = rule_set->num;
  result.num_childnode = gChildCount;
  result.wst_depth = gWstDepth;
  result.avg_depth = (float) gAvgDepth/gChildCount;
  result.num_tree_node = gNumTreeNode;
  result.num_leaf_node = gNumLeafNode;
  result.total_mem_kb = (double)(gNumTreeNode*8 + gNumLeafNode*4)/1024;

/*
  printf("\n\n>>RESULTS:");
  printf("\n>>number of children: %d", gChildCount);
  printf("\n>>worst case tree depth: %d", gWstDepth);
  printf("\n>>average tree depth: %f", (float) gAvgDepth/gChildCount);
  printf("\n>>number of tree nodes:%d", gNumTreeNode);
  printf("\n>>number of leaf nodes:%d", gNumLeafNode);
  printf("\n>>total memory: %f(KB)", (double)(gNumTreeNode*8 + gNumLeafNode*4)/1024);
  printf("\n>>preprocessing time: %ld(ms)", 1000*(gEndTime.tv_sec - gStartTime.tv_sec)
			+ (gEndTime.tv_usec - gStartTime.tv_usec)/1000);
  printf("\n\n>>SUCCESS in building HyperSplit tree :-)\n\n");
*/
}



/*
 * ===  FUNCTION  ======================================================================
 *         Name:  Compare
 *  Description:  for qsort
 *     Comments:  who can make it better?
 * =====================================================================================
 */
int SegPointCompare (const void * a, const void * b)
{
  if ( *(unsigned int*)a < *(unsigned int*)b )
      return -1;
  else if ( *(unsigned int*)a == *(unsigned int*)b )
      return 0;
  else
      return 1;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  BuildHSTree
 *  Description:  building hyper-splitting tree via recursion
 * =====================================================================================
 */
int HyperSplitTrie::BuildHSTree (rule_set_t* ruleset, hs_node_t* currNode, unsigned int depth)
{
	// AR: initiate node with NULL pointers
	currNode->ruleset = NULL;
	currNode->child[0] = NULL;
	currNode->child[1] = NULL;

	// <AR> Add the maximum priority for each node
	int max_priority = 0x7fffffff;
	for (unsigned int i=0; i<ruleset->num; ++i) {
		max_priority = std::min(max_priority, (int)ruleset->ruleList[i].pri);
	}
	currNode->max_priority = max_priority;

    /*Update Leaf node*/
    if(ruleset->num <= binth)
    {
		currNode->d2s = 0;
		currNode->depth = depth;
		currNode->thresh = 0;
		currNode->child[0] = NULL;
		currNode->child[1] = NULL;
		currNode->ruleset = ruleset;

        //printf("\n>>LEAF-NODE: matching rule %d", currNode->ruleset->ruleList[0].pri);

		gChildCount ++;
		gNumLeafNode ++;
		if (gNumLeafNode % 1000000 == 0)
			printf(".");
			//printf("\n>>#%8dM leaf-node generated", gNumLeafNode/1000000);
		if (gWstDepth < depth)
			gWstDepth = depth;
		gAvgDepth += depth;
		return	SUCCESS;

    }


	/* generate segments for input filtset */
	unsigned int	dim, num, pos;
	unsigned int	maxDiffSegPts = 1;	/* maximum different segment points */
	unsigned int	d2s = 0;		/* dimension to split (with max diffseg) */
	unsigned int   thresh;
	unsigned int	range[2][2];	/* sub-space ranges for child-nodes */
	unsigned int	*segPoints[DIM];
	unsigned int	*segPointsInfo[DIM];
	unsigned int	*tempSegPoints;
	unsigned int	*tempRuleNumList;
	float			hightAvg, hightAll;
	rule_set_t		*childRuleSet = nullptr;

	// ==== Alon Extension Start ====
	// Bypass compiler warning for uninitialized variables
	for (int i=0; i<2; ++i) {
		for (int j=0; j<2; ++j) {
			range[i][j]=0;
		}
	}
	thresh=0;
	// ==== Alon Extension End   ====

#ifdef	DEBUG
	/*if (depth > 10)	exit(0);*/
	printf("\n\n>>BuildHSTree at depth=%d", depth);
	printf("\n>>Current Rules:");
	for (num = 0; num < ruleset->num; num++) {
		printf ("\n>>%5dth Rule:", ruleset->ruleList[num].pri);
		for (dim = 0; dim < DIM; dim++) {
			printf (" [%-8x, %-8x]", ruleset->ruleList[num].range[dim][0], ruleset->ruleList[num].range[dim][1]);
		}
	}
#endif /* DEBUG */

	/*Generate Segment Points from Rules*/
	for (dim = 0; dim < DIM; dim ++) {
		/* N rules have 2*N segPoints */
		segPoints[dim] = (unsigned int*) malloc ( 2 * ruleset->num * sizeof(unsigned int));
		segPointsInfo[dim] = (unsigned int*) malloc ( 2 * ruleset->num * sizeof(unsigned int));
		for (num = 0; num < ruleset->num; num ++) {
			segPoints[dim][2*num] = ruleset->ruleList[num].range[dim][0];
			segPoints[dim][2*num + 1] = ruleset->ruleList[num].range[dim][1];
		}
	}
	/*Sort the Segment Points*/
	for(dim = 0; dim < DIM; dim ++) {
		qsort(segPoints[dim], 2*ruleset->num, sizeof(unsigned int), SegPointCompare);
	}

	/*Compress the Segment Points, and select the dimension to split (d2s)*/
	tempSegPoints  = (unsigned int*) malloc(2 * ruleset->num * sizeof(unsigned int));
	hightAvg = 2*ruleset->num + 1;
	for (dim = 0; dim < DIM; dim ++) {
		unsigned int	i, j;
		unsigned int	*hightList;
		unsigned int	diffSegPts = 1; /* at least there are one different segment point */
		tempSegPoints[0] = segPoints[dim][0];
		for (num = 1; num < 2*ruleset->num; num ++) {
			if (segPoints[dim][num] != tempSegPoints[diffSegPts-1]) {
				tempSegPoints[diffSegPts] = segPoints[dim][num];
				diffSegPts ++;
			}
		}
		/*Span the segment points which is both start and end of some rules*/
		pos = 0;
		for (num = 0; num < diffSegPts; num ++) {
			unsigned int	i;
			int ifStart = 0;
			int	ifEnd	= 0;
			segPoints[dim][pos] = tempSegPoints[num];
			for (i = 0; i < ruleset->num; i ++) {
				if (ruleset->ruleList[i].range[dim][0] == tempSegPoints[num]) {
					/*printf ("\n>>rule[%d] range[0]=%x", i, ruleset->ruleList[i].range[dim][0]);*/
					/*this segment point is a start point*/
					ifStart = 1;
					break;
				}
			}
			for (i = 0; i < ruleset->num; i ++) {
				if (ruleset->ruleList[i].range[dim][1] == tempSegPoints[num]) {
					/*printf ("\n>>rule[%d] range[1]=%x", i, ruleset->ruleList[i].range[dim][1]);*/
					/* this segment point is an end point */
					ifEnd = 1;
					break;
				}
			}
			if (ifStart && ifEnd) {
				segPointsInfo[dim][pos] = 0;
				pos ++;
				segPoints[dim][pos] = tempSegPoints[num];
				segPointsInfo[dim][pos] = 1;
				pos ++;
			}
			else if (ifStart) {
				segPointsInfo[dim][pos] = 0;
				pos ++;
			}
			else {
				segPointsInfo[dim][pos] = 1;
				pos ++;
			}

		}

		/* now pos is the total number of points in the spanned segment point list */

		if (depth == 0) {
			gNumNonOverlappings[dim] = pos;
			gNumTotalNonOverlappings *= (unsigned long long) pos;
		}

#ifdef	DEBUG
		printf("\n>>dim[%d] segs: ", dim);
		for (num = 0; num < pos; num++) {
			/*if (!(num % 10))	printf("\n");*/
			printf ("%x(%u)	", segPoints[dim][num], segPointsInfo[dim][num]);
		}
#endif /* DEBUG */

		if (pos >= 3) {
			hightAll = 0;
			hightList = (unsigned int *) malloc(pos * sizeof(unsigned int));
			for (i = 0; i < pos-1; i++) {
				hightList[i] = 0;
				for (j = 0; j < ruleset->num; j++) {
					if (ruleset->ruleList[j].range[dim][0] <= segPoints[dim][i] \
							&& ruleset->ruleList[j].range[dim][1] >= segPoints[dim][i+1]) {
						hightList[i]++;
						hightAll++;
					}
				}
			}

			if (hightAvg > hightAll/(pos-1)) {	/* possible choice for d2s, pos-1 is the number of segs */
				float hightSum = 0;

				/* select current dimension */
				d2s = dim;
				hightAvg = hightAll/(pos-1);

				/* the first segment MUST belong to the left child */
				hightSum += hightList[0];
				for (num = 1; num < pos-1; num++) {  /* pos-1 >= 2; seg# = num */
					if (segPointsInfo[d2s][num] == 0)
						thresh = segPoints[d2s][num] - 1;
					else
						thresh = segPoints[d2s][num];

					if (hightSum > hightAll/2) {
						break;
					}
					hightSum += hightList[num];
				}
				/*printf("\n>>d2s=%u thresh=%x\n", d2s, thresh);*/
				range[0][0] = segPoints[d2s][0];
				range[0][1] = thresh;
				range[1][0] = thresh + 1;
				range[1][1] = segPoints[d2s][pos-1];
			}
			/* print segment list of each dim */
#ifdef	DEBUG
			printf("\n>>hightAvg=%f, hightAll=%f, segs=%d", hightAll/(pos-1), hightAll, pos-1);
			for (num = 0; num < pos-1; num++) {
				printf ("\nseg%5d[%8x, %8x](%u)	",
						num, segPoints[dim][num], segPoints[dim][num+1], hightList[num]);
			}
#endif /* DEBUG */
			free(hightList);
		} /* pos >=3 */

		//printf("\n>>298: maxDiffSegPts = %d    pos = %d", maxDiffSegPts,pos);

		if (maxDiffSegPts < pos) {
			maxDiffSegPts = pos;
			//printf("\n>>300: maxDiffSegPts = %d", maxDiffSegPts);
		}
	}
	free(tempSegPoints);


	/*Update Leaf node*/
	if (maxDiffSegPts <= 2) {
		currNode->d2s = 0;
		currNode->depth = depth;
		currNode->thresh = 0;
		currNode->child[0] = NULL;
		currNode->child[1] = NULL;
		currNode->ruleset = ruleset;

		for (dim = 0; dim < DIM; dim ++) {
			free(segPoints[dim]);
			free(segPointsInfo[dim]);
		}

		//printf("\n>>LEAF-NODE: matching rule %d", ruleset->ruleList[0].pri);

		gChildCount ++;
		gNumLeafNode ++;
		/*printf("\n>>#%8dM leaf-node generated", gNumLeafNode/1000000);*/
		if (gWstDepth < depth)
			gWstDepth = depth;
		gAvgDepth += depth;
		return	SUCCESS;
	}


	/*Update currNode*/
	/*Binary split along d2s*/


#ifdef DEBUG
	/* split info */
	printf("\n>>d2s=%u; thresh=0x%8x, range0=[%8x, %8x], range1=[%8x, %8x]",
			d2s, thresh, range[0][0], range[0][1], range[1][0], range[1][1]);
#endif /* DEBUG */


	if (range[1][0] > range[1][1]) {
		printf("\n>>maxDiffSegPts=%d  range[1][0]=%x  range[1][1]=%x",
				maxDiffSegPts, range[1][0], range[1][1]);
		printf("\n>>fuck\n"); exit(0);
	}


	for (dim = 0; dim < DIM; dim ++) {
		free(segPoints[dim]);
		free(segPointsInfo[dim]);
	}

	gNumTreeNode ++;
	currNode->d2s = d2s;
	currNode->depth = depth;
	currNode->thresh = thresh;
	currNode->child[0] = (hs_node_t *) malloc(sizeof(hs_node_t));
	//printf("\n>>line 354 ----------- currNode->d2s = %u\n", currNode->d2s);

	/*Generate left child rule list*/
	tempRuleNumList = (unsigned int*) malloc(ruleset->num * sizeof(unsigned int)); /* need to be freed */
	pos = 0;
	for (num = 0; num < ruleset->num; num++) {
		if (ruleset->ruleList[num].range[d2s][0] <= range[0][1]
		&&	ruleset->ruleList[num].range[d2s][1] >= range[0][0]) {
			tempRuleNumList[pos] = num;
			pos++;
		}
	}
	childRuleSet = (rule_set_t*) malloc(sizeof(rule_set_t));
	childRuleSet->num = pos;
	childRuleSet->ruleList = (rule_t*) malloc( childRuleSet->num * sizeof(rule_t) );
	for (num = 0; num < childRuleSet->num; num++) {
		childRuleSet->ruleList[num] = ruleset->ruleList[tempRuleNumList[num]];
		/* in d2s dim, the search space needs to be trimmed off */
		if (childRuleSet->ruleList[num].range[d2s][0] < range[0][0])
			childRuleSet->ruleList[num].range[d2s][0] = range[0][0];
		if (childRuleSet->ruleList[num].range[d2s][1] > range[0][1])
			childRuleSet->ruleList[num].range[d2s][1] = range[0][1];
	}
	free(tempRuleNumList);

	BuildHSTree(childRuleSet, currNode->child[0], depth+1);


/*#ifndef	LOOKUP
	free(currNode->child[0]);
	free(childRuleSet->ruleList);
	free(childRuleSet);
#endif*/

	/*Generate right child rule list*/
	currNode->child[1] = (hs_node_t *) malloc(sizeof(hs_node_t));
	tempRuleNumList = (unsigned int*) malloc(ruleset->num * sizeof(unsigned int)); /* need to be free */
	pos = 0;
	for (num = 0; num < ruleset->num; num++) {  // low<=end && high>=start
		if (ruleset->ruleList[num].range[d2s][0] <= range[1][1]
		&&	ruleset->ruleList[num].range[d2s][1] >= range[1][0]) {
			tempRuleNumList[pos] = num;
			pos++;
		}
	}

	childRuleSet = (rule_set_t*) malloc(sizeof(rule_set_t));
	childRuleSet->num = pos;
	childRuleSet->ruleList = (rule_t*) malloc( childRuleSet->num * sizeof(rule_t) );
	for (num = 0; num < childRuleSet->num; num++) {
		childRuleSet->ruleList[num] = ruleset->ruleList[tempRuleNumList[num]];
		/* in d2s dim, the search space needs to be trimmed off */
		if (childRuleSet->ruleList[num].range[d2s][0] < range[1][0])
			childRuleSet->ruleList[num].range[d2s][0] = range[1][0];
		if (childRuleSet->ruleList[num].range[d2s][1] > range[1][1])
			childRuleSet->ruleList[num].range[d2s][1] = range[1][1];
	}

	free(tempRuleNumList);

	BuildHSTree(childRuleSet, currNode->child[1], depth+1);
/*#ifndef	LOOKUP
	free(currNode->child[1]);
	free(childRuleSet->ruleList);
	free(childRuleSet);
#endif*/

	return	SUCCESS;
}

/*
 * ===  FUNCTION  ======================================================================
 *         Name:  LookupHSTtree
 *  Description:  test the hyper-split-tree with give 4-tuple packet
 * =====================================================================================
 */
int LookupHSTree(hs_node_t* node, const uint32_t* header, int priority) {

	int cover = 1;
	int match = 0;
	unsigned int i,k;

	while (node->child[0] != NULL) {
		// <AR> priority optimization
		if ( (priority>=0) && (priority < node->max_priority) ) return priority;

		if (header[node->d2s] <= (node->thresh)){
			node = node->child[0];
		}
		else{
			node = node->child[1];
		}
	}

	if(node != NULL){
		// <AR> priority optimization
		if ( (priority>=0) && (priority < node->max_priority) ) return priority;

		for(i = 0; i < node->ruleset->num; i++){
			cover = 1;
			for(k = 0; k < DIM; k++){
				if(node->ruleset->ruleList[i].range[k][0] > header[k] ||
				    node->ruleset->ruleList[i].range[k][1] < header[k]){
					cover = 0;
					break;
				}
			}
			if(cover == 1){
				match = 1;
				break;
			}
		}
	}


	if(match == 1){
		//printf("\n>>Matched Rule %d\n", node->ruleset->ruleList[i].pri);
		return std::min((uint32_t)priority, (uint32_t)node->ruleset->ruleList[i].pri);
	}else{
		return priority;
	}

}



/**
 * @brief Export hs_node_t array to byte array
 * @param node_array The nodes to pack
 * @returns The packed object
 * @note  <AR> Modified by Alon Rashelbach
 */
ObjectPacker hstrie_pack(const hs_node_t* node_array) {

	ObjectPacker output;

	// Do DFS of the node_array and build the output byte array
	queue<const hs_node_t*> node_queue;
	uint32_t counter = 0;

	// Pack until no available nodes
	if (node_array) {
		node_queue.push(node_array);
		while (!node_queue.empty()) {
			// Get next node
			auto current = node_queue.front();
			node_queue.pop();
			++counter;

			// Pack the node's variables
			output << current->d2s;
			output << current->depth;
			output << current->thresh;

			uint32_t rule_num = (current->ruleset != nullptr) ? current->ruleset->num : 0;
			output << rule_num;

			// Pack the rule-set, if exists
			if (current->ruleset != nullptr) {
				for (uint32_t r=0; r<current->ruleset->num; ++r) {
					// Pack the rule
					rule_t* rule = &current->ruleset->ruleList[r];
					output << rule->pri;
					// Pack the fields
					for (uint32_t d=0; d<DIM; ++d) {
						output << rule->range[d][0];
						output << rule->range[d][1];
					}
				}
			}

			// Pack the children
			for (uint32_t i=0; i<2; ++i) {
				uint32_t value;
				if (current->child[i]) {
					// Add child to queue, write its index
					value = counter + node_queue.size();
					node_queue.push(current->child[i]);
				} else {
					// No child, write 0xffffffff
					value = 0xffffffff;
				}
				output << value;
			}
		}
	}

	// Write total node count in the beginning
	output.insert(counter);

	// Pack object
	return output;
}

/**
 * @brief Unpacks byte array to hs node array
* @note  <AR> Created by Alon Rashelbach
 */
hs_node_t* HyperSplitTrie_unpack(ObjectReader& reader) {

	// Get number of nodes
	uint32_t num_of_nodes;
	reader >> num_of_nodes;

	if (num_of_nodes == 0) {
		return nullptr;
	}

	hs_node_t* output = new hs_node_t[num_of_nodes];

	// Read all nodes
	for (uint32_t i=0; i<num_of_nodes; ++i) {

		// Get current node
		hs_node_t* current = &output[i];
		reader >> current->d2s;
		reader >> current->depth;
		reader >> current->thresh;
		current->ruleset = nullptr;
		current->max_priority = -1;

		// Unpack the rule-set
		uint32_t rule_num;
		reader >> rule_num;

		if (rule_num > 0) {

			// Allocate rules
			current->ruleset = new rule_set_t;
			current->ruleset->num = rule_num;
			current->ruleset->ruleList = new rule_t[rule_num];

			for (uint32_t r=0; r<rule_num; ++r) {
				// Unpack the rule
				rule_t* rule = &current->ruleset->ruleList[r];
				reader >> rule->pri;
				// Pack the fields
				for (uint32_t d=0; d<DIM; ++d) {
					reader >> rule->range[d][0];
					reader >> rule->range[d][1];
				}
			}
		}

		// Unpack children
		for (uint32_t j=0; j<2; ++j) {
			uint32_t index;
			reader >> index;
			if (index == 0xffffffff) {
				current->child[j] = nullptr;
			} else {
				current->child[j] = &output[index];
			}
		}
	}

	// Calculate maximum priority using DFS
	stack<hs_node_t*> node_stack;
	node_stack.push(&output[0]);
	while (!node_stack.empty()) {
		hs_node_t* current =node_stack.top();
	   // In case the current node has already priority
	   if (current->max_priority > 0) {
		   node_stack.pop();
		   continue;
	   }

	   bool has_children =
			   (current->child[0] != nullptr) ||
			   (current->child[1] != nullptr);

	   // In case of leaf
	   if (!has_children) {
		   int max_priority = 0x7fffffff;
		   for (unsigned int i=0; i<current->ruleset->num; ++i) {
			   max_priority = std::min(max_priority, (int)current->ruleset->ruleList[i].pri);
		   }
		   current->max_priority = max_priority;
		   node_stack.pop();
	   }
	   // In case the child nodes of this should be processed
	   else if (current->max_priority == -1) {
		   for (uint32_t j=0; j<2; ++j) {
			   if (current->child[j] != nullptr) {
				   node_stack.push(current->child[j]);
			   }
		   }
		   current->max_priority = -2;
	   }
	   // The childs of this were already processed
	   else if (current->max_priority == -2) {
		   int child_0_prio = (current->child[0] != nullptr) ? current->child[0]->max_priority : 0;
		   int child_1_prio = (current->child[1] != nullptr) ? current->child[1]->max_priority : 0;
		   current->max_priority = (child_0_prio < child_1_prio) ? child_0_prio : child_1_prio;
		   node_stack.pop();
	   }
	}
	return output;
}


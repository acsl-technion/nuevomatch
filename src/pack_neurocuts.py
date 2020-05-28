#!/usr/bin/env python3
## MIT License
##
## Copyright (c) 2019 Alon Rashelbach
##
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.


import argparse
import struct
import pickle
import sys
import numpy as np

# Import modules
try:
	from rule_handlers import *
	from object_reader import ObjectReader
except ImportError as e:
	# Cannot find the library
	print('One of the library files was not found (%s). Did you compile the library?' % e)
	exit(1)

# Global variables
rule_dict = {} # Priority to object
node_dict = {} # Id to object


def parse_arguments():
	""" Parse script arguments """
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', required=True, type=str, help='NeuroCuts tree pkl file to load and analyze.')
	parser.add_argument('--output', required=True, type=str, help='Binary classifier filename.')
	parser.add_argument('--neurocuts', required=True, type=str, help='Original NeuroCuts code directory.')

	parser.add_argument('--remainder-indices', required=False, type=str, default=None, help='The remainder indices file for priority translation')

	return parser.parse_args()


def pack_neurocuts_tree(tree, translation_dict=None):
	""" Packs NeuroCuts tree to bytearray

	Args:
		tree: A NeuroCuts tree object
		translation_dict: A dictionary to convert from the priorities in the tree to other priorities (optional)
	"""

	global rule_dict, node_dict

	# Output byte array and pack methods
	output = bytearray()
	pack_uint32 = lambda num: output.extend(struct.pack('<I', int(max(0, min(num, 2**32-1))) ))

	# Map between integers and objects
	all_nodes = []
	all_rules = []

	# Fill object maps
	node_queue = [tree.root]
	while len(node_queue)>0:
		node = node_queue.pop()

		# Index current node
		if node.id not in node_dict.keys():
			node_dict[node.id] = len(all_nodes)
			all_nodes.append(node)

		# Index the node's rules
		for rule in node.rules:
			if rule.priority not in rule_dict.keys():
				rule_dict[rule.priority] = len(all_rules)
				all_rules.append(rule)

		# Update the queue
		for child in node.children:
			node_queue.append(child)

	# Pack build time, if exists
	build_time = 2**32-1
	if hasattr(tree, 'start_training_time') and hasattr(tree, 'stop_training_time'):
		build_time = (tree.stop_training_time - tree.start_training_time) * 1000
	else:
		print('Warning: tree does not have build time attribute')
	pack_uint32(build_time)

	# Store all rules
	pack_uint32(len(all_rules))
	for rule in all_rules:
		output.extend(pack_rule(rule, translation_dict))

	# Store all nodes
	pack_uint32(len(all_nodes))
	for node in all_nodes:
		output.extend(pack_node(node))

	return output


def pack_node(node):
	""" Packs NeuroCuts node into byte array """
	global rule_dict, node_dict

	output = bytearray()
	pack_uint32 = lambda num: output.extend(struct.pack('<I', int(max(0, min(num, 2**32-1))) ))

	# Pack attributes
	pack_uint32(node.id)
	pack_uint32(node.depth)
	pack_uint32(0 if node.is_partition() else 1)

	for i in range(5):
		pack_uint32(node.ranges[i*2+0]) # Range Low
		pack_uint32(node.ranges[i*2+1]) # Range High

	# Pack the node's rules (by index)
	pack_uint32(len(node.rules))
	for rule in node.rules:
		pack_uint32(rule_dict[rule.priority])

	# Pack ids of children
	pack_uint32(len(node.children))
	for child in node.children:
		pack_uint32(node_dict[child.id])

	return output


def pack_rule(rule, translation_dict=None):
	""" Packs NeuroCuts rule into byte array

	Args:
		rule: A rule
		translation_dict: A dictionary to convert from the priorities in the tree to other priorities (optional)
	"""

	output = bytearray()
	pack_uint32 = lambda num: output.extend(struct.pack('<I', int(max(0, min(num, 2**32-1))) ))

	# Extract the rule priority from the translation dictionary
	priority = rule.priority
	if translation_dict is not None:
		if rule.priority not in translation_dict.keys():
			raise ValueError('Priority %d was found in tree but not in dictionary' % rule.priority)
		priority = translation_dict[rule.priority]

	pack_uint32(priority)
	for i in range(5):
		pack_uint32(rule.ranges[i*2+0]) # Range Low
		pack_uint32(rule.ranges[i*2+1]) # Range High

	return output


# =============================
# Script start
# =============================

args=parse_arguments()

# Load the NeuroCuts tree module
try:
	sys.path.append(args.neurocuts)
	import tree
except ImportError as e:
	print('Cannot find tree module in NeuroCuts directory: %s' % e)
	exit(1)

# Read NeuroCuts output format
print('Loading NeuroCuts generated tree...')
with open(args.input, 'rb') as f:
	tree = pickle.load(f)

translation_dict=None

# Check whether the rules' priorities should be renumbered
if args.remainder_indices is not None:
	# Read the rules (original & remainder) from disk
	print('Reading remainder indices...')
	translation_dict={}
	with open(args.remainder_indices, 'rb') as f:
		reader = ObjectReader(f.read())
	total=reader.read_uint32()
	for i in range(total):
		translation_dict[i] = reader.read_uint32()

# Write output binary format
print('Packing tree...')
with open(args.output, 'wb') as f:
	f.write(pack_neurocuts_tree(tree, translation_dict))

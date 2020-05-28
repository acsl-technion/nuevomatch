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


import sys
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt

# Add the "bin" directory to import path
my_path=os.path.dirname(os.path.realpath(__file__))
bin_path=os.path.join(my_path, '..', 'bin')
if not os.path.exists(bin_path):
	print('Bin directory was not found. Did you compile the library?')
	exit(1)
sys.path.append(bin_path)

# Import modules
try:
	from compatible_interval_set import CompatibleIntervalSet
	from rule_handlers import *
	from object_reader import ObjectReader
	from object_packer import ObjectPacker
except ImportError as e:
	# Cannot find the library
	print('One of the library files was not found (%s). Did you compile the library?' % e)
	exit(1)

def parse_arguments():
	""" Parse script arguments
	"""

	parser = argparse.ArgumentParser()

	# Mandatory arguments
	parser.add_argument('-f', '--filter',   required=True, help='Filename of valid Classbench text file to process')

	# Global arguments
	parser.add_argument('-v', '--verbosity', type=int, default=0, help='Verbosity. Higher number == higher verbosity.')

	# Ruleset index control
	parser.add_argument('--indices-load', default=None, type=str, help='Load a ruleset indices file to filter rules from the ruleset')
	parser.add_argument('--indices-save', default=None, type=str, help='Generate a ruleset indices file of the remainder subset')

	# Control Modes
	parser.add_argument('--coverage', action='store_true', help='Print the coverage of the iSets as a function of size and iSet index')
	parser.add_argument('--find', type=int, default=None, help='Find the iSet index that holds the given rule')
	parser.add_argument('--print', type=int, default=None, help='Print an iSet to the screen')
	parser.add_argument('--save', default=None, help='Save packed format of the ruleset to file')
	parser.add_argument('--print-table', action='store_true', default=False, help='Print the rule-set to the screen (as numeric ranges)')
	parser.add_argument('--get-field-attributes', type=int, nargs='+', default=None, help='Analyze a subset of fields and print their attributes')
	parser.add_argument('--get-ruleset-type', action='store_true', default=False, help='Checks whether the ruleset have 5-tuple or OpenFlow format')
	parser.add_argument('--get-ruleset-size', action='store_true', default=False, help='Print the number of rules in the ruleset')
	parser.add_argument('--randomize', action='store_true', default=False, help='Randomize new rule-table based on the statistics of the loaded one. ' +
																'The arguments --size and --num_of_fields are required.')
	parser.add_argument('--save-remainder-as-classbench', default=None, type=str, help='Works with Classbench model only. ' +
																'Saves the remainder classifier as a new Classbench file.')

	# Coverage Mode Knobs
	parser.add_argument('--sample-half', action='store_true', help='(Coverage Mode) Perform random sampling of the rules, \
																	half rules in each iteration. Stop when table size is less than 1K.')
	parser.add_argument('--repeat', type=int, default=1, help='(Coverage Mode) Repeat the process X times')
	parser.add_argument('--print-fields-without-isets', action='store_true', default=False, help='(Coverage Mode) Print indices of fields without iSets')
	parser.add_argument('--size', type=int, default=0, help='(Randomize Mode) Set the number of rules to randomize')
	parser.add_argument('--num-of-fields', type=int, default=0, help='(Randomize Mode) Set the number of 128bit fields to randomize')

	# Control hyper-parameters
	parser.add_argument('--max-subsets',	type=int,   default=6,  help='Hyper-Parameters: Number of maximum allowed subsets')
	parser.add_argument('--min-size',		type=int,   default=64, help='Hyper-Parameters: Number of minimum intervals per subset')
	parser.add_argument('--precision',		type=str,   default='uint32', help='Hyper-Parameters: Experimental. Field precision. Longer fields are split to precision. Values in [uint32, float32]')

	# In case of no arguments, print usage
	if len(sys.argv)<=1:
		parser.print_usage()
		parser.exit()

	return parser.parse_args()

# ================================================ #
# ================ Script Start ================== #
# ================================================ #

args=parse_arguments()

# Print the ruleset type
if args.get_ruleset_type:
	print(RulesetHandler.get_format(args.filter))
	exit(0)

print('Reading Ruleset file...')
rule_handler=RulesetHandler.read(args.filter, verbose=args.verbosity)

# In case of filtering out rule indices
if args.indices_load:
	# Read file
	with open(args.indices_load, 'rb') as f:
		reader = ObjectReader(f.read())
	# Generate indices
	indices=[]
	num_of_indices=reader.read_uint32()
	for _ in range(num_of_indices):
		val = reader.read_uint32()
		# Validate correctness
		if val > len(rule_handler):
			raise ValueError('Indices file has invalid value (%d)' % val)
		indices.append(val)
	# Apply
	rule_handler.apply(indices)

# Print the ruleset size after filtering out
if args.get_ruleset_size:
	print(len(rule_handler))

# In case of randomizing new table
if args.randomize > 0:
	if args.size == 0 or args.num_of_fields == 0:
		raise ValueError('Randomize: one of these arguments is missing: --size, --num-of-fields')
	print('Randomizing new rules...')
	rule_handler.randomize_table(args.size, args.num_of_fields)

# Pack the ruleset into a binary file for fast handling
if args.save is not None:
	print('Saving ruleset in binary format (total %d rules)...' % len(rule_handler))
	with open(args.save, 'wb') as f:
		f.write(bytes(rule_handler))

# Which kind of ruletable to use?
if args.precision == 'uint32': rule_table = rule_handler.get()
elif args.precision == 'float32': rule_table = rule_handler.generate_imprecise_32bit_table()
else: raise ValueError('Precision type is not supported. See help for valid precision values.')

N=rule_table.shape[0]

# In case of printing the rule-set to stdout
if args.print_table:
	for r in range(rule_table.shape[0]):
		for c in range(rule_table.shape[1]):
			print('%.0f\t' % rule_table[r,c], end='')
		print('\n', end='')

# In case of field analyze mode
if args.get_field_attributes is not None:

	F=int(rule_table.shape[1]/2)

	# Check input validity
	indices = args.get_field_attributes
	invalid=[]
	for i, indx in enumerate(indices):
		if indx < 0 or indx >= F:
			print('Warning: ignoring invalid field %d (valid values: 0-%d)' % (indx, F-1))
			invalid.append(i)

	# Remove invalid fields
	invalid.reverse()
	for x in invalid:
		indices.pop(x)

	# Count how many unique start and end values per index
	for indx in indices:
		value = np.unique(rule_table[:, [indx, indx+F]]).shape[0]
		print('Field %d has %d unique values' % (indx, value))

	# Try to generate an iSet from the field alone
	print('Trying to generate at most %d iSets from the given fields alone, each iSet with minimum size of %d rules...' % (args.max_subsets, args.min_size))

	# Indices now hold the field-start values in the rule-table.
	# Add the corresponding field-end values, and the rule priority
	indices.extend([x+F for x in indices])
	indices.append(-1)

	# Generate iSets from the reduced table
	rule_table=rule_table[:, indices]
	compatible_set=CompatibleIntervalSet(rule_table)
	compatible_set.process(args.max_subsets, args.min_size, verbose=args.verbosity)

	# Print results
	if len(compatible_set) == 0:
		print('Could not generate any iSet from the given field')
	else:
		print('Generated %d iSets:' % len(compatible_set))
	for i in range(len(compatible_set)):
		# Get the field of the current iSet
		current_field = compatible_set[i].get_field_index()
		# Translate to original field value
		original_field = sum([indices[x] for x in range(len(indices)) if x==current_field])
		# Print message
		print('iSet %d (field %d) with coverage %f' % (i, original_field, compatible_set[i].coverage()))
	exit(0)



# In case of coverage mode
if args.coverage:

	# Generate classes for random sampling
	sample_table=[0]
	if args.sample_half:
		if N < 1e3: sample_table=[0]
		elif N < 1e4: sample_table=[0, 1e3]
		elif N < 1e5: sample_table=[0, 1e3, 1e4]
		elif N < 2e5: sample_table=[0, 1e3, 1e4, 1e5]
		else: sample_table=[0, 1e3, 1e4, 1e5, 2e5]

	# Store original table. What is the table class? (1K, 10K, etc.)
	org_table = rule_table
	F=int(org_table.shape[1]/2)

	# Repeat the process X times
	for _ in range(args.repeat):

		# Reset table as original table
		rule_table=np.copy(org_table)
		current_sample_table = sample_table[:]

		# Calculate iSets
		print('Generating iSets from a %dx%d rule-table' % (rule_table.shape[0], F))
		compatible_set = CompatibleIntervalSet(rule_table)
		compatible_set.process(args.max_subsets, args.min_size, verbose=args.verbosity)
		print('Done.')

		# Hold fields without iSets
		non_iset_fields=[i for i in range(F)]

		while len(current_sample_table) != 0:

			K = rule_table.shape[0]

			# What is the class of the current table
			if K <= 1e3: size='1k'
			elif K <= 1e4: size='10k'
			elif K <= 1e5: size='100k'
			elif K <= 2e5: size='200k'
			else: size='500k'

			# Print iSet coverage
			for i in range(len(compatible_set)):
				value = compatible_set[i].coverage()
				field_idx = compatible_set[i].get_field_index()
				validation_phases = compatible_set[i].get_validation_phase_length()
				print('Size: %s (%d) iSet: %d (field: %d) Coverage: %f Validation-Phases: %d' %
					(size, K, i, field_idx ,value, validation_phases))
				# Update the list of fields with no iSet
				if field_idx in non_iset_fields:
					non_iset_fields.remove(field_idx)

			# Remove first element
			K = int(current_sample_table.pop())

			# Sample rule table
			print('Sampling rule-table to hold %d rules' % K)
			indices = np.random.randint(0, rule_table.shape[0], K)
			rule_table = rule_table[indices,:]

			# Re-calculate iSets
			print('Generating iSets from a %dx%d rule-table' % (rule_table.shape[0], F))
			compatible_set = CompatibleIntervalSet(rule_table)
			compatible_set.process(args.max_subsets, args.min_size, verbose=args.verbosity)
			print('Done.')

		# Print all fields with no iSet
		if args.print_fields_without_isets:
			for f in non_iset_fields:
				print('Field %d has no iSet' % f)

	exit(0)

# In case of find mode
if args.find is not None:
	print('Rule %d is being searched across all subsets.' % args.find)
	print('Output format: list of tuples (X,Y,S,R) where:')
	print('X - indicates the subset index (-1 for the remainder set)')
	print('Y - the position of rule within that subset')
	print('S - the field index by which the subset was is indexed')
	print('R - the range of the filed by which the subset was is indexed')
	compatible_set = CompatibleIntervalSet(rule_table)
	compatible_set.process(args.max_subsets, args.min_size, verbose=args.verbosity)
	set = compatible_set.rule_query(args.find)
	for item in set:
		print(str(item))

# In case of print mode
if args.print is not None:
	# Calculate iSets
	compatible_set = CompatibleIntervalSet(rule_table)
	compatible_set.process(args.max_subsets, args.min_size, verbose=args.verbosity)
	print(compatible_set[args.print])

# In case the ruleset indices should be saved
if args.indices_save:
	# Calculate iSets
	compatible_set = CompatibleIntervalSet(rule_table)
	compatible_set.process(args.max_subsets, args.min_size, verbose=args.verbosity)
	# Get the remainder indices
	remainder_indices = compatible_set.remainder_indices()
	print('Saving %d remainder indices to file %s...' % (len(remainder_indices), args.indices_save))
	# Generate data
	packer = ObjectPacker()
	packer.append(len(remainder_indices))
	for x in remainder_indices:
		packer.append(int(x))
	# Save file
	with open(args.indices_save, 'wb') as f:
		f.write(bytes(packer))

# In case of saving remainder classifier as Classbench file
if args.save_remainder_as_classbench is not None:
	if RulesetHandler.get_format(args.filter) != 'Classbench':
		raise ValueError('Cannot save remainder set as Classbench file: input ruleset must be of Classbench format')
	# Calculate iSets
	compatible_set = CompatibleIntervalSet(rule_table)
	compatible_set.process(args.max_subsets, args.min_size, verbose=args.verbosity)
	# Get the remainder indices
	remainder_indices = compatible_set.remainder_indices()
	# Read all Classbench lines
	with open(args.filter, 'r') as f:
		lines = f.readlines()
	# Select only the relevant rules
	lines = [lines[i] for i in remainder_indices]
	# Save to output file
	print('Saving remainder set as Classbench file to "%s"...' % args.save_remainder_as_classbench)
	with open(args.save_remainder_as_classbench, 'w') as f:
		f.write(''.join(lines))

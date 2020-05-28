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
import re
import numpy as np
import argparse
import os
from timeit import default_timer as timer

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
	from lookup_cpu import LookupCpu, initialize_library
	from rule_handlers import *
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
	parser.add_argument('-o', '--output',   required=True, help='Filename of output NuevoMatch classifier')

	# Control hyper-parameters
	parser.add_argument('--max-subsets',	type=int,   default=6,  help='Hyper-Parameters: Number of maximum allowed subsets')
	parser.add_argument('--min-size',	   type=int,   default=64, help='Hyper-Parameters: Number of minimum intervals per subset')
	parser.add_argument('--max-error',	  type=int,   default=64, help='Hyper-Parameters: Number of maximum allowed subsets')

	# In case of no arguments, print usage
	if len(sys.argv)<=1:
		parser.print_usage()
		parser.exit()

	return parser.parse_args()


# ================================================ #
# ================ Script Start ================== #
# ================================================ #

args=parse_arguments()

# Initiate RQRMI library without log
initialize_library(False)

print('Reading Ruleset file...')
rule_handler=RulesetHandler.read(args.filter)

# Read the rule table and process with hyper-parameters
print('Creating Compatible Interval Set...')
compatible_set = CompatibleIntervalSet(rule_handler.get())
compatible_set.process(args.max_subsets, args.min_size, verbose=1)

# Packs the output
output = ObjectPacker()

# Statistics
training_time = []
total_size = 0

# Process iSets
for i, iset in enumerate(compatible_set):
	print('Processing iSet %d' % i)
	lookup = LookupCpu(len(iset), verbose=1)
	lookup.extend(iset.get_index())

	# In case the iSet has less rules than the maximum errors allows, skip iSet
	if len(iset) < args.max_error:
		print('Warning: Skipping iSet as it has %d rules < max-error = %d' % (len(iset), args.max_error))
		continue

	# Set hyper parameters for NuevoMatch
	lookup.epochs=[40, 40, 40]
	lookup.samples_per_bucket = 2000
	lookup.retraining_multiplier = 3
	lookup.threshold_submodel_error = args.max_error
	lookup.retraining_times = 3

	# Measure training time
	time_start = timer()
	lookup.compile()

	# Update statistics
	current_train_time =  (timer() - time_start)
	training_time.append(current_train_time)
	total_size += lookup.get_size()
	print('Total RQRMI training time: %d sec' % current_train_time)

	print('Packing iSet')
	with output as packed_iset:
		packed_iset.append(lookup.pack())
		packed_iset.append(iset)

if total_size==0:
	print('Error: Cannot create NuevoMatch for classifier: vo valid iSets!')
	exit(1)


total_training_time = sum(training_time)
print('Total training time: %d secs' % total_training_time)

# Write static information
output.insert(int(total_training_time * 1000)) # Build time in ms
output.insert(int(total_size)) # Size in bytes
output.insert(len(rule_handler))
output.insert(len(compatible_set))

# Write remainder set
remainder = compatible_set.remainder()
output.append(remainder)

# Write output to file
print('Writing output file')
with open(args.output, 'wb') as f:
	f.write(bytes(output))


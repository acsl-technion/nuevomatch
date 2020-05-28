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


import struct
import numpy as np
import sys
from rule_handlers import Ruleset
from object_packer import ObjectPacker

_log_last_length=0
def _log(verbose, msg, delete_last=False):
	""" Prints log messages according to verbosity

	Args:
		verbose: integer. message is printed when not zero
		msg: message to print, newline char should be included
		delete_last: If true, the last message will be deleted from buffer
	"""
	global _log_last_length
	if verbose==0:
		return
	if delete_last:
		sys.stderr.write('\b'*_log_last_length)
	sys.stderr.write(msg)
	sys.stderr.flush()
	_log_last_length=len(msg)


class iSet:
	""" A set of rules which do not intersect each other on field f
	"""

	def __init__(self, rule_table, iset_indices, iset_field, total_rules, validation_indices, verbose=0):
		""" Initiate this

		Args:
			rule_table: a reference to the iSet rule-table (uint32)
			iset_indices: a list of the relevant iSet indices within the rule-table
			iset_field: the field of the iSet
			total_rules: total rules in classifier after iset partitioning
			verbose: verbosity of this
		"""

		F = int(rule_table.shape[1]/2)
		self.F = F
		self.coverage_value = float(iset_indices.shape[0]) / float(total_rules)
		self.index = rule_table[iset_indices].astype(np.float32)
		self.index = self.index[:, [iset_field, iset_field+F]]
		self.indexed_field = int(iset_field)
		self.database =  rule_table[iset_indices, :].astype(np.uint32)
		self.verbose = verbose

		# Fix the case when the last interval has end=start
		if self.index[-1,0] == self.index[-1,1]:
			self.index[-1,1] = np.nextafter(self.index[-1,0], 2*self.index[-1,0], dtype=np.float32)

		# Validate non overlapping ranges in index
		test_vector = np.concatenate([self.index[:, 0], [self.index[-1,1]] ])
		test_vector = (np.roll(test_vector, -1) - test_vector) <= 0
		test_vector = test_vector[:-1]
		if np.logical_or.reduce(test_vector):
			raise ValueError('iSet has overlapping ranges!')

		# At this point we wish to create the validation-matrix per indexed rule
		# We need to partition all indexed rules and validation rules to groups
		# A group holds all rules with the same priorty and projection of the iSet field

		# Relevant field indices
		f_start = self.indexed_field
		f_end = f_start + F

		# Copy rules that participate in the validation phase, as exact representation (uint32)
		# Complexity: O(r*log(r)) for r the number of rules
		db_indices = np.unique(np.concatenate([validation_indices, iset_indices])).astype(np.uint32)
		database = rule_table[db_indices, :].astype(np.uint32)

		# At this point 'database' holds both non-overlapping indexed rules
		# and validation rules with sampe projections
		# We sort the rules by their start value - it is promised that
		# rules with the same start value overlap and have the same priority!
		# Complexity: O(r*log(r)) for r the number of rules
		indices = np.lexsort([database[:, f_start]])
		database = database[indices, :]

		# Divide the database to gropus of rules with the same priority and projection value
		# Each group holds: (priority, start-idx (inclusive), end-idx (exclusive))
		# Complexity: O(r) for r the number of rules
		groups=[]
		last_idx = 0
		last_prio = database[0, -1]
		last_start = database[0, f_start]
		N = database.shape[0]
		for rule_idx in range(1, N):
			current_start = database[rule_idx, f_start]
			# In the the start value of the current rule differs from the
			# start value of the last rule, they must belong to different
			# groups
			if current_start != last_start:
				groups.append((last_prio, last_idx, rule_idx))
				last_idx = rule_idx
				last_start = current_start
				last_prio = database[rule_idx, -1]
		groups.append((database[-1, -1], last_idx, N))

		# Validate that the number of groups equals the number of rules
		if(self.index.shape[0] != len(groups)):
			raise Exception('Number of groups (%d) differ from number of indexed rules (%d)' %
				(len(groups), self.index.shape[0]))

		# Validation phase of rules in this
		# List of validation matrics, each corresponds to a different rule-group (i.e, rules with same priority)
		# Each column is a set of disjunction conditions on the same field
		# Each row is a different validation phase
		# Note: while the binary (& memory) format of NuevoMatch holds matrices with the same size,
		# At this point different matrices may have different sizes
		self.validation_matrices = []
		self.validation_priorities = []

		# For each group (all rules with the same priority)
		for prio, start, stop in groups:
			# Check that all rules in group have the same start-value, end-value, and priority
			start_values = np.unique(database[start:stop, f_start]).shape[0]
			end_values = np.unique(database[start:stop, f_end]).shape[0]
			prio_values = np.unique(database[start:stop, -1]).shape[0]
			if (start_values != 1) or (end_values != 1) or (prio_values != 1):
				_log(self.verbose, 'Number of start-values: %d, end-values: %d, prio-values: %d\n' %
					(start_values, end_values, prio_values))
				raise Exception('Error in group partitioning.')

			# Get the set of valid ranges in each field.
			# Total complexity: O(r*F) for r number of rules, F number of fields
			valid_ranges=[]
			for f in range(F):
				# Initialize as empty set
				values = set()
				# Get all possible values
				# Complexity: O(r) for r number of rules
				possible_values = database[start:stop, [f, f+F]]
				# Filter only unique values
				for r in range(possible_values.shape[0]):
					current = (int(possible_values[r,0]), int(possible_values[r,1]))
					values.add(current)
				valid_ranges.append(values)

			# How many validation phases are required?
			phase_num = np.max([len(x) for x in valid_ranges])
			# Allocate memory for validation phase of current rule
			# Note: The allocation invalidates all fields in validation phase
			# by putting [0xffffffff, 0] as range in all fields (invalid range)
			validation_matrix = np.zeros([phase_num, 2*F], dtype=np.uint32)
			validation_matrix[:, 0:F] = np.iinfo(np.uint32).max
			# Build validation phases
			for i in range(phase_num):
				for f in range(F):
					if len(valid_ranges[f]) != 0:
						current_range = valid_ranges[f].pop()
						validation_matrix[i, f] = current_range[0]
						validation_matrix[i, f+F] = current_range[1]
			self.validation_matrices.append(validation_matrix)
			self.validation_priorities.append(prio)


		assert(len(self.validation_matrices) == self.index.shape[0])

		# Print messages
		_log(self.verbose, 'iSet has %d rules, %d validation phases, and %d columns\n' %
			(self.index.shape[0], self.get_validation_phase_length(), 2*F))


	def __len__(self):
		""" Returns the number of rules in this """
		return self.index.shape[0]


	def __str__(self):
		output = ''
		for i in range(self.index.shape[0]):
			output += '%d:' % i
			for f in range(self.database.shape[1]):
				if f==self.indexed_field: output += '*'
				output += '%d' % self.database[i,f]
				if f==self.indexed_field: output += '*'
				output+=' '
			else:
				output += '\n'
		return output


	def get_index(self):
		""" Return the index of the iSet """
		return self.index

	def get_field_index(self):
		""" Returns the field index this iSet was created by """
		return self.indexed_field

	def get_validation_phase_length(self):
		""" Returns the number of validation phases of this """
		return int(np.max([x.shape[0] for x in self.validation_matrices]))

	def __bytes__(self):
		""" Packs this to byte array object """
		output = ObjectPacker()
		K = self.get_validation_phase_length()
		F = self.F
		with output as iset_packer:
			# Pack iSet version 1
			iset_packer.append(0)
			iset_packer.append(1)

			# Pack number of validation phases, fields, and field index
			iset_packer.append(K)
			iset_packer.append(F*2)
			iset_packer.append(self.indexed_field)

			# Pack Validation phases
			# Pack format: elements are stored row-wise: E00 E01 E02 ... E10 E11 E12 ...
			# Note: Due to SIMD implementation in interval_set.cpp, in runtime,
			# the data is stored in memory in a transposed format (column-wise)
			# More deatils in interval_set.cpp
			for i in range(len(self.validation_matrices)):
				# Get the current matrix
				matrix = self.validation_matrices[i]
				for k in range(K):
					for f in range(2*F):
						# The default value invalidates match by
						# assigning MAX_UINT for range start and 0 for range end
						default_value = int(np.iinfo(np.uint32).max) if (f<F) else 0
						# In case the matrix does not have a value at the current index,
						# Store the default value
						if matrix.shape[0] <= k:
							iset_packer.append(int(default_value))
						else:
							iset_packer.append(int(matrix[k,f]))
				# Pack the rule priority
				iset_packer.append(self.validation_priorities[i])

		# Print messages
		_log(self.verbose, 'iSet total pack size is %d bytes\n' % len(output))
		return bytes(output)


	def coverage(self):
		""" Returns the coverage of this """
		return self.coverage_value


class RemainderSet:
	""" Represents a remainder set """

	def __init__(self, rule_table, indices):
		""" Initiate this

		Args:
			rule_table: a reference to the complete rule table
			indices: a list of the relevant iSet indices within the rule-table
		"""

		self._db = Ruleset(rule_table[indices, :])


	def get_db(self):
		""" Returns the rules of this """
		return self._db.get()


	def __bytes__(self):
		""" Packs this to byte-array """
		return bytes(self._db)


class CompatibleIntervalSet:
	""" Holds Compatible Interval Subsets of a rule table
	"""

	def __init__(self, rule_table):
		""" Initialize new multiset.

		Args:
			rule_table: Matrix of Nx(2*F+1), where N is the number of rules and F is the number of fields.
						Rows' format: [field_0_start, field_1_start, ..., filed_0_end, field_1_end, ..., priority]

		Throws:
			ValueError in case rule_table is not in valid format
		"""

		if type(rule_table) is not np.ndarray:
			raise ValueError('rule_table argument must by Numpy array')

		self.rule_table = rule_table
		self.N = self.rule_table.shape[0]
		self.F = int(rule_table.shape[1]/2)
		self.total_rules_for_coverage=self.N
		self.subsets = []
		self.subset_field = []
		self.remainder_indx = None
		self.isets = []

		# Store extra validation phases for rules based on their priority
		# Item i is a list of all expanded rule indices for iSet i.
		self.extra_validation_phases=[]

		# Used for iterator
		self.iterator = 0


	def _find_compatible_subset_in_field(self, field, available_rules):
		""" Private method. Extract possible compatible subset from the current rule-set

		Args:
			field: The number of field (0 <= f < F) to extract from
			available_rules: An array of indices, the available rules within the rule_table

		Returns:
			A list of indices (from the rule-set) of the compatible subset
		"""

		# Extract the field's intervals from rule table
		intervals = self.rule_table[available_rules, :].astype(np.float32)
		intervals = intervals[:, [field, field+self.F]]
		lengths = intervals[:, 1] - intervals[:, 0]
	
		# Remove negarive lengths (this is possible to split fields with modulo on thier 32bit values
		negative = (lengths < 0)
		intervals = intervals[~negative]
		lengths = lengths[~negative]
		N = intervals.shape[0]

		# In case there are no valid intervals, return an empty list
		if N == 0:
			return np.empty(shape=[0,1])
	
		# Apply the greedy algorithm for finding the largest compatible subset

		# Sort based on interval length (large to small) and finish value (small to large)
		# This kind of sorting perfers short intervals on long ones
		sorted_idx = np.lexsort([ -lengths, intervals[:, 1] ])
		sorted_intervals = intervals[sorted_idx]

		subset=[]
		i = 0
		while i<N:
			# Select the first interval
			current_interval = i
			subset.append(current_interval)
			i+=1
			# Ignore intervals that are non compatible with current_interval
			while (i < N) and (sorted_intervals[i, 0] <= sorted_intervals[current_interval, 1]):
				i+=1

		# At this point, subset consists of an optimal set of intervals
		# Return a vector with the subset's indices
		return sorted_idx[subset]


	def process(self, max_subset_count, min_items_per_subset, verbose=0):
		""" Extract optimal compatible sets from the rule-table

		Args:
			max_subset_count: The maximum number of allowed subsets
			min_items_per_subset: The minimum allowed number of items in a subset
			verbose: Verbosity
		"""

		# Cannot process empty ruleset
		if self.N == 0:
			return

		available_rules = np.arange(self.N)
		N = available_rules.shape[0]

		for i in range(max_subset_count):
			# Stop in case thre are no available rules left
			if available_rules.shape[0]==0:
				break
			# Extract the subset with maximum intervals
			max_subset = np.empty(shape=[0])
			max_field = 0
			for f in range(self.F):
				subset = self._find_compatible_subset_in_field(f, available_rules)
				if subset.shape[0] > max_subset.shape[0]:
					max_subset = subset
					max_field = f
			# In case the current subset is smaller than minimum, stop
			if max_subset.shape[0] < min_items_per_subset:
				break

			# Update the subsets of this
			_log(verbose, 'Generated subset %d with %d rules (field index: %d) \n' % (i, max_subset.shape[0], max_field))
			self.subsets.append(available_rules[max_subset])
			self.subset_field.append(max_field)
			self.extra_validation_phases.append([])

			# Remove the iSet rules from the remaining rule-set
			predicat = np.full(N, True)
			predicat[max_subset] = False

			# Update available rules
			available_rules = available_rules[predicat]
			N = available_rules.shape[0]

			# Get the priorities of the rules in current iSet
			# Note: as priorities are unique per compact rules, two rules with the same
			# priority indicates that they both originate from the same compact rule
			# Complexity: O(r*log(r)) for r the number of rules
			iset_priorities = np.unique(self.rule_table[self.subsets[-1], -1])

			# Get the priorities of the rules in the remainder set.
			# Complexity: O(r*log(r)) for r the number of rules
			remainder_priorities = np.unique(self.rule_table[available_rules, -1])

			# Compute the intersection of priorities
			# Complexity: O(r) for r the number of rules
			intersect_priorities = set(np.intersect1d(iset_priorities, remainder_priorities))

			# Get the range in the field by which the last iSet was partitioned,
			# for all rules with priorities in the intersection.
			# Complexity: O(r) for r the number of rules
			value_prio_tuple_set=set()
			for j in self.subsets[-1]:
				# Set membership: O(1)
				if self.rule_table[j, -1] in intersect_priorities:
					range_start, range_end, value = self.rule_table[j, [max_field, max_field+self.F, -1]]
					value_prio_tuple_set.add((range_start, range_end, value))

			# Create new predicat to remove rules
			predicat = np.full(N, True)
			counter = 0

			# Go over all rules in the remainder set, in case their range&prio exists in
			# the value_prio_tuple_set, remove the rule
			# Complexity: O(r) for r the number of rules
			for j in range(N):
				idx=available_rules[j]
				range_start, range_end, value = self.rule_table[idx, [max_field, max_field+self.F, -1]]
				if (range_start, range_end, value) in value_prio_tuple_set:
					predicat[j] = False
					counter += 1
					# Update extra validation phases
					self.extra_validation_phases[i].append(available_rules[j])

			# Update available rules
			available_rules = available_rules[predicat]
			N = available_rules.shape[0]

 			# Log
			_log(verbose, 'Removed %d expanded-rules from remainder set \n' % counter)

			# Do not continue to next iteration if remainder is less than minimum
			if available_rules.shape[0] < min_items_per_subset:
				break

		self.remainder_indx = available_rules
		_log(verbose, 'Remainder subset with %d rules \n' % available_rules.shape[0])

		indices = np.arange(self.N)
		_log(verbose, 'Extra validation phase covers %d rules \n' % sum([len(x) for x in self.extra_validation_phases]))

		# Update the total rules of this (duplicates might have changed this)
		self.total_rules_for_coverage=sum([x.shape[0] for x in self.subsets]) + self.remainder_indx.shape[0]
		_log(verbose, 'Total size after removing expanded rules: %d\n' % self.total_rules_for_coverage)

		# Build all iSet objects of this
		self.isets = [iSet(self.rule_table, self.subsets[key], self.subset_field[key],
				self.total_rules_for_coverage, self.extra_validation_phases[key], verbose) for key in range(len(self.subsets))]


	def __eq__(self, rhs):
		""" Returns true whether two CompatibleIntervalSet equal """

		if type(self) != type(rhs): return False

		if  (len(self.subsets) != len(rhs.subsets)) or (len(self.subset_field) != len(rhs.subset_field)):
			return False

		for lhs_subset, rhs_subset in zip(self.subset_field, rhs.subset_field):
			if lhs_subset != rhs_subset:
				return False

		for lhs_subset, rhs_subset in zip(self.subsets, rhs.subsets):
			if not np.all(lhs_subset == rhs_subset):
				return False

		if not np.all(self.remainder_indx == rhs.remainder_indx):
			return False

		if not np.all(self.rule_table == rhs.rule_table):
			return False

		return True


	def __len__(self):
		""" Returns the number of iSets """
		return len(self.subsets)


	def __getitem__(self, key):
		""" Get a specific iSet """
		return self.isets[key]


	def __iter__(self):
		""" Get iterator for this """
		return self


	def __next__(self):
		""" Iterator next """
		if self.iterator == len(self):
			raise StopIteration
		self.iterator += 1
		return self[self.iterator-1]


	def remainder(self):
		""" Returns the remainder set of this """
		return RemainderSet(self.rule_table, self.remainder_indx)


	def remainder_indices(self):
		""" Returns the remainder subset indices of the original rules """
		return self.remainder_indx


	def rule_query(self, rule_priority):
		""" Returns the iSet index that contains the rule

		Args:
			rule_idx: The rule index

		Returns: list of tuples (X,Y,S,R) where:
			X - indicates the subset index (-1 for the remainder set)
			Y - the position of rule within that subset
			S - the field index by which the subset was is indexed
			R - the range of the filed by which the subset was is indexed
		"""
		output = []
		# Check all iSets
		for x,s in enumerate(self.subsets):
			prio_array = self.rule_table[s, -1]
			for y,v in enumerate(prio_array):
				if v == rule_priority:
					field = self.subset_field[x]
					field_start, field_end = self.rule_table[s, :][y, [field, field+self.F]]
					output.append((x, y, field, (int(field_start), int(field_end))))
		# Check the remainder subset
		prio_array = self.rule_table[self.remainder_indx, -1]
		for y,v in enumerate(prio_array):
			if v == rule_priority:
				output.append((-1, y, -1, -1))
		return output


	def coverage(self):
		""" Returns the coverage of this """
		return 1 - len(self.remainder_indx) / self.total_rules_for_coverage


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


import json
import re
import math
import struct
import os
import numpy as np
import sys
import traceback
from functools import reduce
from os.path import stat

_log_last_length=0
def _log(verbose, msg, delete_last=False):
	""" Prints log messages according to verbose

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


class Ruleset:
	""" Holds any number of rules, supports common operations """

	def __init__(self, rules):
		self.rule_table = rules


	def __bytes__(self):
		""" Packs this to binary representation """

		# See https://docs.python.org/3/library/struct.html#struct.pack
		# for struck pack format

		output = bytearray()
		# Local methods to pack numbers in little-endian format
		write_uint32 = lambda num: output.extend(struct.pack('<I', num))

		# Write version (version 1)
		write_uint32(0)
		write_uint32(1)

		# Write number of rules and fields
		write_uint32(len(self))
		write_uint32(self.width())

		# Sort rules by their priority
		if self.rule_table.shape[0] > 0:
			self.rule_table = self.rule_table[np.lexsort([self.rule_table[:,-1]]), :]

		for i in range(self.rule_table.shape[0]):
			# Write the rule priority
			write_uint32(self.rule_table[i,-1])
			# Write rule fields
			for f in range(self.width()):
				write_uint32(self.rule_table[i,f])
				write_uint32(self.rule_table[i,f+self.width()])

		return bytes(output)


	def get(self):
		""" Returns the rule-table of this """
		return self.rule_table.astype(np.uint32)


	def __len__(self):
		""" Returns the number of rules in this """
		return self.rule_table.shape[0]

	def __getitem__(self, indices):
		""" Gets a subset of rules of this """
		return self.rule_table[indices]


	def apply(self, indices):
		""" Applies a rulset subset on this """
		self.rule_table = self.rule_table[indices, :]


	def width(self):
		""" Returns the number of 32bit fields in this """
		return math.floor(self.rule_table.shape[1] / 2)


	def __str__(self):
		""" Returns a string representation of this """
		out=''
		for r in range(self.rule_table.shape[0]):
			for c in range(self.rule_table.shape[1]):
				out += '%f\t' % self.rule_table[r,c]
			out += '\n'
		return out



class RulesetHandler(Ruleset):
	""" Hold information of rulesets to be used with iSet partitioning and the validation phase """

	def __init__(self, filename, verbose=0):
		""" Initiate this, bind to ruleset file.

		Args:
			filename: A valid filename to read the ruleset from

		Throws:
			ValueError in case the filename is invalid
		"""

		if not os.path.isfile(filename):
			raise ValueError('Ruleset file %s not found' % filename)

		# Save the filename of this
		self.filename = filename

		# The rule-table is composed of 32bit unsigned integers
		# Any field larger than 32bit is split to 32bit parts
		super().__init__(None)

		# Holds a mapping from original field name to rule-table's columns
		self.field_dict={}


	def split_to_32bit_parts(self, value, num_of_parts):
		""" Splits a number of any bit number of 32bit parts """
		return [value >> 32*x & 0xFFFFFFFF for x in range(num_of_parts)]


	def split_start_end_32bit_parts(self, start ,end, parts):
		""" Splits two long integers (start - end) to list of 32bit parts

		Args:
			start: An integer, start value
			end: An integer, stop value
			parts: Number of parts in the field

		Returns:
			A list of [start, end] values per 32bit field-part
		"""

		# We do not support splitting ranges to more than 1 part
		if parts != 1:
			raise ValueError('Splitting ranges to 32bit parts is not supported. Got range [%d,%d] with %d parts.' % (start, end, parts))
		return [(start, end)]
		

	def calculate_start_end_32bit_parts(self, value, mask, parts):
		""" Split logical field (may have as arbitrary number of bits) to 32bit fields

		Args:
			value: An arbitrary-long integer
			mask: An arbitrary-long mask
			parts: Number of parts in the field

		Returns:
			A list of [start, end] values per 32bit field-part
		"""

		# Process mask to be long as the number of parts
		mask = bin(mask)[2:]
		mask = int('1' * (32*parts-len(mask)) + mask, 2)

		# Split value and maks to 32bit fields
		values=self.split_to_32bit_parts(value, parts)
		mask=self.split_to_32bit_parts(mask, parts)

		# Create ranges (start-end) per field
		output=[]
		for i in range(len(values)):
			start=values[i] & mask[i] & 0xffffffff
			end=start | ~mask[i] & 0xffffffff
			output.append([start, end])
		return output


	def calculate_start_end_long(self, value, mask, bits):
		""" Calculates the start and the end ranges of value-mask pair for (very) long integers

		Args:
			value: An arbitrary-long integer
			mask: An arbitrary-long mask
			bits: Number of bits in field

		Returns: A tuple (start, end) of python integers
		"""

		# Split to 32bit parts
		num_of_parts=int(np.ceil(bits/32))
		parts=self.calculate_start_end_32bit_parts(value, mask ,num_of_parts)
		parts.reverse()

		# Concatenate parts to long integers
		start = reduce(lambda acc, x: acc << 32 | x, [x[0] for x in parts])
		end= reduce(lambda acc, x: acc << 32 | x, [x[1] for x in parts])
		return (start ,end)


	def generate_imprecise_32bit_table(self):
		""" Generates a table with the original fields as float32 values
			Note that due to float32 precision error, only the 23 most
			significant bits in an integer are recorded
		"""

		F = len(self.field_dict.keys())
		N=len(self)
		imprecise_table=np.zeros([N, 2*F+1], dtype=np.float32)
		imprecise_table[:, -1] = np.arange(N)

		for rule_idx in range(self.rule_table.shape[0]):
			for col_idx, field_name in enumerate(self.field_dict.keys()):
				imprecise_table[rule_idx, [col_idx, col_idx+F]] = self.get_exact_field_value(field_name, rule_idx)

		return imprecise_table.astype(np.float32)


	def get_exact_field_value(self, field_name, rule_idx):
		""" Returns a long integer (may be even 128bit long) of the exact value of field in row

		Args:
			field_name: The name of the field
			rule_idx: The row number

		Returns:
			A tuple [low, high] of values relevant for the rule
		"""
		low_value, high_value = 0, 0
		for part_idx in self.field_dict[field_name]:
			low_value |= abs(int(self.rule_table[rule_idx, part_idx]))
			high_value |= abs(int(self.rule_table[rule_idx, part_idx+self.width()]))
			low_value <<= 32
			high_value <<= 32
		return (low_value >> 32, high_value >> 32)


	def collect_statistics(self, rule_table):
		""" Collects statistics on the ruleset fields. The following statistics are extracted:
			1) Min & Max values 2) Field start mean & stddev 3) Field length mean & stddev

		Args:
			rule_table: A table of Nx(2F) where N is the number of rules and F the number of fields

		Returns:
			A numpy table of Vx7, for V the number of fields with valid statistics, where the columns are:
			[field index, start-unique, end-unique, lengths-mean, lengths-stddev, start-mean, start-stddev]
		"""
		F = int(rule_table.shape[1]/2)
		output=np.zeros(shape=[F, 7], dtype=np.float64)
		valid_fields=0

		for f in range(F):
			# Extract information from table
			start_values = rule_table[:,f]
			end_values = rule_table[:,f+F]
			lengths = end_values - start_values

			# Remove empty fields
			empty_rows = (start_values==0) & (end_values==0)
			start_values=start_values[~empty_rows]
			end_values=end_values[~empty_rows]
			lengths=lengths[~empty_rows]

			# Skip fields with no statistics
			if len(lengths) == 0: continue
			valid_fields += 1

			# Get start&end uniqueness
			output[f, 0] = f
			output[f, 1] = np.unique(start_values).shape[0]
			output[f, 2] = np.unique(end_values).shape[0]

			# Calculate statistics of length
			output[f, 3:5] = np.mean(lengths), np.std(lengths)

			# Calculate statistics of start value
			output[f, 5:7] = np.mean(start_values), np.std(start_values)

		return output[0:valid_fields, :]


	def randomize_table(self, size, number_of_fields):
		""" Generate a new table based on the statistics of the current table. New fields are 128bit long.

		Args:
			size: Number of rules in the generated table
			number_of_fields: Number of fields to generate

		Note:
			Modifies the state of this, adds columns to rule-table
		"""

		# Collect the statistics of current table
		statistics = self.collect_statistics(self.rule_table)

		# Allocate memory for new table (each field is 128 bit, thus 4 32bit fields long)
		F = number_of_fields * 4
		N = size
		output=np.zeros(shape=[N,2*F+1], dtype=np.uint32)

		# Populate
		for f in range(number_of_fields):
			# Randomize field statistics
			k = np.random.randint(0, statistics.shape[0])
			start_mean, start_stddev = statistics[k, 5], statistics[k, 6]
			length_mean, length_stddev = statistics[k, 3], statistics[k, 4]

			# Randomize skew, generate rule start values
			skew = np.random.randint(0, 2**62)
			start_values = np.random.normal(start_mean+skew, start_stddev, N)

			# Randomize rule length
			lengths = abs(np.random.normal(length_mean, length_stddev, N))

			# Populate table
			for n in range(N):
				# Get start and end values
				# TODO
				c=f*4
				output[n, c:c+4] = self.split_to_32bit_parts(int(start_values[n]), 4)
				output[n, (c+F):(c+F+4)] = self.split_to_32bit_parts(int(start_values[n] + lengths[n]), 4)

		# Update the rule table of this
		self.rule_table=output


	@staticmethod
	def get_format(filename):
		""" Returns a ruleset file format. """

		# Check whether the file is binary by checking the first 16 bytes contain \0
		# Note: the first byte in the binary format is 0 (from version 1)
		if not os.path.isfile(filename):
			raise ValueError('Ruleset file %s not found' % filename)
		with open(filename, 'rb') as f:
			if 0 in f.read(16): return 'Binary'

		# In case of a text file
		# In case of a textual file
		handler = RulesetHandler(filename)
		with open(filename, 'r') as f:
			l = f.readlines()[0].strip()
		if (l[0] == '@'): return 'Classbench'
		elif (l[0] == '{'): return 'JSON'
		else: return 'Classbench-ng'


	@staticmethod
	def read(filename, verbose=0):
		""" Returns an instance of the relevant handler according to format """

		fmt=RulesetHandler.get_format(filename)
		if fmt=='Classbench': return ClassbenchRuleHandler(filename, verbose)
		elif fmt=='Classbench-ng': return ClassbenchNgRuleHandler(filename, verbose)
		elif fmt=='JSON': return JSONRuleHandler(filename, verbose)
		elif fmt=='Binary': return BinaryRuleHandler(filename, verbose)


class BinaryRuleHandler(RulesetHandler):
	""" Handle binary format of rule-sets """

	def __init__(self, filename, verbose=0):
		super().__init__(filename, verbose)
		with open(filename, 'rb') as f:
			data = f.read()

		# See https://docs.python.org/3/library/struct.html#struct.pack
		# for struck pack format

		# Local methods to unpack numbers in little-endian format
		idx={'x':0}
		def read_uint32():
			idx['x']+=4
			return struct.unpack('<I', data[idx['x']-4:idx['x']])[0]

		# Read version
		N = read_uint32()
		version = 0 if (N != 0) else read_uint32()

		if version == 0:
			# Version 0 supports only 5-tuple fields
			F = 5
		elif version == 1:
			# Version 1 supports any field number
			N = read_uint32()
			F = read_uint32()

		self.rule_table=np.zeros([N, 2*F+1], dtype=np.uint32)
		for r in range(N):
			self.rule_table[r, -1] = read_uint32()
			for f in range(F):
				self.rule_table[r, f] = read_uint32()
				self.rule_table[r, f+F] = read_uint32()



class JSONRuleHandler(RulesetHandler):
	""" Handle rule-set of custom JSON format (using logical conjunction and disjunction) """

	def __init__(self, filename, verbose=0):
		""" Reads a custom JSON rulset file format and return a table of 32bit fields (longer fields are split to several fields) """

		super().__init__(filename, verbose)

		# Parse metadata and update attributes
		with open(filename, 'r') as f:
			filedata = f.readlines()
		data=json.loads('\n'.join(filedata))

		# Store how many bits each field requires
		self.field_bits = {}

		# Build the field dict of this
		column=0
		for item in data["fields"]:
			num_of_parts = math.ceil(int(item["size"])/32)
			self.field_dict[item["name"]] = list(range(column, column+num_of_parts))
			self.field_bits[item["name"]] = int(item["size"])
			column += num_of_parts

		# Number of 32bit fields of this
		F=column

		# Read all expanded-rules as dictionaries
		ruleset=[]
		for rule in data["rules"]:
			current_subset = self.parse_rule(rule["match"])
			for x in range(len(current_subset)):
				 current_subset[x]['prio'] = rule['prio']
			ruleset.extend(current_subset)

		N=len(ruleset)

		# Create a double-precision table of columns for gathering statistics
		K=len(data["fields"])
		stat_table=np.zeros(shape=[N,K*2] , dtype=np.float64)
		for row_idx, rule in enumerate(ruleset):
			for field_idx, field in enumerate(data["fields"]):
				# Skip fields that are not in rule
				if field['name'] not in rule.keys(): continue
				# Extract values
				if rule[field['name']]['type'] == 'value':
					# The current field holds mask-value tuple
					bits=self.field_bits[field['name']]
					value = rule[field['name']]['value']
					mask = rule[field['name']]['mask']
					# Populate table
					start, end = self.calculate_start_end_long(value, mask, bits)
				else:
					# The current field holds start-end tuple
					start = int( rule[field['name']]['start'])
					end = int( rule[field['name']]['end'])
				stat_table[row_idx, field_idx] = start
				stat_table[row_idx, field_idx+K] = end

		# Collect & print statistics
		stat_table = self.collect_statistics(stat_table)
		_log(verbose, 'JSON rules statistics:\n')
		for f in range(stat_table.shape[0]):
			_log(verbose, 'Field %d statistics: start-unique: %f, end-unique: %f, lengths-mean: %f, lengths-stddev: %f, start-mean: %f, start-stddev: %f\n' %
				(stat_table[f,0], stat_table[f,1], stat_table[f,2], stat_table[f,3], stat_table[f,4], stat_table[f,5], stat_table[f,6]))
		_log(verbose, '----------------------\n')

		# Create an empty wildcard table
		self.rule_table=np.zeros(shape=[N,2*F+1], dtype=np.uint32)
		self.rule_table[:, F:2*F-1] = 0xffffffff

		# Build the exact rule-table
		for row_idx, rule in enumerate(ruleset):
			# Set the priority
			self.rule_table[row_idx, -1] = rule['prio']
			for field_name, field_data in rule.items():
				# Skip field in case its the rule's priority
				if field_name =='prio': continue

				parts = len(self.field_dict[field_name])
				# Get start & end range values per 32bit field
				if field_data['type'] == 'value': ranges=self.calculate_start_end_32bit_parts(field_data['value'], field_data['mask'], parts)
				else: ranges=self.split_start_end_32bit_parts(field_data['start'], field_data['end'], parts)

				for i, current in enumerate(ranges):
					col = self.field_dict[field_name][i]
					self.rule_table[row_idx, [col, col+F]] = current

		_log(verbose, 'Final rule-table statistics:\n')
		stat_table = self.collect_statistics(self.rule_table)
		for f in range(stat_table.shape[0]):
			_log(verbose, 'Field %d statistics: start-unique: %f, end-unique: %f, lengths-mean: %f, lengths-stddev: %f, start-mean: %f, start-stddev: %f\n' %
				(stat_table[f,0], stat_table[f,1], stat_table[f,2], stat_table[f,3], stat_table[f,4], stat_table[f,5], stat_table[f,6]))
		_log(verbose, '----------------------\n')


	def convert_rules(self, rule_array):
		""" Converts rules from a representation of <value, mask> to <start, end>.
			Field values are stored as Python integers

		Args:
			rule_array: A list of dicts, each holds the <value, mask> attributes of the classifier's fields
		"""
		output = []
		for rule in rule_array:
			pass

	def logical_disjunction(self, rule):
		""" Perform a logical disjunction on the arguments of the given rule """
		output=[]
		for item in rule['args']:
			output.extend(self.parse_rule(item))
		return output


	def logical_conjunction(self, rule):
		""" Perform a logical conjunction on the arguments of the given rule """
		output=[]
		for item in rule['args']:
			# Do a conjunction on both lists (A+B)(C+D) = AC + BC + AD + BD
			current_list=self.parse_rule(item)
			if len(output) == 0:
				output = current_list
			else:
				new_list = []
				for x in output:
					for y in current_list:
						new_list.append({**x, **y})
						output=new_list
		return output


	def parse_rule(self, rule):
		""" Parses a rule, and returls a list of expanded rules with priority """
		try:
			# In case of a literal
			if 'op' not in rule.keys():
				if 'range' not in rule:
					obj = {
						'type': 'value',
						'value': int(rule['value']),
						'mask': int(rule['mask'])
					}
				else:
					obj = {
						'type': 'range',
						'start': int(rule['range'][0]),
						'end': int(rule['range'][0])
					}
				return [{rule['field']: obj}]
			elif rule['op'] == 'conj': return self.logical_conjunction(rule)
			else: return self.logical_disjunction(rule)
		except Exception as e:
			raise Exception('Error while parsing JSON rule: %s.' % rule)


class ClassbenchRuleHandler(RulesetHandler):
	""" Handles rule-sets written in Classbench format """

	def __init__(self, filename, verbose=0):
		""" Parse the Classbench file """

		super().__init__(filename, verbose=0)

		with open(filename, 'r') as f:
			filedata = f.readlines()

		# Prepare the rule table
		F=0
		N=len(filedata)

		# Parse each line
		for i, line in enumerate(filedata):

			# Used for fast Regex group extraction
			idx={'val': 0}
			def e():
				idx['val']+=1
				return reg_match.group(idx['val'])
			def g():
				return int(e(), 0)

			# Extract an IPv4 range from the regex match
			def get_ipv4_range():
				ip_address = ( (g()<<24) | (g()<<16) | (g()<<8) | (g()<<0) )
				ip_mask_bits=g()
				ip_mask = 0xffffffff << (32-ip_mask_bits) & 0xffffffff
				ip_start = ip_address & ip_mask
				ip_end = ip_start | ~ip_mask
				return [ip_start, ip_end]

			# Extract an IPv6 range from the regex match
			def get_ipv6_range():
				ip_string = re.sub(':', '', e())
				if len(ip_string)==0: ip_string='0'
				ip_address=int(ip_string, 16) << (128 - len(ip_string) * 4)
				ip_mask_bits=g()
				ip_mask=(2**ip_mask_bits)-1 << (128 - ip_mask_bits)
				ip_start = ip_address & ip_mask
				ip_end = ip_start | ~ip_mask
				return [ip_start, ip_end]

			# Try to parse line as IPv4
			reg_match=re.match(r'@([0-9]+)\.([0-9]+)\.([0-9]+)\.([0-9]+)\/([0-9]+)\s*([0-9]+)\.([0-9]+)\.([0-9]+)\.([0-9]+)\/([0-9]+)\s*([0-9]+)\s\:\s([0-9]+)\s*([0-9]+)\s\:\s([0-9]+)\s*(0x[0-9a-fA-F]+)/(0x[0-9a-fA-F]+)', line)
			if reg_match is not None:
				src_ip_start, src_ip_end = get_ipv4_range()
				dst_ip_start, dst_ip_end = get_ipv4_range()
				# Initiate the tables, if there are not yet initialized (only now we know which IP version the file holds)
				if F == 0:
					F=5
					self.rule_table=np.zeros(shape=[N,2*F+1], dtype=np.uint32)
				elif F != 5:
					raise ValueError('Cannot parse Classbench with both ipv6 and ipv4 fields!')

			# Try to parse line as IPv6
			else:
				reg_match = re.match(r'@([0-9a-fA-F:]+)\/([0-9]+)\s*([0-9a-fA-F:]+):\/([0-9]+)\s*([0-9]+)\s\:\s([0-9]+)\s*([0-9]+)\s\:\s([0-9]+)\s*(0x[0-9a-fA-F]+)\/(0x[0-9a-fA-F]+)', line)
				src_ip_start, src_ip_end = get_ipv6_range()
				dst_ip_start, dst_ip_end = get_ipv6_range()
				# Initiate the tables, if there are not yet initialized (only now we know which IP version the file holds)
				if F == 0:
					F=11
					self.rule_table=np.zeros(shape=[N,2*F+1], dtype=np.uint32)
				elif F != 11:
					raise ValueError('Cannot parse Classbench with both ipv6 and ipv4 fields!')

			# Extract all other fields
			src_port_start = g()
			src_port_end = g()
			dst_port_start = g()
			dst_port_end = g()

			protocol = g()
			protocol_mask=g()

			protocol_start = protocol & protocol_mask
			protocol_end = protocol_start | ~protocol_mask & 0xff

			# Update rule table according to IP version
			if F == 5:
				self.rule_table[i,:] = [src_ip_start, dst_ip_start, src_port_start, dst_port_start, protocol_start,
							   			src_ip_end, dst_ip_end, src_port_end, dst_port_end, protocol_end, i]
			elif F == 11:
				self.rule_table[i,0:4] = self.split_to_32bit_parts(src_ip_start, 4) # 4 32bit fields for an IPv6 address
				self.rule_table[i,4:8] = self.split_to_32bit_parts(dst_ip_start, 4) # 4 32bit fields for an IPv6 address
				self.rule_table[i,8:11] = [src_port_start, dst_port_start, protocol_start]
				self.rule_table[i,11:15] = self.split_to_32bit_parts(src_ip_end, 4) # 4 32bit fields for an IPv6 address
				self.rule_table[i,15:19] = self.split_to_32bit_parts(dst_ip_end, 4) # 4 32bit fields for an IPv6 address
				self.rule_table[i,19:22] = [src_port_end, dst_port_end, protocol_end]
				self.rule_table[i,22] = i
			else:
				raise Exception('Number of fields for Classbench ruleset was not initialized.')



class ClassbenchNgRuleHandler(RulesetHandler):
	""" Handle rule set with classbench-ng format """

	def __init__(self, filename, verbose=0):
		super().__init__(filename, verbose=0)

		# The known fields of Classbnech-ng and their sizes in 32bit parts
		field_order=['dl_src', 'dl_dst', 'eth_type', 'in_port', 'nw_src', 'nw_dst', 'nw_proto', 'tp_dst', 'tp_src']
		known_fields={
			'dl_src': 		{'action': 'mac', 'size': 2},
			'dl_dst': 		{'action': 'mac', 'size': 2},
			'eth_type': 	{'action': 'hex', 'size': 1},
			'in_port':		{'action': 'int', 'size': 1},
			'nw_src':		{'action': 'ip4', 'size': 1},
			'nw_dst':		{'action': 'ip4', 'size': 1},
			'nw_proto':		{'action': 'int', 'size': 1},
			'tp_dst':		{'action': 'int', 'size': 1},
			'tp_src':		{'action': 'int', 'size': 1},
		}

		# Initialize the field dict of this
		column=0
		for key in field_order:
			val=known_fields[key]
			self.field_dict[key]=list(range(column, column+val['size']))
			column+=val['size']

		with open(filename, 'r') as f:
			filedata = f.readlines()

		# Craete an empty rule-table
		F = 11
		N = len(filedata)

		# Create an empty wildcard table
		self.rule_table=np.zeros(shape=[N,2*F+1], dtype=np.uint32)
		self.rule_table[:, F:2*F-1] = 0xffffffff
		self.rule_table[:, -1] = np.arange(N)

		for row_index, line in enumerate(filedata):
			# Get all attributes
			attributes = re.sub(r'\s+', '', line).split(',')
			for attr in attributes:
				(field_name, val) = attr.split('=')
				if field_name not in known_fields.keys():
					raise ValueError('Field %s is not recognized' % field_name)

				# Set the value according to the action
				if known_fields[field_name]['action'] == 'mac':   value=self.prase_mac_address(val)
				elif known_fields[field_name]['action'] == 'hex': value=[int(val,0), int(val,0)]
				elif known_fields[field_name]['action'] == 'int': value=[int(val), int(val)]
				elif known_fields[field_name]['action'] == 'ip4': value=self.pase_ip_mask(val)

				# Get the number of 32bit parts
				parts = known_fields[field_name]['size']
				assert(parts == len(self.field_dict[field_name]))

				# The value may be a long integer. Split to 32bit parts
				low_values =self.split_to_32bit_parts(value[0], parts)
				high_values=self.split_to_32bit_parts(value[1], parts)

				# Update the table
				for i, col in enumerate(self.field_dict[field_name]):
					self.rule_table[row_index, [col,col+F]] = [low_values[i], high_values[i]]


	def pase_ip_mask(self, s):
		""" Parses IP/mask string and returns it as 32bit integer range """
		reg_match=re.match(r'([0-9]+)\.([0-9]+)\.([0-9]+)\.([0-9]+)\/([0-9]+)', s)
		if reg_match is None:
			raise ValueError('Cannot parse IP/mask %s' % s)

		# Used for fast Regex group extraction
		idx={'val': 0}
		def g():
			idx['val']+=1
			return int(reg_match.group(idx['val']), 0)

		ip_address = ( (g()<<24) | (g()<<16) | (g()<<8) | (g()<<0) )
		ip_mask_bits=g()
		ip_mask = 0xffffffff << (32-ip_mask_bits) & 0xffffffff

		ip_start = abs(ip_address & ip_mask)
		ip_end = abs(ip_start | ~ip_mask)
		return [ip_start, ip_end]


	def prase_mac_address(self, s):
		""" Parses HEX mac address and return it as an integer """
		parts=s.lower().split(':')
		out = 0
		for p in parts:
			out = (out << 8) | int(p, 16)
		out=abs(out)
		return [out, out]


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
import os
import struct
import random

# Add the "bin" directory to import path
my_path=os.path.dirname(os.path.realpath(__file__))
bin_path=os.path.join(my_path, '..', 'bin')
if not os.path.exists(bin_path):
    print('Bin directory was not found. Did you compile the library?')
    exit(1)
sys.path.append(bin_path)

# Import the RQRMI library
try: from lookup_cpu import LookupCpu
except ModuleNotFoundError as e:
    # Cannot load a prerequisites module
    print(e)
    exit(1)
except ImportError:
    # Cannot find the rqrmi_library module
    print('RQRMI library was not found is bin directory. Did you compile the library?')
    exit(1)

import numpy as np

# ======== Script Start ======== #

# Initialize objects
validation_dict = {}
incorrect_lookup_values = 0

# Callback method for Lookup
def on_lookup_cpu_new_result(key, output, index, found):
    """ Is invoked by Lookup each time there is a new available result

    Args:
        key: The requested key to search
        output: The found value
        index: The pair index
        found: Whether a valid pair was found
    """
    global incorrect_lookup_values
    global lookup
    print('yer', flush=True)
    if found: correct = (validation_dict[key][1] == output)
    else: correct = (key not in validation_dict.keys())
    if not correct: incorrect_lookup_values+=1


# Randomize lognormal database for ranges
num_of_records=int(5e4)
max_range_value=3e9

print('** Generating random lognormal database with %d records' % num_of_records)
database = np.random.lognormal(0, 1, num_of_records)
# Normalize lognormal from 0 to max_range_value
database -= np.min(database)
database /= np.max(database)
database *= max_range_value
# Remove non unique uint32 values
database = np.unique(database.astype(np.float32))
num_of_records = database.shape[0] - 1

# Create the range-value database
print('** Creating range-value pair database')
range_value_pairs = np.empty([num_of_records, 3])
range_value_pairs[:,0] = database[:-1]
range_value_pairs[:,1] = database[1:]

print('** Generating validation set')
to_float32 = lambda num: struct.unpack('<f', struct.pack('<f', num))[0]
for i in range(num_of_records):

    # Validate database correctness
    assert(range_value_pairs[i,0] <= range_value_pairs[i,1])

    # Insert random sample to validation dictionary
    start = np.ceil(range_value_pairs[i,0]).astype(np.float32)
    end = np.ceil(range_value_pairs[i,1]).astype(np.float32)

    while True:
        if end > start:
            key = to_float32(np.random.randint(start, end))
            if (key >= range_value_pairs[i,0] and key < range_value_pairs[i,1]):
                break
        else:
            key = start
            break

    # Append sample to validation set
    validation_dict[key] = (i, range_value_pairs[i,1])


print('** Initializing Lookup')
lookup = LookupCpu(num_of_records)
for i in range(num_of_records):
    lookup.insert(range_value_pairs[i,0], range_value_pairs[i,1])

print('** Training LookupCpu')
lookup.compile()
lookup.start(8, on_lookup_cpu_new_result)

file='lookup.data'
print('*** Saving database to file %s' % file)
with open(file, 'wb') as f:
    f.write(lookup.pack())

print('** Performing infinite lookup. Interrupt to stop.', flush=True)
try:

    num_of_search=0
    sys.stdout.write('')

    while (1):
        # Request random validation sample
        key = random.choice(list(validation_dict.keys()))
        lookup.search(key)
        num_of_search+=1

        # Let Lookup process internal messages
        if num_of_search % 1000 == 0:
            sys.stdout.write('\rNumber of searches: %d. Incorrect lookups: %d. Correctness percent: %.2f' %
                             (num_of_search, incorrect_lookup_values, (1-incorrect_lookup_values/num_of_search)*100))

except KeyboardInterrupt as e:
    print('\nGot keyboard interrupt, exiting script')


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
from asyncore import read

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
    from object_reader import ObjectReader
    from rqrmi_library import RQRMI
except ImportError as e:
    # Cannot find the library
    print('One of the library files was not found (%s). Did you compile the library?' % e)
    exit(1)


def parse_arguments():
    """ Parse script arguments
    """

    parser = argparse.ArgumentParser()

    # Mandatory arguments
    parser.add_argument('-f', '--file', required=True, help='Filename of valid NuevoMatch classifier file')
    parser.add_argument('-o', '--out',  required=True, help='Output filename')

    # Control Modes
    parser.add_argument('--extract', type=int, default=None, help='Extract an RQRMI model from a specific iSet')

    # In case of no arguments, print usage
    if len(sys.argv)<=1:
        parser.print_usage()
        parser.exit()

    return parser.parse_args()

# ================================================ #
# ================ Script Start ================== #
# ================================================ #

args=parse_arguments()

# Read input file
with open(args.file, 'rb') as f:
    reader = ObjectReader(f.read())

# Extract RQRMI models
if args.extract is not None:

    # Extract classifier properties
    print('Reading classifier file...')
    num_of_isets = reader.read_uint32()
    num_of_rules = reader.read_uint32()
    size = reader.read_uint32()
    build_time = reader.read_uint32()

    # Get the desired iSet number
    print('Extracting model from iSet %d...' % args.extract)
    for _ in range(args.extract+1):

        # Extract the packed objects
        iset_data = reader.read_object()

        # Read the packed model
        packed_model = iset_data.read_object()

    # Write the output model
    print('Writing extracted model to file...')
    with open(args.out, 'wb') as f:
        f.write(bytes(packed_model))



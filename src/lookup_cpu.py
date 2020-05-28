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

# Load the RQRMI library
try:
    import rqrmi as rqrmilib
    from rqrmi_library import RQRMI, initialize_library, _pack_numpy_matrix, _numpy_2_native_matrix, _native_matrix_2_numpy, _log
except ImportError:
    raise ImportError('Cannot import RQRMI library. Make sure to compile it using make.')

# Load prerequisites
try:
    import numpy as np
    import tensorflow as tf
    import struct
    import sys
except ImportError as e:
    raise ModuleNotFoundError('Cannot load module %s. Check for prerequisites.' % e.name)


# Suppress TF warning
os.environ['TF_CPP_MIN_LOG_LEVEL']='3'
tf.logging.set_verbosity(tf.logging.ERROR)

# TF work on CPU only (faster for small NN)
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

class LookupCpu:

    def __init__(self, size, verbose=1):

        # Hyper Parameters
        self.submodel_structure=[1, 8, 1]
        self.epochs=[10, 10, 20]
        self.batch_size=32
        self.samples_per_bucket=1500
        self.threshold_submodel_error=64
        self.threshold_bucket_coverage = 0.5
        self.retraining_multiplier = 3
        self.retraining_times = 3

        # Private objects
        self.model = None
        self.lookup = None
        self.database = np.empty([size, 2], dtype=np.float32)
        self.current = 0
        self.verbose = verbose

        # Initialize the RQRMI library, in case not initialized already
        initialize_library(False)


    def insert(self, range_start, range_end):
        """ Insert new range-value pair

        Args:
            range_start: The range start
            range_end: The range end
            value: The corresponding value
        """
        if range_end < range_start:
            raise ValueError('Invalid range [%f, %f]' % (range_start, range_end))
        self.database[self.current, :] = [range_start, range_end]
        self.current += 1


    def extend(self, db):
        """ Extend the range-value pair database with an external database

        Args:
            db: A Numpy matrix (Nx2) where N is the number of new records.
                Each row with the format [range-start, range-end)

        Throws:
            ValueError in case of invalid input
        """

        if type(db) is not np.ndarray:
            raise ValueError('DB argument is invalid')
        if db.shape[1] != 2:
            raise ValueError('DB argument is invalid')

        N = db.shape[0]
        self.database[self.current:self.current+N, :] = db
        self.current += N


    def compile(self):
        """ Compile Lookup: create RQRMI model.

        Throws:
            ValueError in case of overlapping ranges
            RuntimeError in case of internal library error
        """

        # Sort database by range start
        self.database = self.database[np.lexsort([self.database[:,0]]), :]

        # Create an index
        self.index = np.concatenate([self.database[:, 0], [self.database[-1,1]] ])
        num_of_records = self.index.shape[0]

        # Validate non overlapping ranges
        test_vector = (np.roll(self.index, -1) - self.index) <= 0
        test_vector = test_vector[:-1]
        if np.logical_or.reduce(test_vector):
            raise ValueError('Cannot index overlapping ranges')

        # Set stage structure
        if num_of_records < 1e3: stages_strcture = [1, 4]
        elif num_of_records < 1e4: stages_strcture = [1, 4, 16]
        elif num_of_records < 1e5: stages_strcture = [1, 4, 128]
        else: stages_strcture = [1, 8, 256]

        # Train the model
        self.model=RQRMI(self.index, stages_strcture)
        for i,stage in enumerate(stages_strcture):

            self.model.calculate_responsibility(i, verbose=self.verbose)
            submodels_to_train = np.arange(stage).astype(int)
            current_samples_per_bucket = self.samples_per_bucket

            # Repeat the training process
            for t in range(self.retraining_times):

                _log(self.verbose, '** Starting stage training iteration %d (out of %d max) **\n' % (t, self.retraining_times))

                buckets = self.model.generate_buckets(i, submodels_to_train, current_samples_per_bucket, verbose=self.verbose)
                self.model.train_stage(i, self.submodel_structure, buckets, submodels_to_train, None, self.epochs[i], self.batch_size, verbose=self.verbose)
                self.model.calculate_transition_set(i, verbose=self.verbose)

                # Measure error for refinement
                submodels_to_train = []
                error_list=[]
                for j in range(stage):
                    max_error, coverage = self.model.calculate_submodel_error(i, j)
                    # Retrain shallow submodels based on their bucket coverage
                    if i < len(stages_strcture)-1 and (coverage < self.threshold_bucket_coverage):
                        submodels_to_train.append(j)
                        error_list.append(coverage)
                    # Retrain deep submodels based on their maximum error
                    if i == len(stages_strcture)-1 and (max_error > self.threshold_submodel_error):
                        submodels_to_train.append(j)
                        error_list.append(max_error)

                # In case no need to retrain
                if (len(submodels_to_train) == 0): break

                # Print error values
                _log(self.verbose, 'Measured errors/bucket-coverages: %s\n' % \
                    str([ str(x[0]) + ': ' +  str(x[1]) for x in zip(submodels_to_train, error_list) ]))

                # Update number of samples (more samples = less error = more training time)
                current_samples_per_bucket = int(current_samples_per_bucket*self.retraining_multiplier)

        # Pin the RQRMI error values
        error_list = self.model.pin_errors()
        _log(self.verbose, 'RQRMI error list: %s \n' % error_list)


    def start(self, queue_size, callback):
        """ Starts listening to search results

        Args:
            queue_size: Worker thread queue size
            callback: A callback to handle results

        Throws:
            RuntimeError in case of internal library error
        """

        # Set the Lookup object
        self.lookup = rqrmilib.create_lookup_cpu(
            _numpy_2_native_matrix(self.database),
            self.model._get_native_object(),
            queue_size,
            callback)


    def search(self, input):
        """ Search for enclosing range-value pair

        Args:
            input: What to search

        Throws:
            Exception in case this was not started

        Note: blocking
        """
        if self.lookup is None:
            raise Exception('Lookup listener was not started')

        rqrmilib.lookup_search(self.lookup, input)


    def pack(self):
        """ Packs this to byte array """

        output=bytearray()
        write_uint32 = lambda num: output.extend(struct.pack('<I', num))

        def add_object(obj):
            write_uint32(len(obj))
            output.extend(obj)

        # Pack the RQRMI model
        packed_model =self.model.pack()
        _log(self.verbose, 'RQRMI model occupies %d bytes\n' % len(packed_model))
        add_object(packed_model)

        # Pack the database
        add_object( _pack_numpy_matrix(self.database) )

        return output


    def get_size(self):
        """ Returns the RQRMI size of this, in bytes"""
        return len(self.model.pack())








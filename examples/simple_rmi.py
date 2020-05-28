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

# Simple RMI model (Learned Index Structure) example with lognormal record distribution.
# Paper: https://dl.acm.org/citation.cfm?id=3196909
# This implementation uses the RQRMI library to train & evaluate a simple RMI model

import sys
import os

# Add the "bin" directory to import path
my_path=os.path.dirname(os.path.realpath(__file__))
bin_path=os.path.join(my_path, '..', 'bin')
if not os.path.exists(bin_path):
    print('Bin directory was not found. Did you compile the library?')
    exit(1)
sys.path.append(bin_path)

# Import the RQRMI library
try: from rqrmi_library import RQRMI, initialize_library, EPS
except ModuleNotFoundError as e:
    # Cannot load a prerequisites module
    print(e)
    exit(1)
except ImportError:
    # Cannot find the rqrmi_library module
    print('RQRMI library was not found is bin directory. Did you compile the library?')
    exit(1)

# Import dependencies
from timeit import default_timer as timer
import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np

# ======== Script Start ======== #

# Suppress TF warning
os.environ['TF_CPP_MIN_LOG_LEVEL']='3'
tf.logging.set_verbosity(tf.logging.ERROR)

# TF work on CPU only (faster for small NN)
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'

# Set no library debug printing
initialize_library(0)

# Parameters for dataset
num_of_records=int(1e6)
max_uint32=np.iinfo(np.uint32).max

print('** Generating random lognormal database with %d records' % num_of_records)
database = np.random.lognormal(0, 1, num_of_records)
# Normalize lognormal from 0 to max uint32
database -= np.min(database)
database /= np.max(database)
database *= max_uint32
# Remove non unique uint32 values
database = np.unique(database.astype(np.uint32))
num_of_records = database.shape[0]

print('** Generating dataset for training the model. Number of records: %d' % num_of_records)
dataset = np.zeros(shape=[num_of_records, 2]).astype(np.float32)
dataset[:, 0] = database
dataset[:, 1] = np.arange(num_of_records) / num_of_records
dataset_indices = np.arange(num_of_records)

# Parameters for RMI
submodel_structure=[1, 8, 1]
stages_strcture=[1, 4, 16]
epochs=[10, 10, 20]
batch_size=32
samples_per_bucket=int(1.5e3)
optimizer=tf.train.AdamOptimizer(learning_rate=1e-3)

# The model
print('** Initializing RMI model')
model = RQRMI(database, stages_strcture)

# Used for illustrating outputs
fig, ax = plt.subplots(nrows=len(stages_strcture)+1, ncols=2, figsize=(10, 10))

print('** Start training')
traiaing_time_seconds=0

# Random sample from dataset
sampled_indices = np.random.choice(dataset_indices, size=samples_per_bucket)
buckets = [dataset[sampled_indices]]

# Train all stages
for i,stage in enumerate(stages_strcture):

    # Plot the input buckets for the current stage
    print('** Plot input for stage %d' % i)
    ax[i][1].set_xlabel('Input')
    ax[i][1].set_ylabel('Expected output')
    ax[i][1].set_title('Input for stage %d Divided By Buckets' % i)
    for j,b in enumerate(buckets):
        ax[i][1].scatter(b[:,0], b[:,1], s=3)

    print('** Builing stage %d' % i)

    # Measure training time
    time_start = timer()

    submodels_to_train = np.arange(stage).astype(int)
    model.train_stage(i, submodel_structure, buckets, submodels_to_train, optimizer, epochs[i], batch_size, verbose=1)

    print('** Evaluating stage %d' % i)
    output = model.evaluate(database, verbose=1)

    # Split dataset to buckets based on stage's output
    if i < len(stages_strcture)-1:
        print('** Building buckets for stage %d' % (i+1))
        required_size=stages_strcture[i+1]

        # Assign outputs to buckets
        bucket_output=np.minimum(1-EPS, np.maximum(0, output))*required_size
        bucket_output=bucket_output.astype(int)

        # Build the buckets
        buckets=[]
        for x in range(required_size):
            indices=dataset_indices[bucket_output==x]
            # In case the current bucket is empty
            if indices.shape[0] == 0:
                buckets.append(np.empty(shape=[0,2]))
                print('** Bucket %d for stage %d is empty' % (x, i+1))
            # Else, add new bucket
            else:
                sampled_indices=np.random.choice(indices, size=samples_per_bucket)
                buckets.append(dataset[sampled_indices])
                # Print bucket
                print('** Bucket %d for stage %d holds records from %d to %d with %d samples' %
                      (x, i+1, indices[0], indices[-1], samples_per_bucket))


    # Measure time
    traiaing_time_seconds += (timer() - time_start)

    # Plot figure for current stage
    print('** Plotting figure for stage %d' % i)
    ax[i][0].plot(dataset[:,0], dataset[:,1])
    ax[i][0].plot(dataset[:,0], output)
    ax[i][0].legend(('Database', 'RMI output'))
    ax[i][0].set_xlabel('Input')
    ax[i][0].set_ylabel('Output')
    ax[i][0].set_title('Database vs. RMI outputs for stage %d' % i)


print('** Total training time: %.f seconds' % traiaing_time_seconds)

# Calculate misprediction error per record
print('** Calculating misprediction error per record')
error=np.abs(dataset[:,1] - output)
error_percent=np.ceil(error*100)
ax[-1][0].plot(dataset[:,0], error_percent)
ax[-1][0].set_xlabel('Record')
ax[-1][0].set_ylabel('Misprediction Error (percent)')
ax[-1][0].set_title('Misprediction error per record')

max_error = np.max(error)
print('** Max misprediction error: %d%% (%d records)' % (np.max(error_percent), max_error))

# Save the figures to file
figure_name='rmi_output.png'
print('** Saving output figure to %s...' % os.path.join(os.getcwd(), figure_name))
plt.tight_layout()
fig.savefig(figure_name)

# Save the model to file
model_file='rmi_model.model'
print('** Saved model to %s...' % model_file)
with open(model_file, 'wb') as f:
    f.write(model.pack())

print('** Done. It can be loaded with the RQRMI benchmark utility or with Python code')

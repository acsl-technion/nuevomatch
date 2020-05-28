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

# Simple RQRMI model example with lognormal record distribution.

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
try: from rqrmi_library import RQRMI, initialize_library
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
initialize_library(False)

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

print('** Generating dataset for plotting. Number of records: %d' % num_of_records)
dataset = np.arange(num_of_records) / num_of_records
dataset_indices = np.arange(num_of_records)

# Parameters for RQRMI
submodel_structure=[1, 8, 1]
stages_strcture=[1, 4, 16]

epochs=[10, 10, 20]
batch_size=32
optimizer=tf.train.AdamOptimizer(learning_rate=1e-3)

samples_per_bucket=int(1.5e3)

threshold_submodel_error=0.005
threshold_bucket_coverage = 0.5
retraining_multiplier = 1.5
retraining_times = 3

# The model
print('** Initializing RQRMI model')
model = RQRMI(database, stages_strcture)

# Used for illustrating outputs
fig, ax = plt.subplots(nrows=len(stages_strcture)+1, ncols=2, figsize=(10, 10))

print('** Start training')
traiaing_time_seconds=0

# Train all stages
for i,stage in enumerate(stages_strcture):

    print('** Building stage %d' % i)

    # Measure training time
    time_start = timer()

    print('** Calculating responsibilities')
    model.calculate_responsibility(i, verbose=1)

    submodels_to_train = np.arange(stage).astype(int)
    buckets_to_plot = [None for _ in submodels_to_train]
    current_samples_per_bucket = samples_per_bucket

    # Repeat the training process
    for _ in range(retraining_times):

        buckets = model.generate_buckets(i, submodels_to_train, current_samples_per_bucket, verbose=1)
        model.train_stage(i, submodel_structure, buckets, submodels_to_train, optimizer, epochs[i], batch_size, verbose=1)
        model.calculate_transition_set(i, verbose=1)

        # Update the buckets to plot
        for j,x in enumerate(submodels_to_train):
            buckets_to_plot[x] = buckets[j]

        # Measure error for refinement
        submodels_to_train = []
        for j in range(stage):
            max_error, coverage = model.calculate_submodel_error(i, j)
            # Retrain shallow submodels based on their bucket coverage
            if i < len(stages_strcture)-1 and (coverage < threshold_bucket_coverage or coverage > 1):
                submodels_to_train.append(j)
            # Retrain deep submodels based on their maximum error
            if i == len(stages_strcture)-1 and (max_error/num_of_records > threshold_submodel_error):
                submodels_to_train.append(j)

        # In case no need to retrain
        if (len(submodels_to_train) == 0): break

        # Update number of samples (more samples = less error = more training time)
        current_samples_per_bucket = int(current_samples_per_bucket*retraining_multiplier)

    traiaing_time_seconds += (timer() - time_start)

    # Plot the input buckets for the current stage
    print('** Plot input for stage %d' % i)
    ax[i][1].set_xlabel('Input')
    ax[i][1].set_ylabel('Expected output')
    ax[i][1].set_title('Input for stage %d Divided By Buckets' % i)
    for j,b in enumerate(buckets_to_plot):
        ax[i][1].scatter(b[:,0], b[:,1], s=3)

    print('** Evaluating stage %d (with %d samples)' % (i, num_of_records))
    output = model.evaluate(database, verbose=1)

    # Plot figure for current stage
    print('** Plotting figure for stage %d' % i)
    ax[i][0].plot(database, dataset)
    ax[i][0].plot(database, output)
    ax[i][0].legend(('Database', 'RQRMI output'))
    ax[i][0].set_xlabel('Input')
    ax[i][0].set_ylabel('Output')
    ax[i][0].set_title('Database vs. RQRMI outputs for stage %d' % i)

print('** Total training time: %.f seconds' % traiaing_time_seconds)

# Calculate misprediction error per record
print('** Calculating misprediction error per submodel')
submodel_error=[]
for j in range(stages_strcture[-1]):
    max_error, _ = model.calculate_submodel_error(len(stages_strcture)-1, j)
    submodel_error.append(100*max_error/num_of_records)

# Plot error chart
ax[-1][0].plot(np.arange(stages_strcture[-1]), submodel_error)
ax[-1][0].set_xlabel('Submodel')
ax[-1][0].set_ylabel('Misprediction Error (%%)')
ax[-1][0].set_title('Misprediction Error per Submodel')

max_error=np.max(submodel_error)
print('** Max prediction error: %.2f%% (%d records)' % (max_error, max_error/100*num_of_records))

# Save the figures to file
figure_name='rqrmi_output.png'
print('** Saving output figure to %s...' % os.path.join(os.getcwd(), figure_name))
plt.tight_layout()
fig.savefig(figure_name)

# Save the model to file
model_file='rqrmi_model.model'
print('** Saved model to %s...' % model_file)
with open(model_file, 'wb') as f:
    f.write(model.pack())

print('** Done. It can be loaded with the RQRMI benchmark utility or with Python code')

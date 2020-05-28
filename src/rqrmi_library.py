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

# Load the RQRMI library
try:
    import rqrmi as rqrmilib
except ImportError:
    raise ImportError('Cannot import RQRMI library. Make sure to compile it using make.')

# Load prerequisites
try:
    import numpy as np
    import tensorflow as tf
    import struct
    import sys
    import datetime
except ImportError as e:
    raise ModuleNotFoundError('Cannot load module %s. Check for prerequisites.' % e.name)

# Epsilon
EPS=np.finfo(np.float32).eps

_library_initialzied=False
_log_last_length=0

def _log(verbose, msg, delete_last=False):
    """ Prints log messages according to verbosity to stderr

    Args:
        verbose: integer. message is printed when not zero.
        msg: message to print, newline char should be included
        delete_last: If true, the last message will be deleted from buffer
    """
    global _log_last_length
    if verbose<=0:
        return
    if delete_last:
        sys.stderr.write('\b'*_log_last_length)
    sys.stderr.write(msg)
    sys.stderr.flush()
    _log_last_length=len(msg)


def initialize_library(verbose):
    """ Sets the native library verbosity

        Args:
            verbose: Set to 1 for redirecting native library logs to stderr. For full debug
                     information, compile the RQRMI library with DEBUG flag.
    """

    global _library_initialzied
    if _library_initialzied: return
    if verbose:
        print('*** For maximum debug information on the RQRMI library, compile library with debug flag ***')
        rqrmilib.initialize(1)
    else:
        rqrmilib.initialize(0)
    _library_initialzied=True


def _nn_initial_values(structure):
    """ Returns a list of TF tensors with initial value to place in a NN

    Args:
        structure: A list of integers, the NN layer structure (including input & output layers)
    """

    # Use Xavier uniform initializer
    initializer=tf.glorot_uniform_initializer()

    output=[]
    last_width=None

    # Add biases & weights per layer
    for l in structure:
        output.append(tf.zeros(shape=[l])) # layer l biases
        if last_width is not None: # Exclude weights from layer 0
            output.append(initializer(shape=[last_width, l])) # layer l weights
        last_width=l

    return output


def _nn_create(values, data_in=None):
    """ Builds the infrastructure of fully connected NN using TF variables

    Args:
        values: A list of biases & weights to populate the NN
        data_in: (Optional) Use external input Tensor (do not generate input placeholder)

    Returns:
        A list of TF variables. The first is the input, the last is the NN output.
        All TF variables in between are biases & weights of every layer.

    """

    variables=[]
    num_of_layers=int(len(values)/2)+1 # Two values per layer (biases, weights) except the input layer

    # In case no external input, set first layer as placeholders
    if data_in is None:
        input_width = values[0].shape[0]
        data_in=tf.placeholder(dtype=tf.float32, shape=[None, input_width])

    variables.append(data_in)

    # For each layer
    for i in range(num_of_layers):
        # Input layer
        if i==0:
            b=tf.Variable(dtype=tf.float32, initial_value=values[0])
            layer=tf.add(data_in, b)
            variables.append(b)

        # All other layers
        else:
            b=tf.Variable(dtype=tf.float32, initial_value=values[2*i-1])
            w=tf.Variable(dtype=tf.float32, initial_value=values[2*i-0])
            variables.extend([b, w])
            layer=tf.matmul(layer,w)+b

            # Add ReLU after each hidden layer
            if i<num_of_layers-1: layer=tf.nn.relu(layer)
            # Output layer
            else: variables.append(layer)

    return variables


def _nn_train_net(structure, training_set, optimizer, epochs, batch_size, verbose=1):
    """ Creates and trains fully connected NN with the input parameters using supervised training.

    Args:
        structure: A list of integers, the NN layer structure (including input & output layers)
        training_set: A Nx2 matrix (of uint32) where N is the number of samples.
                      The first column holds the inputs, the last column holds the expected outputs.
        optimizer: TF optimizer to use for training
        epochs: Number of epochs to train
        batch_size: Batch-size for training
        verbose: Verbosity

    Note: The training session uses TF feed_dict mechanism (faster for smaller datasets)

    Returns:
        A list of tuples with the following Items:

        Tuple: statistical data:
            - loss over epoch (list of floats)

        Tuple: the necessary data for inference:
            - input_data mean (named mu)
            - input_data stddev (named sig)
            - output normalization factor (named fact) (output = NN_out * factor + min)
            - minimum output (named omin) (output = NN_out * factor + min)
            - a list of the NN biases & weights
    """

    # Preallocate output variables
    net=[]
    loss_stat=[]

    # Collect training time statistics
    time_start = datetime.datetime.now()

    # Convert input data to be in float 32 (compression)
    # Note: the conversion makes information loss
    norm_data=training_set.copy().astype(np.float32)

    # Calculate input data statistics
    mu=np.mean(norm_data[:,0])
    sig=np.std(norm_data[:,0])
    if sig==0: sig=1
    norm_data[:,0]=(norm_data[:,0]-mu)/sig

    #Calculate the output factor of the net
    min_out=np.min(norm_data[:,-1])
    max_out=np.max(norm_data[:,-1])
    output_factor=(max_out-min_out)
    if output_factor == 0: output_factor=1
    # Normalize the net output to be in [0, 1]
    norm_data[:,-1]=(norm_data[:,-1]-min_out)/output_factor

    # New TF graph
    with tf.Graph().as_default():
        # Request initializer for net structure
        values=_nn_initial_values(structure)
        # Build net
        net_variables=_nn_create(values)
        data_in=net_variables[0]
        data_out=net_variables[-1]
        # Build optimizer
        desired_out=tf.placeholder(dtype=tf.float32, shape=[None, structure[0]])
        loss=tf.losses.mean_squared_error(data_out, desired_out)
        opt=optimizer.minimize(loss)

        # Perform training
        with tf.Session() as sess:
            sess.run(tf.global_variables_initializer())
            for epoch in range(epochs):
                avg_cost=0.0
                total_batch=np.ceil(norm_data.shape[0]/batch_size).astype(int)

                # Split samples to batches
                batches = np.array_split(norm_data, total_batch)

                # For each batch
                for i in range(total_batch):
                    _, c = sess.run([opt, loss],
                                    feed_dict={
                                        data_in:     batches[i][:, 0:1],
                                        desired_out: batches[i][:, 1:2]
                                    })

                    avg_cost += c / total_batch

                # Update statistics
                loss_stat.append(avg_cost)
                _log(verbose, 'Epoch %04d/%d: cost=%.9f' % (epoch+1, epochs, avg_cost), epoch>0)

            # Get the trained net's internal values
            for v in net_variables[1:-1]:
                net.append(sess.run(v))

    # Measure net training time
    time_end = datetime.datetime.now()
    time_diff = (time_end - time_start).total_seconds()

    _log(verbose, ' (mu: %.3f, sig: %.3f, fac: %.3f, omin: %.3f, training time: %.3f sec) \n' % (mu, sig, output_factor, min_out, time_diff))

    # Return all data
    return ( (loss_stat), (mu, sig, output_factor, min_out, net) )


def _pack_numpy_matrix(numpy_mat):
    """ Convert numpy matrix to byte-array

    Args:
        numpy_mat: A numpy vector/matrix. Max 2-dimensions

    Returns:
        Byte array representation of the matrix
    """

    if type(numpy_mat) is not np.ndarray:
        raise ValueError('Input is not a valid Numpy matrix')

    # Pack the vector to byte array
    packed_vec=bytearray()

    # See https://docs.python.org/3/library/struct.html#struct.pack
    # for struck pack format

    # Local methods to pack numbers in little-endian format
    write_uint32 = lambda num: packed_vec.extend(struct.pack('<I', num))
    write_float32 = lambda num: packed_vec.extend(struct.pack('<f', num))

    if len(numpy_mat.shape)==1:
        num_of_rows = numpy_mat.shape[0]
        num_of_cols = 1
        numpy_mat = numpy_mat.reshape([num_of_rows, num_of_cols])
    elif len(numpy_mat.shape)==2:
        num_of_rows = numpy_mat.shape[0]
        num_of_cols = numpy_mat.shape[1]
    else:
        raise ValueError('Input must have at most two dimensions. Got %d.' % len(numpy_mat.shape))

    # Write number of rows / columns
    write_uint32(num_of_rows)
    write_uint32(num_of_cols)

    # Pack all matrix values
    for i in range(num_of_rows):
        for j in range(num_of_cols):
            write_float32(numpy_mat[i, j].astype(np.float32))

    return packed_vec


def _numpy_2_native_matrix(numpy_mat):
    """ Converts numpy matrix to rqrmilib native matrix.

    Args:
        numpy_mat: A numpy vector/matrix. Max 2-dimensions

    Returns:
        An rqrmilib native matrix representation

    Throws:
        ValueError in case the input is invalid
        RuntimeError in case of library error
    """

    # Create native matrix object
    packed_vec = _pack_numpy_matrix(numpy_mat)
    return rqrmilib.create_matrix(packed_vec)


def _native_matrix_2_numpy(mat):
    """ Converts an rqrmilib native matrix to Numpy matrix

    Args:
        mat: An rqrmi native matrix

    Returns:
        A numpy matrix

    Throws:
        ValueError in case the input is invalid
        RuntimeError in case of library error
    """

    if 'RQRMI matrix' not in str(mat):
        raise ValueError('Input is not valid rqrmi matrix object')
    return np.array(rqrmilib.matrix_to_list(mat))


class RQRMI:

    def __init__(self, records, stage_width_list):
        """ Initializes an RQRMI data structure
        Args:
            records: Nx1 Numpy matrix. The key values / range start values.
            stage_width_list: A list of integers, the width of each stage

        Throws:
            ValueError in case the records are invalid
            RuntimeError In case of library error
        """

        # Initialize the rqrmilib in case not initialized
        global _library_initialzied
        if not _library_initialzied:
            initialize_library(False)

        self.records = _numpy_2_native_matrix(records)
        self.num_of_records = records.shape[0]
        self.probe = rqrmilib.create_probe(self.records, stage_width_list)
        self.trained_rqrmi = []
        self.stage_width_list = stage_width_list

        # Calculate the model input domain
        self.input_domain_min = np.min(records.astype(np.float64))
        self.input_domain_max = np.max(records.astype(np.float64))

        # Truncate the domain to float32
        if self.input_domain_min < np.finfo(np.float32).min:  self.input_domain_min = np.finfo(np.float32).min
        if self.input_domain_max > np.finfo(np.float32).max:  self.input_domain_max = np.finfo(np.float32).max
        self.input_domain_min=float(self.input_domain_min)
        self.input_domain_max=float(self.input_domain_max)

        # The error of the submodels in last stage
        self.error_list = [0 for _ in range(stage_width_list[-1])]

        # Reduce number of unnecessary operations in case the state of the
        # RQRMI model was not changed
        self.rqrmi_state_changed=True
        self.packed_rqrmi=None
        self.native_object=None


    def __setitem__(self, key, item):
        """ Override an existing stage with another
        Throws: KeyError in case of invalid key
        """
        if key>=len(self.trained_rqrmi):
            raise KeyError('Stage index invalid')
        self.trained_rqrmi[key]=item
        self.rqrmi_state_changed=True


    def __getitem__(self, key):
        """ Get an existing stage
        Throws: KeyError in case of invalid key
        """
        if key>=len(self.trained_rqrmi):
            raise KeyError('Stage index invalid')
        return self.trained_rqrmi[key]


    def __len__(self):
        """ Returns the number of stages in this """
        return len(self.trained_rqrmi)


    def append(self, item):
        """ Adds new stage to the RQRMI model """
        self.trained_rqrmi.append(item)
        self.rqrmi_state_changed=True


    def pack(self):
        """ Pack this to byte array """
        self._update_state()
        return self.packed_rqrmi


    def calculate_responsibility(self, stage_idx, verbose=0):
        """ Calculates the responsibility of an RQRMI stage

        Args:
            stage_idx: The stage index
            verbose: Verbosity
        """

        if stage_idx == 0:
            _log(verbose, 'R<0,0> is [%d, %d]\n' % (self.input_domain_min, self.input_domain_max))
            return

        responsibilities = rqrmilib.calculate_responsibility(self._get_native_object(), self.probe, stage_idx)

        # Print responsibilities
        if verbose > 0:
            for j, responsibility in enumerate(responsibilities):
                _log(verbose, 'R<%u,%u> is ' % (stage_idx, j))
                if len(responsibility) == 0: _log(verbose, 'empty')
                for interval in responsibility:
                    _log(verbose, '[%f, %f]' % (interval[0], interval[1]))
                _log(verbose, '\n')


    def generate_buckets(self, stage_idx, bucket_indices, num_of_samples, verbose=0):
        """ Generates a list of buckets each with dataset for training a different submodel.

        Args:
            stage_idx: The stage index
            bucket_indices: The buckets to generate
            num_of_samples: The number of samples to generate
            verbose: Verbosity

        Returns: A list of buckets for training the specified submodels.
        """

        # Handle invalid stage indices
        if stage_idx > len(self.trained_rqrmi): raise KeyError('Stage index invalid')
        elif stage_idx < 0: stage_idx = stage_idx % (len(self.trained_rqrmi)+1)

        # Generate output buckets
        N = len(bucket_indices)
        output_buckets = [None for _ in range(N)]

        for i, submodel_idx in enumerate(bucket_indices):

            # Generate samples for the current model.
            # Note: max_record parameter in create_dataset is exclusive
            # Random: false, Shuffle: true
            bucket = rqrmilib.create_dataset(self.probe, False, True, stage_idx, submodel_idx, num_of_samples)
            bucket = _native_matrix_2_numpy(bucket)

            if bucket.shape[0] == 0:
                output_buckets[i] = np.empty([0,2])
                _log(verbose, 'Bucket for submodel <%d,%d> is empty\n' % (stage_idx, submodel_idx))
                continue

            # Convert the uint32 matrix to float64 - without information loss
            # Normalize the outputs between 0 and 1
            bucket = bucket.astype(np.float64)
            bucket[:, 1] /= self.num_of_records
            output_buckets[i] = bucket

            # Debug info
            if verbose>0:
                min_bucket_idx, max_bucket_idx = np.argmin(bucket[:,0]), np.argmax(bucket[:,0])
                _log(verbose, 'Bucket for submodel <%d,%d> was generated. requested samples: %d, got: %d. min(x, y) = (%d, %f), max(x, y) = (%d, %f)\n' %
                     (stage_idx, submodel_idx, num_of_samples, bucket.shape[0],
                      bucket[min_bucket_idx, 0], bucket[min_bucket_idx, 1],
                      bucket[max_bucket_idx, 0], bucket[max_bucket_idx, 1]))

        return output_buckets


    def train_stage(self, stage_idx, structure, buckets, submodel_indices, optimizer, epochs, batch_size, verbose=0):
        """ Trains a single stage of an RQRMI model

        Args:
            stage_idx: The index of the stage to train
            structure: A list of integers, the stage's NN layer structure (including input & output)
            buckets: A list of buckets with datasets to train the submodels (float32 or float64)
            submodel_indices: Which submodels to train within stage
            optimizer: TF optimizer to use for training.
            epochs: Number of epochs to train.
            batch_size: Batch-size for training.
            verbose: Verbosity
        """

        _log(verbose, 'Training RQRMI stage %d with structure %s \n' % (stage_idx, structure))

        # Set default optimizer
        if optimizer is None:
            optimizer = tf.train.AdamOptimizer(learning_rate=1e-3)

        # In case of new stage, create a placeholder for nets
        if stage_idx == len(self):
            self.append([None for _ in range(self.stage_width_list[stage_idx])])

        net_list = self[stage_idx]

        # For each model in the stage
        for i, submodel_idx in enumerate(submodel_indices):

            # Skip submodels with no input data
            if buckets[i].shape[0]==0:
                _log(verbose, 'Skipping submodel <%d,%d>, no inputs\n' % (stage_idx, submodel_idx))
                continue

            _log(verbose, 'Training submodel <%d,%d> (dataset size: %d): ' % (stage_idx, submodel_idx, buckets[i].shape[0]))

            # Train the current net
            _, net = _nn_train_net(structure, buckets[i], optimizer, epochs, batch_size, verbose)
            net_list[submodel_idx] = net

        # Update this
        self[stage_idx] = net_list


    def calculate_transition_set(self, stage_idx, verbose=0):
        """ Calculates the transition set U_i for a given stage s_i

        Args:
            stage_idx: The required stage index
            verbose: Verbosity of this
        Throws:
            RuntimeError in case of library error
        """

        U_i = rqrmilib.calculate_transition_set(self._get_native_object(), self.probe, stage_idx)
        # Do not print transition set of last stage (long as hell)
        if stage_idx < len(self.stage_width_list)-1:
            _log(verbose, 'U_%d: %s\n' %(stage_idx, U_i))


    def evaluate(self, input_data, verbose=0):
        """ Evaluate this on the input data and returns the output

        Args:
            input_data: Numpy vector of inputs to the model
            verbose: Verbosity of this

        Returns:
            A vector of outputs (normalized to be in [0,1) )

        Throws:
            ValueError in case the input is invalid
            RuntimeError in case of library error
        """

        input_mat = _numpy_2_native_matrix(input_data)
        (total_time, avg_time), output = rqrmilib.evaluate_model(self._get_native_object(), input_mat)
        _log(verbose, 'Evaluation done. Total time taken: %.2f us. Average time per input: %.2f us\n' % (total_time, avg_time))

        # Convert to Numpy
        output = np.array(output)

        # Debug outputs
        if verbose>0:
            min_input_idx, max_input_idx = np.argmin(input_data), np.argmax(input_data)
            _log(verbose, 'Evaluation outputs: min input: %f -> (out: %f), max input: %f -> (out: %f)\n' %
                 (input_data[min_input_idx], output[min_input_idx], input_data[max_input_idx], output[max_input_idx]))

        return output


    def calculate_submodel_error(self, stage_idx, submodel_idx):
        """ Calculate the maximum error per error of this

        Args:
            stage_idx: The required stage index
            submodel_idx: The required submodel index

        Returns:
            A tuple of (maximum error, bucket coverage)

        Throws:
            RuntimeError in case of library error
        """
        return rqrmilib.calculate_submodel_error(self._get_native_object(), self.probe, stage_idx, submodel_idx)


    def pin_errors(self):
        """ Pins the errors of the submodel in the last stage

        Returns:
            A list with the submodels errors

        Throws:
            RuntimeError in case of library error
        """
        for m in range(self.stage_width_list[-1]):
            error, _ = rqrmilib.calculate_submodel_error(self._get_native_object(), self.probe, len(self)-1, m)
            if error < 0: error = 0
            self.error_list[m] = int(error)
        self.rqrmi_state_changed = True
        return self.error_list


    def _get_native_object(self):
        """ Returns a native RQRMI object representation """
        if self.rqrmi_state_changed:
            self._update_state()
            self.native_object = rqrmilib.create_model(self.packed_rqrmi)
        return self.native_object


    def _update_state(self):
        """ Generates the RQRMI model internal variables """
        if not self.rqrmi_state_changed:
            return

        self.rqrmi_state_changed = False
        self.packed_rqrmi=bytearray()

        # Returns empty byte-array in case the model is empty
        if len(self.trained_rqrmi) == 0:
            return

        # See https://docs.python.org/3/library/struct.html#struct.pack
        # for struck pack format

        # Local methods to pack numbers in little-endian format
        write_uint8 = lambda num: self.packed_rqrmi.extend(struct.pack('<B', num))
        write_uint32 = lambda num: self.packed_rqrmi.extend(struct.pack('<I', num))
        write_float32 = lambda num: self.packed_rqrmi.extend(struct.pack('<f', num))

        # Write input doamin for RQRMI model
        write_float32(self.input_domain_min)
        write_float32(self.input_domain_max)

        num_of_stages = len(self.trained_rqrmi)
        write_uint32(num_of_stages)

        # Write the stage data
        for s in range(num_of_stages):

            stage_nets=self.trained_rqrmi[s]
            num_of_models = len(stage_nets)

            # Write stage's properties
            write_uint32(num_of_models)

            for m in range(num_of_models):
                # Write the model version
                if stage_nets[m] is None:
                    # 0: empty model
                    write_uint8(0)
                    continue

                write_uint8(1)

                mu, sig, fac, omin, net_vals = stage_nets[m]
                num_of_layers=int(len(net_vals)/2)+1 # Two values per layer (biases, weights) except the input layer

                # Write the model properties
                write_float32(mu)
                write_float32(sig)
                write_float32(fac)
                write_float32(omin)
                write_uint32(num_of_layers)

                # Write layer widths
                for l in range(num_of_layers):
                    layer=net_vals[2*l]
                    if layer.ndim==1: layer_width=layer.shape[0]
                    else: layer_width=layer.shape[1]
                    write_uint32(layer_width)

                # Write weights & biases
                write_float32(net_vals[0].flatten()) # First layer bias
                write_float32(1) # First layer weight is one (always)
                for l in range(1, num_of_layers):
                    w=net_vals[2*l].flatten()
                    b=net_vals[2*l-1].flatten()
                    for x in b: write_float32(x)   # Biases
                    for x in w: write_float32(x)   # Weights

        # Write the maximum error of each last stage submodel
        for e in self.error_list:
            write_uint32(e)


    def unpack(self, buff, verbose=0):
        """ Updates the stages of this from an external byte array

        Args:
            buff: Byte array
            verbose: Verbosity
        """


        # See https://docs.python.org/3/library/struct.html#struct.pack
        # for struck pack format

        # Local methods to unpack numbers in little-endian format
        idx={'x':0}

        def read_uint8():
            idx['x']+=1
            return struct.unpack('<B', buf[idx['x']-1:idx['x']])[0]
        def read_uint32():
            idx['x']+=4
            return struct.unpack('<I', buf[idx['x']-4:idx['x']])[0]
        def read_float32():
            idx['x']+=4
            return struct.unpack('<f', buf[idx['x']-4:idx['x']])[0]

        # Return empty model in case the byte-array contains no information
        if len(buf) == 0:
            return None

        # Read global stddev and mean (not used in RQRMI version 1.1)
        _=read_float32()
        _=read_float32()

        num_of_stages=read_uint32()
        _log(verbose, 'Num of stages: %d' % num_of_stages)

        # Preallocate array
        trained_rqrmi=[None for _ in range(num_of_stages)]

        for s in range(num_of_stages):

            # Read the current stage
            num_of_models=read_uint32()

            _log(verbose, '\nStage %d num of models: %d' % (s, num_of_models))

            # Preallocate net_list
            net_list=[None for _ in range(num_of_models)]

            for m in range(num_of_models):
                # Read version
                version=read_uint8()
                if version==0:
                    _log(verbose, '\nSkipping model <%d,%d>: model not compiled' % (s, m))
                    continue
                elif version!=2:
                    _log(verbose, '\nUnsupported version for model <%d,%d>' % (s, m))
                    continue

                _log(verbose, '\nLoading model <%d, %d>: ' % (s,m))

                # Read model parameters
                mu=read_float32()
                sig=read_float32()
                fac=read_float32()
                omin=read_float32()
                num_of_layers=read_uint32()
                _log(verbose, 'layers: %d, ' % num_of_layers)

                # Preallocate net values
                net_values=[None for _ in range(2*num_of_layers-1)]

                # Read network structure
                structure=[None for _ in range(num_of_layers)]
                for l in range(num_of_layers):
                    structure[l]=read_uint32()

                # Layer 0 bias
                net_values[0]=np.empty(structure[0])

                # Preallocate all other layers
                for l in range(1, num_of_layers):
                    net_values[2*l-1]=np.empty(structure[l]) # Layer bias
                    net_values[2*l-0]=np.empty([structure[l-1], structure[l]]) # Layer weights

                _log(verbose, 'structure: [%s]' % ','.join([str(x) for x in structure]))

                # Read values of first layer
                net_values[0][0]=read_float32()
                _=read_float32() # First layer weight is one (always)

                # Read values
                for l in range(1, num_of_layers):
                    # Read bias
                    for i in range(structure[l]):
                        net_values[2*l-1][i]=read_float32()
                    # Read weights
                    for y in range(structure[l-1]):
                        for x in range(structure[l]):
                            net_values[2*l][y,x]=read_float32()

                # Update stage's net list
                net_list[m]=(mu, sig, fac, omin, net_values)

            # Update output with stage
            trained_rqrmi[s] = net_list

        # Read the maximum error of each last stage submodel
        self.error_list = []
        for e in range(len(self.trained_rqrmi[-1])):
            self.error_list.append(read_uint32())

        _log(verbose, '\n')
        self.trained_rqrmi = trained_rqrmi


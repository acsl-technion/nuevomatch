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

class ObjectPacker:
    """ Packs objects into byte array format which is recognized by the RQRMI library """

    def __init__(self):
        # Object queue
        self._object_queue = [bytearray()]

        # Append to buffer methods
        self._push_back = lambda arr: self._object_queue[-1].extend(arr)

        # See https://docs.python.org/3/library/struct.html#struct.pack
        # for struck pack format
        # Methods to pack numbers in little-endian format
        self._write_uint8 = lambda num: struct.pack('<B', num)
        self._write_uint32 = lambda num: struct.pack('<I', num)
        self._write_float32 = lambda num: struct.pack('<f', num)


    def _push_front(self, arr):
        """ Private method. Extend last object in the beginning """
        self._object_queue[-1] = bytearray(arr + self._object_queue[-1])


    def __bytes__(self):
        """ Returns a byte-array representation of this

        Throws:
            Exception in case this has opened sub-objects
        """
        if len(self._object_queue) != 1:
            raise Exception('Multiple sub-object exist. Did you forget to close a sub-object?')
        return bytes(self._object_queue[0])


    def __len__(self):
        """ Returns the length of this as byte-array

        Throws:
            Exception in case this has opened sub-objects
        """
        if len(self._object_queue) != 1:
            raise Exception('Multiple sub-object exist. Did you forget to close a sub-object?')

        return len(self._object_queue[0])


    def __enter__(self):
        """ Opens a sub-object. All following writes will be to a packed buffer """
        self._object_queue.append(bytearray())
        return self


    def __exit__(self, type, value, traceback):
        """ Close the last sub-object. Packs last object to previous one """
        last_object = self._object_queue.pop()
        self._push_back(self._write_uint32(len(last_object)))
        self._push_back(last_object)


    def append(self, buffer):
        """ Append byte-array buffer to the object

        Args:
            buffer: An object that can be converted to bytes, an integer, a float, or None

        Throws:
            ValueError in case of invalid input
        """

        if buffer is None:
            self._push_back(self._write_uint32(0))
        elif (type(buffer) is int) or (type(buffer) is np.uint32):
            self._push_back(self._write_uint32(buffer))
        elif (type(buffer) is float) or (type(buffer) is np.float32):
            self._push_back(self._write_float32(buffer))
        else:
            try:
                self._push_back(bytes(buffer))
            except:
                raise ValueError('ObjectPacker does not support type %s' % type(buffer))


    def insert(self, buffer):
        """ Insert byte-array buffer to the object (at the beginning)

        Args:
            buffer: An object that can be converted to bytes, an integer, a float, or None

        Throws:
            ValueError in case of invalid input
        """
        if buffer is None:
            self._push_front(self._write_uint32(0))
        elif (type(buffer) is int) or (type(buffer) is np.uint32):
            self._push_front(self._write_uint32(buffer))
        elif (type(buffer) is float) or (type(buffer) is np.float32):
            self._push_front(self._write_float32(buffer))
        else:
            try:
                self._push_front(bytes(buffer))
            except:
                raise ValueError('ObjectPacker does not support type %s' % type(buffer))

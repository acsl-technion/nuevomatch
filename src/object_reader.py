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

class ObjectReader:
    """ Packs objects into byte array format which is recognized by the RQRMI library """

    def __init__(self, buffer):

        # See https://docs.python.org/3/library/struct.html#struct.pack
        # for struck pack format

        self.buffer = buffer
        self.idx=0


    def read_uint8(self):
        self.idx+=1
        return struct.unpack('<B', self.buffer[self.idx-1:self.idx])[0]


    def read_uint32(self):
        self.idx+=4
        return struct.unpack('<I', self.buffer[self.idx-4:self.idx])[0]


    def read_float32(self):
        self.idx+=4
        return struct.unpack('<f', self.buffer[self.idx-4:self.idx])[0]


    def extract(self, num_of_bytes):
        """ Extract chunk of bytes from front of reader """
        self.idx+=num_of_bytes
        return self.buffer[self.idx-num_of_bytes:self.idx]


    def __len__(self):
        """ Returns the number of bytes left to read in this """
        return len(self.buffer) - self.idx


    def __bytes__(self):
        """ Returns a byte-array representation of this
        """
        return self.buffer


    def read_object(self):
        """ Reads the next packed object inside this """
        num_of_bytes = self.read_uint32()
        buffer = self.extract(num_of_bytes)
        return ObjectReader(buffer)


"""This module implements classes for a FORTRAN binary (unformatted) 
file parser"""

# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import io
import struct


class RecordError(Exception):
    """Exception for failed parsing operations"""
    pass


class ForBinFile(object):
    """A FORTRAN binary (unformatted) file parser

    Arguments:

            * **fname**: name of the binary file to open.
            * **endianness**: a flag indicating the endian convention for the
                      file, either '>' (big-endian) or '<' (small-endian).
    """

    def __init__(self, fname, endianness):

        try:
            self.rawfile = io.open(fname, 'rb')
        except IOError:
            raise ValueError("File " + fname + " not found")

        if endianness in ('>', '<'):
            self.endianness = endianness
        else:
            raise ValueError("Invalid endianness value")

    def rewind(self):
        """Move the pointer back to the beginning of the file"""

        self.rawfile.seek(0)

    def backspace(self):
        """Move the pointer back of one record"""

        prev_i = self.rawfile.tell()

        try:
            # Read the previous record's length
            self.rawfile.seek(prev_i-4)
            reclen = self.rawfile.read(4)
            try:
                reclen = struct.unpack(self.endianness + 'i', reclen)[0]
            except struct.error:
                raise RecordError(
                    "Invalid record length at beginning of record")
            # Go back by reclen spaces + 8 for opening and ending sequences
            self.rawfile.seek(prev_i-reclen-8)
        except RecordError as RErr:
            self.rawfile.seek(prev_i)
            raise RecordError(RErr)

    def read_record(self, rectype):
        """Read a record of type 'rectype'. Returns a tuple for multiple 
        values found.
        If failed rewinds to the reading point it was at the beginning.

        Arguments:

        * **rectype**: a type definition, following the conventions of the 
                       *struct* module
        """

        # We store the position in file reading to restore it in case of a
        # failed read
        prev_i = self.rawfile.tell()

        # We also check that the type makes sense
        try:
            type_size = struct.calcsize(rectype)
        except struct.error:
            raise ValueError("Invalid type passed to read_record")

        try:

            # First we need to check the record length. This comes from
            # reading the next 4 bytes in the file and interpreting
            # them as an integer
            reclen = self.rawfile.read(4)
            try:
                reclen = struct.unpack(self.endianness + 'i', reclen)[0]
            except struct.error:
                raise RecordError(
                    "Invalid record length at beginning of record")

            # Then we need to calculate the size of the bytestring to extract
            # The last i is for the record lenght at the end
            recfmt = self.endianness + rectype*int(reclen/type_size) + 'i'
            recsize = struct.calcsize(recfmt)

            # Now we extract it
            recraw = self.rawfile.read(recsize)
            rec = struct.unpack(recfmt, recraw)

            if rec[-1] != reclen:
                raise RecordError("Invalid record length at end of record")

        except RecordError as RErr:
            self.rawfile.seek(prev_i)
            raise RecordError(RErr)

        # Return the record itself

        return rec[:-1]

    def read_string_record(self):
        """Read a character record and return it joined into a single stripped 
        string."""

        recstr = self.read_record('c')
        recstr = ''.join(recstr).strip()

        return recstr

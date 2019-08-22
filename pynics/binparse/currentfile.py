# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
#from __future__ import unicode_literals

import sys
from struct import unpack, error as StructError
from pynics.binparse import ForBinFile, RecordError


class CurrentError(Exception):
    pass


class CurrentFile(object):
    """ A wrapper class to load and parse a .castep_bin file, storing its 
        contents into data structures

            Arguments:

            * **fname**: path of the current.dat file to open
    """

    def __init__(self, fname):
        """ Open and parse fname, storing its contents into data structures.
        """

        self.binfile = ForBinFile(fname, '>')

        # First things first, read chi

        try:
            self.chi = self.binfile.read_record('d')
        except RecordError as RErr:
            raise CurrentError(
                "Invalid record found in file {0}".format(fname))

        # Then go on to reading the actual current

        # Current grid will be stored here
        self.current = []
        # This is a convenient function to convert double to complex
        d_to_c = lambda l: [l[i]+1.0j*l[i+1] for i in range(0, len(l), 2)]
        # And this one groups together the data in 3x3 matrices
        c_to_3x3 = lambda l: [[[l[int(len(l)/9*(j+i*3)+k)] for i in range(3)]
                               for j in range(3)]
                              for k in range(0,
                                             int(len(l)/9))]

        while True:
            try:
                # Avoid stripping spaces since we still need to unpack this
                current_col = b''.join(self.binfile.read_record('c'))
                # First, read the first two values as integers
                nx, ny = unpack('>ii', current_col[:8])
                # Then check if they are present as indices of current
                if nx > len(self.current):
                    # Append until nx reached
                    self.current += [[] for x in range(nx-len(self.current))]
                # Same goes for ny
                if ny > len(self.current[nx-1]):
                    # Append until ny reached
                    self.current[
                        nx-1] += [None
                                  for y in range(ny-len(self.current[nx-1]))]
                # And now save
                try:
                    current_d_col = unpack(
                        '>' + 'd'*int(len(current_col)/8-1), current_col[8:])
                except StructError:
                    raise CurrentError(
                        "Invalid operation during unpacking current grid in "
                        "file {0}".format(fname))
                try:
                    self.current[nx-1][ny-1] = c_to_3x3(d_to_c(current_d_col))
                except IndexError:
                    # Something wrong in the size of the column
                    raise CurrentError(
                        "Invalid operation during parsing current grid in "
                        "file {0}".format(fname))

            except RecordError as RErr:
                # Ok, file over
                break

        # Finally, grid shape and validation
        self.grid = [0, 0, 0]

        try:
            self.grid[0] = len(self.current)
            self.grid[1] = len(self.current[0])
            self.grid[2] = len(self.current[0][0])
        except IndexError:
            raise CurrentError(
                "Invalid current grid shape found in file {0}".format(fname))

        for l in self.current:
            if len(l) != self.grid[1]:
                raise CurrentError(
                    "Invalid current grid shape found in "
                    "file {0}".format(fname))
            for l2 in l:
                if len(l2) != self.grid[2]:
                    raise CurrentError(
                        "Invalid current grid shape found in "
                        "file {0}".format(fname))

if __name__ == "__main__":

    cf = CurrentFile(sys.argv[1])

    print(cf.chi)
    print(cf.current[0][0][0])
    print(cf.grid)

# Python 2-to-3 compatibility code
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import itertools
import collections
import sys
from pynics.binparse import ForBinFile, RecordError
from pynics.binparse.castep_bin_params import (cbin_parameters_parse,
                                               CastepBinError)
from pynics.binparse.castep_bin_cell import cbin_cell_parse
from pynics.binparse.castep_bin_results import (cbin_elec_parse,
                                                cbin_optional_parse,
                                                cbin_results_parse)


class CastepBinFile(object):

    """ A wrapper class to load and parse a .castep_bin file, storing its
    contents into data structures

            Arguments:

            * **fname**: path of the .castep_bin file to open
            * **tolerant**: optional flag, if set to True any unknown optional 
                            sections will be ignored instead of raising an 
                            exception
    """

    def __init__(self, fname, tolerant=False):
        """ Open and parse fname, storing its contents into data structures.
        If the 'tolerant' flag is set to true, ignore any unknown optional 
        sections instead of raising an exception."""

        try:
            self.binfile = ForBinFile(fname, '>')
        except ValueError as VErr:
            raise ValueError(VErr)

        # Start by checking that the file is in fact a CASTEP_BIN file

        try:
            assert (self.binfile.read_string_record() == 'CASTEP_BIN')
        except AssertionError:
            raise ValueError("File " + fname +
                             " is not written in the CASTEP_BIN format")

        self.castep_version = 0.0
        self.params = self._parameters_restore()
        self.model_cell = self._cell_restore()
        self.current_cell = self._cell_restore()
        self.results = self._results_restore(tolerant)
        #self.results['elec'] = self._elec_quick_restore()
        # self._optional_restore(tolerant)

    def _parameters_restore(self):
        """ Sequentially load parameters from the currently open file """

        params = collections.OrderedDict()

        # Start with reading the header

        header = self.binfile.read_string_record()
        if (header != 'BEGIN_PARAMETERS_DUMP'):
            # Either something went wrong, or this is from too old a file. For
            # now skip
            raise CastepBinError(
                'Castep_bin file is too old to be read with this class')

        # Version number
        self.castep_version = float(self.binfile.read_string_record())

        # Parse everything else
        cbin_parameters_parse(self.binfile, params, self.castep_version)

        return params

    def _cell_restore(self):
        """ Sequentially load cell contents from the currently open file """

        cell = collections.OrderedDict()

        # Start by reading the header
        try:
            header = self.binfile.read_string_record()
            if (header == 'BEGIN_UNIT_CELL'):
                cbin_cell_parse(self.binfile, cell, self.castep_version)
            else:
                self.binfile.backspace()
                # Reader not yet implemented
                raise CastepBinError(
                    'Castep_bin file is too old to be read with this class')
        except RecordError:
            raise CastepBinError('End of file reached while parsing cell')

        return cell

    def _results_restore(self, tolerant=False):
        """Restore all the results

                Arguments:

                * **tolerant**: optional flag, if set to True any unknown optional sections will be ignored instead of raising an exception

        """

        results = collections.OrderedDict()
        try:
            cbin_results_parse(self.binfile, results, self.castep_version,
                               self.params, self.current_cell, tolerant)
        except RecordError:
            raise CastepBinError('End of file reached while parsing results')

        return results


if __name__ == "__main__":

    testfile = CastepBinFile(sys.argv[1], tolerant=True)
    elec = testfile.results['elec']

    # print elec['total_energy']*27.2107

    print(testfile.params)
    print(testfile.current_cell)
    print(testfile.results)

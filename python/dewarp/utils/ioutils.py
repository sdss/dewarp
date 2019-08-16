# encoding: utf-8
#
# @Author:    Adam Mendenhall
# @Date:      August 5, 2019
# @Filename:  ioutils.py
# @License:   BSD 3-Clause
# @Copyright: Adam Mendenhall
#

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

"""
returns an array whose length is a multiple of the length of float_data_indicies
returned array contains float attributes for first instrument, then second, then third...
only contains information for instruments which satisfy which_instruments(name) is True, where name is the name of the instrument as it appears in the file (e.g. 'BOSS' 'Fiducial')
assumes each instrument is on its own line
ignores (only) lines that begin with '#' or emptylines
assumes attribute fields for each instrument are separated (only) by white space
assumes that the fifth field on a line is the name of the instrument
"""
def specific_instrument_entries_from_file(filename_with_extension, which_instruments, float_data_indicies):
    f = open(filename_with_extension, 'r')
    data = []
    for line in f:
        if len(line)>0 and line[0]!='#':
            attributes = line.split()
            if len(attributes)>4:
                if which_instruments(attributes[4]):
                    try:
                        datum = []
                        for i in range(min(len(attributes)-1,len(float_data_indicies))):
                            datum.append(float(attributes[float_data_indicies[i]]))
                        for a in datum:
                            data.append(a)
                    except ValueError:
                        pass
    return data

def fiducial_xys_from_file(filename_with_extension):
    return specific_instrument_entries_from_file(filename_with_extension, lambda a: a.lower()=='fiducial', [2,3])

def fiducial_hexs_from_file(filename_with_extension):
    return specific_instrument_entries_from_file(filename_with_extension, lambda a: a.lower()=='fiducial', [0,1])

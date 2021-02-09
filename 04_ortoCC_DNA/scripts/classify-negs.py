"""
classify-negs.py - parse CIGAR string and classify the errors
Nathan Lubock

REQUIRES Sam 1.4 CIGAR strings (e.g. = for a match and X for a mismatch)
"""

# Ensure Python 2/3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import itertools
import multiprocessing
import sys
import re
from signal import signal, SIGPIPE, SIG_DFL
from collections import defaultdict

# catch broken pipe errors to allow ex) python pyParse.py foo bar | head
# see: https://stackoverflow.com/a/30091579
signal(SIGPIPE, SIG_DFL)

#===============================================================================

def split_tag(my_tag):
    """
    Split the CIGAR string into various operations with some regex. Assumes the
    tag is alphanumeric!
    Idea:
        regex (\d+) = repeats of digits, (\D) = [^0-9] (not numbers)
        (\D+) allows for matching of deletion characters
        ex) 30M1I4M -> [('30','M'), ('1','I'), ('4','M')]
    Note:
        WILL LEAVE OFF STRAGGLING NUMERICS
        ex) 30M1I4M4 -> [('30','M'), ('1','I'), ('4','M')]
    """
    my_split = re.findall(r'(\d+)(\D+)', my_tag)
    #print(my_split)
    return ((int(x[0]), x[1]) for x in my_split)

#-------------------------------------------------------------------------------

def classify(cigar):
    """
    Classify a read based on its CIGAR string
    """
    split = list(split_tag(cigar))
    indels = [x for x in split if x[1] in ['D', 'I']]

    if len(indels) > 0:
        # only track position of the first indel
        pos = sum(x[0] for x in itertools.takewhile(lambda x: x[1] not in ['D', 'I'], split))
        length, err = indels[0]
        if length == 3 and err == 'D' :
            out_class = 'Skip'
        elif length == 3 and err == 'I' :
            out_class = 'Add'
        else:
            out_class = 'Indel'
    else:
        out_class = 'Mismatch'
        pos = 0
    #print(out_class, pos)

    return (out_class, pos)

#-------------------------------------------------------------------------------

def wrap(pack):
    """
    Quick wrapper for multiprocessing
    """
    cigar = pack[-1]
    pos = pack[-2]

    cigar_class, pos_offset = classify(cigar)
    return pack[:-2] + [cigar] + [str(int(pos) + pos_offset), cigar_class]

#===============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Classify CIGAR strings. MUST BE SAM 1.4 VERSION!')
    parser.add_argument('infile',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        nargs='?',
                        help='path to flat file with CIGARS (or stdin if none).' +
                        ' Assumes CIGAR is last column, and left-most position' +
                        ' is the second last column.')
    parser.add_argument('-j',
                        '--proc',
                        dest='proc',
                        type=int,
                        default=1,
                        metavar='N',
                        choices=range(1, multiprocessing.cpu_count()+1),
                        help='number of processors (default=1, max={})'.format(
                            multiprocessing.cpu_count()))
    args = parser.parse_args()

    # lazily collapse whitespace and parse last item
    pool = multiprocessing.Pool(args.proc)
    collapse = (line.split() for line in args.infile)
    out_class = pool.imap_unordered(wrap, collapse, chunksize=10000)
    for line in out_class:
        print('\t'.join(line), file=sys.stdout)

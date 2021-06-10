"""
Simple barcode mapping utility.
Nathan Lubock

Outputs
"""

# Ensure Python 2/3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import itertools
import argparse
import multiprocessing
import sys
import csv
from signal import signal, SIGPIPE, SIG_DFL
from collections import Counter, defaultdict
import numpy as np

# External Deps
import Levenshtein

# catch broken pipe errors to allow ex) python pyParse.py foo bar | head
# see: https://stackoverflow.com/a/30091579
signal(SIGPIPE, SIG_DFL)

#===============================================================================
# recepies/helper functions


def grouper(iterable, size, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * size
    return itertools.izip_longest(fillvalue=fillvalue, *args)

def flatten(iterable):
    "Flatten one level of nesting"
    return itertools.chain.from_iterable(iterable)

def most_common(hashable):
    "Simple alias for grabbing the most common element"
    return Counter(hashable).most_common(1)[0][0]

#===============================================================================
# I/O Functions

def text_reader(text):
    """
    Read in a text file which is just sequences
    """
    seq = itertools.islice(text, None)
    for s in seq:
        yield s.strip()



def fastq_reader(fastq):
    """
    Read in a fastq file laizy and return a generator of the sequence

    Parameters:
    -----------
    fastq :: FileType
        opened file to read

    Yields:
    -------
    generator :: str
        read portion of the fastq file

    Requires:
    ---------
    itertools
    """
    fourth = itertools.islice(fastq, 1, None, 4)
    for seq in fourth:
        yield seq.strip()

def fasta_reader(fasta):
    """
    Read in a fasta file lazily and return a generator of the name and sequence

    Parameters:
    -----------
    fasta :: FileType
        opened file

    Yields:
    -------
    generator :: (name, seq)
        name :: str
            Name of the read taken from the fasta file
        read :: str
            Sequence taken from the fasta file

    Requires:
    ---------
    itertools

    Example:
    --------
    itertools.groupby takes a key function and groups all items into a list
    until that key changes. We can key on lines beginning with >, then grab
    every line until the next record in the fasta. This makes our method robust
    to some fasta formats that have forced line breaks at given characters.

    foo = '>ABC>DEF>GHI'
    [(k, list(g)) for k,g in itertools.groupby(foo, lambda x: x == '>')]
    --> [(True, ['>']), (False, ['A', 'B', 'C']), (True, ['>']), ... ]

    Note:
    -----
    Adapted from: https://www.biostars.org/p/710/#1412
    """
    # ditch the boolean (x[0]) and just keep the header/seq grouping
    fa_iter = (x[1] for x in itertools.groupby(fasta, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        name = next(header)[1:].strip()
        # join all sequence lines to one by iterating until the next group.
        read = "".join(s.strip() for s in next(fa_iter))
        yield name, read

#===============================================================================

def get_idx(seq, bc_start, bc_len):
    """
    Get the slices corresponding to barcode and variant

    Parameters:
    -----------
    seq :: str
        Sequence to slice
    bc_start :: int
        Starting location of barcode (can be negatively indexed)
    bc_len :: int

    Returns:
    --------
    bc_slice :: slice
        location of the barcode
    var_slice :: slice
        location of the variant
    """
    if bc_start <= 0:
        bc_slice = slice(len(seq) + bc_start, len(seq) + bc_start + bc_len)
        var_slice = slice(bc_start)
    else:
        # compensate for 0-based indexing
        bc_slice = slice(bc_start - 1, bc_start + bc_len - 1)
        var_slice = slice(bc_start + bc_len - 1, len(seq))
    return (bc_slice, var_slice)

#-------------------------------------------------------------------------------
def merge_reads(reads):
    """
    Merge list of strings into a consensus sequnce. Pad Ns for longer sequences.

    Parameters:
    -----------
    reads :: [str]
        List of strings corresponding to variants

    Returns:
    --------
    out_seq :: str
        Merged read

    Requires:
    ---------
    itertools
    most_common -> collections.Counter().most_common()[0][0]

    Note:
    -----
    Assumes all reads are aligned at the left-most base (e.g. no indels at the
    very beginning of the sequence. Briefly, we transpose the sequences, grab
    the most common base at each position, and then trim the resulting
    consensus to the most common length (for >2 reads). For two reads, we will
    call an N at any disagreeing bases. Also note that in the case of ties for
    >3 reads, we will report at random (according to Counter) the most common
    base.
    """
    if len(reads) == 1:
        return reads[0]

    trans = itertools.izip_longest(*reads, fillvalue='N')
    if len(reads) == 2:
        consensus = (x[0] if x[0] == x[1] else 'N' for x in trans)
        return ''.join(consensus)
    else:
        raw_seq = ''.join(most_common(x) for x in trans)
        common_len = most_common((len(x) for x in reads))
        return raw_seq[:common_len]

#-------------------------------------------------------------------------------

def mismatch(word, num_mismatches):
    """
    Generate all mismatches of a sequence upto a specified distance

    Parameters:
    -----------
    word :: str
    num_mismatches :: int

    Yields:
    -------
    _ :: str
        generator of sequences at num_mismatches from word

    Requires:
    ---------
    itertools

    Note:
    -----
    Boosted from: http://stackoverflow.com/a/11679968
    Briefly, generate all positions that will be edited. Insert a list of
    substitutions at those positions, then generate all possible combinations
    of that list.
    """
    letters = 'ATGC'
    for dist in range(1, num_mismatches + 1):
        for locs in itertools.combinations(range(len(word)), dist):
            this_word = [[char] for char in word]
            for loc in locs:
                orig_char = word[loc]
                this_word[loc] = [l for l in letters if l != orig_char]
            for poss in itertools.product(*this_word):
                yield ''.join(poss)

#-------------------------------------------------------------------------------

def bc_overlap(bcmap, dist):
    """
    Check if there are any sequences within a given distance to the input

    Parameters:
    -----------
    bcmap :: dict
        dict of sequences to check
    dist :: int
        distance to generate variants

    Returns:
    --------
    collector :: set
        Set of sequences from bcmap that are within dist of eachother

    Requires:
    ---------
    mismatch

    Note:
    -----
    Assumes that checking len(mismatch(seq, dist)) x len(bcmap) is less
    opperations than checking len(bcmap) x len(bcmap). Also uses a dict for
    constant time lookups.
    Briefly, iterate through every seq in bcmap, generate its variants, check
    to see if any of those variants are in bcmap, and then add them to a
    collector.
    """
    collector = set()
    for bc in bcmap:
        bc_var = mismatch(bc, dist)
        for var in bc_var:
            if var in bcmap and var not in collector:
                # assumes implicitly that bc -> var and var -> bc will be added
                collector.add(var)
    return collector

#-------------------------------------------------------------------------------

def bootstrap_dist(lib, num=50000, percentile=1):
    """
    Bootstrap the distribution of Levenshtein distances between library members

    Parameters:
    -----------
    lib :: [str]
        list of library members
    num :: int
        number of pairwise comparisons
    percentile :: float
        percentile of distribution to set cutoff at

    Returns:
    --------
    dist :: int
        Levenshtein distance cutoff

    Requires:
    ---------
    itertools
    Levenshtein - https://pypi.python.org/pypi/python-Levenshtein/0.12.0
    numpy as np
    """
    # we want to do n pairwise comparisons, so multiply n by 2
    # if library is smaller than n*2 then just do the actual pairwise
    size = len(lib)
    if size**2 - size < num*2:
        samps = itertools.permutations(lib, 2)
    else:
        samps = grouper(np.random.choice(lib, num*2), 2)
    dist = np.percentile([Levenshtein.distance(*x) for x in samps], percentile)

    # set a distance floor at 5 for really low cutoff
    if dist < 5:
        return 5
    else:
        return dist

#-------------------------------------------------------------------------------

def intra_bc_dist(pack, dist, debug=False):
    """
    Return the percentage of sequences further than a given distance

    Parameters:
    -----------
    pack :: (str, [str])
        bc :: str
            input barcode
        variants :: [str]
            variants associated with that barcode
    dist :: int
        Levenshtein distance cutoff
    debug :: bool
        switch to output the list of distances

    Returns:
    --------
    bc :: str
        input barcode
    percent :: float
        percent of sequences > dist

    Requires:
    ---------
    collections.Counter
    Levenshtein - https://pypi.python.org/pypi/python-Levenshtein/0.12.0
    """
    # ~ 5000 bc/sec
    bc, variants = pack
    uniqd = Counter(variants).most_common()
    dists = [0]
    if len(uniqd) > 1:
        compare = ((uniqd[0][0], x[0]) for x in uniqd[1:])
        dists.extend([Levenshtein.distance(*x) for x in compare])
    else:
        return (bc, 0.0)
    if debug:
        return dists
    else:
        total_reads = sum(x[1] for x in uniqd)
        out_reads = sum(tup[1] for lev, tup in zip(dists, uniqd) if lev > dist)
        percent = float(out_reads) / total_reads * 100.0
        return (bc, percent)

#===============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Collapse reads into a barcode-variant map. Barcodes are \
                   filtered at three steps. First, we remove any barcodes that \
                   are a Hamming distance of -dN away from eachother. Then we \
                   remove barcodes where -pN percent of reads for that barcode \
                   are > -cN away from the most common read.')
    # parser.add_argument('library',
    #         type=argparse.FileType('r'),
    #         help='path to a *.fasta file of your pefect library')
    parser.add_argument('infile',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        nargs='?',
                        help='path to a *.fastq file of the reads (or stdin if none)')
    parser.add_argument('-v',
                        '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='Verbose logging')
    parser.add_argument('-j',
                        '--proc',
                        dest='proc',
                        type=int,
                        default=1,
                        metavar='N',
                        choices=range(1, multiprocessing.cpu_count()+1),
                        help='number of processors (default=1, max={})'.format(multiprocessing.cpu_count()))
    parser.add_argument('-s',
                        '--bc-start',
                        dest='bc_start',
                        type=int,
                        default=1,
                        metavar='N',
                        help='starting position of barcode (1-indexed). Can also be negative \
                                for relative indexing from the end of a read (default=1)')
    parser.add_argument('-l',
                        '--bc-length',
                        dest='bc_length',
                        type=int,
                        default=20,
                        metavar='N',
                        help='length of barcode (default=20)')
    parser.add_argument('-m',
                        '--min-reads',
                        dest='min_reads',
                        type=int,
                        default=5,
                        metavar='N',
                        help='min reads for each barcode (default=5)')
    parser.add_argument('-d',
                        '--bc-overlap',
                        dest='overlap',
                        type=int,
                        metavar='N',
                        default=1,
                        choices=range(5),
                        help='filter barcodes within a certain distance of eachother \
                            (Default=1). NOTE - computationally intensive ~30k bc/min for \
                            dist of 2.')
    parser.add_argument('-c',
                        '--bc-contam',
                        dest='dist',
                        type=int,
                        metavar='N',
                        default=10,
                        help='set distance cutoff for intra-barcode filtering \
                             (default=10). Used in conjunction with -p.')
    parser.add_argument('-p',
                        '--percent',
                        dest='percent',
                        type=int,
                        metavar='N',
                        default=5,
                        choices=range(101),
                        help='set percentage cutoff for intra-barcode filtering \
                             (default=5). Used in conjunction with -c.')
    parser.add_argument('-b',
                        '--bad-bcs',
                        dest='bad_bcs',
                        type=str,
                        metavar='foo',
                        help='output bad barcodes to this file')
    args = parser.parse_args()

    #---------------------------------------------------------------------------
    # DATA LOADING
    # about 1 min for about 13 million reads
    if args.verbose:
        print('Mapping BCs...', file=sys.stderr)

    bcmap = defaultdict(list)
    for seq in text_reader(args.infile):
        bc_slice, var_slice = get_idx(seq, args.bc_start, args.bc_length)
        bcmap[seq[bc_slice]].append(seq[var_slice])

    # filter the dict in place and add bad bcs to a list
    bad_bcs = []
    count=0
    for k,v in bcmap.iteritems():
        if len(v) < args.min_reads:
            count += 1
            for seq in Counter(v).iteritems():
                bad_bcs.append((k, seq[1], 'min_reads', seq[0]))
    bcmap = defaultdict(list,
                        {k:v for k, v in bcmap.iteritems() if len(v) >= args.min_reads})
    if args.verbose:
        print('Found {} BCs with < {} reads'.format(count, args.min_reads),
              file=sys.stderr)

    # # bootstrap cutoff based on input library (can be done in parallel with ^^^)
    # # takes about 30 secs for a decent lenght library
    # lib = list(fastaReader(args.library))
    # cutoff = bootstrap_dist([x[1] for x in lib])

    #---------------------------------------------------------------------------
    # BARCODE FILTERING

    if args.verbose:
        print('Mapped {} BCs'.format(len(bcmap.keys())), file=sys.stderr)
        print('Finding overlapping BCs at dist={}'.format(args.overlap),
              file=sys.stderr)

    # vvv expensive! Consider distributing if args.overlap > 1
    # ~ 2 mins for 900k bcs (0.06% of bcs) @ dist=1
    overlap_bcs = bc_overlap(bcmap, args.overlap)
    for bc in overlap_bcs:
        popped = Counter(bcmap.pop(bc)).iteritems()
        bad_bcs.extend((bc, x[1], 'overlap', x[0]) for x in popped)

    if args.verbose:
        print('Found {} overlapping BCs'.format(len(overlap_bcs)), file=sys.stderr)
        print('Finding contaminated BCs', file=sys.stderr)

    # vvv also expensive... should consider distributing
    # ~3 mins for 900k bcs (removes ~2% at 15)
    bc_dist = (intra_bc_dist(x, args.dist) for x in bcmap.iteritems())
    contam_bc = [x[0] for x in bc_dist if x[1] > args.percent]
    for bc in contam_bc:
        popped = Counter(bcmap.pop(bc)).iteritems()
        bad_bcs.extend((bc, x[1], 'contam', x[0]) for x in popped)

    if args.verbose:
        print('Found {} contaminated BCs'.format(len(contam_bc)), file=sys.stderr)
        print('Total BCs remaining after filter {}'.format(len(bcmap.keys())),
              file=sys.stderr)

    #---------------------------------------------------------------------------
    # READ MERGING

    def dummy(pack):
        """since we cant pickle lambdas define a dummy function to map over the tups"""
        key, value = pack
        return (key, merge_reads(value), len(value))

    # set up pool and distribute
    pool = multiprocessing.Pool(args.proc)
    writer = csv.writer(sys.stdout, lineterminator='\n')

    # ideally results are written as they come in from pool, unsure if true...
    consensus = pool.imap_unordered(dummy, bcmap.iteritems(), chunksize=10000)
    writer.writerows(consensus)

    pool.close()
    pool.join()

    if args.bad_bcs is not None:
        with open(args.bad_bcs, 'w') as f:
            writer = csv.writer(f, lineterminator='\n')
            writer.writerow(['bc', 'count', 'reason', 'var'])
            writer.writerows(bad_bcs)

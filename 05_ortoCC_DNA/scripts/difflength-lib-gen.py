"""
Simple library generation script.
Nathan Lubock

Takes in:
    -- Fasta of Library Members
    -- List of combinations (otherwise all pairwise)
    -- List of restriction enzyme sites to check (optional)
    -- Number of codon usages (otherwise random per variant)
Outputs:
    -- Library of variants (stdout)
"""


# Ensure Python 2/3 compatibility
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
# from __future__ import unicode_literals

import itertools
import argparse
import re
import sys
from signal import signal, SIGPIPE, SIG_DFL
from numpy.random import choice, RandomState
from numpy import mean

#===============================================================================
# DEFAULTS/GLOBALS

# Reproducible builds
# RandomState(42)

# catch broken pipe errors to allow ex) python pyParse.py foo bar | head
# see: https://stackoverflow.com/a/30091579
signal(SIGPIPE, SIG_DFL)

# default codon usage dictionary (from sri, but see:
# http://www.openwetware.org/wiki/Escherichia_coli/Codon_usage)
CODON_USAGE = {
    '*': [('TAA', 0.64), ('TGA', 0.29), ('TAG', 0.07)],
    'A': [('GCT', 0.16), ('GCC', 0.27), ('GCA', 0.21), ('GCG', 0.36)],
    'C': [('TGT', 0.44), ('TGC', 0.56)],
    'D': [('GAT', 0.63), ('GAC', 0.37)],
    'E': [('GAA', 0.69), ('GAG', 0.31)],
    'F': [('TTT', 0.57), ('TTC', 0.43)],
    'G': [('GGT', 0.34), ('GGC', 0.41), ('GGA', 0.10), ('GGG', 0.15)],
    'H': [('CAT', 0.57), ('CAC', 0.43)],
    'I': [('ATT', 0.51), ('ATC', 0.42), ('ATA', 0.07)],
    'K': [('AAA', 0.76), ('AAG', 0.24)],
    'L': [('TTA', 0.13), ('TTG', 0.13), ('CTT', 0.10), ('CTC', 0.10), ('CTA', 0.04), ('CTG', 0.50)],
    'M': [('ATG', 1.00)],
    'N': [('AAT', 0.45), ('AAC', 0.55)],
    'P': [('CCT', 0.16), ('CCC', 0.12), ('CCA', 0.19), ('CCG', 0.53)],
    'Q': [('CAA', 0.35), ('CAG', 0.65)],
    'R': [('CGT', 0.38), ('CGC', 0.40), ('CGA', 0.06), ('CGG', 0.10), ('AGA', 0.04), ('AGG', 0.02)],
    'S': [('TCT', 0.15), ('TCC', 0.15), ('TCA', 0.12), ('TCG', 0.15), ('AGT', 0.15), ('AGC', 0.28)],
    'T': [('ACT', 0.16), ('ACC', 0.44), ('ACA', 0.13), ('ACG', 0.27)],
    'V': [('GTT', 0.26), ('GTC', 0.22), ('GTA', 0.15), ('GTG', 0.37)],
    'W': [('TGG', 1.00)],
    'Y': [('TAT', 0.57), ('TAC', 0.43)]
}

spacerx = "TCCAGAAGAGCACGGCGGGGTCTTCCGGAGTTGCCGACGACGAAAGCGACTTTGGGTTCTGTCTGTTGCACGGTTTAGCAGAAGGTTTGAGGAATAGGTTAAATTGAGTGGTTTAATAACGGTATGTCTGGGATTAAAGTGTAGTATAGTGTGATTATCGGAGACGGTTTTAAGACACGAGTTCCCAAAATCAAGCGGGGTCATTACAACGGTTGGATGGGCAGACGGG"
spacery = "CCTTCATTCACTCGGCAGCTCTGTAATAGGGACTAAAAAAGTGATGATAATCATGAGTGCCGCGTTATCAGTCGTATGCCTTCTCGAGTTCCGTCCAGTTAAGCGTGACAGTCCCAGTGTACCCACAAACCGTGATGGCTGTGCTTGGAGTCAATCGCAAGTAGGGTGGTGTCGGAACAGAGCGGTCTTACGGCCCCTCATCGCCTAGAATTACCTACTACGGTCGACCATACCTTCGATTATCGCGGCCACTCTCGCATTAGTCGGCAGAGGTGGTTGTGTTGCTTGAAGACCCGGCGGCACAGCAGTGGA"

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

def rev_comp(seq):
    pairs = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    #return ''.join(reversed([pairs[x] for x in seq]))
    return ''.join([pairs[x] for x in seq][::-1]) # ~2x speedup

def restricted_float(x):
    """
    See: http://stackoverflow.com/a/12117065
    """
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("{!r} not in range [0.0, 1.0]".format(x))
    return x

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
# user functions

def norm_codons(codon_dict, cutoff = 0.15):
    """
    Drop rare codons (according to cutoff) and renormalize percentages

    Parameters:
    -----------
    codon_dict :: dict([(str, float), ...])
        Key - one letter code, Val - list of codons and their frequency
    cutoff :: float
        Drop codons beneath this cutoff

    Returns:
    --------
    out_dict :: dict([(str, float), ...])
    """
    out_dict = {}
    for key in codon_dict:
        # filter codons, then norm
        common_codons = [x for x in codon_dict[key] if x[1] > cutoff]
        freq_sum = sum(x[1] for x in common_codons)
        out_dict[key] = [(k, v / freq_sum) for k, v in common_codons]
    return out_dict

#-------------------------------------------------------------------------------

def rand_codon(amino, codon_dict):
    """
    Randomly pull a codon for a given amino acid weighted on it's frequency

    Requires:
    ---------
    numpy.random.choice
    """
    codons, freqs = zip(*codon_dict[amino])
    out_codon = choice(codons, p=freqs)
    return out_codon

#-------------------------------------------------------------------------------

def test_re(seq, re_list):
    """
    Check a sequence for restriction enzyme sites
    """
    checks = (x in seq for x in re_list)
    return sum(checks)

def generate_seq(peptide, codon_dict, primer, spacer, re_list, flip=False, revcomp=False):
    # iterate forever until you find a good match
    loops = 0
    while True:
        if revcomp:
            seq = rev_comp(''.join(rand_codon(x, codon_dict) for x in peptide))
        else:
            seq = ''.join(rand_codon(x, codon_dict) for x in peptide)

        # test with spacer but do not keep for output (since its common to x&y)
        if flip:
            test_seq = ''.join((spacer, seq, primer))
            out_seq = ''.join((seq, primer))
        else:
            test_seq = ''.join((primer, seq, spacer))
            out_seq = ''.join((primer, seq))

        loops += 1
        if test_re(test_seq, re_list) == 0:
            break
        if loops > 50:
            break

    return out_seq, loops


#===============================================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate a codon-optimized DNA library from peptides')
    parser.add_argument('peptides',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        nargs='?',
                        help='path to a *.fastA file of the peptides (or stdin if none)')
    parser.add_argument('-p',
                        '--pairs',
                        dest='pairs',
                        type=argparse.FileType('r'),
                        default=None,
                        metavar='foo.csv',
                        help='CSV of explicit pairwise peptides to make')
    parser.add_argument('-s',
                        '--re-sites',
                        dest='re_sites',
                        type=argparse.FileType('r'),
                        default=None,
                        metavar='foo.txt',
                        help='DNA sequence of restriction enzyme sites to check (5 -> 3) MUST INCLUDE REVERSE COMPLIMENT FOR NON-PALINDROMIC SEQUENCES!')
    parser.add_argument('-v',
                        '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='Verbose logging')
    parser.add_argument('-c',
                        '--cutoff',
                        dest='cutoff',
                        type=restricted_float,
                        default=0,
                        metavar='N',
                        help='drop codons that occur with a frequency < N')
    # parser.add_argument('-j',
    #                     '--proc',
    #                     dest='proc',
    #                     type=int,
    #                     default=1,
    #                     metavar='N',
    #                     choices=range(1, multiprocessing.cpu_count()+1),
    #                     help='number of processors (default=1, max={})'.format(multiprocessing.cpu_count()))
    parser.add_argument('-f',
                        '--forward',
                        dest='f_primer',
                        type=str,
                        default='',
                        help='forward primer sequence (5 -> 3)')
    parser.add_argument('--outerxfwdlist',
                        dest='fx_list',
                        type=str,
                        nargs='*',
                        help='outer forward x primer list')
    parser.add_argument('-r',
                        '--reverse',
                        dest='r_primer',
                        type=str,
                        default='',
                        help='reverse primer sequence (5 -> 3)')
    parser.add_argument('--outerxrevlist',
                        dest='rx_list',
                        type=str,
                        nargs='*',
                        help='out reverse x primerlist')
    parser.add_argument('--spacer',
                        dest='spacer',
                        type=str,
                        default='',
                        help='spacer sequence (5 -> 3)')
    parser.add_argument('-x',
                        '--revcomp-x',
                        dest='rc_x',
                        action='store_true',
                        help='Reverse compliment X')
    parser.add_argument('-y',
                        '--revcomp-y',
                        dest='rc_y',
                        action='store_true',
                        help='Reverse compliment y')
    # parser.add_argument('-o',
    #                     '--other_x',
    #                     dest='rev_x')
    # parser.add_argument('-q',
    #                     '--other_y',
    #                     dest='fwd_y')
    parser.add_argument('--yBsaI',
                        dest='yBsaI',
                        type=str,
                        help='should contain the major reverse primer and the BsaI site')
    parser.add_argument('--outeryfwdlist',
                        dest='fy_list',
                        type=str,
                        nargs='*',
                        help='outer forward y primer list')
    parser.add_argument('--outeryrevlist',
                        dest='ry_list',
                        type=str,
                        nargs='*',
                        help='out reverse y primerlist')
    args = parser.parse_args()

#===============================================================================

# load everything
peptides = {key:val.upper() for key, val in fasta_reader(args.peptides)}
spacer = args.spacer.upper()

if args.pairs is not None:
    pairs = [x.strip().split(',') for x in args.pairs.readlines()]
    if args.verbose:
        print('Loaded {} explicit pairs'.format(len(pairs)), file=sys.stderr)
        #print("pairs ", pairs,file=sys.stderr)
else:
    pairs = itertools.product(*itertools.repeat(peptides.keys(), 2))
    if args.verbose:
        print('Generating all pairwise members', file=sys.stderr)

if args.re_sites is not None:
    re_sites = [x.strip().upper() for x in args.re_sites.readlines()]
    if args.verbose:
        print('Will check for: {}'.format(re_sites), file=sys.stderr)
else:
    re_sites = []


#-------------------------------------------------------------------------------

if args.verbose:
        print('Dropping codons that are less frequent than {}'.format(args.cutoff), file=sys.stderr)
codon_dict = norm_codons(CODON_USAGE, args.cutoff)

# test up to n-mer into the spacer region
mer = 7

# main loop
if args.verbose:
    print('Generating variants...', file=sys.stderr)
x_acc = []
y_acc = []
variants = 0
print_mult = [2 << i for i in range(24)]
for codonnum in range(1,4):
    for name in peptides:
        #print("Pairs", x_name,y_name, file=sys.stderr)
        x_dna, x_loops = generate_seq(peptides[name], codon_dict, args.f_primer.upper(), spacer[:mer], re_sites, flip=False, revcomp=args.rc_x)
        x_acc.append(x_loops)
        y_dna, y_loops = generate_seq(peptides[name], codon_dict, args.r_primer.upper(), spacer[-mer:], re_sites, flip=True, revcomp=args.rc_y)
        y_acc.append(y_loops)
        for i in range(len(args.ry_list)):
            ofwdx = args.fx_list[i]
            #print("this is ofwdx ", ofwdx," and this is i ", i, " and this is args.fx_list ", args.fx_list)
            orevx = args.rx_list[i]
            ofwdy = args.fy_list[i]
            orevy = args.ry_list[i]
            x_dna_pool = str(ofwdx) + x_dna

            # combine together and output fasta
            variants += 1
            spacer=spacerx[:215-len(x_dna_pool)]
            print('>{}_x,{}'.format(''.join((name,'-c',str(codonnum),'-',str(i +1))), ''.join((x_dna_pool, spacer,orevx))), file=sys.stdout)
            spacer = spacery[len(spacery)-(180-len(y_dna)):]

            print('>{}_y,{}'.format(''.join((name,'-c',str(codonnum),'-',str(i+1))), ''.join((ofwdy, spacer, y_dna,args.yBsaI, orevy))), file=sys.stdout)
            if args.verbose and variants in print_mult:
                print('Number of Variants: {}'.format(variants), file=sys.stderr)


            # seqlen = len(x_dna + y_dna)
            # if seqlen < 263:
            #     space =  SPACR[1:(263 - seqlen)]
            #     spacer = "tccAGAAGAGC" + space + "GCAGTGga"
            #     print('>{}+{},{}'.format(x_name, y_name, ''.join((x_dna, spacer, y_dna))), file=sys.stdout)
            #     if args.verbose and variants in print_mult:
            #         print('Number of Variants: {}'.format(variants), file=sys.stderr)


if args.verbose:
    print('Mean loops for X - {}'.format(mean(x_acc)), file=sys.stderr)
    print('Mean loops for Y - {}'.format(mean(y_acc)), file=sys.stderr)


'''
@author: Kimberly Insigne
kiminsigne@gmail.com

This script demultiplexes a FASTQ file of single reads into separate FASTQ files 
based on their index. If an index read is 1bp away from an index in the given 
reference file, it will be included in that sample FASTQ file.

Input:
reads_file	: FASTQ file of single-end reads
index_reads_file : FASTQ file of index read, where index is first N bp of read
index_file	: tab-separated text file, no headers, first column is sample name, 
second column is index
index_length : integer, length of index 

Optional:
output	: name of output directory, default is current directory
rev		: if enabled, sequences in index file are reverse complement 
relative to index reads file

Output:
Separate FASTQ files, one for each sample
'''

import argparse
import gzip
from itertools import islice
from helpful_utils import reverse_complement, mutate_1_bp


if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('reads_file', help='.fastq file')
	parser.add_argument('index_reads_file', 
		help='.fastq file, where the index is the first nbp of the read')
	parser.add_argument('index_file', 
		help='tab separated text file where first column is sample name \
												and second column is index')
	parser.add_argument('index_length', type=int, 
		help="length of index, assuming it's the first n bases")
	parser.add_argument('--output', nargs='?', const='./', default='./', 
		help='Name of output directory, default is current directory')
	parser.add_argument('-rev', dest="rev", action='store_true', 
		help='If sequences in index_file are reverse complement\
		relative to index reads file')
	parser.add_argument('-c','--compressed', dest="compressed", action='store_true',
		help="Flag for if the reads_files are in .gz form")
	args = parser.parse_args()	

	output = args.output

	# read in index file

	index_file = open(args.index_file, 'rU')
	indices = {}
	for line in index_file:
		name, index = line.strip().split('\t')
		# convert index to uppercase
		index = index.upper()
		if args.rev:
			index = reverse_complement(index)
		indices[index] = name
		# generate all 1bp mutants and add to look-up so we don't have to 
		# check each index
		for x in mutate_1_bp(index):
			indices[x] = name

	# because we added 1bp mutant indices there will be lots of repeated 
	# sample entries in the dictionary, so make this a set
	samples = set(indices.values())
	# for each sample, create a dictionary where keys are sample names and 
	# values are opened files for writing
	index_to_handle = {sample : open(output+sample+'.fastq', 'w') for sample in samples}
	index_to_handle_c = {sample : gzip.open(output+sample+'.fastq.gz', 'w') for sample in samples}
	bad_index = open('bad_index.fastq', 'w')

	if args.compressed:
		reads_file = gzip.open(args.reads_file)
		index_reads_file = gzip.open(args.index_reads_file)
	else:
		reads_file = open(args.reads_file)
		index_reads_file = open(args.index_reads_file)

	count = 0

	while True:
		reads_info = list(islice(reads_file, 4))
		index_info = list(islice(index_reads_file, 4))
		if not reads_info:
			break

		if count % 10000000 == 0:
			print count, '...'

		index = index_info[1].strip()[:args.index_length]

		if index in indices and args.compressed:
			# grab the sample name for this index, then grab the open file 
			# for this sample
			filehandle = index_to_handle_c[indices[index]]
		elif index in indices and not compressed:
			filehandle = index_to_handle[indices[index]]
		else:
			filehandle = bad_index
		
		filehandle.writelines(reads_info)

		count += 1

	# close files
	for x in index_to_handle:
		index_to_handle[x].close()

	reads_file.close()
	index_reads_file.close()
	bad_index.close()



#! python
# -*- coding: utf-8 -*-
"""
bzipscore-all.py

wrapper to score all the interaction pairs in a given fasta file
"""

def get_ids_from_fasta(fasta_file):
    "Returns all the ids (names) from a fasta file."
    from Bio import SeqIO

    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    ids = []
    for fasta in fasta_sequences:
        ids.append(fasta.id)
    
    return ids   


def get_pairs(ids):
    """
    Returns all possible (symetric) pairs from a list of ids.
    
    For example get_pairs([A B C]) returns:
      AA, BB, CC, AB, AC, CB""" 
    import itertools
    return list(itertools.combinations_with_replacement(ids, 2))
    
def get_interaction_string(ids):
    """
    Returns a string of all symetric pairs, comma separated. 
    This string is sutable for writing in a interactions file.
    """
    pairs = get_pairs(ids)
    s = ""
    for p in pairs:
      s = s + ",".join(p)+'\n'
    return s        
            
            
import string
import random

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
   return ''.join(random.choice(chars) for _ in range(size))

         

def write_interaction_file(ids, file_name=None):
    if file_name is None:
        #get temporary name, so nothing is overwritten
        file_name = id_generator(16) + '.int'

    #int_str = get_interaction_string(ids)
    # 'wb -- don't convert new lines lines'
    pairs = get_pairs(ids)
    with open(file_name, 'wb') as f:
         for p in pairs:
            f.write(",".join(p)+'\n')
        
    return file_name


def do_bzip_scoring(fasta_file, sfunc='complete', script = "bzipscore.pl", out_name=None):
    """Writes the interaction file and runs the bzipscore. The result is printed to stdout."""
    import os
    import subprocess
    if not os.path.isfile(fasta_file):
        raise NameError("File "+fasta_file+" does not exist!")
    
    ids = get_ids_from_fasta(fasta_file)
    #print ids
    interactions_file = write_interaction_file(ids)
    
    
    if not os.path.isfile(script):
        script_dir = os.path.dirname(__file__)
        script = os.path.join(script_dir, script)
        assert os.path.isfile(script)

    if out_name is None:
        out_name = os.path.splitext(fasta_file)[0]+'.out'

    out_str = subprocess.check_output(['perl',script, '--sfunc',sfunc, interactions_file, fasta_file ])  
    

    
    with open(out_name, "w") as text_file:
        text_file.write(out_str)

    print "Saved interaction matrix to: "+out_name  
    os.remove(interactions_file)    
    #return out_str


if __name__ == "__main__":
    import argparse
    import os
    parser = argparse.ArgumentParser(
                        description="""Wrapper to score all the interaction pairs in a given fasta file.
                        VERSION: 0.1
                        WARNING: All sequences should start at the f position!""")
    parser.add_argument('fasta_file', type=str,
                       help='path to fasta file (file with sequences)')
    parser.add_argument('--sfunc', choices=['complete','rfe','vinson_ce','fong_svm','bcipa'], default='complete',
                       help='function used for scoring the interactions (default is complete)')
    parser.add_argument('-o', '--out-name',  default=None,
                       help='Outname. Default is inputfile name with .out extension')

    a = parser.parse_args()
    #print args.fasta_file
    #print args.sfunc


    os.system("dos2unix '{}'".format(a.fasta_file) )
    do_bzip_scoring(a.fasta_file, sfunc=a.sfunc,out_name=a.out_name);


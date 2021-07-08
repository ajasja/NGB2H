import os
import hepran.utils as u
import hepran.bzipscore as bz
import hepran
import hepran.bzipscore as bz
from hepran.bzipscore import *
import hepran.utils as u
import hepran.registers as r
import numpy as np
from glob import glob
import pandas as pd
from StringIO import StringIO




def final_name(df_row):
    if u.is_str(df_row):
        name = df_row
    else:
        name = df_row["!full_name"]
    name = os.path.normpath(name)
    p = name.split(os.sep)
    fname = p[-2]+"-"+p[-1]
    return fname
    
#final_name(sdf.iloc[0])    

def open_excel(df_row, in_dir_prefix=""):
    
    from subprocess import call

    if u.is_str(df_row):
        name = df_row
    else:
        name = df_row["!full_name"]
    name = in_dir_prefix + u.replace_extension(name, '.xlsx')
    print(name)
    return call(["open", name])
    #!open $name


def open_path(df_row, in_dir_prefix=""):
    
    from subprocess import call

    if u.is_str(df_row):
        name = df_row
    else:
        name = df_row["!full_name"]
    name = in_dir_prefix + os.path.dirname(name)
    print(name)
    return call(["open", name])
    #!open $name

def copy_to_out(df_row, dirname="!01_OUT", exts=[".set", ".xlsx", ".fasta", ".set_info", ".png", ".bin"], in_dir_prefx="" ,open_in_excel=False):
    from shutil import copyfile
    if u.is_str(df_row):
        name = df_row
    else:
        name = df_row["!full_name"]
   
    
    fname = final_name(name)
    fname = os.path.join(dirname, fname)
    u.ensure_dir(dirname)
    
    name = in_dir_prefx+name

    print(name +"->\n"+ fname+"\n")           
    for ext in exts:
        dest_fname = u.replace_extension(name, ext)
        target_fname = u.replace_extension(fname,ext)
        if os.path.exists(dest_fname):
            copyfile(dest_fname, target_fname)  
        else:
            print("Warning: file does not exists: "+dest_fname)
    
    if open_in_excel:
        open_excel(df_row)

def filter_df(df, peptides_type, set_type):
    """Return dataset filtered acording to above criteria"""
    
    if set_type is None:
        sdf=df[(df.peptides_type==peptides_type)]
    else:        
        sdf=df[(df.peptides_type==peptides_type) & 
           (df.type==set_type)]
    sdf[sdf.N_pairs==sdf.N_pairs.max()]
    sdf = sdf.sort_values(by="N_pairs N_heterodimers total_IN_mismatches total_electrostatic_mismatches orthogonality_average".split(),
                         ascending=[False, False, True, True, False]
                         )
    return sdf
    
display_cols = "finalname basename binding_cutoff N_pairs N_heterodimers N_homodimers\
                 total_IN_mismatches total_electrostatic_mismatches orthogonality_average orthogonality_gap".split()  


def get_out_file(set_name, bcf_name, out_dir="!OUT_bcf"):
    name = os.path.basename(set_name)
    name = os.path.join(out_dir, name)
    name = u.replace_extension(name,bcf_name+'.set')
    #print(name)
    return name
    
#get_out_file('!OUT\\4h-1or2N-not_first_only_B07_bc-8.60_nc-7.60-hetero.01.set','bN')== \
#             '!OUT_bcf\\4h-1or2N-not_first_only_B07_bc-8.60_nc-7.60-hetero.01.bN.set'
    
    


def convert_set_to_bcf(set_name, bcf_name, out_dir="!OUT_bcf", bcf_dict=None, fasta=None):
    """Give backgrounds to a set"""
    ofile = get_out_file(set_name, bcf_name, out_dir)
    a_set = u.load_set(set_name)
    a_IDs = bz.get_ids_from_pairs(a_set)
    
    #local fasta file has priority
    if os.path.isfile(u.replace_extension(set_name, 'fasta')):
        ffasta = u.load_fasta(u.replace_extension(set_name, 'fasta'))
    else:        
        ffasta = bz.filter_fasta(fasta, a_IDs)

    rfasta = u.rename_fasta(ffasta, bcf_name)
    rfasta = r.merge_background_fasta(bcf_dict[bcf_name], rfasta, inplace=True, strip_whitespace=True)
    r_set = u.rename_set(a_set, bcf_name)
    
    u.ensure_dir(out_dir)
    u.save_set(r_set, ofile)
    u.save_dict_to_fasta(rfasta, u.replace_extension(ofile, '.fasta'))        
    
convert_set = convert_set_to_bcf

    
#convert_set(sets[0], "bN")



def get_aliases(final_fasta):
    """Returns the diffrent names for the same sequence"""
    #make dataframe
    dff = pd.DataFrame(data=dict(ID=final_fasta.keys(), seq=final_fasta.values()))
    #check for duplicates
    dff.seq = dff.seq.str.replace("-", "")
    dff.seq = dff.seq.str.replace(" ", "")
    dff['dup']=dff.duplicated(subset='seq', keep=False)
    #group duplicates by sequence and return list
    gb = dff.query('dup==True').groupby(by='seq')
    return gb['ID'].apply(list)


def get_alias_dict(aliases):
    """Takes a dataframe of lists for the same seq and create a dictionary that can be used in normalizing"""
    alias_dict = {}
    for alist in aliases.values:
        first = alist[0]
        others = alist[1:]
        for o in others:
            alias_dict[o] = first

    return alias_dict


def make_aliased_pair(pair, alias_dict, sort=True):
    """Makes an aliased pair that is also sorted"""
    p1 = pair[0]
    p2 = pair[1]
    if p1 in alias_dict:
        p1 = alias_dict[p1]
    if p2 in alias_dict:
        p2 = alias_dict[p2]

    result = (p1, p2)
    if sort:
        result = sorted(result)

    return tuple(result)

assert make_aliased_pair(('P3', 'P15'), alias_dict=dict(P3='P2', P15='P1')) == ('P1', 'P2')


def zip_pairs(pairs, alias_dict):
    zipped_pairs = [make_aliased_pair(p, alias_dict) for p in pairs]
    return list(collections.OrderedDict.fromkeys(zipped_pairs))


def zip_fasta(fasta, alias_dict):
    zipped_fasta = fasta.copy()
    for key in zipped_fasta.keys():
        if (key in alias_dict) and (alias_dict[key] in zipped_fasta):
            del zipped_fasta[key]
    return zipped_fasta


def flip_pairs(pairs):
    flipped_pairs = []
    for p1, p2 in pairs:
        flipped_pairs.append((p1,p2))
        if p1!=p2:
            flipped_pairs.append((p2,p1))
    return flipped_pairs
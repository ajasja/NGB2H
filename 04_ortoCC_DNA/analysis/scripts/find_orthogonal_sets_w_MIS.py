#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 10:48:20 2019

@author: Cliff Boldridge
This file takes in interaction data where two proteins have a measured level of interaction
it maps this data to different sets of proteins that were tried in an all-against-all fashion. Representing the data as
a graph it then finds the largest orthogonal set of interactions for different classifying levels of on/off. Each
maximum independent set is the largest set of orthogonal interactions at that classifying level.
"""

import os
import re
import numpy as np
import networkx as nx


DESIGN_DIR1 = "../designs/180419 18k constructs/OUT_SETS"
DESIGN_DIR2 = "../designs/170317 R8000 all sets !OUT"
DATAFILE1 = "../data/190627_R18k_avg.csv"
DELLIST1 = ['18k_flipped_zipped_final_pairs.fasta', 'flipped_final_pairs.fasta', 'aa_ge_crosspairs.fasta', 'NIEK_crosspairs.fasta']
DATAFILE2 = "../data/200306medR8000.csv"
OUTFILE1 = "210602_R18k_max_ortho_sets_gap0.csv"
OUTFILE2 = "210602_R8000_max_ortho_sets_gap0.csv"
OUTFILE3 = "210602_R18k_max_ortho_sets_gap0.5.csv"
OUTFILE4 = "210602_R8000_max_ortho_sets_gap0.5.csv"
OUTFILE5 = "210602_R18k_max_ortho_sets_gap1.0.csv"
OUTFILE6 = "210602_R8000_max_ortho_sets_gap1.0.csv"
OUTFILE7 = "210602_R18k_max_ortho_sets_gap1RMS.csv"
OUTFILE8 = "210602_R8000_max_ortho_sets_gap1RMS.csv"
OUTFILE9 = "210518_R18k_max_ortho_sets_gap_enforced.csv"

def read_data(file):
    """Takes in interaction data and returns it as a list of lists
    Parameters:
        file: csv file containing data on interactions. should be of the form X_peptide, Y_peptide, interaction score

    Returns:
        list of lists where inner lists are [X_peptide, Y_peptide, interaction score]"""
    protlist = []
    with open(file, 'r') as f:
        f.readline()
        for line in f:
            line = line.strip()
            X_peptide, Y_peptide, medRD = line.split(',')
            protlist.append([X_peptide.upper(), Y_peptide.upper(), np.log(float(medRD) + 0.0001)]) # 0.0001 is well below our threshold of detection
    return protlist


def read_designed_fastas(folder, dellist):
    """Finds all fastas files stored in a folder or subfolder and creates a dictionary of that fasta's name to the
    proteins inside it which are presumed to be tested for interactions in an all-against-all manner. Deletes those
    files that aren't intended to be mapped all-on-all
    Parameters:
        folder: a path to a folder containing fasta files in it or subfolders
        dellist: a list of strings containing the names of fastas that are not to be analyzed
    Returns:
        fasta_lists: a dictionary of fastas filenames to proteins in the fasta files"""
    fasta_lists = {}
    for root, dirs, files in os.walk(folder):
        for name in files:
            curfile = os.path.join(root, name)
            if re.search("fasta", curfile):
                with open(curfile, 'r') as f:
                    count = 0
                    protein_list = []
                    for line in f:
                        if count % 2 == 0:
                            protein_list.append(line[1:].strip())
                        count += 1
                    fasta_lists[name] = protein_list
    for i in dellist:
        del fasta_lists[i]
    return fasta_lists


def convert_fastas_to_R8000_group(fastasdic):
    """The R8000 data is divided into several groups that share backbones. This labels the fasta for the groups
    Parameters:
        fastasdic: dictionary of fasta file names to proteins contained in that file
    Returns:
        a new dictionary where each protein is stamped with a group name
    """
    groups = ["bA", "bN", "bH", "bS", "bP"]
    newdic = {}
    for fasta in fastasdic:
        for groupname in groups:
            newlist = []
            for prot in fastasdic[fasta]:
                newlist.append(prot + "-" + groupname)
            newdic[fasta+"-"+groupname] = newlist
    return newdic


def write_nodeset(filename, grouplist):
    """Takes the maximum orthogonal sets and writes it to a csv file
    Parameters;
        filename: Name of the file to write to
        grouplist: a list of lists of the form [fasta name, classifying level of on/off interaction, difference between
        minimum on target interaction strength and maximum off target interaction strength, the list of proteins
        in this set.
    Returns:
        None
    """
    with open(filename, 'w') as f:
        f.write("name,classifier, minpos, maxneg, lnogap,peptides\n")
        for group in grouplist:
            line = ''
            for cnt, col in enumerate(group):
                line += str(group[cnt]) + ','
            line = line[:-1] + '\n'
            f.write(line)
        f.close()


def select_set(data, group):
    """
    Takes a list of protein pairs in form of [[X_peptide_1,Y_peptide_1,RNA/DNA_1],....]
    and selects those that are in group [peptide1, peptide2, peptide3....]
    Parameters:
        data: a list of lists containing all the data from an experiment [X_protein, Y_protein, interaction strength]
        group: a single all-against-all experiment
    Returns:
        thisset: a list of lists from data where all X_peptide, and Y_peptides appeared in group
    """
    thisset = []
    group = [i.upper() for i in group]
    for peppair in data:
        if peppair[0] in group and peppair[1] in group:
            thisset.append(peppair)
    return thisset


def classify(peplist, classlevel, gap=0):
    """Takes a list of interactions and breaks them out into interactions or non-interactions based on whether their
    interaction strength is greater or less than the classifying level
    Parameters:
        peplist: a list of lists containing [X_protein, Y_protein, interaction strength]
        classlevel: a float that represents what level should be considered an interaction/non-interaction
    Returns:
        a dictionary of interactions (1) and non interactions (0) to lists of lists of interactions"""
    newpeps = {'1': [], '0': []}
    goodpeps = {'1': [], '0': []}
    badpeps = set()
    badpeps2 = {'1': []}
    for pep in peplist:
        if float(pep[2]) > classlevel + gap/2:
            newpeps['1'].append([pep[0], pep[1], pep[2]])
        elif float(pep[2]) < classlevel - gap/2:
            newpeps['0'].append([pep[0], pep[1], pep[2]])
        else:
            badpeps.add(pep[0])
            badpeps.add(pep[1])
            badpeps2['1'].append(pep)
    g = build_normal_graph(badpeps2)
    badnodes = [node for node in nx.nodes_with_selfloops(g)]
    g.remove_nodes_from(badnodes)
    indiv =  maxindset(g)
    badpeps = badpeps.difference(set(indiv))
    for pep in newpeps['1']:
        if pep[0] not in badpeps and pep[1] not in badpeps:
            goodpeps['1'].append(pep)
    for pep in newpeps['0']:
        if pep[0] not in badpeps and pep[1] not in badpeps:
            goodpeps['0'].append(pep)
    return goodpeps


def build_graph(classifiedpeplist):
    """Takes a dictionary of classified interactions and turns it into a graph where each node is a protein and each
    edge is an interaction. Builds the line graph of this (so each node represents an interaction and edges represent
    shared proteins between the nodes. Then add the edges that are one neighbor away since they wouldn't be acceptable
    for orthogonal proteins
    Parameters:
        classifiedpeplist: a dictionary of '1' interactions or '0' non-interactions to lists of lists of protein pairs
        and their interaction strength.
    Returns:
        a networkx graph where the nodes are interacting proteins and the edges represent a shared protein between
        interactions or a neighbor that is shared between interactions"""
    g = nx.Graph()
    for peppair in classifiedpeplist['1']:
        g.add_edge(peppair[0], peppair[1])

    h = nx.line_graph(g)
    #for selfloopnode in nx.nodes_with_selfloops(g):
    #    h.add_node((selfloopnode, selfloopnode))
    #    for nbr in g[selfloopnode]:
    #        if nbr != selfloopnode:
    #            h.add_edge((selfloopnode, selfloopnode), (nbr, selfloopnode))
    hedgelist = []
    for node in h.nodes:
        for nbr in h[node]:
            for nnbr in h[nbr]:
                if nnbr != node:
                    hedgelist.append((node, nnbr))
    h.add_edges_from(hedgelist)
    hnodelist = []
    for selfloopnode in nx.nodes_with_selfloops(g):
        for nbr in h[(selfloopnode, selfloopnode)]:
            if nbr[0] == selfloopnode or nbr[1] == selfloopnode:
                hnodelist.append(nbr)
    h.remove_nodes_from(hnodelist)
    return h

def build_normal_graph(classifiedpeps):
    g = nx.Graph()
    for peppair in classifiedpeps['1']:
        g.add_edge(peppair[0], peppair[1])
    return g


def maxindset(graph):
    """Finds the maximum independent set for a given graph. It does this by solving the equivalent problem of taking
    the complement of the graph and then solving the maximal clique problem which networkx has implemented a nice fast
    algorithm.
    Parameters:
        graph: a networkx graph
    Returns:
        a list of nodes contained in the maximum independent set in the graph presumed to be of the form
         (protein1, protein2)"""
    gprime = nx.complement(graph)
    # This is a rather dumb thing that must be done to remove proteins that are part of homodimers and heterodimers
    # Takes the inverse graph and finds homodimers that share a 2-step neighbor. Cuts that neighbor out of the clique
    #badedges = []
    #for node in nx.nodes(gprime):
    #    for nbr in nx.neighbors(gprime, node):
    #        if nbr == node:
    #            for otnbr in nx.neighbors(gprime, nbr):
    #                if (otnbr == node) and otnbr != nbr:
    #                    badedges.append((otnbr, node))
    #gprime.remove_edges_from(badedges)
    cliqiter = nx.find_cliques(gprime)
    maxcliq = max(cliqiter, key=len, default=[])
    return maxcliq


def parsenodes(maxcliq):
    """Takes a list of nodes from and just finds the unique proteins present
    Parameters:
        maxcliq: a list of proteins in the maximum independant set of interactions in the form (protein1, protein2)
    Returns:
        a list of each of the proteins in the maximal clique"""
    nodes = set()
    for node in maxcliq:
        xpep, ypep = node[0], node[1]
        nodes.add(xpep)
        nodes.add(ypep)
    return list(nodes)


def orthodistance(classdic, orthonodes):
    """Calculates the orthogonality gap for a given group of proteins. Takes those orthogonal nodes and asks what is the
    weakest interaction among the interacting set and what is the strongest interaction in the non-interacting set.
    The difference is the orthogonality gap
    Parameters:
        classdic: dictionary of classified interactions to a list of lists of [X_protein, Y_protein, interaction score]
        orthonodes: A list of those nodes that were orthogonal at a given classification level
    Returns:
        a float that is the orthogonality gap"""
    maxneg, minpos = -1, 100
    for peppair in classdic['1']:
        if peppair[0] and peppair[1] in orthonodes:
            if float(peppair[2]) < float(minpos):
                minpos = float(peppair[2])
    for peppair in classdic['0']:
        if peppair[0] and peppair[1] in orthonodes:
            try:
                if float(peppair[2]) > float(maxneg):
                    maxneg = float(peppair[2])
            except NameError:
                maxneg = float(peppair[2])
    return  float(minpos), float(maxneg), float(minpos) - float(maxneg)

  
def wrapper(data, designs, filename, gap):
    """Main function for finding all the orthogonal interactions in an experiment. Takes a data file, and some designs
      then for all the proteins of a design classifies them, builds them into a graph to find the maximal independent set,
      finds the MIS, calculates the orthogonality gap of it and writes them to file"""
    allmaxes = []
    for fasta in designs:
        oneorthoset = designs[fasta]
        setdata = select_set(data, oneorthoset)
        for classy in [np.log(x*0.1) for x in reversed(range(10, 100))]:
            setdataclone = setdata[:]
            classeddata = classify(setdataclone, classy, gap)
            builtgraph = build_graph(classeddata)
            optnodes = maxindset(builtgraph)
            nodelist = parsenodes(optnodes)
            print("fasta: ", fasta, " classifier: ", classy, " the nodes found ", nodelist)
            minpos, maxneg, lnorthogap = 'bad', 'bad','bad'
            if classeddata['1']:
                minpos, maxneg, lnorthogap = orthodistance(classeddata, nodelist)
            allmaxes.append([fasta, classy, minpos, maxneg, lnorthogap, nodelist])
    write_nodeset(filename, allmaxes)

def large_forced_gaps(data, designs, filename):
    """Main function for finding all the orthogonal interactions in an experiment. Takes a data file, and some designs
      then for all the proteins of a design classifies them, builds them into a graph to find the maximal independent set,
      finds the MIS, calculates the orthogonality gap of it and writes them to file"""
    allmaxes = []
    for gap in [x*0.1  for x in range(0,22)]:
        for fasta in designs:
            oneorthoset = designs[fasta]
            setdata = select_set(data, oneorthoset)
            for classy in [np.log(x*0.1) for x in reversed(range(10, 100))]:
                setdataclone = setdata[:]
                classeddata = classify(setdataclone, classy, gap)
                builtgraph = build_graph(classeddata)
                optnodes = maxindset(builtgraph)
                nodelist = parsenodes(optnodes)
                print("fasta: ", fasta, " classifier: ", classy, " the nodes found ", nodelist)
                minpos, maxneg, lnorthogap = 'bad', 'bad','bad'
                if classeddata['1']:
                    minpos, maxneg, lnorthogap = orthodistance(classeddata, nodelist)
                allmaxes.append([fasta, classy, minpos, maxneg, gap, nodelist])
    write_nodeset(filename, allmaxes)

if __name__ == '__main__':
    R18data = read_data(DATAFILE1)
    fastas1 = read_designed_fastas(DESIGN_DIR1, DELLIST1)
    wrapper(R18data, fastas1, OUTFILE1, 0)
    wrapper(R18data, fastas1, OUTFILE3, 0.5)
    wrapper(R18data, fastas1, OUTFILE5, 1.0)
    wrapper(R18data, fastas1, OUTFILE7, 0.6823)
    large_forced_gaps(R18data, fastas1, OUTFILE9)
    R8000data = read_data(DATAFILE2)
    fastas2 = read_designed_fastas(DESIGN_DIR2, [])
    fastas2 = convert_fastas_to_R8000_group(fastas2)
    wrapper(R8000data, fastas2, OUTFILE2, 0)
    wrapper(R8000data, fastas2, OUTFILE4, 0.5)
    wrapper(R8000data, fastas2, OUTFILE6, 1.0)
    wrapper(R8000data, fastas2, OUTFILE8, 0.7017)


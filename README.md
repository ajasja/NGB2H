# Next generation bacterial two hybrid (NGB2H)
This repository host the code accompanying the article "A Multiplexed Bacterial Two-Hybrid for Rapid Characterization of Protein-Protein Interactions and Iterative Protein Design"

Protein-protein interactions (PPIs) are required for most biological functions as well as applications ranging from drug design to synthetic cell circuits.  As myriad biological functions involve protein-protein interactions (PPIs), engineered PPIs are crucial for applications ranging from drug design to synthetic cell circuits. However, engineering an arbitraryany PPI is still challenging, and designing sets of orthogonal PPIs that do not cross-interact (a necessity for protein nanostructures, synthetic signaling networks and gene circuit design) is even harder. The main issues are inaccurate predictions of interactions and limited ability to assay large numbers of PPIs. Here we address both problems.  First, we developed a method called the Next-Generation Bacterial Two-Hybrid (NGB2H), which combines gene synthesis, a bacterial two-hybrid, and a high-throughput next-generation sequencing readout, allowing rapid exploration of interactions of programmed protein libraries in a quantitative and scalable way. After rigorously validating it, we used the NGB2H system to design, build, and test large sets of orthogonal synthetic coiled-coils. In an iterative set of experiments, we assayed thousands of PPIs, used the datasets to improve the accuracy of coiled-coil scoring algorithms and then built the largest set of orthogonal PPIs identified to date. 

# Contents

The code is divided into several subdomain problems:

- **01_ortoCC_score**: Finding orthogonal sets and fast scoring of Coiled-coil interactions.
- **02_ortoCC_design**: Designing orthogonal sets
- **03_iCipa**: Creating a new scoring function (iCipa)
- **04_ortoCC_design_iCipa**: Designing orthogonal sets using iCipa
- **05_ortoCC_DNA**: scripts used in making the DNA library.
- **06_set_visualizations**: visualization of set orthogonality and heptade alignment
- **07_supplementary_files**: supplementary tables and plasmid sequences 

Each folder contains instructions for installation and usage.

To clone the repository (and initialize submodules) do:
```bash
    git clone --recurse-submodules https://github.com/ajasja/NGB2H
```

# Requirements

The installation time is 30-60 min, depending on previous familiarity with python. 

## Hardware Requirements

Hardware requirments vary based on the size of the orthogonal set being created. 4096*4096 interactions can be scored on a desktop computer with 6 cores @3.3Ghz and 16GB of RAM in a few seconds.


## Software Requirements
The package uses python 3 and Jupyter notebooks. It was tested using python 3.6 and 3.7. We recommend using the conda package manager. Exact versions of all need packages are given in the conda environment file.   



# Citation

When using this package or data please cite [A Multiplexed Bacterial Two-Hybrid for Rapid Characterization of Protein-Protein Interactions and Iterative Protein Design | bioRxiv](https://www.biorxiv.org/content/10.1101/2020.11.12.377184v1)

# Next generation bacterail two hybrid (NGB2H)
This repository host the code accompanying the article "A Multiplexed Bacterial Two-Hybrid for Rapid Characterization of Protein-Protein Interactions and Iterative Protein Design"

Predictiong orthoognality and specificity in 

#Contents

The code is divided into several subdomain problems:

1) 01_ortoCC_score: Finding orthogonal sets and fast scoring of Coild-coil interactions.
2) 02_ortoCC_design: Designing orthogonal sets
3) 03_iCipa: Creating a new scoring function (iCipa)
4) 04_ortoCC_DNA: scripts used in making the DNA library. 

Each folder contains instructions for installation and usage.

#Requirements
##Hardware Requirements

Hardware requirments vary based on the size of the orthogonal set being created. 4096*4096 interactions can be scored on a desktop computer with 4 cores @3.3Ghz and 16GB of RAM in ??30 min??


##Software Requirements
The package uses python 3 and Jupyter notebooks. It was tested using python 3.6 and 3.7. We recommend using the conda package manager. Exact versions of all need packages are given in the conda environment file.   

#Citation

When using this package or data please cite [A Multiplexed Bacterial Two-Hybrid for Rapid Characterization of Protein-Protein Interactions and Iterative Protein Design | bioRxiv](https://www.biorxiv.org/content/10.1101/2020.11.12.377184v1)

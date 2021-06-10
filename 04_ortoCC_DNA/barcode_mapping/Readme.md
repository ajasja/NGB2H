This directory and subdirectories contain the material necessary to take NGS reads of proteins and barcodes and uniquely deidentify each barcode and corresponding protein pair. 

# Requirements

BBSuite tools which can be obtained from https://jgi.doe.gov/data-and-tools/bbtools/.

Starcode which can be obtained from https://github.com/gui11aume/starcode

Python 3. 

Also make sure to configure paths in the Makefile.

Included are reference files in the form name, DNA sequence. 
Raw fastqs need to be obtained from the short read archive at:

# Usage

Run with

```bash
make all
```
To obtain those sequences with insertions and deletions
```bash
make noncodes
```

# Output
Outputs a comma separated text file in the form X-hybrid protein,  Y-hybrid protein, DNA barcode.



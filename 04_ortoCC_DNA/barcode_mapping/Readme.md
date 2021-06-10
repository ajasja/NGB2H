This directory and subdirectories contain the material necessary to take NGS reads of proteins and barcodes and uniquely deidentify each barcode and corresponding protein pair. 

# Requirements

BBSuite tools which can be obtained from https://jgi.doe.gov/data-and-tools/bbtools/.
Starcode which can be obtained from https://github.com/gui11aume/starcode
Python 3. 
Also make sure to configure paths in the Makefile.

# Usage

Run with

```bash
make all

Requires reference files in the form name, DNA sequence and outputs files in the form X-hybrid protein,  Y-hybrid protein, barcode.


Here's where we take all of our NGS reads of DNA barcodes and condense it down to just barcode counts. This is done on an independent basis for each replicate and nucleic acid fraction.

# Requires
Starcode, which can be obtained https://github.com/gui11aume/starcode

awk

Input fastq.gz files of the raw sequencing reads, available at:
# Usage

In bash 
```
make all
```

# Output
A comma seperated file for each unique condition with form barcode, count.

#! /bin/sh

for f in */*.set
do
  echo "$f"
  bzipscore-all.py "${f%.set}.fasta"
  draw-matrix.py "$f"
done
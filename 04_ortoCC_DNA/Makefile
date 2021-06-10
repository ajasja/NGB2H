#===============================================================================
# SHELL
SHELL := /bin/bash

# PATHS
DATA := data
SCRIPTS := ../scripts

# PROGRAMS
BBDUK := bbduk.sh
BBMERGE := bbmerge.sh
BBMAP := bbmap.sh
STARCODE := starcode

# VARS
THREADS := 20 # how many threads to run the BB* portion of pipeline

#===============================================================================

# RECIPIES
all: map xy noncodes
map: $(addprefix pipeline/, $(addsuffix .map.csv, RXYb-c1 R8000-c1  R18k-c1))
xy: $(addprefix pipeline/, $(addsuffix .x-y.txt, RXYb-c1 R8000-c1  R18k-c1))
noncodes: $(addprefix pipeline/, $(addsuffix .x.negs.txt, RXYb-c1 R8000-c1  R18k-c1)) $(addprefix pipeline/, $(addsuffix .y.negs.txt, RXYb-c1 R8000-c1  R18k-c1))

clean:
	rm -f pipeline/*

.PRECIOUS: $(addprefix pipeline/, %.map.csv %.merge.fastq %.filter.fastq %.x.txt %.y.txt %.map.trim.fasta %.x.map.sam %.y.map.sam %.scmerge.txt )

#===============================================================================
# BARCODE MAPPING

pipeline/phiX.fasta:
	curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&amp;id=NC_001422.1&amp;rettype=fasta&amp;retmode=text" >> $@

pipeline/%.filter.fastq: pipeline/phiX.fasta $(DATA)/%_R1_001.fastq.gz $(DATA)/%_R2_001.fastq.gz
	@echo Filtering - $(word 2, $^)
	@$(BBDUK) ref=$< \
	    in1=$(word 2, $^) \
	    in2=$(word 3, $^) \
	    out=$@ \
	    stats=$(@:.fastq=.stats.txt) \
	    pigz=t \
	    k=27 \
	    hdist=1 \
	    ktrim=f \
	    maxns=-1 \
	    overwrite=t \
	    threads=$(THREADS) \
	    -Xmx12g 2> $(@:.fastq=.err)

# merge reads
pipeline/%.merge.fastq: pipeline/%.filter.fastq
	@echo "Merging - $<"
	@$(BBMERGE) in=$< \
	    outm=$@ \
	    outu=$(@:.merge.fastq=.unmerged.fastq) \
	    maxloose=t \
	    interleaved=t \
	    threads=$(THREADS) \
	    2> $(@:.fastq=.err)


# starcode to combine sequencing errors into a single read
pipeline/%.scmerge.txt: pipeline/%.merge.fastq
	@echo "Using starcode to remove sequencing errors - $<"
	@$(STARCODE) -d 3 -t 20 -r 3 -i $<   | \
		awk '{print $$1}' > $@


# concat reads together and map barcodes
pipeline/%.map.csv: pipeline/%.scmerge.txt
	@echo "Barcode mapping - $<"
	@python $(SCRIPTS)/bc-map.py $< \
	    -v \
	    -j15 \
	    -s-20 \
	    --min-reads 1 \
	    -b $(@:.map.csv=.bad-bcs.txt) > $@ \
	    2> $(@:.csv=.err)

#===============================================================================
# GENERATE X - Y map

# appends codon group to end of name
# leaves X/Y in comment section of fasta header
# 


trimmed_fastas/RXYb-c1_x-oligos.fasta: ref_files/R20_x-oligos.csv
	@echo "Converting R20 x oligos to fasta"
	@tr '[atgc]' '[ATGC]' < $< | \
		awk -F, '{if($$1 ~ /^.*-C1/) {print $$1" \n"gensub(/.*AGCGAAACCGTGCGTTTA(.*)AGAAGAGC.*/, "\\1", "g", $$2)} \
		else if ($$1 ~ /^.*-C2/) {print $$1" \n"gensub(/.*TGTCCCAGGTCGCAGTTA(.*)AGAAGAGC.*/, "\\1", "g", $$2)} \
		else if ($$1 ~ /^.*-C3/) print $$1" \n"gensub(/.*TCGCGGAGTTGAGGTTTA(.*)AGAAGAGC.*/, "\\1", "g", $$2)}' > $@
	@echo "All done with R20 x fasta"


trimmed_fastas/RXYb-c1_y-oligos.fasta: ref_files/R20_y-oligos.csv
	@echo "Converting R20 y oligos to fasta"
	@tr '[atgc]' '[ATGC]' < $< | \
	   awk -F, '{print $$1" \n"gensub(/.*GCAGTG(.*)GAGACC.*/, "\\1", "g", $$2)}' > $@
	@echo "All done with R20 y fasta"

trimmed_fastas/R8000-c1_x-oligos.fasta: ref_files/R8000_Design_Full.csv
	@echo "Converting R8000 x oligos to fasta"
	@tr '[atgc' '[ATGC]' < $< | \
		awk -F, '{print $$1" \n"gensub(/.*AAGAGGGACGCAGCATTA(.*)AGAAGAGC.*/, "\\1", "g", $$2)}' > $@
	@echo "All done with R8000 x fasta"

trimmed_fastas/R8000-c1_y-oligos.fasta: ref_files/R8000_Design_Full.csv
	@echo "Converting R8000 y oligos to fasta"
	@tr '[atgc]' '[ATGC]' < $< | \
	   awk -F, '{print $$1" \n"gensub(/.*AGAGCAGT(.*)GAGACC.*/, "\\1", "g", $$2)}' > $@
	@echo "All done with R8000 y fasta"

trimmed_fastas/R18k-c1_x-oligos.fasta: ref_files/all_18k_3c
	@echo "Converting R18k x oligos to fasta"
	@tr '[atgc]' '[ATGC]' < $< | \
		awk -F, '{print $$1" \n"gensub(/.*GCTGGAGGCGAGGTTA(.*)AGAAGAGC.*/, "\\1", "g", $$2)}' > $@
	@echo "All done with R18k x fasta"

trimmed_fastas/R18k-c1_y-oligos.fasta: ref_files/all_18k_3c
	@echo "Converting R18k y oligos to fasta"
	@tr '[atgc]' '[ATGC]' < $< | \
	   awk -F, '{print $$1" \n"gensub(/.*AGAGCAGT(.*)GAGACC.*/, "\\1", "g", $$2)}' > $@
	@echo "All done with R18k y fasta"




pipeline/%.x.txt: pipeline/%.map.csv trimmed_fastas/%_x-oligos.fasta
	@echo "Mapping X variants for $<"
	@awk -F, '{print ">"$$1"\n"substr($$2, 1, length($$2) / 2)}' $< | \
	    bbmap.sh ref=stdin.fasta \
	    in=$(word 2, $^) \
	    outm=$(@:.txt=.sam) \
	    nodisk=t \
	    noheader=t \
	    semiperfectmode=t \
	    maxindel=500 \
	    ambiguous=all \
	    secondary=t \
	    ssao=t \
	    maxsites=1000000 \
	    overwrite=t \
	    -Xmx16g \
	    threads=$(THREADS) 2> $(@:.txt=.err)
	@awk '{print $$3, $$1}' $(@:.txt=.sam) | sort > $@



pipeline/%.map.trim.fasta: pipeline/%.map.csv
	awk -F, '{print ">"$$1"_"$$3"\n"$$2}' $< >$@


# Intermediate Alignment
# ----------------------------------------------
pipeline/%.x.map.sam: pipeline/%.map.trim.fasta trimmed_fastas/%_x-oligos.fasta
	bbmap.sh \
            in=$< \
            ref=$(lastword $^) \
            nodisk=t \
            noheader=t \
            k=8 \
            vslow=t \
            maxindel=500 \
            overwrite=t \
            threads=$(THREADS) \
            outm=$@ \
            outu=$(@:.sam=.no-map.sam) \
            -Xmx16g \
            2> $(@:.sam=.err)




pipeline/%.y.map.sam: pipeline/%.map.trim.fasta trimmed_fastas/%_y-oligos.fasta
	bbmap.sh \
            in=$< \
            ref=$(lastword $^) \
            nodisk=t \
            noheader=t \
            k=8 \
            vslow=t \
            maxindel=500 \
            overwrite=t \
            threads=$(THREADS) \
            outm=$@ \
            outu=$(@:.sam=.no-map.sam) \
            -Xmx16g \
            2> $(@:.sam=.err)


# Classify Negative Controls (via DNA sequence):
# -------------------------------------
# Parse CIGAR string for relevant information
pipeline/%.x.negs.txt: pipeline/%.x.map.sam
	awk '{print $$1, $$4, $$6}' $< | \
            python $(SCRIPTS)/classify-negs.py -j15 - > $@


pipeline/%.y.negs.txt: pipeline/%.y.map.sam
	awk '{print $$1, $$4, $$6}' $< | \
            python $(SCRIPTS)/classify-negs.py -j15 - > $@




pipeline/%.y.txt: pipeline/%.map.csv trimmed_fastas/%_y-oligos.fasta
	@echo "Mapping Y variants for $<"
	@awk -F, '{print ">"$$1"\n"substr($$2, length($$2) / 3)}' $< | \
	    bbmap.sh ref=stdin.fasta \
	    in=$(word 2, $^) \
	    outm=$(@:.txt=.sam) \
	    nodisk=t \
	    noheader=t \
	    semiperfectmode=t \
	    maxindel=500 \
	    ambiguous=all \
	    secondary=t \
	    ssao=t \
	    maxsites=1000000 \
	    overwrite=t \
	    -Xmx16g \
	    threads=$(THREADS) 2> $(@:.txt=.err)
	@awk '{print $$3, $$1}' $(@:.txt=.sam) | sort > $@ 

pipeline/%.x-y.txt: pipeline/%.map.csv pipeline/%.x.txt pipeline/%.y.txt
	@echo "Merging X and Y for $(filter-out $<, $^)"
	@echo "Barcode X_peptide Y_peptide Reads" > $@
	@join -j1 -a1 -a2 -e "NA" -o 0,1.2,2.2 $(filter-out $< ,$^) | \
	    join -j1 -a1 -e "NA" -o 0,1.2,1.3,2.2 - <(awk -F, '{print $$1, $$3}' $< | sort -k1,1) >> $@	

#-------------------------------------------------------------------------------


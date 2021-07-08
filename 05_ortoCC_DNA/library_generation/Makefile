SCRIPTS := ../scripts
CUTOFF ?= 0.05

clean:
	rm -f all-8k.lib.fasta all-8k.lib.err 18k.lib.fasta  18k.lib.err mason.lib.fasta  mason.lib.err


# Note we do not include the full reverse primer as it has a restriction site
# we are searching for. Thus we include all but one base of the site (so that
# the last codon does not produce the site) and paste the rest of the constant
# sequence on later
all-8k.lib.fasta: all-8k.fasta all-8k.pairs re-sites_R8000.txt
	python $(SCRIPTS)/lib-gen.py \
	    --verbose \
	    $< \
	    --pairs $(word 2, $^) \
	    --re-sites $(word 3, $^) \
	    --spacer 'GTCAgaagagcagtgAC' \
	    --forward 'AAGAGGGACGCAGCATTA' \
	    --reverse 'CCAGAga' \
	    --revcomp-x \
	    --cutoff $(CUTOFF) \
	    2> $(@:.fasta=.err) \
	    | awk '{if(NR % 2 == 0) print $$0"GACCTTCAGCGCGACCAGCC"; else print $$0}' \
	    > $@

18k.lib.fasta: 18491_remove.fasta 18491_flipped_zipped_final_pairs.pairs re-sites.txt
	python $(SCRIPTS)/lib-gen.py \
	    --verbose \
	    $< \
	    --pairs $(word 2, $^) \
	    --re-sites $(word 3, $^) \
	    --spacer 'gtcAGAAGAGCCTCGCAGTGac' \
	    --forward 'GCTGGAGGCGAGGTTA' \
	    --reverse 'ggagAGA' \
	    --revcomp-x \
	    --cutoff $(CUTOFF) \
	    2> $(@:.fasta=.err) \
	    | awk '{if(NR % 1 == 0) print $$0"GACCAGAGTACCGCCGT"; else print $$0}' \
	    > $@

mason.lib.fasta: mason.fasta re-sites.txt
	python $(SCRIPTS)/difflength-lib-gen.py \
	    --verbose \
	    $< \
	    --re-sites $(word 2, $^) \
	    --spacer 'gtcAGAAGAGCTTCTCGCAGTGac' \
	    --forward 'GAGTCTGAGCGGCGTTTA' \
	    --reverse 'ggagAGA' \
	    --yBsaI 'GACCGTCAGGCGAGCTAGGC' \
	    --revcomp-x \
	    --cutoff $(CUTOFF) \
	    --outerxfwdlist ACTGGTGCGTCGTCT AGCGAAACCGTGCGT ACCGGTTTCCACGCA GGGTTCGAGCGGGAG TGGGCGCCAAGAACC \
		--outerxrevlist GGCGAACACTTCCCG GGGACGCGGGAGAGT GTTCGCGCTTGGAGG CCTCTCTGCGCGCTA CTTGGGCGCAGTGGG \
		--outeryfwdlist ACTCGACGGCCTCTG CCCGGATCCCTTGCT GCTTCCGCCCGTAGA CGAGCCGTGGTTCCT ACGCCGAATCCCACA \
		--outeryrevlist CTTCAACGCGCGTGT GAGGCCAATGACGCG CCACCAGCGCGAATG AACACCCTCCCAGCG CCCGACGTTCCAGCA \
	    2> $(@:.fasta=.err) \
	    > $@

#===============================================================================

# RECIPIES
all: chops sccount
chops: $(addprefix bcpipeline/, $(addsuffix .chops, R18k-cDNA-1 R18k-cDNA-2 R18k-pDNA-1 R18k-pDNA-2 \
R20_RNA_0h R20_RNA_0.5h R20_RNA_1h R20_RNA_2h R20_RNA_4h R20_DNA_0h R20_DNA_0.5h R20_DNA_1h R20_DNA_2h R20_DNA_4h R8000-cDNA-1 R8000-cDNA-2 R8000-pDNA-1 R8000-pDNA-2 \
Mason_new-1_DNA-1 Mason_new-1_DNA-2 Mason_new-1_RNA-1 Mason_new-1_RNA-2 Mason_old-1_DNA-1 Mason_old-1_DNA-2 Mason_old-1_RNA-1 Mason_old-1_RNA-2 Mason_old-1_DNA-30 Mason_old-1_RNA-30 \
Mason_old-2_DNA-1 Mason_old-2_DNA-2 Mason_old-2_DNA-30 Mason_old-2_RNA-1 Mason_old-2_RNA-2 Mason_old-2_RNA-30 Mason-1_DNA-1 Mason-1_DNA-2 Mason-1_RNA-1 Mason-1_RNA-2 Mason-2_DNA-1 \
Mason-2_DNA-2 Mason-2_RNA-1 Mason-2_RNA-2))
sccount: $(addprefix bcpipeline/, $(addsuffix .sccount, R18k-cDNA-1 R18k-cDNA-2 R18k-pDNA-1 R18k-pDNA-2 \
R20_RNA_0h R20_RNA_0.5h R20_RNA_1h R20_RNA_2h R20_RNA_4h R20_DNA_0h R20_DNA_0.5h R20_DNA_1h R20_DNA_2h R20_DNA_4h R8000-cDNA-1 R8000-cDNA-2 R8000-pDNA-1 R8000-pDNA-2 \
Mason_new-1_DNA-1 Mason_new-1_DNA-2 Mason_new-1_RNA-1 Mason_new-1_RNA-2 Mason_old-1_DNA-1 Mason_old-1_DNA-2 Mason_old-1_RNA-1 Mason_old-1_RNA-2 Mason_old-1_DNA-30 Mason_old-1_RNA-30 \
Mason_old-2_DNA-1 Mason_old-2_DNA-2 Mason_old-2_DNA-30 Mason_old-2_RNA-1 Mason_old-2_RNA-2 Mason_old-2_RNA-30 Mason-1_DNA-1 Mason-1_DNA-2 Mason-1_RNA-1 Mason-1_RNA-2 Mason-2_DNA-1 \
Mason-2_DNA-2 Mason-2_RNA-1 Mason-2_RNA-2))


clean:
	rm -f bcpipeline/*

#==============================================================================
# CHOPPING
bcpipeline/%.bcs: data/%.fastq.gz
	@echo "Extracting barcodes for " $(word 1, $^)
	@zcat $< | awk 'NR%4==2' > $@



bcpipeline/%.chops: bcpipeline/%.bcs
	@echo "Chopping barcodes down to 20bp for" $(word 1, $^)
	@cut -c 1-20 $< | \
	sort -T ./ --parallel=20 | \
	uniq -c | \
	awk '{print $$2,$$1}' > $@ 

bcpipeline/%.sccount: bcpipeline/%.chops
	@echo "Using starcode to remove sequencing errors for $<"
	@starcode -d 1 -t 20 -r 3 -i $< -o $@





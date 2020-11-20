NTHREADS:=4

all: $(addprefix data/interim/seq/, $(addsuffix /read_umi_counts_clean.tsv.gz, $(shell ls -d data/raw/seq/all/*/* | xargs -n1 basename)))

# Symlink fastq files into one directory 
# and organise them by batch (plate) and sample (demultiplexed).
# Also create whitelists of (sample) barcodes.

# NOTE: This has to be adapted to the particular file naming scheme

SAMPLES := $(shell for i in $$(ls data/raw/seq/batch_1/*txt.gz); do \
	j=$$(basename $$i); \
    sample_id=$$(echo $$j | sed -E 's/^(.+)-1-1_.*_lane1(.+)_.*$$/\1-\2/'); \
    sid=$$(echo $$j | sed -E 's/^(.+)-1-1_.*_lane1(.+)_.*$$/\2/'); \
    batch_id=$$(echo $$j | sed -E 's/^(.+)-1-1_.*_lane1(.+)_.*$$/\1/'); \
    batch_name=$$(basename $$(dirname $$i)); \
    mkdir -p data/raw/seq/all/$$batch_id/$$sample_id; \
    ln -sh ../../../$$batch_name/$$j data/raw/seq/all/$$batch_id/$$sample_id/$${sample_id}_1.txt.gz 2> /dev/null;  \
    awk -v sid=$$sid '$$1 == sid {print $$2}' data/raw/tables/$${batch_name}_whitelist.txt > data/raw/seq/all/$$batch_id/$$sample_id/barcodes.txt;\
    done)

# UMI-tools based pipeline
# See more here: https://github.com/CGATOxford/UMI-tools/blob/master/doc/QUICK_START.md

# 1. Extract UMIs

# Barcode pattern is 15 bases for the sample barcode,
# then 15 bases for the UMI
BC_PATTERN:=$(shell printf 'C%.0s' {1..15})$(shell printf 'N%.0s' {1..15})

# We only work with the read 1
data/interim/seq/%/extracted.fastq.gz: data/raw/seq/all/*/%
	# Generate barcodes to be used and corrected
	# from know whitelist barcodes
	mkdir -p data/interim/seq/$*
	@echo "Creating barcodes whitelist $*..."
	awk '{print $$1"\t"1"\t"}' $</barcodes.txt > $(shell dirname $@)/whitelist.txt
	@echo "Extracting UMIs for $*..."
	umi_tools extract --stdin=$</$(shell basename $<)_1.txt.gz \
					  --bc-pattern=$(BC_PATTERN) \
					  --log=log/umi-extract/$*-extracted.log \
					  --stdout=$@ \
					  --filter-cell-barcode \
					  --error-correct-cell \
					  --whitelist=$(shell dirname $@)/whitelist.txt


# 2. Map reads to amplicons

# Mapped file is sorted by default
data/interim/seq/%/mapped.bam: data/interim/seq/%/extracted.fastq.gz
	@echo "Mapping reads for $*..."
	# This might be problematic for large samples when using not enough memory:
	# bwa mem -t $(NTHREADS) $(REF_FASTA) $< | samtools sort -@$(NTHREADS) -o $@ -
	bwa mem -t $(NTHREADS) $(REF_FASTA) $< | samtools view -bS -o $@
	samtools sort -o $@.sorted $@
	mv $@.sorted $@



# Index BAM files
data/interim/seq/%/mapped.bam.bai: data/interim/seq/%/mapped.bam
	@echo "Indexing $* BAM file..."
	samtools index $<



# 2.1. Clean mapped reads by subsetting known viral regions

data/interim/seq/%/mapped_sub.bam: data/interim/seq/%/mapped.bam data/interim/seq/%/mapped.bam.bai
	@echo "Subsetting %* BAM file..."
	samtools view -H $< > $(basename $@).sam

	# E: 140-200 for shorter reads, 70-200 for longer reads (121bp)
	samtools view $< "E_gene_amplicon:150-200" >> $(basename $@).sam
	# If UMI blacklists to be used:
	# 	samtools view $< "E_gene_amplicon:140-200" | \
	# 		grep -v -F -f $(shell dirname $@)/umi_blacklist.txt >> $(basename $@).sam

	# nsp14: 150-200 for shorter reads, 70-200 for longer reads (121bp)
	samtools view $< "nsp14_amplicon:150-200" >> $(basename $@).sam

	# RP
	samtools view $< "RP_amplicon:10-140"     >> $(basename $@).sam
	samtools view -S -b $(basename $@).sam > $@
	rm $(basename $@).sam


# Index BAM files
data/interim/seq/%/mapped_sub.bam.bai: data/interim/seq/%/mapped_sub.bam
	@echo "Indexing $* BAM file..."
	samtools index $<


# 3. Count UMIs per amplicon

# READS: For all cleaned reads
data/interim/seq/%/read_umi_counts_clean.tsv.gz: data/interim/seq/%/mapped_sub.bam data/interim/seq/%/mapped_sub.bam.bai
	@echo "Computing read counts per amplicon for $*..."
	samtools view -F 4 $< | \
		awk -F'\t' '{split($$1,s,"_")}; {print s[2], s[3], $$3}' | \
		sort | uniq -c | \
		awk 'BEGIN{OFS="\t"}; {print $$4, $$2, $$3, $$1}' | sort | \
		awk 'BEGIN{print "gene\tcell\tumi\tcount"}; {print}' | \
		gzip -9 \
		> $@



# Create genome index

REFDIR:=data/ref
REF_FASTA:=$(REFDIR)/amplicons.fa
REF_INDEX:=$(REFDIR)/amplicons.bwt

$(REF_INDEX): $(REF_FASTA)
	bwa index $<


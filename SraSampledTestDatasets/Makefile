SHELL:=/bin/bash

all: tables figures
	@echo Hallo

report:
	R -e "rmarkdown::render('sra_sampled_datasets-analysis.R', output_file='test-datasets-analysis.html')"

.PHONY: tables figures

tables: TODO

figures: TODO

krakenhll_results: ../../Databases/std ../../Databases/nt
	./run_krakenhll.sh

fastq:
	wget CCB_FTP/fastqs TODO

tax-sra-file1.txt:
	ls fastq-samples1/*.1.fq | sed 's#fastq-samples1/taxid##;s/reads.1.fq//' | tr . '\t' > tax-sra-file1.txt

truth_sets: fastq tax-sra-file1.txt
	mkdir -p truth_sets.tmp
	#for F in fastq/*1.fq.gz; do \
	#	T=truth_sets.tmp/`basename $$F .1.fq.gz`_TRUTH.txt; 
	for F in fastq/*1.fq; do \
		T=truth_sets.tmp/`basename $$F .1.fq`_TRUTH.txt; \
	   	join -1 2 -2 1 -t$$'\t' <( sort -k2,2 tax-sra-file1.txt ) <( cat $$F | sed -n '1~4p' |  sed 's/.//;s/\..*//' | uniq | sort) > $$T; \
	done
	mv truth_sets.tmp truth_sets


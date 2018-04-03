
THREADS=10
KRAKENHLL_DIR=../install
KRAKENHLL_SRC_DIR=../krakenhll

all: 
	@echo "Specify one of the targets std, nt, or nt-helminths."

nt: nt/database.kdb

std: $(KRAKENHLL_DIR)/krakenhll-build std/download.complete
	/usr/bin/time -v -o logs/std-build-time.log $(KRAKENHLL_DIR)/krakenhll-build --db std --threads $(THREADS) --build --taxids-for-genomes --taxids-for-sequences


std/download.complete: $(KRAKENHLL_DIR)/krakenhll-download
	time $(KRAKENHLL_DIR)/krakenhll-download --db std taxonomy contaminants
	## Execute download of microbial sequences twice, as files fail sometimes
	for i in 1 2; do \
		time $(KRAKENHLL_DIR)/krakenhll-download --db std --threads $(THREADS) --dust --include-viral-neighbors refseq/viral/Any refseq/archaea refseq/bacteria; \
	done
	time $(KRAKENHLL_DIR)/krakenhll-download --db std --threads $(THREADS) --fna rna,genomic refseq/vertebrate_mammalian/Chromosome/taxid9606
	touch std/download.complete

nt/database.kdb: $(KRAKENHLL_SRC_DIR)/indices/Makefile $(KRAKENHLL_SRC_DIR)/indices/microbial-taxa.tsv
	make -f $(KRAKENHLL_SRC_DIR)/indices/Makefile TAXONLIST=$(KRAKENHLL_SRC_DIR)/indices/microbial-taxa.tsv nt

std-stats: std/database.report.tsv
	awk '$$7 == 2 || $$7 == 10239 || $$7 == 9606 || $$7 == 32630 || $$7 == 2157' std/database.report.tsv

nt-stats: nt/database.report.tsv $(KRAKENHLL_SRC_DIR)/indices/microbial-taxa.tsv
	join -1 7 -2 1 <( sort -k 7,7 nt/database.report.tsv ) <( sort -k 1,1 $(KRAKENHLL_SRC_DIR)/indices/microbial-taxa.tsv )

nt-helminths: $(KRAKENHLL_SRC_DIR)/indices/Makefile $(KRAKENHLL_SRC_DIR)/indices/microbial-taxa-incl-helminths.tsv
	make -f $(KRAKENHLL_SRC_DIR)/indices/Makefile TAXONLIST=$(KRAKENHLL_SRC_DIR)/indices/microbial-taxa-incl-helminths.tsv nt





all: tables figures
	@echo Hallo

report:
	R -e "rmarkdown::render('test-datasets-analysis.R', output_file='test-datasets-analysis.html')"

.PHONY: tables figures
tables: 
	cp -v tables_formatted/test_datasets-f1-recall-for-all-tools-bio.xlsx ../Tables/Table2-F1_Recall-AllTools-Bio.xlsx
	cp -v tables_formatted/test_datasets-f1-recall-for-all-tools-sim.xlsx ../Tables/SupplTable4-F1_Recall-AllTools-Sim.xlsx
	cp -v tables_formatted/test_datasets-f1-recall-for-krakenuniq.xlsx ../Tables/Table1-F1_Recall-KrakenUniq.xlsx
	cp -v tables_formatted/test_datasets-description-plus-timings.xlsx ../Tables/SupplTable3-TestDatasetsDescriptionAndTimings.xlsx

figures:
	cp -v figures/scatter-hmp-vs-hc1.pdf ../Figures/SupplFigure4-ReadsVsKmers-HmpAndHC1.pdf


truth_sets:
	wget https://ftp-private.ncbi.nlm.nih.gov/nist-immsa/IMMSA/truth_sets.zip
	unzip truth_sets.zip -x "__MACOSX/*"
	mv truth_sets.zip archive

truth_sets/fixed1: truth_sets truth_sets/patches1
	mkdir -p truth_sets/fixed1/{species,genus}
	for SG in species genus; do \
		for TRUTH in truth_sets/$$SG/*_TRUTH.txt; do \
			BN=`basename $$TRUTH`; \
			NTRUTH=truth_sets/fixed1/$$SG/$$BN; \
			PATCH=truth_sets/patches1/$$SG/$${BN/_TRUTH.txt/_TRUTH_fix.patch}; \
			[[ -f $$PATCH ]] && patch $$TRUTH $$PATCH -o $$NTRUTH || cp $$TRUTH $$NTRUTH; \
		done \
	done
	
truth_sets/fixed2: truth_sets/fixed1 truth_sets/patches2
	mkdir -p truth_sets/fixed2/{species,genus}
	for SG in species genus; do \
		for TRUTH in truth_sets/fixed1/$$SG/*_TRUTH.txt; do \
			BN=`basename $$TRUTH`; \
			NTRUTH=truth_sets/fixed2/$$SG/$$BN; \
			PATCH=truth_sets/patches2/$$SG/$${BN/_TRUTH.txt/_TRUTH_fix.patch}; \
			[[ -f $$PATCH ]] && patch $$TRUTH $$PATCH -o $$NTRUTH || cp $$TRUTH $$NTRUTH; \
		done \
	done

tool_results_reformatted/genus:
	wget https://pbtech-vc.med.cornell.edu/git/mason-lab/benchmarking_metagenomic_classifiers/raw/master/genus_results_formatted.tar.gz
	tar xvvf genus_results_formatted.tar.gz
	mv genus_results_formatted.tar.gz archive

tool_results_reformatted/species:
	wget https://pbtech-vc.med.cornell.edu/git/mason-lab/benchmarking_metagenomic_classifiers/raw/master/species_results_formatted.tar.gz
	tar xvvf species_results_formatted.tar.gz
	mv species_results_formatted.tar.gz archive

krakenuniq_results: ../../Databases/std ../../Databases/nt ../../Databases/orig
	./run_krakenuniq.sh

fastq:
	wget -r -l1 --no-parent --directory-prefix fastq.tmp --no-directories -A "*.gz" -R "P00*" -R "*NegControl*" https://ftp-private.ncbi.nlm.nih.gov/nist-immsa/IMMSA
	mv fastq.tmp fastq

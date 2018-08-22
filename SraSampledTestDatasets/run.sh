#!/bin/bash

## This folder contains scripts for and results from gathering FASTQs from SRA for 
## bacterial assemblies in RefSeq, and classifiying them with Kraken

## Write BioSample XMLs to biosample-xmls/
./get-biosample-xmls.pl biosample-xmls/ ../library/bacteria/assembly_summary_filtered.txt

## Get SRA IDs from BioSample XMLs and get more information on them to bactSra_{urls.tab.txt,info.tab.txt}
grep -rhom1 '<Id db="SRA">[^<]*' biosample-xmls | sed 's/<.*>//' | ./getSraRunsFromAccIds.py --outStub bactSra

## Download FASTQ files and map them with Kraken
perl dl-fastqs-and-map.pl biosample-xmls/ bactSra_info.tab.txt ../library/bacteria/assembly_summary_filtered.txt

## Refactor ...
join -t$'\t' -1 1 -2 2 <(cut -f 6,7 ../library/bacteria/assembly_summary_filtered.txt | sort -k 1,1 | uniq)  <(cut -f2,4,7,9  bactSra_info.tab.txt | grep -v dbgap | grep ILLUMINA | sed 's#[^\t]*/##' | sort -k 2,2 -k 1,1 | uniq -f 1) > tax-sra-file.txt
## TODO: Only use big files
## cut -f4,7,9,13,18  bactSra_info.tab.txt | awk -F$'\t' '$5 > 1000000'

## TODO: Check that number of reads (reported in bactSra_info.tab) is right

## Sample from the files
perl dl-fastqs-and-map.pl -S fastq-samples1 biosample-xmls/ bactSra_info.tab.txt ../library/bacteria/assembly_summary_filtered.txt
for S in 1; do
	cat fastq-samples$S/taxid*reads.1.fq > sampled-res/sample$S.1.fq
	cat fastq-samples$S/taxid*reads.2.fq > sampled-res/sample$S.2.fq
	krakenuniq --paired --threads 20 --fastq --db ../../dbs/refseq-oct2017-k31 --report-file sampled-res/sample$S.report.tsv --output off sampled-res/sample$S.1.fq sampled-res/sample$S.2.fq
	krakenuniq --paired --threads 20 --fastq --db ../../dbs/refseq-oct2017-k31 --report-file sampled-res/sample$S.report.tsv --output off sampled-res/sample$S.1.fq sampled-res/sample$S.2.fq
done
ls fastq-samples1/*reads.1.fq | grep -o 'taxid[0-9]*' | sed 's/taxid//' > sampled-res/true-positives.txt

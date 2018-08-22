#!/bin/bash

set -xeu

THREADS=10

DB=${1:-NA}
DB_PATH="../Databases/$DB"

if [[ ! -f "$DB_PATH/datbase.kdb" ]]; then
  echo "Usage: `basename $0` DB" >&2
  echo "DB specifies a database present in the directory ../Databases" >&2
  exit 1
fi

  ## Preload database
  DB=`basename $DB_PATH`
  echo $DB_PATH $DB
  kraken --preload --db $DB_PATH --fasta <(printf ">A\nA")
  HR=krakenuniq_results/$DB
  mkdir -p ${HR}/{report,output,log}
  for R in ../fastq/*.gz; do
	PARAM="--db ../../Databases/$DB $FAQ --threads $THREADS --gzip"
	RR=${R%.gz}        ## Remove gz ending
	RR=${RR#../fastq/} ## Remove ../fastq/ from the beginning

	## Handle paired and single-end files separately
	if [[ "$R" == *"-R1"* ]] || [[ "$R" == *"_R1_"* ]]; then
	  R1=$R
	  R2=${R/R1/R2}
	  RR=${RR/R1/paired}
	  if [[ -f "$R2" ]]; then
		echo "ERROR: Got paired file $R1, but I cannot find $R2!" >&2;
		exit 1;
	  fi
	  PARAM="$PARAM --paired $R1 $R2"
	elif [[ "$R" == *"-R2"* ]] || [[ "$R" == *"_R2_"* ]]; then
	  continue;
	else
	  PARAM="$PARAM $R"
	fi
	[[ -s $HR/report/$RR.krakenuniq.report.tsv ]] || /usr/bin/time -vo $HR/log/$RR.krakenuniq.tlog ../../install/krakenuniq --report-file $HR/report/$RR.krakenuniq.report.tsv --output $HR/output/$RR.tsv $PARAM 2> $HR/log/$RR.krakenuniq.log
	if [[ "$DB" == "std" ]]; then
	  KR=kraken_results/$DB
	  mkdir -p ${KR}/{report,output,log}
	  [[ -s $KR/output/$RR.tsv ]] || /usr/bin/time -vo $KR/log/$RR.kraken.tlog kraken --output $KR/output/$RR.tsv $PARAM 2> $KR/log/$RR.kraken.log
	  [[ -s $KR/report/$RR.report.tsv ]] || /usr/bin/time -vo $KR/log/$RR.kraken-report.tlog kraken-report  $PARAM $KR/output/$RR.tsv $KR/report/$RR.report.tsv
	fi
  done

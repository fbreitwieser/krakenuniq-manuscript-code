#!/bin/env perl
use strict;
use warnings;
use File::Fetch;
use Parallel::ForkManager;
use File::Basename;
use Getopt::Std;
use Term::ANSIColor;
use Time::HiRes;

my %options=(t=>10, m=>5);
my $db = "../../dbs/refseq-oct2017-k31";
my $db_euk = "../../dbs-kraken/refseq-euk-k31/";

my $USAGE = basename($0)." OPTS <biosample xml dir> <sra info txt> <assembly summary file>

Downloads FASTQ files for assemblies.

OPTS
	-t NUM  Number of threads, default $options{t}.
	-s      Run Kraken against standard database.
	-e      Run Kraken against eukaryote database.
	-E      Run Kraken against both eukaryote and standard database
	-S DIR  Sample reads from the FASTQs into DIR!
	-m NUM  Sample up to 10^NUM reads, logarithmically distrubted. Default $options{m}.

Standard database:  $db
Eukaryote database: $db_euk
";

getopts("ht:seES:", \%options);

sub system_l {
  my $start_time = [Time::HiRes::gettimeofday()];
  print STDERR colored(join(" ", @_)."\n", "blue");
  system(@_);
  my $timing = (Time::HiRes::tv_interval($start_time));
  my ($s,$m,$h) = ($timing, 0, 0);
  if ($s > 60) { $m = floor($s/60); $s -= 60*$m; }
  if ($m > 60) { $h = floor($m/60); $m -= 60*$h; }
  $s = "$m:$s" if $m>0;
  $s = "$h:$s" if $h>0;

  print STDERR colored("Took $s\n", "blue");
}

#my $pm = new Parallel::ForkManager($options{"t"});

defined $ARGV[1] || die $USAGE;
my $BDIR=shift;
-d $BDIR || mkdir($BDIR);
my $sra_info_f=shift;

my ($OF1, $OF2);
if (defined $options{"S"}) {
	-d $options{"S"} || system_l("mkdir -p $options{S}");
	open($OF1, ">", ($options{"S"})."/sampled_reads.1.fq");
	open($OF2, ">", ($options{"S"})."/sampled_reads.2.fq");
}

my %taxids_w_sras;
my %sampled_species;
my @all_sras;
while (<>) {
	next if /^#/;
	my ($assembly_accession,$bioproject,$biosample,$wgs_master,$refseq_category,$taxid,$species_taxid,$organism_name,$infraspecific_name,$isolate,$version_status,$assembly_level,$release_type,$genome_rep,$seq_rel_date,$asm_name,$submitter,$gbrs_paired_asm,$paired_asm_comp,$ftp_path,$excluded_from_refseq,$relation_to_type_material) = split(/\t/);

	next unless defined $biosample;
	my $biosample_entryf = "$BDIR/taxid$taxid-biosample$biosample.xml";
	if ((defined($options{"S"}) && !defined $sampled_species{$species_taxid}) || (!defined $options{"S"} && !defined $taxids_w_sras{$taxid})) {
		#my $pid = $pm->start and next;
		if (!-s $biosample_entryf) {
			print STDERR "No BioSample XML for $biosample / $organism_name";
			next;
			#system_l("curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=$biosample&retmode=xml' > $biosample_entryf 2> /dev/null");
		}
		my $SRA = `grep -o '<Id db="SRA">[^<]*' $biosample_entryf | sed 's/.*>//'`; chomp $SRA;
		if (length($SRA) > 3) {
			my $has_illumina = `grep -m1 "ILLUMINA.*$SRA" $sra_info_f | cut -f 2`; chomp $has_illumina;
			if ($has_illumina ne "") {
				my $SRR = basename($has_illumina);
				$SRR =~ s/.sra//;
				print STDERR "$SRR ($organism_name, taxid $species_taxid): ";
				#system_l("wget $has_illumina");
				my @files = glob("fastq/$SRR*.fastq.gz");
				if (@files) {
					print STDERR colored("FASTQ", "green"),"\n";
				} else {
    				#$pm->finish && next;
					print STDERR colored("no FASTQ", "red"),"\n";
					#system_l("fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip $SRR");
					next;
					@files = glob("fastq/$SRR*.fastq.gz");
				}

				my $basename = "taxid$species_taxid.$SRR";
				if (-f "report/$SRR.report.tsv" && !-f "report-s/$basename.report.tsv") {
					system_l("cp -v report/$SRR.report.tsv report-s/$basename.report.tsv");
				}
				
				if (@files == 3) {
					# Remove single-ended file
					@files = grep(!/_pass.fastq.gz/, @files);
				}
				next if @files > 2 || @files < 1;
				next if @files != 2; ## temp - only take data from paired experiments

				my $kraken_cmd = "krakenuniq --threads $options{t} --fastq --gz";
				$kraken_cmd .= " --paired" if (@files == 2);
			
				if (defined $options{"s"} && !-s "report-s/$basename.report.tsv") {
					system_l("$kraken_cmd --db $db --report-file report-s/$basename.report.tsv --output kraken/$basename.kraken.tsv @files");
				}

				if (defined $options{"e"} && !-s "report-euk/$basename.report.tsv") {
					system_l("$kraken_cmd --db $db_euk --report-file report-euk/$basename.report.tsv --output kraken-euk/$basename.kraken.tsv @files");
				}

				if (defined $options{"E"} && !-s "report-euk/$basename.report.tsv") {
					system_l("krakenuniq --db $db --db $db_euk --paired --report-file report-euk2/$basename.report.tsv --fastq --gz --threads 2 --output kraken-euk2/$basename.kraken.tsv @files");
				}
				if (-s "report-s/$basename.report.tsv" && defined $options{"S"}) {
					my ($perc, $reads) = split(/\t/, `awk '\$7 == $species_taxid' report-s/$basename.report.tsv`);
					if (defined $reads && $reads > 1000000 && $perc > 75) {
						$sampled_species{$species_taxid} = 1;
						my $n_reads_to_sample = int(10**(rand(4)));
						my $current_reads = 0;
						print "$perc\t$reads Sampling $n_reads_to_sample\n";
						open(my $IF1, "-|", "gunzip -c $files[0]");
						open(my $IF2, "-|", "gunzip -c $files[1]");
						#if ($individual_files_per_taxon) {
						#	open(my $OF1, ">", ($options{"S"})."/$basename.${n_reads_to_sample}reads.1.fq");
						#	open(my $OF2, ">", ($options{"S"})."/$basename.${n_reads_to_sample}reads.2.fq");
						#}
						my $prob = $n_reads_to_sample / 1000000; # Choose every read with this chance
						while (!eof($IF1) && !eof($IF2) && $n_reads_to_sample > $current_reads) {
							if (rand() < $prob) {
								++$current_reads;
								for (my $i=0; $i<4; ++$i) {
									print $OF1 scalar readline($IF1);
									print $OF2 scalar readline($IF2);
								}
							} else {
								for (my $i=0; $i<4; ++$i) {
									<$IF1>;
									<$IF2>;
								}
							}
						}
						close($IF1);
						close($IF2);
						#if ($individual_files_per_taxon) {
						#	close($OF1);
						#	close($OF2);
						#}
					}
				}

				$taxids_w_sras{$taxid} = 1;
			}
			#push @all_sras, $SRA;
		}
    	#$pm->finish; # Terminates the child process
	}
}
#$pm->wait_on_children;

if (defined $options{"S"}) {
	close($OF1);
	close($OF2);
}

#system_l("./getSraRunsFromAccIds.py --identifiers ".(join(",",@all_sras))." --outStub bacteria_sras");

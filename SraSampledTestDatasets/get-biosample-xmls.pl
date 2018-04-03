#!/bin/env perl
use strict;
use warnings;
use File::Fetch;
use Parallel::ForkManager;
use File::Basename;

my $pm = new Parallel::ForkManager(10);
my $USAGE = basename($0)." <biosample xml dir> <assembly summary file>

Downloads BioSample XML files for all assembly projects to <biosample xml dir>
";
defined $ARGV[0] || die $USAGE;
my $BDIR=shift;
-d $BDIR || mkdir($BDIR);

my %taxids_w_sras;
my @all_sras;
my $i=0;
while (<>) {
	next if /^#/;
	my ($assembly_accession,$bioproject,$biosample,$wgs_master,$refseq_category,$taxid,$species_taxid,$organism_name,$infraspecific_name,$isolate,$version_status,$assembly_level,$release_type,$genome_rep,$seq_rel_date,$asm_name,$submitter,$gbrs_paired_asm,$paired_asm_comp,$ftp_path,$excluded_from_refseq,$relation_to_type_material) = split(/\t/);

	next unless defined $biosample;
	my $biosample_entryf = "$BDIR/taxid$taxid-biosample$biosample.xml";
	if (!defined $taxids_w_sras{$taxid}) {
		++$i;
		if (! -s $biosample_entryf) {
			print STDERR "\r$i: Getting BioSample $biosample for $organism_name";
			my $pid = $pm->start and next;
			system("curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=$biosample&retmode=xml' > $biosample_entryf 2> /dev/null");
    		$pm->finish; # Terminates the child process
		} else {
			print STDERR "\r$i: Got BioSample $biosample for $organism_name";
		}
		#my $SRA = `grep -o 'Id db="SRA".[^<]*' $biosample_entryf | sed 's/.*>//'`; chomp $SRA;
		#if (length($SRA) > 3) {
		#	print "$SRA\n";
		#}
	}
}
$pm->wait_all_children;
print STDERR "\r\nGot BioSample XMLs for $i assemblies.\n";

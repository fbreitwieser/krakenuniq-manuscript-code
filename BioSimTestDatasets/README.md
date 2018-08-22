# McIntyre et al. datasets analysis

This folder contains scripts for the analysis of the 35 datasets compiled by McIntyre et al. 

### Truth sets
- `truth_sets/{species,genus}` contains truth information from https://ftp-private.ncbi.nlm.nih.gov/nist-immsa/IMMSA/truth_sets.zip
- `truth_sets/patches1/` contains patches that fix mistakes in three truth sets (applied to old and new results)
- `truth_sets/fixed1/` contains truth sets with patches 1 applied
- `truth_sets/patches2/` contains patches for the truth sets to reflect NCBI taxonomy changes since the original publication
- `truth_sets/fixed2` contains truth sets with patches 1 and 2 applied

### Tool results
- `tool_results_reformatted/` contains genus and species-level results from https://pbtech-vc.med.cornell.edu/git/mason-lab/benchmarking_metagenomic_classifiers

There are results for BLAST (filtered with Megan), Clark, Clark Spaced, Diamond (filtered with Megan), Gottcha, Kraken, LMAT, MetaFlow, MetaPhlAn, NBC and PhyloSift.

### Kraken and KrakenUniq results
- `krakenuniq_results/{std,nt,orig}/` contains report and log files for KrakenUniq against three different databases. `orig` is the database used in the original publication.
-`kraken_results/std/log` contains log files for the runs of Kraken on test datasets, used for comparison with KrakenUniq

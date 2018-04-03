#'--- 
#'title: "SRA sampled datasets analysis"
#'author: "Florian Breitwieser"
#'date: March 2018
#'output:
#'  html_document:
#'    style: yeti
#'    toc: true
#'    toc_float: true
#'---

#'<style>
#'pre { overflow-x: auto; }
#'pre code { word-wrap: normal; white-space: pre; }
#'</style>

source("../BioSimTestDatasets/functions.R")

RANKS=c("genus", "species")
DBS <- c("orig", "std", "nt")

truth_set <- read.delim("tax-sra-file1.txt", col.names = c("species_taxid", "SRR", "n_reads"), stringsAsFactors = FALSE)

bla <- list(species=list(), genus=list())
for (frac in c("frac0.03", "frac0.06", "frac0.15", "frac0.3", "frac0.6", "full")) { 
    report <- read.delim(sprintf("krakenhll_results/std/sample1.%s.krakenhll.report.tsv", frac), stringsAsFactors = FALSE)
                                        #truth_set <- read.delim(sprintf("truth_sets/species/sample1.%s_TRUTH.txt", frac), col.names = c("SRR", "species_taxid", "species_taxid", "n_reads"), stringsAsFactors = FALSE)
    report$TP <- report$taxID %in% truth_set$species_taxid
    report$depth <- nchar(sub("[^ ].*", "", report$taxName))/2
    tp_depth = 100
    for (i in rev(seq_len(nrow(report)))) {
        if (report[i, "TP"]) {
            tp_depth <- report[i, "depth"]
        } else if (report[i, "depth"] < tp_depth) {
            report[i, "TP"] <- TRUE
            tp_depth <- report[i, "depth"]
        }
    }
    sum(report$TP[report$rank=="species"])
    report <- report[!report$taxID %in% c(9606, 9605), ] ## Ignore human identifications
    report <- report[order(report$kmers, decreasing = TRUE),]
    rres <- data.frame()
    for (sg in c("genus", "species")) {
        rs <- subset(report, rank==sg)
        print(plot_reads_vs_kmers_gg(rs, min.reads = 10, add.f1 = TRUE, title=paste(frac, sg)))
        
        rrs <- rs[order(-rs$reads),]
        stats_r <- get_f1_and_recall(rrs$TP, max_fdr = 0.05, sum(rrs$TP))
        stats_r['f1_thresh'] <- rrs[stats_r['which.max_f1'],"reads"]
        stats_r['recall_thresh'] <- rrs[stats_r['which.max_recall'],"reads"]
        read_thresh <- rrs$reads[stats_r['which.max_recall']]
        rrs <- rs[order(-rs$kmers),]
        stats_k <- get_f1_and_recall(rrs$TP, max_fdr = 0.05, sum(rrs$TP))
        stats_k['f1_thresh'] <- rrs[stats_k['which.max_f1'],"kmers"]
        stats_k['recall_thresh'] <- rrs[stats_k['which.max_recall'],"kmers"]
        kmer_thresh <- rrs$kmers[stats_k['which.max_recall']]

        res <- data.frame(max=c(stats_r[1:2], stats_k[1:2]), which=c(stats_r[3:4], stats_k[3:4]), thresh=c(stats_r[5:6], stats_k[5:6]))
        rownames(res) <- c("F1 reads", "Recall reads", "F1 k-mers", "Recall k-mers")
        bla[[sg]][[frac]] <- res
    }
    rs <- report
    
}
bla[["species"]] %>% sapply(function(x) x["Recall k-mers", c("max", "thresh")])
bla[["genus"]] %>% sapply(function(x) x["Recall k-mers", c("max", "thresh")])

sdf <- data.frame(recall=bla[["species"]] %>% sapply(function(x) x["Recall k-mers", "max"]),
                  thresh=bla[["species"]] %>% sapply(function(x) x["Recall k-mers", "thresh"]),
                  rank="species")
sdf$frac <- c(0.03, 0.06, 0.15, 0.3, 0.6, 1)
sdf$nr <- c(10^6, 2*10^6, 5*10^6, 10^7, 2*10^7, 34.3*10^6)

gdf <- data.frame(recall=bla[["genus"]] %>% sapply(function(x) x["Recall k-mers", "max"]),
                  thresh=bla[["genus"]] %>% sapply(function(x) x["Recall k-mers", "thresh"]),
                  rank="genus")
gdf$frac <- c(0.03, 0.06, 0.15, 0.3, 0.6, 1)
gdf$nr <- c(10^6, 2*10^6, 5*10^6, 10^7, 2*10^7, 34.3*10^6)

34000000/(sdf$thresh/sdf$frac)
(34300000*sdf$frac)/sdf$thresh
sdf$thresh/(34.3*sdf$frac)

ggplot(rbind(sdf,gdf)) + geom_point(aes(x=frac, y=thresh, size=recall)) + facet_wrap(~rank) + geom_abline(xintercept=0, slope=3500/34.3 * 10^3)

ggplot(rbind(sdf, gdf)) + geom_point(aes(x=nr, y=thresh, size=recall, col=recall), alpha=.8) +
    facet_wrap(~rank) + geom_abline(xintercept=0, slope=2000/10^6, alpha=.8) +
    scale_color_gradient(high = "#132B43", low = "#56B1F7") +
    xlab("Sequencing depth") + ylab("Unique k-mer threshold") +
    cowplot::background_grid(major="xy") + scale_x_continuous(breaks=c(10^6, 10*10^6,20*10^6, 30*10^6))
ggsave("figures/sra_sampled_dataset-kmer_thesholds.pdf", width=6, height=2.5, dpi=600, scale=1.2)



frac <- "frac0.3"
report <- read.delim(sprintf("krakenhll_results/std/sample1.%s.krakenhll.report.tsv", frac), stringsAsFactors = FALSE)
report$TP <- report$taxID %in% truth_set$species_taxid
report$depth <- nchar(sub("[^ ].*", "", report$taxName))/2
report <- report[!report$taxID %in% c(9606, 9605), ] ## Ignore human identifications
report <- report[order(report$kmers, decreasing = TRUE),]
rs <- subset(report, rank=="species")
a <- plot_reads_vs_kmers_gg(rs, min.reads = 10, add.f1 = TRUE, dens.x = TRUE, dens.y = TRUE, show_thres=TRUE)
ggdraw(a)
ggsave("figures/sra_sampled_dataset-1mio-reads_vs_kmers.pdf", a, width=5, height=5, dpi=600, scale=1.2)


a <- plot_reads_vs_kmers_gg(rs, min.reads = 10, add.f1 = TRUE)
ggdraw(a + geom_vline(xintercept=log10(1000)))

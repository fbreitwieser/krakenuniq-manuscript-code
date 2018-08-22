#'--- 
#'title: "McIntyre et al. test dataDBS analysis"
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

source("functions.R")
options(width=10000) 

RANKS=c("genus", "species")
DBS <- c("orig", "std", "nt")

## Read truth for species and genus levels
truth_sets <- slapply(RANKS, function(rank_) read_truth_sets(sprintf("truth_sets/fixed1/%s", rank_), ext = "_TRUTH.txt"))
truth_sets_fixed <- slapply(RANKS, function(rank_) {
                              tf2 = read_truth_sets(sprintf("truth_sets/fixed2/%s", rank_), ext = "_TRUTH.txt");
                              list(orig=truth_sets[[rank_]],
                                   std=tf2,
                                   nt=tf2)
})

# Change in Carma_eval_carma (i.e. move of [Haemophilus] parasuis to genus Glaesserella) only applies to nt database
truth_sets_fixed[["genus"]][["std"]][["Carma_eval_carma"]] <- truth_sets[["genus"]][["Carma_eval_carma"]]

## Read KrakenUniq results for all databases
krr <- slapply(RANKS, function(rank_) slapply(DBS, function(db_)  
         read_mason_reports(report_dir = sprintf("krakenuniq_results/%s/report",db_), 
                            rank_ = rank_, truth_set_ = truth_sets_fixed[[rank_]][[db_]])))

## database std does has an older taxonomy H. parasuis has been moved to Glaesserella 

message("Got truth mappings for ",sum(names(krr[[1]][[1]]) %in% names(truth_sets_fixed[[1]][[1]])), " datasets.")
truth_sets_fixed <- slapply(RANKS, function(rank_) add_detected(truth_sets_fixed[[rank_]], krr[[rank_]]))

report_names <- names(krr[[1]][[1]])
dataset_description <- read_dataset_description("mcintyre-mason-datasets.csv", report_names)
biological_sets <- report_names[dataset_description[report_names, "Data.Type"] == "Biological"]

## Calculate the F1 scores for KrakenUniq on all databases
max.f1.scores <- slapply(RANKS, function(rank_) slapply(DBS, function(db_)
  get_max_f1_scores(krr[[rank_]][[db_]], truth_sets_fixed[[rank_]][[db_]])))

max.recall <- slapply(RANKS, function(rank_) slapply(DBS, function(db_)
  get_max_recall(krr[[rank_]][[db_]], truth_sets_fixed[[rank_]][[db_]], max.fdr = 0.05)))
#ge <- function(x,a,b) ifelse(x[,a]==x[,b], "eq", ifelse(x[,a]>x[,b], a, b))
#lapply(max.recall, function(x) sapply(x, ge, a="reads", b="kmers")) 

dataset_description <- dataset_description[report_names, ]
report_names_bio <- dataset_description$Full_Name[dataset_description$Data.Type=="Biological"]
report_names_sim <- dataset_description$Full_Name[dataset_description$Data.Type=="Simulated"]

## Extend dataset description to contain F1 scores of KrakenUniq
#biological_sets <- grep("ABRF|")
#dataset_description <- set_max_f1_scores(dataset_description, max.f1.scores)

#slapply(RANKS, function(rank_) slapply(DBS, function(db_) slapply(names(krr[[rank_]][[db_]]), 
#  function(n) subset(truth_sets_fixed[[rank_]][[db_]][[n]], !taxID %in% krr[[rank_]][[db_]][[n]]$taxID, "name",drop=T)) %>% unlist %>% table %>% sort))

programs <- list.files("tool_results_reformatted/species", pattern = "ABRF_MGRG_10ng_[A-Z]") %>%
    sub("ABRF_MGRG_10ng_", "", .) %>%
    sub(".txt$", "", .) %>% grep("DiamondMeganKrakenFiltered", ., value = T, invert=TRUE)

## Note: Truth is not 1 in all cases due to changes in the truth sets
( bio_program_f1s_recalls <- f_tool_results_f1_recall(report_names_bio) )
( sim_program_f1s_recalls <- f_tool_results_f1_recall(report_names_sim) )

write.csv(bio_program_f1s_recalls %>% arrange.matrix("avg", decreasing=T), "tables/test_datasets-f1-recall-for-all-tools-bio.csv")
write.csv(sim_program_f1s_recalls %>% arrange.matrix("avg", decreasing=T), "tables/test_datasets-f1-recall-for-all-tools-sim.csv")
( program_f1s_recalls <- f_tool_results_f1_recall(report_names) )
write.csv(program_f1s_recalls %>% arrange.matrix("avg", decreasing=T), "tables/test_datasets-f1-recall-for-all-tools.csv")

f1rc_bs <- cbind(bio_program_f1s_recalls %>% add_to_colname("bio"), sim_program_f1s_recalls %>% add_to_colname("sim"))
f1rc_bs <- cbind(f1rc_bs, avg = rowMeans(f1rc_bs[,c("bio.avg", "sim.avg")]))
write.csv(f1rc_bs %>% arrange.matrix("avg", decreasing=T), "tables/test_datasets-f1-recall-for-all-tools-bioandsim.csv")

df_f1rc <- rbind(summarize_stat(max.recall, "Recall"), summarize_stat(max.f1.scores, "F1"))
df_f1rc$Statistic <- factor(df_f1rc$Statistic, levels=c("Recall", "F1"))
(krakenuniq_stat <- (df_f1rc %>% 
    dplyr::select(Data.Type, Rank, Statistic, dplyr::everything())) %>%
    dplyr::arrange(Data.Type, Rank, Statistic))
write.csv(krakenuniq_stat, "tables/test_datasets-f1-recall-for-krakenuniq.csv")

pdf(sprintf("figures/test-datasets-krakenuniq-kmers-reads.pdf"),
    width=6,height=6, pointsize = 8,
    title="Reads vs k-mers on test datasets of McIntyre et al.")
for (rn in report_names) {
  message(rn)
  df_all <- data.frame()
  df_cnt <- data.frame()
  for (rank_ in RANKS) {
    for (db_ in DBS) {
      ts <- truth_sets_fixed[[rank_]][[db_]][[rn]]
      rs <- krr[[rank_]][[db_]][[rn]]
      rs$rank <- rank_
      rs$db <- db_
      df_cnt <- rbind(df_cnt, data.frame(rank=rank_, db=db_, 
                                         text=sprintf("detected %s/%s\nreads f1: %.2f\nkmers f1: %.2f", 
                                                      sum(ts$taxID %in% rs$taxID), 
                                                      nrow(ts),
                                                      max.f1.scores[[rank_]][[db_]][rn,"reads"],
                                                      max.f1.scores[[rank_]][[db_]][rn,"kmers"])))
      df_all <- rbind(df_all, rs)
    }
  }
  print(plot_reads_vs_kmers_gg(df_all, title=sprintf("%s (%s)",rn, dataset_description[rn,"Data.Type"])) + 
    facet_grid(rank~db) + geom_text(data=df_cnt, aes(label=text), x=Inf, y=-Inf, hjust=1, vjust=-0.2), col="gray")
}
dev.off()

## Print suppl figure 2: HMP_even_illum and Huttenhower_HC1
df_all <- data.frame()
df_cnt <- data.frame()
db_ = "nt"
for (rank_ in RANKS) {
  for (rn in c("HMP_even_illum_SRR172902", "Huttenhower_HC1")) {
      ts <- truth_sets_fixed[[rank_]][[db_]][[rn]]
      rs <- krr[[rank_]][[db_]][[rn]]
      rs$rank <- rank_
      rs$db <- db_
      df_cnt <- data.frame(rank=rank_, db=db_, rn=rn,
                           text=sprintf("detected %s/%s\nreads f1: %.2f\nkmers f1: %.2f", 
                                        sum(ts$taxID %in% rs$taxID), 
                                        nrow(ts),
                                        max.f1.scores[[rank_]][[db_]][rn,"reads"],
                                        max.f1.scores[[rank_]][[db_]][rn,"kmers"])) %>% rbind(df_cnt, .)
      df_all <- rbind(df_all, rs %>% mutate(rn = rn))
   }
}

print(plot_reads_vs_kmers_gg(df_all) + 
      facet_wrap(rank~rn, scales="free") + geom_text(data=df_cnt, aes(label=text), x=Inf, y=-Inf, hjust=1, vjust=-0.2), col="gray")
ggsave("figures/scatter-hmp-vs-hc1.pdf", width=5, height=5, dpi = 600, scale=1.2)


## Plots for individual datasets
f1ss <- max.f1.scores[["species"]][["std"]]
krakhll_f1s <- data.frame(Rank="Species", f1ss,Data.Type = dataset_description[rownames(f1ss), "Data.Type"])
ggplot(krakhll_f1s) + 
  geom_jitter(aes(x=reads,y=kmers, col=Data.Type), pch=21, size=4, alpha=.8, width=.001, height=.001) + 
  geom_abline(slope = 1, intercept=0, alpha=.2) + 
  #scale_x_continuous(limits = c(0.83,1)) + scale_y_continuous(limits = c(0.83,1)) +
  scale_color_hue("") +
  theme(legend.position = c(0.5,.2))
ggsave(filename = "figures/test_datasets-datasets-kmer-vs-reads-std-species.pdf", width = 3, height=3, scale = 1.2)

## Get logs (Suppl table 3)

zlogs <- read_mason_logs("kraken_results/std/log", "krakenuniq_results/std/log" )
write.csv(parse_log_results(dataset_description, zlogs), "tables/test_datasets-desc-plus-timings.csv")

#########################################################


req <- c("ggplot2", "cowplot", "magrittr", "lubridate", "plyr")
for (rr in req) {
    if (!require(rr, character.only=TRUE, quietly=TRUE)) {
        install.packages(rr)
        require(rr, character.only=TRUE, quietly=TRUE)
    }
}


fix_names <- function(nn, ext = "") {
  nn <- sub(paste0(".fa?s?t?[aq]",ext,"$"), "", basename(nn))
  
  ## Change names to match with dataset description
  sel_mgrg <- startsWith(nn,"ABRF_MGRG")
  nn[sel_mgrg] <-  sub("^([^_]*_[^_]*_[^_]*)_.*", "\\1", nn[sel_mgrg])
  nn[startsWith(nn,"BioPool_BioPool")] <- "QiagenFX_Assay_BioPool"
  nn <- sub("MGRG_nanopore", "MGRGnanopore", nn)
  nn[startsWith(nn,"P00")] <- paste0("NYCSubway","_",nn[startsWith(nn,"P00")])
  
  nn
}

read_mason_reports <- function(report_dir, rank_, truth_set_, main_dir = "krakenuniq_results") {
  
  kr <- list.files(report_dir,pattern = "*.krakenuniq.report.tsv", full.names = TRUE)
  krr <- lapply(kr, function(x) {
    read.delim(x, strip.white=TRUE, stringsAsFactors=F) %>% subset(rank == rank_) %>% arrange(-reads)
  })
  
  names(krr) <- fix_names(basename(kr), ".krakenuniq.report.tsv")
  
  krr["MGRGnanopore_R6_2d_pass"] <- NULL ## Remove datasets w/o description
  krr["MGRGnanopore_R6_2d_fail"] <- NULL ## Remove datasets w/o description
  krr["UnAmbiguouslyMapped_ds.frankengenome"] <- NULL ## Remove datasets w/o description
  krr["UnAmbiguouslyMapped_ds.frankengenome.mix"] <- NULL ## Remove datasets w/o description
  krr["JGI_SRR033547"] <- NULL ## Exclude b/c only 100 reads
  krr["Mavromatis_simHC"] <- NULL ## Exclude dataset from 2007
  krr["Mavromatis_simLC"] <- NULL ## Exclude dataset from 2007
  krr["Mavromatis_simMC"] <- NULL ## Exclude dataset from 2007
  
  for (rn in intersect(names(krr), names(truth_set_))) {
      krr[[rn]] <- subset(krr[[rn]], rank == rank_)
      krr[[rn]]$TP <- krr[[rn]]$taxID %in% truth_set_[[rn]]$taxID
      krr[[rn]] <- subset(krr[[rn]], taxName != "uncultured bacterium")
      krr[[rn]] <- subset(krr[[rn]], taxName != "synthetic construct")
      krr[[rn]] <- subset(krr[[rn]], taxName != "Homo sapiens")
      krr[[rn]] <- subset(krr[[rn]], taxName != "Homo")
      krr[[rn]] <- subset(krr[[rn]], taxName != "Mus musculus")
      #krr[[rn]] <- subset(krr[[rn]], taxName != "Pseudomonas antarctica")
      
      #truth_set_[[rn]]$Detected <- truth_set_[[rn]]$taxID %in% krr[[rn]]$taxID
      #dataset_description[rn, "NDetected"] <- sum(krr[[rn]]$TP, na.rm=T)
      #dataset_description[rn, "No..of.Genomes"] <- nrow(truth_sets[[rn]])
  }
  sel <- sapply(krr, function(x) "TP" %in% colnames(x))
  krr <- krr[sel] #Only keep datasets with truth set
  
  krr
}

tts <- c(s=1, m=60, h=60*60, d=60*60*60)

time_to_sec <- function(s) {
    sapply(strsplit(s, ":"), function(ss) {
        tt <- rev(as.numeric(ss))
        stopifnot(all(!is.na(tt)))
        sum(tt * tts[seq_along(tt)])
    })
}

read_mason_logs <- function(krak_log_dir, krakhll_log_dir) {
  krak_logs <- list.files(path=krak_log_dir, pattern = "*.kraken.log", full.names = T)
  res <- do.call(rbind,lapply(krak_logs, function(krak_log_f) {
    krak_log <- strsplit(readLines(krak_log_f),"[ ()]")
    krak_tlog <- strsplit(readLines(sub("kraken.log$","kraken.tlog", krak_log_f)),": ")
    
    krakhll_log_f<- sub("kraken.log$","krakenuniq.log", krak_log_f) %>% sub(krak_log_dir, krakhll_log_dir, ., fixed=T)
    krakhll_tlog_f<- sub("kraken.log$","krakenuniq.tlog", krak_log_f) %>% sub(krak_log_dir, krakhll_log_dir, ., fixed=T)
    krakhll_log <- strsplit(readLines(krakhll_log_f),"[ ()]")
    krakhll_tlog <- strsplit(readLines(krakhll_tlog_f),": ")
    
    krak_rep_tlog <- strsplit(readLines(sub("kraken.log$","kraken-report.tlog", krak_log_f)),": ")
    
    data.frame(n_seq = as.numeric(krak_log[[1]][1]),
               krak_mbp = as.numeric(krak_log[[1]][13]),
               krakhll_mbp = as.numeric(krakhll_log[[5]][13]),
               diff_mbp = NA,
               krak_walltime = krak_tlog[[5]][[2]],
               krak_rep_walltime = krak_rep_tlog[[5]][2],
               krakhll_walltime = krakhll_tlog[[5]][[2]],
               diff_walltime = NA,
               krak_mrss_gb = as.numeric(krak_tlog[[10]][[2]]) / 1024 / 1024,
               krakhll_mrss_gb = as.numeric(krakhll_tlog[[10]][[2]]) / 1024 / 1024,
               diff_mrss_gb = NA,
               stringsAsFactors = F)
  }))
  res$diff_mbp <- res$krakhll_mbp / res$krak_mbp - 1
  res$diff_walltime <- time_to_sec(res$krakhll_walltime) / (time_to_sec(res$krak_walltime) + time_to_sec(res$krak_rep_walltime)) - 1
  res$diff_mrss_gb <- res$krakhll_mrss_gb / res$krak_mrss_gb
  rownames(res) <- fix_names(krak_logs, ".kraken.log")
  res
}

read_dataset_description <- function(fname = "mcintyre-mason-datasets.csv", report_names = NULL) {
  dataset_description <- read.csv(fname, stringsAsFactors = F)
  dataset_description$Read.Length <- round(dataset_description$Read.Length, 1)
  
  # Make name compatible with the report/fastq filenames
  dataset_description$Name <- sub(".fa?s?t?[aq].gz$", "", dataset_description$Dataset.File)
  sel_mgrg <- startsWith(dataset_description$Name,"MGRG")
  dataset_description$Name[sel_mgrg] <-  sub("^([^_]*_[^_]*)_.*", "\\1", dataset_description$Name[sel_mgrg])
  dataset_description$Name[startsWith(dataset_description$Name,"QiagenFX_Assay_BioPool")] <- "Assay_BioPool"
  dataset_description$Name <- sub("-R1/2$", "-paired", dataset_description$Name)
  dataset_description$Name <- sub("^SL", "NegControl_SL", dataset_description$Name)
  dataset_description <- dataset_description[dataset_description$Name != "", ]
  
  if(!is.null(report_names)) {
    if(!all(sub("^[^_]*_","", report_names) %in% dataset_description$Name)) {
      print(report_names[!sub("^[^_]*_","", report_names) %in% dataset_description$Name])
      print(sort(setdiff(sub("^[^_]*_","", report_names), dataset_description$Name)))
    #sort(setdiff(dataset_description$Name, sub("^[^_]*_","", report_names)))
      stop("Didn't find description for all reports")
    }
    message("Got dataset descriptions for all datasets.")
    
    dataset_description$Full_Name <- sapply(dataset_description$Name, function(n) report_names[sub("^[^_]*_","",report_names) == n]  %>% ifelse(length(.) == 0, n, .))
    rownames(dataset_description) <- dataset_description$Full_Name
    dataset_description$Project <- sapply(dataset_description$Name, function(n) sub("_.*","",report_names)[sub("^[^_]*_","",report_names) == n]  %>% ifelse(length(.) == 0, NA, .))
    
  }
  dataset_description
}

get_tool_results_f1_recall <- function(res_dir, 
                                       truth_sets,
                                       programs,
                                       max.fdr=0.05) {
  
  stopifnot(dir.exists(res_dir))
  
  datasets <- list.files(res_dir, pattern = "TRUTH") %>%
    sub(paste0("_TRUTH.txt"), "", .)
  
  datasets <- datasets[fix_names(datasets) %in% names(truth_sets)]
    slapply(programs, function(program) {
        res <- t(sapply(datasets, function(rn) {
        fname <- sprintf("%s/%s_%s.txt", res_dir, rn, program)
        program_res <- read.delim(fname, header=F, col.names=c("taxid","n_reads","perc","rank","name"), stringsAsFactors = FALSE)
        if (length(program_res) == 0 || nrow(program_res) == 0)
          return(c(f1_score=0,recall=0))
        
        program_res$TP = program_res$taxid %in% truth_sets[[fix_names(rn)]]$taxID
        TPs <- program_res[order(program_res[,"n_reads"], decreasing = T), "TP"]
        stats <- get_f1_and_recall(TPs, max.fdr, nrow(truth_sets[[fix_names(rn)]]))
        c(f1_score = stats['max_f1'], recall = stats["max_recall"])
      }))
      rownames(res) <- fix_names(datasets)
      res
    })
}

get_f1_and_recall <- function(rsi, max_fdr, total_p) {
  sensitivity = sapply(seq_along(rsi), function(j) sum(rsi[1:j])) / total_p
  specificity = sapply(seq_along(rsi), function(j) sum(rsi[1:j])/j)
  fdr <- sapply(seq_along(rsi), function(j) sum(!rsi[1:j])/j)
  recall <- sapply(seq_along(rsi), function(j) sum(rsi[1:j])) / total_p
  if (any(fdr < max_fdr, na.rm=T)) {
    which.max_recall <- which.max(recall[fdr < max_fdr])[1]
    max_recall <- recall[which.max_recall]
  } else {
    which.max_recall <- -1
    max_recall <- 0
  }
  
  f1_scores <- 2*sensitivity*specificity/ (sensitivity+specificity)
  which.max_f1 <- which.max(f1_scores)[1]
  max_f1 <- f1_scores[which.max_f1]
  c(max_f1 = max_f1, max_recall= max_recall, which.max_f1 = which.max_f1, which.max_recall = which.max_recall)
}

summarize_stat <- function(f1rc, name) {
  res <- ddply(dataset_description, "Data.Type", function(dd) {
      lapply(f1rc, function(x) { 
      res = sapply(x, function(y) apply(y[dd$Full_Name,], 2, mean))
      rbind(res, diff=res[2,]/res[1,]-1) %>% flatten
      }) %>% rbindc("Rank")
  })
  res$Statistic <- name
  res
}

slapply <- function(n, ...) setNames(lapply(n, ...), n) 
cbindl <- function(x) { 
  do.call(cbind, lapply(names(x), function(y) {
    if (is.null(dim(x[[y]]))) {
      cn <- names(x[[y]])
      x[[y]] <- matrix(x[[y]], nrow=1)
    } else {
      cn <- colnames(x[[y]])
    }
    colnames(x[[y]]) <- sprintf("%s.%s", y, cn)
    x[[y]]
  }))
}

## Order a matrix by a column
arrange.matrix <- function(m, col, ...)
  m[order(m[,col], ...),]

## Prefix colnames
add_to_colname <- function(m, prefix) {
    colnames(m) <- sprintf("%s.%s", prefix, colnames(m))
    m
}


flatten <- function(x) { 
  cn <- colnames(x)
  cr <- rownames(x)
  setNames(as.numeric(x), expand.grid(cr, cn) %>% apply(1, paste, collapse="."))
}


rbindl <- function(x) { 
  if (is.null(dim(x[[1]]))) {
    res <- do.call(rbind, x)
    rownames(res) <- names(x)
    return(res)
  } else {
  do.call(rbind, lapply(names(x), function(y) {
    rownames(x[[y]]) <- sprintf("%s.%s", y, rownames(x[[y]]))
    x[[y]]
  }))
  }
}

rbindc <- function(x, n, r=NULL) {
  if (is.null(r)) {
    res <- do.call(rbind, lapply(x, function(y) {
      if (is.null(dim(y))) {
        cn <- names(y)
        y <- data.frame(matrix(y, nrow=1))
        colnames(y) <- cn
      } else {
        y <- data.frame(y)
      }
      y
    }))
    res[[n]] <- names(x)
    return(res)
  } else {
    do.call(rbind, lapply(names(x), function(y) {
      if (is.null(dim(x[[y]]))) {
        x[[y]] <- data.frame(val=x[[y]]) 
      } else {
        x[[y]] <- data.frame(x[[y]])
      }
      x[[y]][[n]] <- y
      x[[y]][[r]] <- rownames(x[[y]])
      rownames(x[[y]]) <- NULL
      x[[y]]
    }))
  }
}


read_truth_sets <- function(truth_dir, ext = "_TRUTH.txt") {
  pattern <- sprintf("*%s", ext)
  truth <- lapply(list.files(path=truth_dir, full.names = T, pattern = pattern), function(f) {
    res <- read.delim(f, header=F, col.names = c("taxID","ab1","ab2","rank","name"), stringsAsFactors = FALSE)
    res[order(res$name),]
  })
  names(truth) <- 
    list.files(path=truth_dir, pattern = pattern) %>% sub(ext, "", .)
  names(truth)[names(truth) == "BioPool_BioPool"] <- "QiagenFX_Assay_BioPool"
  
  truth
}

log10.axis <- function(side, at, ...) {
  at.minor <- log10(outer(1:9, 10^(min(at):max(at))))
  lab <- sapply(at, function(i) as.expression(bquote(10^ .(i))))
  axis(side=side, at=at.minor, labels=NA, tcl=par("tcl")*0.5, ...)
  axis(side=side, at=at, labels=lab, ...)
}

get_max_f1_scores <- function(krr, truth_sets) {
  fpx <- lapply(names(krr), function(rn) {
    rs <- subset(krr[[rn]], reads >= 10)
    res <- data.frame(j=seq_len(nrow(rs)))
    total_taxa <- nrow(truth_sets[[rn]])
    
    for (metric in c("reads","kmers")) {
      rsi <- rs[order(rs[,metric], decreasing = T), "TP"]
      res[,sprintf("sensitivity.%s",metric)] <- sapply(res$j, function(j) sum(rsi[1:j])/total_taxa)
      res[,sprintf("specificity.%s",metric)] <- sapply(res$j, function(j) sum(rsi[1:j])/j)
      res[,sprintf("f1.%s",metric)] <- 2*res[,sprintf("sensitivity.%s",metric)]*res[,sprintf("specificity.%s",metric)]/
        (res[,sprintf("sensitivity.%s",metric)]+res[,sprintf("specificity.%s",metric)])
    }
    res
  })
  
  max.f1.scores <- t(sapply(fpx,function(j) {
    c(reads = max(j$f1.reads,na.rm=T),
      kmers = max(j$f1.kmers,na.rm=T))
  }))
  rownames(max.f1.scores) <- names(krr)
  max.f1.scores  
}

get_max_recall <- function(krr, truth_sets, max.fdr = 0.05) {
  t(sapply(names(krr), function(rn) {
    rs <- subset(krr[[rn]])
    res <- data.frame(j=seq_len(nrow(rs)))
    total_taxa <- nrow(truth_sets[[rn]])
    
    sapply(c("reads","kmers"), function(metric) {
      rsi <- rs[order(rs[,metric], decreasing = T), "TP"]
      fdr <- sapply(seq_along(rsi), function(j) sum(!rsi[1:j])/j)
      recall <- sapply(seq_along(rsi), function(j) sum(rsi[1:j])/total_taxa)
      if (!any(fdr < max.fdr)) {
        return(NA)
      }
      #round(max(recall[fdr < 0.05]),2)
      max(recall[fdr < max.fdr])
    })
  }))
}

set_max_f1_scores <- function(dataset_description, max.f1.scores_) {
  for (db_ in names(max.f1.scores_)) {
    for (rank_ in names(max.f1.scores_[[db_]])) {
      max.f1.scores <- max.f1.scores_[[db_]][[rank_]]
      f1_reads <- sprintf("F1.Reads.%s.%s", db_, rank_)
      f1_kmers <- sprintf("F1.kmers.%s.%s", db_, rank_)
      dataset_description[, f1_reads] <- NA
      dataset_description[, f1_kmers] <- NA
      for (rn in rownames(max.f1.scores)) {
        dataset_description[rn, f1_reads] = round(max.f1.scores[rn, "max.f1.reads"],3)
        dataset_description[rn, f1_kmers] = round(max.f1.scores[rn, "max.f1.kmers"],3)
      }
    }
  }
  #  print(list(All=table(dataset_description$Better_F1),
  #          Sim=table(dataset_description[dataset_description$Data.Type == "Simulated", "Better_F1"]),
  #          Bio=c(table(dataset_description[dataset_description$Data.Type == "Biological", "Better_F1"]))))
  
  dataset_description
}

plot_reads_vs_kmers_gg <- function(rss, dsd=NULL, add.title=T, title=NULL, desc.wrap=50, dens.x=FALSE, dens.y=FALSE, min.reads = 10,
                                   add.f1 = FALSE, show_thresh = FALSE,
                                   cols = c("Data.Type", "Description", "Publication.Source", "Total.No..of.Reads",
                                            "No..of.Genomes", "F1.Reads / Kmers", "Missing"), legend.position = "none") {
  if (!is.null(min.reads)) {
    rs <- subset(rss, reads >= min.reads)
  } else {
      rs <- rss
   }
  scatter <- ggplot(rs, aes(x=log10(reads), y=log10(kmers))) + 
    geom_point(aes(col=TP), size = 2, alpha=.5) + 
    geom_point(aes(col=TP), size = 2, pch=1, alpha=.5) + 
    scale_x_continuous("Reads", breaks=seq(0,10), labels=scales::math_format(10^.x)) + 
    scale_y_continuous("Unique k-mers", breaks=seq(0,10), labels=scales::math_format(10^.x)) + 
    annotation_logticks(sides="lb", scaled = T, size = 0.4, alpha = 0.5) + cowplot::background_grid() +
    theme(legend.position=legend.position) + 
    scale_color_manual("",values=c("TRUE"="#D55E00", "FALSE"="#000000"))
  #+ theme(plot.margin=unit(c(0,0,1,1),"points"))
  if (add.title && !is.null(dsd)) {
    scatter <- scatter + ggtitle(sprintf("%s: %s", dsd$Project, dsd$Name))
  }
  
  if (add.f1) {
    rrs <- rss[order(-rs$reads),]
    stats_r <- get_f1_and_recall(rrs$TP, max_fdr = 0.05, sum(rrs$TP))
    stats_r['f1_thresh'] <- rrs[stats_r['which.max_f1'],"reads"]
    stats_r['recall_thresh'] <- rrs[stats_r['which.max_recall'],"reads"]
    read_thresh <- rrs$reads[stats_r['which.max_recall']]
    rrs <- rss[order(-rs$kmers),]
    stats_k <- get_f1_and_recall(rrs$TP, max_fdr = 0.05, sum(rrs$TP))
    stats_k['f1_thresh'] <- rrs[stats_k['which.max_f1'],"kmers"]
    stats_k['recall_thresh'] <- rrs[stats_k['which.max_recall'],"kmers"]
    kmer_thresh <- rrs$kmers[stats_k['which.max_recall']]
    print(sum(rss$TP))

    res <- data.frame(max=c(stats_r[1:2], stats_k[1:2]), which=c(stats_r[3:4], stats_k[3:4]), thresh=c(stats_r[5:6], stats_k[5:6]))
    rownames(res) <- c("F1 reads", "Recall reads", "F1 k-mers", "Recall k-mers")
    print(res)
    
    scatter <- 
      scatter + geom_text(label=sprintf("F1 reads %.2f, kmers %.2f\nRecall reads %.2f, kmers %.2f",
                                        stats_r['max_f1'], stats_k['max_f1'],stats_r['max_recall'], stats_k['max_recall']), 
                          x=Inf, y=-Inf, hjust=1, vjust=-0.5, alpha = .8, col="black") #+
                                        #geom_rug(data=data.frame(x=log10(stats_r['f1_thresh']),y=log10(stats_k['f1_thres'])), aes(x=x, y=y), sides="b")

    if (show_thresh) {
        #scatter <- scatter +
        #    geom_vline(xintercept=stats_r['f1_thresh'] %>% log10, linetype=2, alpha=.8) +
        #    geom_vline(xintercept=stats_r['recall_thresh'] %>% log10, linetype=3, alpha=.8)
        scatter <- scatter +
            geom_hline(yintercept=stats_k['f1_thresh'] %>% log10, color="#D55E00", linetype=2, alpha=.8) +
            geom_hline(yintercept=stats_k['recall_thresh'] %>% log10, color="#D55E00", linetype=3, alpha=.8)
    }
  }
  
  if (!is.null(title)) {
    scatter <- scatter + ggtitle(title)
  }
  #marginal density of y - plot on the right
  if (isTRUE(dens.y)) {
    dens_plot_y <- scatter %>% axis_canvas("y", coord_flip = T) + 
      geom_density(aes(log10(kmers), fill=TP), alpha=.5, data=rs) + coord_flip() +
    scale_fill_manual("",values=c("TRUE"="#D55E00", "FALSE"="#000000"))
  }
  
  if (isTRUE(dens.x)) {
    dens_plot_x <- scatter %>% axis_canvas("x") + 
      geom_density(aes(log10(reads), fill=TP), alpha=.5, data=rs) +
    scale_fill_manual("",values=c("TRUE"="#D55E00", "FALSE"="#000000"))
    if (add.f1) {
      #dens_plot_x <- dens_plot_x + geom_vline(xintercept = log10(stats_r['f1_thresh']))
    }
    scatter <- insert_xaxis_grob(scatter, dens_plot_x, position="top")
  }
  
  if (isTRUE(dens.y)) 
    scatter <- insert_yaxis_grob(scatter, dens_plot_y, position="right")
  
  if (is.null(dsd) || nrow(dsd) == 0) {
    return(scatter)
  } else {
    dsd$No..of.Genomes = sprintf("%s/%s",dsd$NDetected,dsd$No..of.Genomes)
    dsd$`F1.Reads / Kmers` = sprintf("%.2f / %.2f", dsd$Max.F1.Reads, dsd$Max.F1.Kmers)
  dsd <- dsd[,cols[cols %in% colnames(dsd)]]
  dsd$Description <- paste(strwrap(dsd$Description, desc.wrap), collapse="\n")
  if ("Missing" %in% colnames(dsd)) {
    dsd$Missing <- paste(strwrap(dsd$Missing, desc.wrap), collapse="\n")
  }
  colnames(dsd) <- gsub(".", " ", colnames(dsd), fixed = T)
  colnames(dsd) <- gsub("  ", ". ", colnames(dsd), fixed = T)
  cowplot::plot_grid(scatter, 
                     tableGrob(t(dsd), cols=NULL, theme = ttheme_default(9)))
  }
}

plot_reads_vs_kmers_r <- function(rs, rn) {
  ## base R plots
  plot(log10(rs$reads), log10(rs$kmers), col=ifelse(rs$TP, "red", "black"), main=rn, xaxt="n", yaxt="n", 
       pch=ifelse(rs$taxID==9606,2, ifelse(rs$taxID==32630,3,1)),
       xlab="Reads", ylab="k-mers")
  log10.axis(1, at=seq(0,10,1))
  log10.axis(2, at=seq(0,10,1))
  legend("topleft", pch=1, col=c("red","black"), legend=c("TP","FP"))
  # identify(log10(rs$reads), log10(rs$kmers), rs$taxName)
  for (rr in c("kmers")) {
    dens_t <- density(log10(rs[[rr]])[rs$TP])
    dens_f <- density(log10(rs[[rr]])[!rs$TP])
    plot(dens_t, col="red", main = "",
         xlim=c(0,max(c(dens_t$x,dens_f$x))), 
         ylim=c(0, max(c(dens_t$y,dens_f$y))),
         xlab="k-mers", xaxt="n")
    lines(dens_f, col="black")
    log10.axis(1, at=seq(0,10,1))
    
  }
}



read_report <- function(fname, report_taxid, rank_) {
  rep <- read.delim(fname, stringsAsFactors = F, strip.white = T) %>% subset(rank == rank_) %>% arrange(-reads)
  within(rep, {
    TP <- taxID %in% report_taxid
    species <- ifelse(TP, "true taxon", ifelse(taxID==32630, "artificial seq", ifelse(taxID==9606, "human", "other bacteria")))
    readsp <- reads/sum(reads)
    kmersp <- kmers/sum(kmers)
  })
}

add_detected <- function(truth_sets_, krr_) {
  for (db_ in names(krr_)) {
    for (rn in names(krr_[[db_]])) {
      truth_sets_[[db_]][[rn]][["detected"]] <- truth_sets_[[db_]][[rn]]$taxID %in% krr_[[db_]][[rn]]$taxID
    }
  }
  truth_sets_
}

fix_truth_sets <- function(truth_sets, rank_) {

    # UnAmbiguouslyMapped_ds.nycsm: Enterobacter => Klebsiella
    # Huttenhower_LC1: Prosthecochloris => Chlorobium
    #                  Prosthecochloris vibrioformis => Chlorobium phaeovibrioides", stringsAsFactors = F)
    # Huttenhower_HC1, Huttenhower_LC2: Peptoclostridium => Clostridioides
    # Carma_eval_carma: Haemophilus => Glaesserella
    # UnAmbiguouslyMapped_ds.7: Propionibacterium => Cutibacterium
}

parse_log_results <- function(dataset_description, logs) {
  do.call(rbind,lapply(seq_along(rownames(dataset_description)), function(rn) {
    dsd <- dataset_description[rn, c("Data.Type","Description","Publication.Source","Total.No..of.Reads")]
    dsd[["kraken mbpm"]] <- logs[rn, "krak_mbp"]
    dsd[["krakenuniq mbpm"]] <- logs[rn, "krakhll_mbp"]
    dsd[["diff mbpm"]] <- logs[rn, "diff_mbp"] 
    dsd[["kraken walltime"]] <- logs[rn, "krak_walltime"]
    dsd[["kraken-report walltime"]] <- logs[rn, "krak_rep_walltime"]
    dsd[["krakenuniq walltime"]] <- logs[rn, "krakhll_walltime"]
    dsd[["diff walltime"]] <- logs[rn, "diff_walltime"]    
    dsd[["kraken mem [gb]"]] <- logs[rn, "krak_mrss_gb"]
    dsd[["krakenuniq mem [gb]"]] <- logs[rn, "krakhll_mrss_gb"]
    dsd[["diff mem [gb]"]] <- logs[rn, "diff_mrss_gb"]
    dsd
  }))
}

read_result_sets <- function(SET, truth_sets, rank_) {
  message(SET)
  postfix = switch(SET, nt = ".nt", comp=".comp", chrom=".chrom", extra=".extra", std=".std")
  nn <- function(n,ext) { paste0(n,postfix,ext) }
  
  #truth_sets[["BMI_bmi_reads"]] <- subset(truth_sets[["BMI_bmi_reads"]], name != "Homo sapiens")
  
  
  
  ## Was there before, now missing. E.,g. Shigella boydii and Pyrobaculum arsenaticum (suppressed in RefSeq)
  
  #dsd_set$No..of.Genomes <- as.numeric(dsd_set$No..of.Genomes)
  #dsd_set$Genomes <- " <= 10"
  #dsd_set$Genomes[dsd_set$No..of.Genomes>10] <- "< 100"
  #dsd_set$Genomes[dsd_set$No..of.Genomes>=100] <- ">= 100"
  
  
  for (rn in names(krr)) {
    dataset_description[rn, "Max F1 Reads"] = round(max.f1.scores[rn,"max.f1.reads"],3)
    dataset_description[rn, "Max F1 Kmers"] = round(max.f1.scores[rn,"max.f1.kmers"],3)
  }
  
  list(dataset_description=dataset_description,
       dataset_timings=dsd_set,
       truth_sets=truth_sets,
       reports=krr,
       program_f1s=dsd_program_f1s)
}


read_reports <- function(report_dir, rank_, no_duplicate_species = FALSE) {
  report_files <- list.files(report_dir, full.names = TRUE)
  if (length(report_files) == 0) {
    warning("No report files")
    return()
  }
  report_taxids <- as.numeric(sub("\\..*","", sub(".*taxid","",report_files)))
  if (isTRUE(no_duplicate_species)) {
    sel <- duplicated(report_taxids)
    report_files <- report_files[!sel]
    report_taxids <- report_taxids[!sel]
  }
  kra <- lapply(report_files, read.delim, stringsAsFactors = F)
  for (i in 1:length(kra)) {
    kra[[i]] <- subset(kra[[i]], rank == rank_)
    kra[[i]] <- within(kra[[i]], {
      taxName <- sub("^  *","", taxName)
      TP <- taxID %in% report_taxids[i]
      species <- ifelse(TP, "true taxon", ifelse(taxID==32630, "artificial seq", ifelse(taxID==9606, "human", "other bacteria")))
      readsp <- reads/sum(reads)
      kmersp <- kmers/sum(kmers)
    })
  }
  names(kra) <- sub(".report.tsv", "", sub(".*/","",report_files))
  kra <- kra[sapply(kra,nrow) > 5]
  kra
}

plot_report <- function(krai, main=NULL, id=FALSE, col=ifelse(krai$TP,"red","black"), 
                        f=plot, x="reads",y="kmers", log="xy", ma=FALSE,  ...) {
  if (is.character(f) && isTRUE(f == "density")) {
    if (ma) {
      dens_t <- density(log10(krai[krai$TP,"kmers"]/krai[krai$TP,"reads"]))
      dens_f <- density(log10(krai[!krai$TP,"kmers"]/krai[!krai$TP,"reads"]))
      plot(dens_t, col="red", main = ifelse(is.null(main), a, main),
           xlim=c(0,max(c(dens_t$x,dens_f$x))), ylim=c(0, max(c(dens_t$y,dens_f$y))))
      lines(dens_f, col="black") 
    } else {
      for (a in c(x,y)) {
        dens_t <- density(log10(krai[krai$TP,a]))
        dens_f <- density(log10(krai[!krai$TP,a]))
        plot(dens_t, col="red", main = ifelse(is.null(main), a, main),
             xlim=c(0,max(c(dens_t$x,dens_f$x))), ylim=c(0, max(c(dens_t$y,dens_f$y))))
        lines(dens_f, col="black") 
      }
    }
  } else {
    rat <- krai[[y]]/krai[[x]]
    if (ma) {
      plot(krai[[x]], rat, col = col, bty="n",
           pch = sapply(as.character(krai$taxID),switch,"9606"=3,"32630"=2,1),
           log=log, main = main, xlab=x, ylab=paste(y,"/",x), ...) 
      abline(h=median(rat[krai$TP]), col="red")
      abline(h=median(rat[!krai$TP]))
    } else {
      f(krai[[x]],krai[[y]],col = col, bty="n",
        pch = sapply(as.character(krai$taxID),switch,"9606"=3,"32630"=2,1),
        log=log, main = main , xlab=x, ylab=y, ...)
      abline(0, mean(rat[krai$TP]), col="red")
      abline(0, mean(rat[!krai$TP]))
    }
  }
  if (isTRUE(id)) {
    krai[identify(krai$reads, krai$kmers, krai$taxName),]  
  }
}

gplot_report <- function(krai, name, id=FALSE, col=ifelse(krai$TP,"red","black"), 
                         f=plot, x="reads",y="kmers", log="xy", ...) {
  if (is.character(f) && isTRUE(f == "density")) {
    dens_t <- density(log10(krai[krai$TP,y]))
    dens_f <- density(log10(krai[!krai$TP,y]))
    plot(dens_t, col="red", main = name,
         xlim=c(0,max(c(dens_t$x,dens_f$x))), ylim=c(0, max(c(dens_t$y,dens_f$y))))
    lines(dens_f, col="black") 
  } else {
    ggplot(krai) + ggplot(aes_text(x=x, y=y,col="TP"),  
                          #pch = sapply(as.character(krai$taxID),switch,"9606"=3,"32630"=2,1),
                          log=log, main = name , xlab=x, ylab=y, ...)
    abline(0,1)
  }
  if (isTRUE(id)) {
    krai[identify(krai$reads, krai$kmers, krai$taxName),]  
  }
}


gplot_kra <- function(krad, col1=F, wrap=T, dens=T, legend=F, limits=T, solid=F, smooth=F, alpha = .5, rug=FALSE, size=1) {
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  if (!is.logical(col1)) {
    rep_text <- subset(krad, taxName %in% col1 | species != "other bacteria")
    #for (cc in col1)
    #  krad$species[krad$taxName == cc] <- cc
  }
  g <- ggplot(subset(krad, reads>100),
              aes(x=log10(reads),y=log10(kmers)))
  if (!solid) g <- g + scale_shape_discrete(solid=F)
  if (!is.logical(col1)) {
    g <- g + geom_point(aes(col=species), alpha=alpha) + 
      geom_text_repel(aes(label=taxName),rep_text, 
                      # Add extra padding around each text label.
                      box.padding = 0.5,
                      # Add extra padding around each data point.
                      point.padding = .6,
                      # Color of the line segments.
                      segment.color = '#cccccc',
                      # Width of the line segments.
                      segment.size = 0.5,
                      # Draw an arrow from the label to the data point.
                      arrow = arrow(length = unit(0.01, 'npc'))#,
                      #nudge_y = ifelse(log10(rep_text$kmers) < 5, .3, 0)
      )
  } else if (isTRUE(col1)) {
    g <- g + geom_point(aes(col=kmers<10*reads&kmers<1000), alpha=alpha)
  } else {
    #g <- g + geom_point(aes(col=species, shape=kmers<reads&kmers<1000), alpha=alpha)
    g <- g + geom_point(aes(col=species), alpha=alpha, size=size) + 
      scale_color_manual("",values=c("true species"="#D55E00", "other bacteria"="#000000", "human"="#E69F00", "artificial seq"="#56B4E9"), 
                         breaks=c("true species", "other bacteria", "human", "artificial seq"))
  }
  if(rug) g <- g+geom_rug(aes(col=species))
  if (dens) g <- g + geom_density_2d(aes(group=species, col=species), alpha=.5)
  if (smooth) g <- g + geom_smooth(alpha=.25)
  if (wrap) g <- g + facet_wrap(~species, nrow=1)
  if (!legend) g <- g + theme(legend.position = "none")
  if (limits) g <- g + 
    scale_x_continuous(limits = c(1.8, 8), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 7.1), expand = c(0, 0)) 
  g+
    background_grid(major = "xy", minor = "none")
} 

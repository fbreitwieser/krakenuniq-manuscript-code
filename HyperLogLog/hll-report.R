#'---
#' title: HyperLogLog performance testing
#' author: Florian Breitwieser
#' date: March 2018
#'---

library(ggplot2)
library(plyr)
library(cowplot)

#' ## Results on random data / comparison of Flajolet, Heule and Ertl
mean_error <- ldply(c(10,12,14,16,18), function(P) {
  res <- read.delim(sprintf("random/hll-p%s.txt",P))
  res1 <- ldply(c("heule","flajolet","ertl"), function(M) {
    res2 <- as.data.frame(do.call(rbind,
                                  tapply(res[[paste0("rel_error_",M)]],
                                         factor(res$observed), quantile, 
                                         probs = c(0.005,0.025, 0.341, 0.5, 0.682, 0.975, 0.995))))
    
    colnames(res2) <- sub("%","",colnames(res2))
    colnames(res2) <- paste0("p",colnames(res2))
    res2$method <- M
    res2$cardinality <- as.numeric(rownames(res2))
    res2
  })
  res1$p = paste("precision",P)
  res1
})
mean_error$method <- factor(mean_error$method, 
                            levels = c("flajolet","heule", "ertl"), 
                            labels = c("Flajolet","Heule","Ertl++"))

## annotation log-ticks!! http://ggplot2.tidyverse.org/reference/annotation_logticks.html
head(mean_error)

ggplot(subset(mean_error, cardinality > 800 & cardinality < 10000 & method != "ertl" & p == "precision 14")) + 
  geom_line(aes(x=cardinality, y=abs(p50), color=method)) 

plot_hll <- function(x){ 
  ggplot(x, alpha=.25) + 
    geom_ribbon(aes(x=cardinality, ymin=p2.5, ymax=p97.5), color="yellow", fill="yellow", alpha=.25) +
    geom_ribbon(aes(x=cardinality, ymin=p34.1, ymax=p68.2), color="orange", fill="orange", alpha=.25) +
    geom_line(aes(x=cardinality, y=p50), alpha=.75) +
    annotation_logticks(sides="tb", size = .4) + xlab("Cardinality") + ylab("Relative error") +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      expand = c(0,0)) +
    scale_y_continuous(labels=scales::percent) +
    theme_bw()
}
plot_hll(subset(mean_error, cardinality>=100)) + facet_grid(p~method)
ggsave("figures/hll-method-comparison-ertl_heule_flajolet.pdf", width = 6, height = 6)

plot_hll(subset(mean_error, cardinality>=100 & method !="Flajolet"))+ facet_grid(p~method, scales="free")
ggsave("figures/hll-method-comparison-ertl_heule.pdf", width = 6, height = 6)

#' ## Results on random sampling of one million k-mers from Kraken database

## Read HLL results - skip headers, starting with 'precision', since the concatenate files have multiple
res <- read.delim("hll-on-1m-db-kmers.txt", 
                  colClasses = c("factor","numeric","numeric","numeric"), 
                  header=F, comment.char="p")
colnames(res) <- c("precision", "true_count", "estimate", "ertl_estimate")
res$rel_error <- (res$true_count-res$ertl_estimate)/res$true_count

#/*
#' Don't inclued scatterplot
ggplot(subset(res,true_count >= 10 & true_count < 100000000),
       aes(log10(true_count),rel_error,color=factor(precision), fill=factor(precision))) + 
  geom_jitter(width=.1) + geom_smooth()
#*/
mean_summary <- function(x, mean_f = mean, sd_f = function(x) sqrt(var(x)/length(x))) {  
  x <- na.omit(x)
  data.frame(y = mean_f(x), ymin = mean_f(x) - sd_f(x), ymax = mean_f(x) + sd_f(x))
}
quant_summary <- function(x, mean_f = mean, lower=0.025, upper=0.975) {  
  x <- na.omit(x)
  data.frame(y = mean_f(x), ymin = quantile(x, probs=lower), ymax = quantile(x, probs=upper))
}
 
mean_se <- mean_summary
mean_sd <- function(x) mean_summary(x, mean_f=mean, sd_f=sd)
median_mad <- function(x) mean_summary(x, mean_f=median, sd_f=mad)

## Calculate minor tick marks - currently not displayed
get_log10_minorticks <- function(x) {
  log10(do.call(c, lapply(seq(from=round(x[1]),to=round(x[2])), 
                    function(y) seq(from=10^y, to=10^(y+1), by=10^y))))
}

std_x_scale <- scale_x_continuous(labels = function(x) format(10^x, scientific = F, big.mark = ","), expand=c(0,0))
log10_x_scale <- scale_x_continuous(labels = scales::math_format(10^.x), expand=c(0,0))
std_y_scale <- scale_y_continuous(breaks=function(x) seq(from=round(x[1],2),to=round(x[2],2),by=.01))
percent_y_scale <- scale_y_continuous(breaks=function(x) seq(from=round(x[1],2),to=round(x[2],2),by=.01), labels=scales::percent)

g <- ggplot(subset(res,true_count >= 100 & true_count <= 10000000),
       aes(log10(true_count),rel_error,color=factor(precision),fill=factor(precision))) + #geom_jitter(width=.1,alpha=.2) + 
  cowplot::background_grid(major="xy", minor="none") +
  scale_color_hue("precision") + scale_fill_hue("precision") +
  theme(legend.position="right") + annotation_logticks(sides="b") +
  xlab("True cardinality") + ylab("Relative error")

#' ### 95% boundaries of relative error
g95 <- g + stat_summary(fun.data=quant_summary, geom="ribbon", alpha=0.25)
#gpos95 <- gpos + stat_summary(fun.data=quant_summary, geom="ribbon", alpha=0.25)
#g95 + std_x_scale + std_y_scale
#ggsave("figures/relative-error-95quant.std_xy.pdf", width = 5, height = 3, scale = 1)

#g95 + log10_x_scale + std_y_scale
#ggsave("figures/relative-error-95quant.log10x.std_y.pdf", width = 5, height = 3, scale = 1)

g95 + std_x_scale + percent_y_scale
ggsave("figures/relative-error-95quant.std_x.percent_y.pdf", width = 5, height = 3, scale = 1)

g95 + log10_x_scale + percent_y_scale
ggsave("figures/relative-error-95quant.log10_x.percent_y.pdf", width = 5, height = 3, scale = 1)

#' ### SD boundaries of relative error
gsd <- g + stat_summary(fun.data=mean_sd, geom="ribbon", alpha=0.25) 
#gpossd <- gpos + stat_summary(fun.data=mean_sd, geom="ribbon", alpha=0.25) 

#gsd + std_x_scale + std_y_scale
#ggsave("figures/relative-error-sd.std_xy.pdf", width = 5, height = 3, scale = 1)

#gsd + log10_x_scale + std_y_scale
#ggsave("figures/relative-error-sd.log10x.std_y.pdf", width = 5, height = 3, scale = 1)

gsd + std_x_scale + percent_y_scale
ggsave("figures/relative-error-sd.std_x.percent_y.pdf", width = 5, height = 3, scale = 1)
ggsave("figures/relative-error-sd.std_x.percent_y.svg", width = 5, height = 3, scale = 1)

gsd + log10_x_scale + percent_y_scale
ggsave("figures/relative-error-sd.log10_x.percent_y.pdf", width = 5, height = 3, scale = 1)

#gpossd + std_x_scale + percent_y_scale
#gpossd + log10_x_scale + percent_y_scale



rep <- pavian::read_reports(list.files("report", full.names = T))

head(rep[[1]])
merged <- pavian::merge_reports(rep, numeric_col = c("cladeReads","kmers"))

merged2 <- pavian::merge_reports2(rep, numeric_cols = c("cladeReads","kmers"))

head(merged)

sel <- merged$Name %in% c('Escherichia coli', 'Cutibacterium acnes', 'synthetic construct', 'Delftia', 'Enterobacteria phage phiX174 sensu lato')
merged[sel,]
library(xlsx)
write.xlsx(merged[sel, ], "common-contaminants.xlsx", showNA = FALSE)
str(merged2)

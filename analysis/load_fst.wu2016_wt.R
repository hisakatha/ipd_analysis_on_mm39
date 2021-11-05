library(data.table)
library(fst)
library(doParallel)

fai.col.names <- c("name", "length", "offset", "linebases", "linewidth")
reference_metadata <- fread("/glusterfs/hisakatha/ensembl_mm39/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.fai", col.names = fai.col.names)

chrs <- reference_metadata[, name]

cat("Start loading fst data\n")
wu2016_wt_data_list <- foreach(chr = chrs) %do% {
    read_fst(sprintf("load_and_save.wu2016_wt.%s.fst", chr), as.data.table = TRUE)
}
names(wu2016_wt_data_list) <- chrs
cat("Loaded wu2016_wt_data_list from load_and_save.wu2016_wt.*.fst\n")

library(data.table)
library(doParallel)
kinetics_deep_region_to_bed <- function (data_list, chrs, coverage_thres, out) {
    file.remove(out)
    foreach(chr = chrs) %do% {
        cat(sprintf("Start processing: %s\n", chr))
        fwrite(data_list[[chr]][coverage >= coverage_thres, .(refName, tpl - 1, tpl, "deep_kinetics", coverage, ifelse(strand == 0, "+", "-"))] , sep = "\t", file = out, col.names = FALSE, append = TRUE)
    }
}

coverage_thres <- 25

source("load_fst.wu2016_ko.R")
kinetics_deep_region_to_bed(wu2016_ko_data_list, chrs, coverage_thres, "deep_kinetics_region.wu2016_ko.bed")

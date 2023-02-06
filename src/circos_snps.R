library(data.table)
library(tidyverse)

mhv1_variants <- fread("MhV1_variants.tsv", skip=3713)
variant_sites <- mhv1_variants[,c(2)]
variant_sites$label <- paste("main")
variant_sites$end <- variant_sites$POS
variant_sites$color <- paste("color=black")
variant_sites_circos <- variant_sites[,c(2,1,3,4)]
write_tsv(variant_sites_circos, "output/circos_snps.txt", col_names = F)

library(data.table)
library(tidyverse)

het_sites <- fread("output/heterozygosity_multiallelic/filter.vcf", skip=3712)
##only contains variants that passed filtering

##sum number of variants per contig
no_sites <- het_sites %>% count(`#CHROM`)
##subset for HI-C or viral
hic_viral_table <- fread("data/hic_viral_contig_lengths.csv")
hic_viral_sites <- subset(no_sites, `#CHROM` %in% hic_viral_table$contig_id)

##length of contigs
variants_lengths <- merge(hic_viral_table, hic_viral_sites, by.x="contig_id", by.y="#CHROM")
variants_lengths$sites_per_kb <- (variants_lengths$n/variants_lengths$contig_length)*1000

##significant difference - with higher heterozygosity in viral contigs (4.29/kb vs 3.38/kb in Hi-C)
t.test(sites_per_kb ~ ID, data=variants_lengths)

summary <- aggregate(cbind(contig_length, n) ~ ID, data=variants_lengths, FUN=sum)
summary$sites_per_kb <- (summary$n/summary$contig_length)*1000

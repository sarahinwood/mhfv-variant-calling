library(data.table)
library(tidyverse)

##########
## Hi-C ##
##########

##read block info
hic_blocks <- fread("output/hapcut2/hic_block_headers.out")
hic_block_table <- hic_blocks[,c(3,5,7,9,11)]
setnames(hic_block_table, old=c("V3", "V5", "V7", "V9", "V11"), new=c("offset", "len","phased", "span", "fragments"))
##mean block length (bp)
mean(hic_block_table$span)
##len: SNV span of block
##SPAN: base pair span of block
hic_block_table$ID <- paste("hic")

##are blocks on all contigs
hic_contig_ids <- fread("output/hapcut2/hic_block_contig_ids.out", header=F)
unique(hic_contig_ids$V1)

###########
## Viral ##
###########

##read block info
viral_blocks <- fread("output/hapcut2/viral_block_headers.out")
viral_block_table <- viral_blocks[,c(3,5,7,9,11)]
setnames(viral_block_table, old=c("V3", "V5", "V7", "V9", "V11"), new=c("offset", "len","phased", "span", "fragments"))
##mean block length (bp)
mean(viral_block_table$span)
viral_block_table$ID <- paste("viral")

##are blocks on all contigs
viral_contig_ids <- fread("output/hapcut2/viral_block_contig_ids.out", header=F)
unique(viral_contig_ids$V1)
##haplotype blocks on all viral contigs - with viral a mean 1313 bp and hic 1083 bp in length

##############
## Analysis ##
##############

all_blocks <- full_join(hic_block_table, viral_block_table)
t.test(span ~ ID, data=all_blocks)

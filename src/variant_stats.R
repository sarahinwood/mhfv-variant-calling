library(data.table)

##https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants

##can't use VQSR for filtering as requires other data sets - use hard filtering instead
#Trinity already filtered for:
# FS < 30 - is the Phred-scaled probability that there is strand bias at that site, 0 = little to no bias, GATK rec.s filtering at least 60
# QD > 2.0 - quality by depth, quality score normalised by read depth, GATK rec.s remove those 2.0 or less
#Then I added
# DP > 10 - depth above 10
##other filters that could be used:
# SOR > 3 - measurement for strand bias, but it takes into account the ratios of reads that cover both alleles
# MQ > 40 - mapping quality, good = ~60
# ReadPosRankSum - measures the relative position of the reference vs the alternative allele within reads. E.g., SNPs towards the ends of reads may be more likely a sequencing error, < -8.0
# DP = sum of read depth over all samples, Variants that have a high mean depth across all samples are indicative of either repetitive regions or paralogs. Remove those with twice mean

##some columns have . rather than value - not sure why these still kept - shouldn't they fail?
variant_stats <- fread("output/gatk/Mh_MhV1_genotyped_stats.txt", header=F, na.strings=".")

setnames(variant_stats, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7"),
         new=c("FS", "SOR", "MQRankSum", "ReadPosRankSum", "QD", "MQ", "DP"))

##remove rows with NA - removes 425 rows
variant_stats_nona <- na.omit(variant_stats)

########
## FS ##
########

plot(density(variant_stats_nona$FS))
sum(variant_stats_nona$QD>30)

########
## QD ## trinity > 2, all above
########

plot(density(variant_stats_nona$QD))
sum(variant_stats_nona$QD<2)

########
## DP ## me > 10, still 1022309 below 10, not on GATK page though - better to use QD **
########

sum(variant_stats_nona$DP<10)

## Depth - remove lower than 10 and above 2*mean
plot(density(variant_stats_nona$DP))

##mean of 28.5
mean(variant_stats_nona$DP)
##increases to 49.6 when less than 10 removed
dp_10 <- subset(variant_stats_nona, variant_stats_nona$DP>10)
mean(dp_10$DP)
##would reduce to 142,256
a <- subset(dp_10, dp_10$DP>100)
length(a$DP)

#########
## SOR ## GATK recc. > 3 ***
#########

plot(density(variant_stats_nona$SOR))
sum(variant_stats_nona$SOR>3)
##only 9624 above 3

#########
## MQ ## GATK recc. > 40, all are
#########

plot(density(variant_stats_nona$MQ))
sum(variant_stats_nona$MQ>40)/length(variant_stats_nona$MQ)

###############
## MQRankSum ## GATK recc. between -&+ 12.5 - all pass
###############

plot(density(variant_stats_nona$MQRankSum))

####################
## ReadPosRankSum ## GATK recc. between -&+ 8.0 - all pass
####################

plot(density(variant_stats_nona$ReadPosRankSum))

#### could still re-filter for depth>10 and SOR, all other filters pass
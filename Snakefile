
#!/usr/bin/env python3

# containers 
minimap_container = 'docker://staphb/minimap2:2.24'
bcftools_container = 'docker://staphb/bcftools:1.16'
hapcut2_container = 'docker://dovetailg/hapcut2:1.0' # v1.3.1
gatk_container = 'docker://broadinstitute/gatk:4.3.0.0'
bwa_container = 'docker://staphb/bwa:0.7.17'

#########
# RULES #
#########

rule target:
    input:
        'output/bwa/sorted_marked_dups_NO_MhV1.bam'
        'output/gatk/Mh_MhV1_genotyped_stats.txt',
        'output/02_filtering/filtered_snps_stats.txt',
        'output/02_filtering/MhV1.vcf.gz',
        'output/02_filtering/Mh_HiC.vcf.gz'

rule Mh_VCF:
    input:
        'output/02_filtering/final_snp_variants.vcf.gz'
    output:
        'output/02_filtering/Mh_HiC.vcf.gz'
    log:
        'output/logs/Mh_HiC_VCF.log'
    singularity:
        bcftools_container
    shell:
        'bcftools view -t "^scaffold_1_MhV1" {input} -o {output} -Oz'

rule MhV1_VCF:
    input:
        'output/02_filtering/final_snp_variants.vcf.gz'
    output:
        'output/02_filtering/MhV1.vcf.gz'
    log:
        'output/logs/MhV1_VCF.log'
    singularity:
        bcftools_container
    shell:
        'bcftools view -t "scaffold_1_MhV1" {input} -o {output} -Oz'

# AC==0 removes all sites where no alternative alleles called for any samples - not variants
# AC==AN removes all sites where only alternative allele called - not true variant (assembly error)
# SnpGap removes variants close to indels - harder to call with certainty
# -m2 -M2 -v snps keeps only biallelic SNPs - minimum and maximum no. alleles is 2 and 10 as pooled sample, indels removed
rule filter_for_variants:
    input:
        'output/02_filtering/filtered_snps_pass.vcf.gz'
    output:
        'output/02_filtering/final_snp_variants.vcf.gz'
    log:
        'output/logs/filter_for_variants.log'
    singularity:
        bcftools_container
    shell:
        'bcftools filter -e "AC==0 || AC==AN" --SnpGap 10 {input} | bcftools view -q 0.1 -m2 -M10 -v snps -O z -o {output} 2> {log}'

rule filtered_variant_stats:
    input:
        'output/02_filtering/filtered_snps_pass.vcf.gz'
    output:
        'output/02_filtering/filtered_snps_stats.txt'
    log:
        'output/logs/filtered_variant_stats.log'
    singularity:
        bcftools_container
    shell:
        'bcftools query {input} -f "%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n" > {output} 2> {log}'

# then select only variants that passed
rule select_pass_merged:
    input:
        'output/02_filtering/filtered_snps.vcf.gz'
    output:
        'output/02_filtering/filtered_snps_pass.vcf.gz'
    singularity:
        bcftools_container
    shell:
        'bcftools view -f "PASS" {input} -o {output} -Oz'

# filter with GATK recc.s
rule filter_snps:
    input:
        'output/02_filtering/merged_snps_only.vcf.gz'
    output:
        'output/02_filtering/filtered_snps.vcf.gz'
    log:
        'output/logs/filter_snps.log'
    singularity:
        gatk_container
    shell:
        'gatk VariantFiltration '
        '-V {input} '
        '-filter "QD < 2.0" --filter-name "QD2" '
        '-filter "QUAL < 30.0" --filter-name "QUAL30" '
        '-filter "SOR > 3.0" --filter-name "SOR3" '
        '-filter "FS > 60.0" --filter-name "FS60" '
        '-filter "MQ < 40.0" --filter-name "MQ40" '
        '-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" '
        '-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" '
        '-O {output} '
        '2> {log}'

rule variant_stats:
    input:
        'output/02_filtering/merged_snps_only.vcf.gz'
    output:
        'output/02_filtering/variant_stats.txt'
    log:
        'output/logs/variant_stats.log'
    singularity:
        bcftools_container
    shell:
        'bcftools query {input} -f "%FS\t%SOR\t%MQRankSum\t%ReadPosRankSum\t%QD\t%MQ\t%DP\n" > {output} 2> {log}'

rule only_snps:
    input:
        'output/gatk/Mh_MhV1_genotyped.vcf.gz'
    output:
        'output/02_filtering/merged_snps_only.vcf.gz'
    log:
        'output/logs/only_snps.log'
    singularity:
        gatk_container
    shell:
        'gatk SelectVariants '
        '-V {input} '
        '-select-type SNP '
        '-O {output} '
        '2> {log}'

# gatk site says it can handle any ploidy/a mix of plodies without specifying
    # should be okay but consider this could yield different results than performing Mh and MhV1 separately
rule genotype_gvcfs_combined:
    input:
        genome = 'data/Mh_MhV1_genomes.fa',
        vcf = 'output/gatk/Mh_MhV1_combined.g.vcf.gz'
    output:
        'output/gatk/Mh_MhV1_genotyped.vcf.gz'
    log:
        'output/logs/genotype_gvcfs.log'
    singularity:
        gatk_container
    threads:
        40
    shell:
        ' gatk --java-options "-Xmx4g" GenotypeGVCFs '
        '-R {input.genome} '
        '-V {input.vcf} '
        '-O {output} '
        '2> {log}'

rule combine_gvcfs:
    input:
        genome = 'data/Mh_MhV1_genomes.fa',
        diploid = 'output/gatk/gatk_diploid.g.vcf.gz',
        viral = 'output/gatk/gatk_viral.g.vcf.gz'
    output:
        'output/gatk/Mh_MhV1_combined.g.vcf.gz'
    log:
        'output/logs/combine_gvcfs.log'
    singularity:
        gatk_container
    threads:
        40
    shell:
        'gatk CombineGVCFs '
        '-R {input.genome} '
        '--variant {input.diploid} '
        '--variant {input.viral} '
        '-O {output} '
        '2> {log}'

rule gatk_diploid:
    input:
        genome = 'data/Mh_MhV1_genomes.fa',
        genome_dict = 'data/Mh_MhV1_genomes.dict',
        genome_index = 'data/Mh_MhV1_genomes.fa.fai',
        bam = 'output/bwa/sorted_marked_dups_NO_MhV1.bam',
        bam_index = 'output/bwa/sorted_marked_dups_NO_MhV1.bam.bai'
    output:
        gvcf = 'output/gatk/gatk_diploid.g.vcf.gz'
    threads:
        40
    log:
        'output/logs/gatk_diploid.log'
    singularity:
        gatk_container
    shell:
        'gatk --java-options "-Xmx4g" HaplotypeCaller '
        '-R {input.genome} '
        '-I {input.bam} '
        '--ploidy 10 ' # ploidy is ploidy per indiv (2) * number samples in pool (5)
        '--native-pair-hmm-threads {threads} '
        '-ERC GVCF '
        '-O {output.gvcf} '
        '2> {log}'

rule gatk_viral:
    input:
        genome = 'data/Mh_MhV1_genomes.fa',
        genome_dict = 'data/Mh_MhV1_genomes.dict',
        genome_index = 'data/Mh_MhV1_genomes.fa.fai',
        bam = 'output/bwa/sorted_marked_dups_rgs.bam',
        bam_index = 'output/bwa/sorted_marked_dups_rgs.bam.bai'
    output:
        gvcf = 'output/gatk/gatk_viral.g.vcf.gz'
    threads:
        40
    log:
        'output/logs/gatk_viral.log'
    singularity:
        gatk_container
    shell:
        'gatk --java-options "-Xmx4g" HaplotypeCaller '
        '-R {input.genome} '
        '-I {input.bam} '
        '-L scaffold_1_MhV1 '
        '--ploidy 5 ' # ploidy is ploidy per indiv (1) * number samples in pool (5)
        '--native-pair-hmm-threads {threads} '
        '-ERC GVCF '
        '-O {output.gvcf} '
        '2> {log}'

# index bam files for haplotype caller
rule samtools_index_NO_MhV1:
    input:
        bam = 'output/bwa/sorted_marked_dups_NO_MhV1.bam'
    output:
        index = 'output/bwa/sorted_marked_dups_NO_MhV1.bam.bai'
    log:
        'output/logs/samtools_index_NO_MhV1.log'
    threads:
        20
    shell:
        'samtools index '
        '{input.bam} '
        '2> {log}'

# bam without MhV1 for diploid variant calling
rule remove_MhV1_bam:
    input:
        bam = 'output/bwa/sorted_marked_dups_rgs.bam',
        index = 'output/bwa/sorted_marked_dups_rgs.bam.bai'
    output:
        not_MhV1_contigs = temp('output/bwa/not_MhV1_contigs.txt'),
        bam = 'output/bwa/sorted_marked_dups_NO_MhV1.bam'
    log:
        'output/logs/remove_MhV1_bam.log'
    shell:
        'samtools idxstat {input.bam} | cut -f1 | grep -E -v scaffold_1_MhV1 > {output.not_MhV1_contigs} || exit 1 ; '
        'cat {output.not_MhV1_contigs} | tr "\n" " " | xargs samtools view -bh {input.bam} > {output.bam} 2> {log}'

rule samtools_index:
    input:
        bam = 'output/bwa/sorted_marked_dups_rgs.bam'
    output:
        index = 'output/bwa/sorted_marked_dups_rgs.bam.bai'
    log:
        'output/logs/samtools_index.log'
    threads:
        20
    shell:
        'samtools index '
        '{input.bam} '
        '2> {log}'

# add readgroups
rule gatk_readgroups:
    input:
        'output/bwa/sorted_marked_dups.bam'
    output:
        'output/bwa/sorted_marked_dups_rgs.bam'
    log:
        'output/logs/gatk_readgroups.log'
    singularity:
        gatk_container
    shell:
        'gatk --java-options "-Xmx8g" AddOrReplaceReadGroups I={input} O={output} RGID=genomelib1 RGLB=genomelib1 RGPU=genomelib1 RGSM=genomelib1 RGPL=ILLUMINA 2> {log}'

# sort sam
rule gatk_sortsam:
    input:
        'output/bwa/marked_dups.bam'
    output:
        'output/bwa/sorted_marked_dups.bam'
    log:
        'output/logs/gatk_sortsam.log'
    singularity:
        gatk_container
    shell:
        'gatk --java-options "-Xmx8g" SortSam '
        'I={input} '
        'O={output} '
        'SORT_ORDER=coordinate '
        '2> {log}'

# mark duplicate reads
rule gatk_mark_dups:
    input:
        bam = 'output/bwa/sorted.bam'
    output:
        md_bam = 'output/bwa/marked_dups.bam',
        metrics = 'output/bwa/mark_dups_metrics.txt'
    log:
        'output/logs/gatk_mark_dups.log'
    singularity:
        gatk_container
    threads:
        20
    shell:
        'gatk --java-options "-Xmx8g" MarkDuplicates '
        'I={input.bam} '
        'ASSUME_SORT_ORDER=coordinate '
        'O={output.md_bam} '
        'M={output.metrics} '
        '2> {log}'

# map reads
rule sam_to_bam:
    input:
        sam = 'output/bwa/bwa_mem.sam'
    output:
        sorted_bam = 'output/bwa/sorted.bam'
    threads:
        20
    log:
        'output/logs/sam_to_bam.log'
    shell:
        'samtools sort '
        '{input.sam} '
        '-o {output.sorted_bam} '
        '--threads {threads} '
        '2> {log}'

rule bwa_mem:
    input:
        index = 'output/bwa/index.bwt',
        r1 = 'data/illumina_reads_trimmed/Mh_trimr1.fq.gz',
        r2 = 'data/illumina_reads_trimmed/Mh_trimr2.fq.gz'
    output:
        sam = 'output/bwa/bwa_mem.sam'
    params:
        index_dir = 'output/bwa/index'
    threads:
        50
    log:
        'output/logs/bwa_mem.log'
    singularity:
        bwa_container
    shell:
        'bwa mem '
        '-t {threads} '
        '{params.index_dir} '
        '{input.r1} {input.r2} '
        '-p '
        '> {output.sam} '
        '2> {log}'

rule bwa_index:
    input:
        genome = 'data/Mh_MhV1_genomes.fa'
    output:
        index = 'output/bwa/index.bwt'
    params:
        outdir = 'output/bwa/index'
    threads:
        20
    log:
        'output/logs/bwa_index.log'
    singularity:
        bwa_container
    shell:
        'bwa index '
        '{input.genome} '
        '-p {params.outdir} '
        '2> {log} '

# gatk prep fasta for reference
rule gatk_seq_dict:
    input:
        'data/Mh_MhV1_genomes.fa'
    output:
        'data/Mh_MhV1_genomes.dict'
    shell:
        'gatk CreateSequenceDictionary -R {input}'

rule faidx:
    input:
        fasta = "data/Mh_MhV1_genomes.fa"
    output:
        fai = "data/Mh_MhV1_genomes.fa.fai"
    log:
        "output/logs/faidx.log"
    shell:
        "samtools faidx {input.fasta} -o {output.fai} 2> {log} "

#### old code ####

rule minimap2_illumina:
    input:
        genome = 'data/Mh_MhV1_genomes.fa',
        r1 = 'data/illumina_reads_trimmed/Mh_trimr1.fq.gz',
        r2 = 'data/illumina_reads_trimmed/Mh_trimr2.fq.gz'
    output:
        sam = temp('output/minimap2/minimap2.sam')
    threads:
        20
    singularity:
        minimap_container
    log:
        'output/logs/minimap2_illumina.log'
    shell:
        'minimap2 '
        '-ax sr '
        '{input.genome} '
        '{input.r1} {input.r2} '
        '-t {threads} '
        '> {output} '
        '2> {log}'

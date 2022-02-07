
#!/usr/bin/env python3

#########
# RULES #
#########

rule target:
    input:
        expand("output/hapcut2/{type}_block_contig_ids.out", type=["viral", "hic"])

##########################
## burke heterozygosity ##
##########################

rule  grep_block_headers:
    input:
        "output/hapcut2/{type}_blocks.out"
    output:
        headers = "output/hapcut2/{type}_block_headers.out",
        ids = "output/hapcut2/{type}_block_contig_ids.out"
    shell:
        'egrep "BLOCK" {input} > {output.headers} || exit 1 ; '
        'egrep -o "scaffold_[0-9]+" {input} > {output.ids} || exit 1 ; '
        'wait'

rule grep_hic_blocks:
    input:
        ids = "data/Mh_hic_contig_ids.txt",
        blocks = "output/hapcut2/haplotypes.out"
    output:
        "output/hapcut2/hic_blocks.out"
    shell:
        'egrep -B 1 -wf {input.ids} {input.blocks} > {output}'

rule grep_viral_blocks:
    input:
        ids = "data/Mh_DNA_virus_contigs_reduced.txt",
        blocks = "output/hapcut2/haplotypes.out"
    output:
        "output/hapcut2/viral_blocks.out"
    shell:
        'egrep -B 1 -wf {input.ids} {input.blocks} > {output}'

rule hapcut2:
    input:
        fragments = "output/hapcut2/extracthairs_fragments.out",
        vcf = "output/heterozygosity_multiallelic/filter.vcf"
    output:
        "output/hapcut2/haplotypes.out"
    log:
        "output/logs/hapcut2.log"
    shell:
        'bin/HapCUT2-1.3.3/build/HAPCUT2 '
        '--fragments {input.fragments} '
        '--VCF {input.vcf} '
        '--output {output} '
        '2> {log}'

rule extracthairs:
    input:
        bam = "/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-depth-gc/output/samtools/Mh/sorted.bam",
        vcf = "output/heterozygosity_multiallelic/filter.vcf"
    output:
        "output/hapcut2/extracthairs_fragments.out"
    log:
        "output/logs/extracthairs.log"
    shell:
        'bin/HapCUT2-1.3.3/build/extractHAIRS '
        '--bam {input.bam} '
        '--VCF {input.vcf} '
        '--out {output} '
        '2> {log}'

rule bcftools_filter:
    input:
        vcf = "output/heterozygosity_multiallelic/call_multiallelic.vcf.gz"
    output:
        vcf = "output/heterozygosity_multiallelic/filter.vcf"
    log:
        'output/logs/bcftools_filter.log'
    threads:
        20
    shell:
        'bcftools filter {input.vcf} -o {output.vcf} -O v -i"%QUAL>20 && AC/AN=0.5" --threads {threads} 2> {log}'

rule bcftools_call_multiallelic:
    input:
        vcf = "output/heterozygosity/mpileup.vcf.gz"
    output:
        vcf = "output/heterozygosity_multiallelic/call_multiallelic.vcf.gz"
    log:
        'output/logs/bcftools_call.log'
    threads:
        20
    shell:
        'bcftools call {input.vcf} -o {output.vcf} -O z -m --threads {threads} 2> {log}'

rule bcftools_mpileup:
    input:
        bam = "/Volumes/archive/deardenlab/sarahinwood/mh_projects/microctonus-depth-gc/output/samtools/Mh/sorted.bam",
        fasta = "data/Mh.hic_assembly.fa",
        fai = "data/Mh.hic_assembly.fa.fai"
    output:
        vcf = "output/heterozygosity/mpileup.vcf.gz"
    log:
        'output/logs/mpileup.log'
    shell:
        'bcftools mpileup {input.bam} -f {input.fasta} -d 50, -C 50 -o {output.vcf} -O z --threads {threads} 2> {log}'

rule faidx:
    input:
        fasta = "data/Mh.hic_assembly.fa"
    output:
        fai = "data/Mh.hic_assembly.fa.fai"
    log:
        "output/logs/faidx.log"
    shell:
        "samtools faidx {input.fasta} -o {output.fai} 2> {log} "

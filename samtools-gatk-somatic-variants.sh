# Identify somatic variants in a human leukemia sample
# using healthy T cells as germline control

# written by Jennifer Reid, 07-Mar-2021

#####################

cd /path/to/directory
# samtools v1.2, bcftools v1.2, gatk v4.2 installed and in $PATH

# download hg19 as reference genome
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/analysisSet/
# hg19.p13.plusMT.no_alt_analysis_set.fa.gz  2020-03-09 10:21  823M  
gunzip hg19.fa.gz
more hg19.fa | grep ">" | wc -l
#	85
#	85 accounts for many "chrUn" lines

# hg19 chr labels are prefixed with "chr", but neither bam nor  gnomad chr's are (labelled "1", etc)
more hg19.fa | grep ">"
#	>chr1  AC:CM000663.1  gi:224384768  LN:249250621  rl:Chromosome  	
#	M5:d949d5ece55a7f11abb62abfe947815d  AS:GRCh37.p13
# remove "chr" from hg19.fa
awk '{gsub(/^chr/,""); print}' hg19.fa > hg19_no_chr.fa
more hg19_no_chr.fa | grep ">"
# 	>1  AC:CM000663.1  gi:224384768  LN:249250621  rl:Chromosome  	
#	M5:d949d5ece55a7f11abb62abfe947815d  AS:GRCh37.p13

# create hg19 index and dict files
samtools faidx hg19_no_chr.fa
gatk CreateSequenceDictionary -R hg19_no_chr.fa

#####################

# leukemia sample:
# AML patient FACS-purified blasts --> leuk.bam
samtools view leuk.bam | more
samtools view leuk.bam | tail -10
# header
samtools view -H leuk.bam
# aligned using hg19 ref genome, bwa index
# on line: UR:http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_0.7.6a_ind/genome

# sort and index, takes ~30 min
nohup samtools sort leuk.bam leuk.sorted &
cat leuk.sorted*.bam > leuk.sorted.bam
samtools view leuk.sorted.bam | head
nohup samtools index leuk.sorted.bam &

samtools view -H leuk.sorted.bam | head
#	@SQ	SN:1	LN:249250621 ...
#	"SN:1" indicates chr 1 without using "chr"
nohup samtools index leuk.sorted.bam &

# summary:
samtools flagstat leuk.sorted.bam
#	1383866 + 0 in total (QC-passed reads + QC-failed reads)
#	11802 + 0 secondary
#	0 + 0 supplementary
#	252647 + 0 duplicates
#	1381670 + 0 mapped (99.84%:nan%)
#	1372064 + 0 paired in sequencing
#	686437 + 0 read1
#	685627 + 0 read2
#	1320716 + 0 properly paired (96.26%:nan%)
#	1367672 + 0 with itself and mate mapped
#	2196 + 0 singletons (0.16%:nan%)
#	14209 + 0 with mate mapped to a different chr
#	6545 + 0 with mate mapped to a different chr (mapQ>=5)

#####################

# healthy 
# AML patient FACS-purified T cells --> healthy.bam
samtools view healthy.bam | more
samtools view healthy.bam | tail
# header
samtools view -H healthy.bam

# sort and index, takes ~30 min
nohup samtools sort healthy.bam healthy.sorted &
cat healthy.sorted*.bam > healthy.sorted.bam
samtools view healthy.sorted.bam | head
nohup samtools index healthy.sorted.bam &
samtools view -H healthy.sorted.bam | head

# summary:
samtools flagstat healthy.sorted_chr.bam
#	1383838 + 0 in total (QC-passed reads + QC-failed reads)
#	11678 + 0 secondary
#	0 + 0 supplementary
#	261338 + 0 duplicates
#	1381163 + 0 mapped (99.81%:nan%)
#	1372160 + 0 paired in sequencing
#	685865 + 0 read1
#	686295 + 0 read2
#	1330018 + 0 properly paired (96.93%:nan%)
#	1366810 + 0 with itself and mate mapped
#	2675 + 0 singletons (0.19%:nan%)
#	15546 + 0 with mate mapped to a different chr
#	6563 + 0 with mate mapped to a different chr (mapQ>=5)

#####################

# Variant calling
samtools mpileup -g -f hg19_no_chr.fa healthy.sorted.bam leuk.sorted.bam > healthy_leuk.bcf	
bcftools view healthy_leuk.bcf | more

# call variants
bcftools call -v -m -O v -o healthy_leuk_var.vcf healthy_leuk.bcf
bcftools view healthy_leuk_var.vcf | more

# how many variants called?
bcftools view healthy_leuk_var.vcf | grep -v “^#” | wc -l
#   10448

# filter variants by quality score
bcftools view healthy_leuk_var.vcf | grep -v “^#” | cut -f6 | sort -k1n -u
bcftools view -i '%QUAL>=20' healthy_leuk_var.vcf > healthy_leuk_var_qual.vcf
bcftools view healthy_leuk_var_qual.vcf | grep -v “^#” | cut -f6 | sort -k1n -u
# quality is now >20, up to 486: 9691 variants remain

# Add back the "chr" to label chromosomes:
bcftools view healthy_leuk_var.vcf

# count number of variants on chr1
bcftools view healthy_leuk_var.vcf | cut -f1 | grep "chr1" | wc -l
# 10345

# View VCF results in IGV

#####################
# GATK somatic variant discovery Mutect2 (tumor vs. normal)
#####################

# preparation:

# get sample name:
samtools view -H leuk.sorted.bam | grep '@RG'
#	@RG     ID:189195       PL:illumina     PU:H33YCALXX.1  LB:A62931       
#	SM:AML 16626 Leukemic blasts    CN:BCCAGSC
# Name: 'AML 16626 Leukemic blasts'

samtools view -H healthy.sorted.bam | grep '@RG'
#	@RG	ID:188924	PL:illumina	PU:H32J2ALXX.1	LB:A62936	
#	SM:AML 16626 Healthy T cells	CN:BCCAGSC
# Name: 'AML 16626 Healthy T cells'

# -L: only end of chromosome 1 p-arm (chr1:1-98,661,704) to reduce comp time

# acquire reference vcf with allele frequencies
# download chr1 genome VCF and TBI (index) from:
# https://gnomad.broadinstitute.org/downloads#v2-variants 
# rename: gnomad.genomes.r2.1.1.sites.1.vcf.bgz --> gnomad-chr1.vcf.bgz
# rename: gnomad.genomes.r2.1.1.sites.1.vcf.bgz.tbi --> gnomad-chr1.vcf.bgz.tbi

# example of data column:
bcftools view -H gnomad-chr1.vcf.bgz | head -1
#	1	10067	rs1489251879	T	...

#####################

# only process chr1
gatk Mutect2 -R hg19_no_chr.fa -I leuk.sorted.bam -I healthy.sorted.bam -tumor 'AML 16626 Leukemic blasts'\
     -normal 'AML 16626 Healthy T cells' --germline-resource gnomad-chr1.vcf.bgz --af-of-alleles-not-in-resource 0.0000025\
     -O 1_somatic_m2.vcf.gz -bamout 2_tumor_normal.bam -L 1:1-249,250,621

# calculate contamination
gatk GetPileupSummaries -R hg19_no_chr.fa -I leuk.sorted.bam -V gnomad-chr1.vcf.bgz\
     -O leuk_getpileupsummaries.table -L 1:1-249,250,621

gatk CalculateContamination -I leuk_getpileupsummaries.table\
     -O leuk_calculatecontamination.table

gatk FilterMutectCalls -R hg19_no_chr.fa -V 1_somatic_m2.vcf.gz -L 1:1-249,250,621\
     --contamination-table leuk_calculatecontamination.table -O somatic_filtered.vcf.gz

gunzip somatic_filtered.vcf.gz
more somatic_filtered.vcf
# first variant record:
# ##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  AML 16626 Healthy T cells       AML 16626 Leukemic blasts
# 1       13813   .       T       G       .       contamination;haplotype;map_qual;normal_artifact        AS_FilterStatus=map_qual,contamination;AS_SB_TABLE=81,124|6,2;DP=217;ECNT=2;GERMQ=93;MBQ=41,41;MFRL=422,467;MMQ=25,24;MPOS=59;NALOD=-6.406e+00;NLOD#=14.49;POPAF=3.04;TLOD=10.40 GT:AD:AF:DP:F1R2:F2R1:PGT:PID:PS:SB     0|#0:89,3:0.042:92:49,2:36,1:0|1:13813_T_G:13813:32,57,3,0       0|1:116,5:0.048:121:56,2:52,3:0|1:13813_T_G:13813:49,67,3,2

# View somatic_filtered.vcf in IGV

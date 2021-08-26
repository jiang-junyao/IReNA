ATAC-seq preprocessing pipeline
================
Jiang junyao  GIBH, CAS

-   [ATAC-seq preprocessing pipeline](#atac-seq-preprocessing-pipeline)
    -   [Step 1: convert sra file to fastq
        file](#step-1-convert-sra-file-to-fastq-file)
    -   [Step 2: Quality control](#step-2-quality-control)
    -   [Step 3: Mapping](#step-3-mapping)
    -   [Step 4: Filter low-quality reads (MAPQ score &lt; 10), then
        remove duplicate with
        samtools](#step-4-filter-low-quality-reads-mapq-score--10-then-remove-duplicate-with-samtools)
    -   [Step 5: Call peaks with macs2](#step-5-call-peaks-with-macs2)
    -   [Step 6: Merge all peaks and generate gtf
        file](#step-6-merge-all-peaks-and-generate-gtf-file)
    -   [Step 7: Calculate raw counts withn
        htseq](#step-7-calculate-raw-counts-withn-htseq)
    -   [Step 8: Combine all samples, and calculate differential peaks
        according to
        edgeR](#step-8-combine-all-samples-and-calculate-differential-peaks-according-to-edger)
    -   [Step 9: Merge bam files](#step-9-merge-bam-files)
    -   [Step 10: Calculate footprints of merged bam files with
        HINT](#step-10-calculate-footprints-of-merged-bam-files-with-hint)

# ATAC-seq preprocessing pipeline

This pipeline is developed for ATAC-seq data pre-processing of
[IReNA](https://github.com/jiang-junyao/IReNA). Result of step4
(sample1\_rmdup.bam.bam)， step8 (differential\_peaks.bed) and step10
(footprints.bed) are input of IReNA. Test data used in
[IReNA](https://github.com/jiang-junyao/IReNA) can be download from
<https://www.ncbi.nlm.nih.gov/biosample?Db=biosample&DbFrom=bioproject&Cmd=Link&LinkName=bioproject_biosample&LinkReadableName=BioSample&ordinalpos=1&IdsFromResult=357084>.

## Step 1: convert sra file to fastq file

If you ATAC-seq data derive from [SRA
database](https://www.ncbi.nlm.nih.gov/sra), you need to use fastq-dump
function in [sra toolkit](https://github.com/ncbi/sra-tools) to convert
sra file to fastq file.

    ### For one sample
    fastq-dump --gzip --split-3 SRR5099531.1
    ### For multiple samples, use for loop to perform this operation on all files
    for (( i = 1; i < 7; i++)); do fastq-dump --gzip --split-3 SRR509953"$i".1 ; done;

## Step 2: Quality control

First, use
[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to
check reads quality.

    fastqc SRR509953*.fastq

Then trim reads. [fastp](https://github.com/OpenGene/fastp) is a user
friendly software which can remove adapator and low quality reads
automatically.

    ### For one sample
    fastp -i SRR5099531.1_1.fastq.gz -I SRR5099531.1_2.fastq.gz -o sample1_1.fastq.gz -O sample1_2.fastq.gz
    ### For multiple samples, use for loop to perform this operation on all files
    for (( i = 1; i < 7; i++)); do fastp -i SRR509953"$i".1_1.fastq.gz -I SRR509953"$i".1_2.fastq.gz -o sample"$i"_1.fastq.gz -O sample"$i"_2.fastq.gz ; done;

If you have adapator sequence, you can use
[cutadapt](https://cutadapt.readthedocs.io/en/stable/) to remove
adaptor.

    ### For one sample
    cutadapt -m 5 -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -o sample1_1.fastq.gz -p sample1_2.fastq.gz SRR5099531.1_1.fastq.gz SRR5099531.1_2.fastq.gz > sample1_log.txt
    ### For multiple samples, use for loop to perform this operation on all files
    for (( i = 1; i < 7; i++)); cutadapt -m 5 -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -A CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -o sample"$i"_1.fastq.gz -p sample"$i"_2.fastq.gz SRR509953"$i".1_1.fastq.gz SRR509953"$i".1_2.fastq.gz > sample"$i"_log.txt ; done;

## Step 3: Mapping

Before mapping, you need to create bowtie2 index with reference genome
which can be download from
[USCS](https://hgdownload.soe.ucsc.edu/downloads.html) or
[NCBI](https://www.ncbi.nlm.nih.gov/refseq/).

    bowtie2-build hg38/genome/hg38.fa hg38/bt2index/hg38_index
    ### For one sample
    bowtie2 -p 30 -x hg38/bt2index/hg38_index -1 sample1_1.fastq.gz -2 sample1_2.fastq.gz -S sample1.sam
    ### For multiple samples, use for loop to perform this operation on all files
    for (( i = 1; i < 7; i++)); do bowtie2 -p 30 -x hg38/bt2index/hg38_index -1 sample"$i"_1.fastq.gz -2 sample"$i"_2.fastq.gz -S sample"$i".sam ; done;

## Step 4: Filter low-quality reads (MAPQ score &lt; 10), then remove duplicate with [samtools](http://www.htslib.org/)

    ### For one sample
    samtools view -f 0x2 -q 10 -bh sample1.sam > sample1.bam
    samtools sort -@ 5 sample1.bam -o sample1_sorted.bam
    samtools rmdup sample1_sorted.bam sample1_rmdup.bam
    ### For multiple samples, use for loop to perform this operation on all files
    for (( i = 1; i < 7; i++)); do samtools view -f 0x2 -q 10 -bh sample"$i".sam > sample"$i"q10.bam ; done;
    for (( i = 1; i < 7; i++)); do samtools sort -@ 5 sample"$i"q10.bam > sample"$i"sorted.bam ; done;
    for (( i = 1; i < 7; i++)); do samtools rmdup sample"$i"sorted.bam sample"$i"rmdup.bam ; done;

## Step 5: Call peaks with [macs2](https://github.com/macs3-project/MACS)

    ### For one sample
    macs2 callpeak -t sample1_rmdup.bam --nomodel --shift -100 --extsize 200 -g mm -n sample1 > macs2_log.txt
    ### For multiple samples, use for loop to perform this operation on all files
    for (( i = 1; i < 7; i++)); do macs2 callpeak -t sample"$i"rmdup.bam --nomodel --shift -100 --extsize 200 -g hs -n sample"$i" > macs2_log.txt ; done;

## Step 6: Merge all peaks and generate gtf file

We merge all peaks through merge\_peak function in IReNA and generate
peaks file of bed format.

    ### R code
    library(IReNA)
    sample1_peaks <- read.delim("sample1_peaks.narrowPeak", header=FALSE)
    sample2_peaks <- read.delim("sample2_peaks.narrowPeak", header=FALSE)
    sample3_peaks <- read.delim("sample3_peaks.narrowPeak", header=FALSE)
    sample4_peaks <- read.delim("sample4_peaks.narrowPeak", header=FALSE)
    sample5_peaks <- read.delim("sample5_peaks.narrowPeak", header=FALSE)
    sample6_peaks <- read.delim("sample6_peaks.narrowPeak", header=FALSE)
    all_peak <- rbind(sample1_peaks, sample2_peaks, sample3_peaks, sample4_peaks, sample5_peaks, sample6_peaks)
    peaks_merged_gtf <- generate_gtf(all_peak)
    write.table(peaks_merged_gtf, 'peaks_merged.gtf', quote = F, row.names = F, col.names = F， sep = '\t')

## Step 7: Calculate raw counts withn [htseq](https://htseq.readthedocs.io/en/master/)

If sample has repeat, combine those repeat for further analysis. In our
test data, sample1rmdup.bam and sample2\_rmdup.bam are repeat of SSC
patient1 sample, sample3\_rmdup.bam and sample4\_rmdup.bam are repeat of
SSC patient2 sample, sample5\_rmdup.bam and sample6\_rmdup.bam are
repeat of esc sample.

    samtools merge SSC_patient1.bam sample1_rmdup.bam sample2_rmdup.bam
    samtools merge SSC_patient2.bam sample3_rmdup.bam sample4_rmdup.bam
    samtools merge esc.bam sample5_rmdup.bam sample6_rmdup.bam

Calculate raw counts

    htseq-count -r pos --idattr gene_id --stranded no -a 10 -f bam SSC_patient1.bam peaks_merged.gtf > SSC_patient1_Counts.txt
    htseq-count -r pos --idattr gene_id --stranded no -a 10 -f bam SSC_patient2.bam peaks_merged.gtf > SSC_patient2_Counts.txt
    htseq-count -r pos --idattr gene_id --stranded no -a 10 -f bam esc.bam peaks_merged.gtf > esc_Counts.txt

## Step 8: Combine all samples, and calculate differential peaks according to [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

edgeR needs more than two samples to make grouping to identify
differential peaks, so if your sample numbers are less than 2, just use
all peaks.

``` r
### R code
### If you have more than two samples, run the following codes
library(IReNA)
SSC_patient1_Counts <- read.delim("SSC_patient1_Counts.txt", header=FALSE)
SSC_patient2_Counts <- read.delim("SSC_patient2_Counts.txt", header=FALSE)
esc_Counts_Counts <- read.delim("esc_Counts.txt", header=FALSE)
peaks_gtf <- read.delim("peak_merged.gtf", header=FALSE)
count_all <- cbind(SSC_patient1_Counts, SSC_patient2_Counts[,2], esc_Counts_Counts[,2])
merged_count <- merge_sort_count(count_all, peaks_gtf)
group1 <- c(1,1,2) ### set group for three sample
differenetial_peaks1 <- diff_peaks(count_all[,4:ncol(diff_peaks)], group1) ### identify differential peaks
differenetial_peaks2 <- count_all[differenetial_peaks1$FDR<0.05,] ### filter peaks whose FDR is more than 0.05
write.table(differential_peaks[, c(1,2,3)], 'differential_peaks.bed', quote = F, row.name = F, col.names = F, sep = '\t')
### If you only have one sample, run the following codes
SSC_patient1_Counts <- read.delim("SSC_patient1_Counts.txt", header=FALSE)
peaks_gtf <- read.delim("peak_merged.gtf", header=FALSE)
merged_count <- merge_sort_count(SSC_patient1_Counts, peaks_gtf)
write.table(merged_count[, c(1,2,3)], 'peaks.bed', quote = F, row.name = F, col.names = F, sep = '\t')
```

## Step 9: Merge bam files

Merge all bam files for footprints calling

    samtools merge merged_all.bam sample1_rmdup.bam sample2_rmdup.bam sample3_rmdup.bam sample4_rmdup.bam sample5_rmdup.bam sample6_rmdup.bam

## Step 10: Calculate footprints of merged bam files with [HINT](http://www.regulatory-genomics.org/hint/introduction/)

Use hint to call footprints

    rgt-hint footprinting --atac-seq --paired-end --organism=hg38 merged_all.bam differential_peaks.bed

Use tag-count score &gt; 80th percentile as threshold to identify
footprints with high quality

``` r
### R code
footprints <- read.table('footprints.bed',sep='\t',header = T)
footprints_80th <- footprints[footprints$V5 > quantile(footprints$V5, 0.8),]
footprints_80th <- footprints_80th[,c(1,2,3,5)]
write.table(footprints_80th, 'filtered_footprints.bed', quote=F, sep = '\t', row.names = F, col.names = F)
```

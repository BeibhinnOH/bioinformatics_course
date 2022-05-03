#!/bin/bash


## 1.0: Pre Alignment Quality Control 

# 1.1: Perform basic quality assessment of untrimmed fastq files using FASTQC.
cd /home/ubuntu/ngs_course/dna_seq_pipeline/data/untrimmed_fastq
fastqc *.fastq.gz

# 1.2: Move fastqc results to results directory
mkdir ~/ngs_course/dna_seq_pipeline/results/fastqc_untrimmed_reads
cd ~/ngs_course/dna_seq_pipeline/data/untrimmed_fastq
mv *fastqc* /home/ubuntu/ngs_course/dna_seq_pipeline/results/fastqc_untrimmed_reads/

# 1.3: Unzip FASTQC zip files and generate summary of FASTQC results
cd /home/ubuntu/ngs_course/dna_seq_pipeline/results/fastqc_untrimmed_reads
for zip in *.zip; do unzip $zip; done
cat */summary.txt > ~/ngs_course/dna_seq_pipeline/logs/fastqc_untrimmed_summaries.txt

# 1.4: Trim fastq data using Trimmomatic to improve read quality.
cd ~/ngs_course/dna_seq_pipeline/data/untrimmed_fastq
trimmomatic PE  \
  -threads 4 \
  -phred33 \
  ~/ngs_course/dna_seq_pipeline/data/untrimmed_fastq/NGS0001.R1.fastq.gz /home/ubuntu/ngs_course/dna_seq_pipeline/data/untrimmed_fastq/NGS0001.R2.fastq.gz \
  -baseout ~/ngs_course/dna_seq_pipeline/data/trimmed_fastq/NGS0001_trimmed_R \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:50

# Remove data to clear space (we are only looking at paired end data)
rm NGS0001_trimmed_R_1U
rm NGS0001_trimmed_R_2U

# 1.5: Perform basic quality assessment of trimmed fastq files using FASTQC.
fastqc NGS0001_trimmed_R_1P NGS0001_trimmed_R_2P

# 1.6: Move FastQC results on paired trimmed data to results directory 
mkdir ~/ngs_course/dna_seq_pipeline/results/fastqc_trimmed_reads
cd ~/ngs_course/dna_seq_pipeline/data/trimmed_fastq
mv *fastqc* ~/ngs_course/dna_seq_pipeline/results/fastqc_trimmed_reads/

# 1.7: Unzip FastQC zip files and save record of FastQC summaries 
cd ~/ngs_course/dna_seq_pipeline/results/fastqc_trimmed_reads
for zip in *.zip; do unzip $zip; done
cat */summary.txt > ~/ngs_course/dna_seq_pipeline/logs/fastqc_trimmed_summaries.txt





## 2.0 Align paired trimmed fastq files using BWA MEM and reference genome hg19 

# 2.1 Generate BWA index files 
bwa index ~/ngs_course/dna_seq_pipeline/data/reference/hg19.fa.gz

# 2.2 Make an aligned_data folder and run BWA MEM with the following read group information: Read group identifier (ID): HWI-D0011.50.H7AP8ADXX.1.NGS0001 -Read group identifier (SM): NGS0001 -Read group identifier (PL): ILLUMINA -Read group identifier (LB): nextera-NGS001-blood -Read group identifier (PU): HWI-D00119 -Read group identifier (DT): 2022-03-23
mkdir ~/ngs_course/dna_seq_pipeline/data/aligned_data
bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001-blood\tDT:2022-03-23\tPU:HWI-D00119' -I 250,50  ~/ngs_course/dna_seq_pipeline/data/reference/hg19.fa.gz ~/ngs_course/dna_seq_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_1P ~/ngs_course/dna_seq_pipeline/data/trimmed_fastq/NGS0001_trimmed_R_2P > ~/ngs_course/dna_seq_pipeline/data/aligned_data/NGS0001.sam

# 2.3 Convert the sam file into bam format, sort it and generate an index using samtools 
cd ~/ngs_course/dna_seq_pipeline/data/aligned_data
samtools view -h -b NGS0001.sam > NGS0001.bam 
samtools sort NGS0001.bam > NGS0001_sorted.bam 

# 2.4 Generate a .bai index file
samtools index NGS0001_sorted.bam 



## 3.0 Perform duplicate marking 

# 3.1. Picard tools is used to mark duplicates. Picard tools examines aligned records in the .bam dataset to locate duplicate molecules. All records are then written to the output file with the duplicate records flagged. 
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt
samtools index NGS0001_sorted_marked.bam

# 3.2 Quality Filter the duplicate marked BAM file. Reads will be filtered according to the following criteria: Minimum MAPQ quality score : 20 -Filter on bitwise flag: yes a. Skip alignments with any of these flag bits set i. The read is unmapped ii. The alignment or this read is not primary iii. The read fails platform/vendor quality checks iv. The read is a PCR or optical duplicate.
$ samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
$ samtools index NGS0001_sorted_filtered.bam



## 4.0 Generate standard alignment statistics 

# 4.1 Use samtools flagstat to calculate and print statistics for NGS001_sorted_filtered.bam 
samtools flagstat NGS0001_sorted_filtered.bam > flagstat_output.txt
mv /home/ubuntu/ngs_course/dna_seq_pipeline/data/aligned_data/flagstat_output.txt /home/ubuntu/ngs_course/dna_seq_pipeline/results

# 4.2 Use samtools idxstats to generate alignment statistics per chromosome
samtools idxstats NGS0001_sorted_filtered.bam > idxstats_output.txt
mv /home/ubuntu/ngs_course/dna_seq_pipeline/data/aligned_data/idxstats_output.txt /home/ubuntu/ngs_course/dna_seq_pipeline/results

# 4.3 Use Picard insert_size_metrics to determine the distribution of insert sizes. To use picard Java is required. 
cd ~
java -version 

./gradlew shadowJar

java -jar build/libs/picard.jar

java -jar /home/ubuntu/picard/build/libs/picard.jar CollectInsertSizeMetrics \
      I= NGS0001_sorted_filtered.bam \
      O=insert_size_metrics.txt \
      H=insert_size_histogram.pdf \
      M=0.5

mv /home/ubuntu/ngs_course/dna_seq_pipeline/data/aligned_data/insert_size_metrics.txt /home/ubuntu/ngs_course/dna_seq_pipeline/results

# 4.4. Use bedtools to determine Depth of Coverage 

# Calculate depth of coverage for all regions in the .bam file using bedtools 
bedtools genomecov -ibam NGS0001_sorted_filtered.bam -bga -split > CoverageTotal.bedgraph.txt
mv /home/ubuntu/ngs_course/dna_seq_pipeline/data/aligned_data/CoverageTotal.bedgraph /home/ubuntu/ngs_course/dna_seq_pipeline/results



## 5.0 Variant Calling 

# Call Variants using Freebayes restricting the analysis to the regions in the bed file provided in assignment. 
# Freebayes determines the most-likely combination of genotypes for a population at each position in the reference. It produces a variant call file (VCF) 

# 5.1 Unzip hg19.fa.gz file 
zcat ~/ngs_course/dna_seq_pipeline/data/reference/hg19.fa.gz > ~/ngs_course/dna_seq_pipeline/data/reference/hg19.fa

# 5.2 Produce a .fai index file for FASTA files
samtools faidx ~/ngs_course/dna_seq_pipeline/data/reference/hg19.fa

# 5.3 Produce the VCF file using Freebayes 
freebayes --bam /home/ubuntu/ngs_course/dna_seq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference /home/ubuntu/ngs_course/dna_seq_pipeline/data/reference/hg19.fa --vcf /home/ubuntu/ngs_course/dna_seq_pipeline/results/NGS0001.vcf

# 5.4 Zip the NGS0001.vcf file 
bgzip /home/ubuntu/ngs_course/dna_seq_pipeline/results/NGS0001.vcf

# 5.5 Index the VCF with tabix
tabix -p vcf /home/ubuntu/ngs_course/dna_seq_pipeline/results/NGS0001.vcf.gz




## 6.0 Quality Filter variants 


# 6.1 Apply the freebayes hard filter for human diploid sequencing using vcffilter: QUAL > 1: removes horrible sites QUAL / AO > 10 : additional contribution of each obs should be 10 log units (~ Q10 per read) SAF > 0 & SAR > 0 : reads on both strands RPR > 1 & RPL > 1 : at least two reads “balanced” to each side of the site
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" /home/ubuntu/ngs_course/dna_seq_pipeline/results/NGS0001.vcf.gz > /home/ubuntu/ngs_course/dna_seq_pipeline/results/NGS0001_filtered.vcf

# 6.2 Filter the vcf file using bedtools for the regions in annotation.bed file provided in this assignment
bedtools intersect -header -wa -a /home/ubuntu/ngs_course/dna_seq_pipeline/results/NGS0001_filtered.vcf -b /home/ubuntu/ngs_course/dna_seq_pipeline/data/annotation.bed 
 	> /home/ubuntu/ngs_course/dna_seq_pipeline/results/NGS0001_filtered_annotation.vcf

# 6.3 Zip the file
bgzip /home/ubuntu/ngs_course/dna_seq_pipeline/results/NGS0001_filtered_annotation.vcf

# 6.4 Index the file with tabix
tabix -p vcf /home/ubuntu/ngs_course/dna_seq_pipeline/results/NGS0001_filtered_annotation.vcf.gz



## 7.0 Annotate variants with Annovar 

# ANNOVAR is used to annotate variants with respect to genes, databases of normal variantion and pathogenicity predictors. 
# Variant filters can be applied to the VCF file in excel to generate a list of candidate genes. 
# ANNOVAR was downloaded from Kai Wang (Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data, Nucleic Acids Research, 38:e164, 2010).

# 7.1 Download Annovar databases that are used for annotation
chmod +x annotate_variation.pl
cd annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/


# 7.2 VCF to Annovar input format 
./convert2annovar.pl -format vcf4 /home/ubuntu/ngs_course/dna_seq_pipeline/results/NGS0001_filtered_annotation.vcf.gz > /home/ubuntu/ngs_course/dna_seq_pipeline/results/NGS0001_filtered_annotation.avinput

# 7.3 Run Annovar table function -- csv output
./table_annovar.pl /home/ubuntu/ngs_course/dna_seq_pipeline/results/NGS0001_filtered_annotation.avinput humandb/ -buildver hg19 -out /home/ubuntu/ngs_course/dna_seq_pipeline/results/NGS0001_filtered_annotation -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout


# The output will be in csv format
# Download it via FileZilla and open it with Office Excel or any other speadsheet software.
# Use excel filtering to prioritise variants. 

# Note: pipeline.sh does not include installation of the tools, reference (hg19), or fastq files download. 
# This has been completed in project set up stage and detailed in the assignment.

# More work on the RNAseq data with BenE and BenF

Since our last meeting Xue has trimmed the sequence data from X. laevis and X. tropicalis using trimmomatic.  We now want to do the following:
- [ ] Set up a 'vanilla' trinity assembly for each transcriptome
- [ ] Download XL and XT genomes using wget
- [ ] Make a blast database out of each genome.  This will be used to assess quanity of the transcriptome assemblies from trinity.


# Summary of chat with BenF, BenE, and Xue

## Issues related to the technical accuracy of the assembly program
* Develop a metric for assembly fidelity based on BLAST (or BLAT) results for each transcriptome query.  This might, for example, include an assessment of whether or not all exons map to one chromosome, and whether they are in the proper order.
* Assess the impact of rate of evolution on the assembly fidelity
* Directions for the future: assess performance of other software

## Issues related to the biology of our organisms (i.e. diploid versus tetraploid
* For 'real' transcripts, do we see more alternate splicing in the tetraploid than the diploid?

# Trinity assembly

First we combined the trimmed forward and reverse reads like this:
```
zcat BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_001_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_002_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_003_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_004_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_005_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_006_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_007_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_001_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_002_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_003_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_004_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_005_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_006_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_007_paired.fastq.gz > BJE3909_paired_R1.fq.gz
```
```
zcat BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_001_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_002_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_003_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_004_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_005_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_006_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_007_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_001_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_002_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_003_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_004_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_005_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_006_paired.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_007_paired.fastq.gz > BJE3909_paired_R2.fq.gz
```

```
zcat BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_001_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_002_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_003_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_004_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_005_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_006_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R1_007_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_001_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_002_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_003_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_004_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_005_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_006_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R1_007_single.fastq.gz > BJE3909_single_R1.fq.gz
```
```
zcat BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_001_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_002_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_003_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_004_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_005_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_006_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L004_R2_007_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_001_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_002_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_003_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_004_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_005_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_006_single.fastq.gz BenEvansBJE3909cDNA_Library_GTGAAA_L005_R2_007_single.fastq.gz > BJE3909_single_R2.fq.gz
```

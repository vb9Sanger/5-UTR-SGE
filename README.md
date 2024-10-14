# 5-UTR-SGE-Pilot

## HDR Oligo Design 

**STEP ONE:** design HDR oligo libraries using [Mutator.VB.py](https://github.com/vb9Sanger/5-UTR/blob/main/Mutator.VB.py)

### Requirements:

This script requires the path to an **input txt file** containing the following tab-delimited columns:
  1. library identifier
  2. Forward TWIST primer (target-specific)
  3. mutated region (target-specific, containing PAM edits and Hap1 background variants)
  4. Reverse TWIST primer (target-specific)
  5. *pre-targeton sequence classed as '5' intronic'
  6. *pre-targeton sequence classed as 'splice site acceptor'
  7. *pre-targeton sequence classed as '5'UTR'
  8. *pre-targeton sequence classed as 'coding'
  9. *pre-targeton sequence classed as 'splice site donor'
  10. *pre-targeton sequence classed as '3' intronic'

  *columns 5-10 contain the corresponsing positions of bases (e.g. 42,228). If a given target does not contain a certain sequence type, this column is assigned a value of 0.

### Running the script:

Run: 
```bash
python3 Mutator.VB.py
```

### Output:

The output will be a txt file corresponding to the unique library identifier, target and pool number. 

### Notes:

The oligos contained in the output txt file need to be de-duplicated and supplemented with additional oligos not coded for in this script (for example, custom clinical indels, uorf-disrupting variants, back-up gRNA oligos)


## HDR Library QC

**STEP ONE:** Retreive fastq files from iRODS 

### Requirements:

[irods_to_lustre_bystudy_id.sh](https://github.com/vb9Sanger/5-UTR/blob/main/irods_to_lustre_bystudy_id.sh)

[irods_to_lustre.sh](https://github.com/vb9Sanger/5-UTR/blob/main/irods_to_lustre.sh)

### Running the script:

First, edit [irods_to_lustre_bystudy_id.sh](https://github.com/vb9Sanger/5-UTR/blob/main/irods_to_lustre_bystudy_id.sh) to contain the correct study ID and path to working directory. 

Then, run: 
```bash
./irods_to_lustre_bystudy_id.sh
```
This will generate a **samples.tsv** file in the output directory, under the folder 'metadata'. 

Select all samples of interest from this file, and copy this over to a folder called 'input'.

Then, run: 
```bash
./irods_to_lustre.sh
```
### Output:

This will generate merged cram files and fastqfiles. 

Output fastqfiles will be used for QC analysis. 


**STEP TWO:** Randomly subsample fastq files to 1M reads 

### Requirements:

raw fastq file 

[subsample.sh](https://github.com/vb9Sanger/5-UTR/blob/main/subsample.sh)

### Running the script:

Run: 
```bash
./subsample.sh input_fastq.qz 
```
### Output:
 
This will create a subdirectory called 'subsampled' containing new fastq files with 1M randomly subsampled reads.


**STEP THREE:** Generate count data 

### Requirements:

[demultiplex_UTR.txt](https://github.com/vb9Sanger/5-UTR/blob/main/demultiplex_UTR.txt) manifest file containing correct sample name(s), working directory, file location, sample specifications and primer sequences (ensure file type is ASCII text using dos2unix) 
[demultiplex_demultiplex_trim.sh](https://github.com/vb9Sanger/5-UTR/blob/main/demultiplex_demultiplex_trim.sh)
[demultiplex.sh](https://github.com/vb9Sanger/5-UTR/blob/main/demultiplex.sh)
[counting_extract_count.sh](https://github.com/vb9Sanger/5-UTR/blob/main/counting_extract_count.sh)
[count.sh](https://github.com/vb9Sanger/5-UTR/blob/main/count.sh)

### See:

Hong Kee's [sge-fastq-to-count](https://gitlab.internal.sanger.ac.uk/hk5/sge-fastq-to-count/-/tree/main)

### Running the script:

Step 1: Run the bsub demultiplex and trimming.

```bash
bash demultiplex.sh
```
Step 2: Run the bsub extract counting. Use cutadapt to remove primer and awk to count

```bash
bash count.sh
```
### Output:

This will generate a 'count' folder containing all_count files needed for QC 


**STEP FOUR:** Assess proportion of reads that map to deigned oligos 

### Requirements:

all_count file(s) generated in STEP THREE
[process_counts.py](https://github.com/vb9Sanger/5-UTR/blob/main/process_counts.py)
unique_trimmed.txt file(s) containing trimmed (F and R TWIST primers removed) deisgned oligos corresponding to a single library

### Running the script:

Run: 
```bash
./process_counts.py unique_trimmed.txt all_count.txt > run1.txt & 
```
### Output:

1. a 'final_counts' folder with FINAL_count.txt files containing all reads that exactly match designed oligos 
2. an 'error_counts' folder with ERROR_count.txt files containing all reads (that passed QC) and did NOT match any designed oligos 

**STEP FIVE:** Conduct independent QC Analysis assessing:

1. Total read counts
2. Proportion of subsampled reads that passed STEP THREE (Accepted reads)
3. Read length distribution (of subsampled reads)
4. Missing library sequences
5. Proportion of subsampled reass that reads that mapped to designed oligos (Mapped reads)
6. Genomic Coverage 














  

  

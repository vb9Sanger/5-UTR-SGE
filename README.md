# 5'UTR SGE Pilot

## HDR Oligo Design 

### **STEP ONE:** design HDR oligo libraries 

### Requirements:

* [Mutator.VB.py](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/Mutator.VB.py)
* an **input txt file** containing the following tab-delimited columns:
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

The oligos contained in the output txt file need to be de-duplicated and supplemented with additional oligos not coded for in this script (for example, custom clinical indels, uORF-disrupting variants, back-up gRNA oligos), prior to submission for TWIST oligo synthesis. 

--- 
--- 

## Guide RNA Cloning 

### Requirements:

* [gDNA cloning protocol](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/gDNA_cloning.md)

--- 

## WT Cloning 

### Requirements: 

* [WT cloning protocol](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/WT_cloning.md)
* [gel purification protocol](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/gel_purification.md) 

--- 

## HDR Library Assembly 

### Requirements: 

* [HDR Library Assembly protocol](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/HDR_library_assembly.md)
* [gel purification protocol](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/gel_purification.md)
* [PCR purification protocol](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/PCR_purification.md)

--- 

## HDR Library Sequencing Prep (NovaSeq)

### Requirements: 

* [HDR Library Sequencing Prep protocol](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/HDR_library_sequencing_prep.md)
* [PCR purification protocol](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/PCR_purification.md)
* [bead purification protocol](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/bead_purification.md)

### Notes: 
* samples are then submitted for KAPA Library Quantification and NovaSeq
  
--- 
--- 
## HDR Library QC

### **STEP ONE:** Retreive fastq files from iRODS 

### Requirements:

* [irods_to_lustre_bystudy_id.sh](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/irods_to_lustre_bystudy_id.sh)
* [irods_to_lustre.sh](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/irods_to_lustre.sh)

### Running the script:

First, edit [irods_to_lustre_bystudy_id.sh](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/irods_to_lustre_bystudy_id.sh) to contain the correct study ID and path to working directory. 

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

This will generate merged cram files and fastq files. 

Output fastq files will be used for QC analysis. 

--- 
### **STEP TWO:** Randomly subsample fastq files to 1M reads 

### Requirements:

* raw fastq file 
* [subsample.sh](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/subsample.sh)

### Running the script:

Run: 
```bash
./subsample.sh raw_fastq.qz 
```

### Output:
 
This will create a subdirectory called 'subsampled' containing new fastq files with 1M randomly subsampled reads.

--- 
### **STEP THREE:** Generate count data 

### Requirements:

* [demultiplex_UTR.txt](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/demultiplex_UTR.txt) manifest file 
* [demultiplex_demultiplex_trim.sh](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/demultiplex_demultiplex_trim.sh)
* [demultiplex.sh](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/demultiplex.sh)
* [counting_extract_count.sh](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/counting_extract_count.sh)
* [count.sh](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/count.sh)

### See:

Hong Kee's [sge-fastq-to-count](https://gitlab.internal.sanger.ac.uk/hk5/sge-fastq-to-count/-/tree/main) GitHub Repo.

### Running the script:

Step 1: Edit [demultiplex_UTR.txt](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/demultiplex_UTR.txt) manifest file in excel to contain correct sample name, working directory, path to Read1, exon identifier, and correct F/R TWIST primer sequences.  

Then, run:

```bash
bash dos2unix demultiplex_UTR.txt
```

Step 2: Run the bsub demultiplex and trimming.

```bash
bash demultiplex.sh
```

Step 3: Run the bsub extract counting. Use cutadapt to remove primer and awk to count.

```bash
bash count.sh
```

### Output:

This should generate 'trim', 'log', 'extracted', 'tempo' and 'count' folders, where the 'count' folder contains **all_count.txt** files needed for QC analysis. 

--- 
### **STEP FOUR:** Assess proportion of reads that map to deigned oligos 

### Requirements:

* **all_count.txt** file(s) generated in STEP THREE
* [process_counts.py](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/process_counts.py)
* unique_trimmed.txt file(s) containing trimmed (F and R TWIST primers removed) deisgned oligos corresponding to a single library

### Running the script:

Run: 
```bash
./process_counts.py unique_trimmed.txt all_count.txt > run1.txt & 
```

### Output:

1. a 'final_counts' folder with **FINAL_count.txt** files containing all reads that exactly match designed oligos 
2. an 'error_counts' folder with **ERROR_count.txt** files containing all reads (that passed QC) and did NOT match any designed oligos 

### Notes:

**unique_trimmed.txt** files were generated by manually editing de-duplicated txt files submitted for TWIST oligo synthesis. 

--- 
### **STEP FIVE:** Conduct independent QC Analysis assessing:

1. **Total read counts**
2. **Proportion of subsampled reads that passed STEP THREE (Accepted reads)**

Requirement: [accepted_reads.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/accepted_reads.R)

3. **Read length distribution (of subsampled reads)**

Requirement: 
* [length_distribution.py](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/length_distribution.py)
* [length_distribution.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/length_distribution.R)

Run: 
```bash
./length_distribution.py all_count.txt 
```
Output: **length_count.txt**

Then create corresponding histogram in R.

4. **Missing library sequences**
* Pass criterion: less than 1% of expected variants are missing
  
5. **Proportion of subsampled reads that map to designed oligos (Mapped reads)**

Requirement: [mapped_reads.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/mapped_reads.R)
* Pass criterion: more than 40% of accepted reads map to library reads

6. **Genomic Coverage**

Requirement: 
* [positions.py](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/positions.py)
* [genomic_coverage.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/genomic_coverage.R)

Run: 
```bash
./positions.py FINAL_count.txt
```
Output: **FINAL_count_positions.txt** 

Then, calculate log2(count+1) for each variant in excel and create a plot in R to visualise log2(count+1) variant variant position. All variants for which there is no corresponding position (e.g. whole targeton inversions) were assigned a position of 0.  
* distribution of variants should appear relatively tight
  
### Notes: 

* If libraries pass QC, OK to continue to Screening.
--- 
--- 
## SGE Screening 

### **STEP ONE:** 











  

  

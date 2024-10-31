# 5'UTR SGE Pilot

## HDR Oligo Design 

### STEP ONE: design HDR oligo libraries 

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
* The oligos contained in the output txt file need to be de-duplicated and supplemented with additional oligos not coded for in this script (for example, custom clinical indels, uORF-disrupting variants, back-up gRNA oligos), prior to submission for TWIST oligo synthesis. 

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
* Samples are then submitted for KAPA Library Quantification and NovaSeq
  
--- 
--- 
## HDR Library QC

### STEP ONE: Retreive fastq files from iRODS 

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

*Output fastq files will be used for QC analysis.*

--- 
### STEP TWO: Randomly subsample fastq files to 1M reads 

### Requirements:
* raw fastq file 
* [subsample.sh](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/subsample.sh)

### Running the script:
Run: 
```bash
./subsample.sh raw_fastq.qz 
```

### Output:
This will create a subdirectory called **'subsampled'** containing new fastq files with 1M randomly subsampled reads.

--- 
### STEP THREE: Generate count data 

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
This should generate **'trim'**, **'log'**, **'extracted'**, **'tempo'** and **'count'** folders, where the 'count' folder contains **all_count.txt** files needed for QC analysis. 

--- 
### STEP FOUR: Assess proportion of reads that map to deigned oligos 

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
1. A **'final_counts'** folder with **FINAL_count.txt** files containing all reads that exactly match designed oligos 
2. An **'error_counts'** folder with **ERROR_count.txt** files containing all reads (that passed QC) and did NOT match any designed oligos 

### Notes:
* **unique_trimmed.txt** files were generated by manually editing de-duplicated txt files submitted for TWIST oligo synthesis. 

--- 
### STEP FIVE: Conduct independent QC Analysis assessing

1. **Total read counts**
2. **Proportion of subsampled reads that passed STEP THREE (Accepted reads)**

**Requirement:** 
* [accepted_reads.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/accepted_reads.R)

3. **Read length distribution (of subsampled reads)**

**Requirement:** 
* [length_distribution.py](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/length_distribution.py)
* [length_distribution.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/length_distribution.R)

**Run:**
```bash
./length_distribution.py all_count.txt 
```
**Output: length_count.txt**

Then create corresponding histogram in R.

4. **Missing library sequences**
*Pass criterion: less than 1% of expected variants are missing*
  
5. **Proportion of subsampled reads that map to designed oligos (Mapped reads)**

**Requirement:** 
* [mapped_reads.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/mapped_reads.R)
*Pass criterion: more than 40% of accepted reads map to library reads*

6. **Genomic Coverage**

**Requirement:**
* [positions.py](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/positions.py)
* [genomic_coverage.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/genomic_coverage.R)

**Run:**
```bash
./positions.py FINAL_count.txt
```
**Output: FINAL_count_positions.txt** 

Then, calculate log2(count+1) for each variant in excel and create a plot in R to visualise log2(count+1) variant variant position. All variants for which there is no corresponding position (e.g. whole targeton inversions) were assigned a position of 0.  

*The distribution of variants should appear relatively tight*
  
### Notes: 
* If libraries pass QC, OK to continue to Screening.

--- 
--- 
## SGE Screening 

### Background:
This set of protocols follows a No-weekend protocol following the schedule specified below:

| Monday | Tuesday | Wednesday | Thursday| Friday |
| :--- | :--- | :--- | :--- |:--- |
| D0: electroporation | D1: selection | D2: selection| D3: media change | D4: split and harvest |
| D7: split and harvest | | D9: split | | D11: split and harvest|
|D14: split| | D16: split and harvest || D18: split|
|D21: harvest |||||

--- 
### STEP ONE: Electroporation

#### <ins>Day 0: Electroporation</ins>

#### Requirements: 
* [electroporation protocol](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/Electroporation.md)

##### Materials:
1x T75 per replicate

--- 
### STEP TWO: Screening

#### <ins>Day 1-2: Selection</ins>

##### Materials:
Hap1 media (IMDM + 10% FBS)
Blasticidin
Puromycin
Penicillin/Streptomycin (optional)

##### Notes: 
* From this day onward, you can put pen/strep in the media as an extra precaution for preventing contamination (using pen/strep does NOT replace using good sterile technique!). If you do decide to use pen/strep, use it the whole way throughout the screening process. 
* Pen/Strep is used at 1% working concentration (ie 5.5ml for a 550ml bottle of media). 

1.	24-30 hours after the electroporation, change the media to Hap1 media + 10ug/ml blasticidin + 2ug/ml puromycin
2.	Take a picture of each replicate using the EVOS
a.	Find a representative area of cells and focus
b.	Click the “freeze” button either by clicking on screen with the mouse or using the button on the front
c.	Click the “settings” icon on screen with the mouse and make sure Quick Save is turned OFF
d.	Attach a usb to the outlet (right hand side of screen)
e.	Click save (a pop-up will let you name the file)
3.	The next day, repeat steps 1 & 2


#### <ins>Day 3: Media change</ins>

##### Materials:
Hap1 media (IMDM + 10% FBS)
Blasticidin
Penicillin/Streptomycin (optional)

1.	After 48 hours of puromycin selection, change the media to Hap1 media + 10ug/ml blasticidin


#### <ins>Day 4: Split & harvest</ins>

##### Materials:
3x T150 per replicate 
Hap1 media (IMDM + 10% FBS)
Blasticidin
Penicillin/Streptomycin (optional)
TrypLE Express
PBS(-/-)
Countess slides/haemocytometer + Trypan Blue
Preprinted labels for your harvest tubes containing targeton name, screen timepoint, date


1.	Take a picture of the cells using the EVOS (see above for detailed instructions)
2.	Wash cells once with PBS (10ml)
3.	Add 2ml trypLE to each T75 and incubate for 3 minutes at 37C
4.	Add 8ml of Hap1 media to each flask and collect into a falcon tube
5.	Count cells using the countess. Take two counts. If they deviate by 10+%, re-suspend the collection buffer and count again. If the live counts are above 5e6/ml, add more media to the suspension and count again.
6.	For each replicate, 
a.	aliquot 5m cells into one Falcon tube and label it for re-seeding
b.	aliquot 6m cells into another Falcon tube and label it for harvest
c.	If your cell numbers are insufficient, please refer to the end of this protocol for a guide on how to proceed
7.	Centrifuge tubes at 300rcf for 3 minutes

**Re-Seeding:**
8.	For each replicate, prepare 3x T150 flasks with Hap1 media + 10ug/ml blasticidin
9.	Aspirate the supernatant media from the re-seeding pellet and re-suspend it in 3ml Hap1 media
10.	Put 1ml each of the cell suspension into the prepared T150s (~1.7m cells/flask). Gently agitate the flask and transfer them to the incubator

**Harvests:**
11.	Aspirate the supernatant media from the harvest pellet and re-suspend in 2ml PBS
12.	Put 1ml each into two Eppendorf tubes*
13.	Spin the Eppendorfs in the microcentrifuge for 1 min
14.	Aspirate the PBS carefully without disturbing the pellet
15.	Immediately snap-freeze pellets on dry ice, or transfer them to a box in the -80C freezer (for safety, do not store SGE screen pellets in your own -20 freezers)

##### Notes:
* If you are submitting DNA pellets to pipelines, you should use 2ml Eppendorf tubes. If we are processing DNA/library prep in our lab you can use either 1.5ml or 2ml tubes. This applies to all timepoints. 


#### <ins>Day 7: Split & harvest</ins>
*From this day, you can leave out blasticidin from the media*

##### Materials:
2x T150 per replicate 
Hap1 media (IMDM + 10% FBS)
Penicillin/Streptomycin (optional)
TrypLE Express
PBS(-/-)
Countess slides/haemocytometer + Trypan Blue
Preprinted labels for your harvest tubes containing targeton name, screen timepoint, date

1.	Take a picture of the cells using the EVOS (see above for detailed instructions)
2.	Wash cells once with PBS (10-20ml)
3.	Add 3ml trypLE to each T150 and incubate for 3 minutes at 37C
4.	Add 7ml of Hap1 media to each flask and collect into a falcon tube
5.	Count cells using the countess. Take two counts. If they deviate by 10+%, re-suspend the collection buffer and count again. If the live counts are above 5e6/ml, add more media to the suspension and count again.
6.	For each replicate, 
a.	aliquot 5m cells into one Falcon tube and label it for re-seeding
b.	Aliquot 6m cells into another Falcon tube and label it for harvest
7.	Centrifuge tubes at 300rcf for 3 minutes

**Re-Seeding:**
8.	For each replicate, prepare 2x T150 flasks with 30ml Hap1 media each
9.	Aspirate the supernatant media from the re-seeding pellet and re-suspend it in 2ml Hap1 media
10.	Put 1ml each of the cell suspension into the prepared T150s (~2.5m cells/flask). Gently agitate the flask and transfer them to the incubator

**Harvests:**
11.	Aspirate the supernatant media from the harvest pellet and re-suspend in 2ml PBS
12.	Put 1ml each into two Eppendorf tubes
13.	Spin the Eppendorfs in the microcentrifuge for 1 min
14.	Aspirate the PBS carefully without disturbing the pellet
15.	Immediately snap-freeze pellets on dry ice, or transfer them to a box in the -80C freezer (for safety, do not store SGE screen pellets in your own -20 freezers)


#### <ins>Day 9: Split</ins>

##### Materials:
2x T150 per replicate 
Hap1 media (IMDM + 10% FBS)
Penicillin/Streptomycin (optional)
TrypLE Express
PBS(-/-)
Countess slides/haemocytometer + Trypan Blue

1.	Wash cells once with PBS (10-20ml)
2.	Add 3ml trypLE to each T150 and incubate for 3 minutes at 37C
3.	Add 7ml of Hap1 media to each flask and collect into a falcon tube
4.	Count cells using the countess. Take two counts. If they deviate by 10+%, re-suspend the collection buffer and count again. If the live counts are above 5e6/ml, add more media to the suspension and count again.
5.	For each replicate aliquot 5m cells into one Falcon tube
6.	Centrifuge tubes at 300rcf for 3 minutes
7.	For each replicate, prepare 2x T150 flasks with 30ml Hap1 media
8.	Aspirate the supernatant media from the re-seeding pellet and re-suspend it in 2ml Hap1 media
9.	Put 1ml each of the cell suspension into the prepared T150s (~2.5m cells/flask). Gently agitate the flask and transfer them to the incubator


#### <ins>Day 11: Split & harvest</ins>

##### Materials:
3x T150 per replicate 
Hap1 media (IMDM + 10% FBS)
Penicillin/Streptomycin (optional)
TrypLE Express
PBS(-/-)
Countess slides/haemocytometer + Trypan Blue
Preprinted labels for your harvest tubes containing targeton name, screen timepoint, date

1.	Wash cells once with PBS (10-20ml)
2.	Add 3ml trypLE to each T150 and incubate for 3 minutes at 37C
3.	Add 7ml of Hap1 media to each flask and collect into a falcon tube
4.	Count cells using the countess. Take two counts. If they deviate by 10+%, re-suspend the collection buffer and count again. If the live counts are above 5e6/ml, add more media to the suspension and count again.
5.	For each replicate, 
a.	aliquot 5m cells into one Falcon tube and label it for re-seeding
b.	Aliquot 6m cells into another Falcon tube and label it for harvest
6.	Centrifuge tubes at 300rcf for 3 minutes

**Re-Seeding**
7.	For each replicate, prepare 3x T150 flasks with 30ml Hap1 media each
8.	Aspirate the supernatant media from the re-seeding pellet and re-suspend it in 3ml Hap1 media
9.	Put 1ml each of the cell suspension into the prepared T150s (~1.7m cells/flask). Gently agitate the flask and transfer them to the incubator

**Harvests:**
10.	Aspirate the supernatant media from the harvest pellet and re-suspend in 2ml PBS
11.	Put 1ml each into two Eppendorf tubes
12.	Spin the Eppendorfs in the microcentrifuge for 1 min
13.	Aspirate the PBS carefully without disturbing the pellet
14.	Immediately snap-freeze pellets on dry ice, or transfer them to a box in the -80C freezer (for safety, do not store SGE screen pellets in your own -20 freezers)


#### <ins>Day 14: Split</ins>
*Repeat the steps from Day 9*

##### Materials:
2x T150 per replicate 


#### <ins>Day 16: Split & harvest</ins>
*Repeat the steps from Day 7*

##### Materials:
2x T150 per replicate 


#### <ins>Day 18: Split</ins>
*Repeat the steps from Day 9*

##### Materials:
3x T150 per replicate 


#### <ins>Day 21: Harvest</ins>
*Repeat the steps from Day 7 (no splitting)* 

*Screen done!*

##### Notes:
* If you are extending your screen follow the Steps from Day 16. Split again on day 23 (repeating steps from Day 9) and 25 (as day 18), and harvest on day 28.








  

  

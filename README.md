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

#### <ins>Day 1-21: Screening</ins>

#### Requirements: 
* [SGE screening protocol](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/SGE_screening.md) 

--- 
## gDNA Extraction 

### Requirements:
* [gDNA extraction protocol](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/gDNA_extraction.md)

### Notes:
* Where applicable, do not extract DNA from all pellets corresponding to a single replicate and timepoint at once

--- 

## NGS Library Prep (gDNA) 

### Requirements:
* [NGS library prep - gDNA protocol](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/NGS_library_prep_gDNA.md)

--- 
--- 

## Data Analysis 

### Background:
Downstream analysis of SGE screening data is performed using the established Sanger **sge-metapipeline**, which wraps sample retrieval (iRODS to lustre), quantification (QUANTS) and QC (MAVEQC) into a single pipeline.

* Full pipeline repo: [sge_production_screen_qc](https://gitlab.internal.sanger.ac.uk/team302/sge/sge_production_screen_qc/-/tree/feature/run46503_Jan2023?ref_type=heads)
* QUANTS output docs: [QUANTS output.md](https://github.com/cancerit/QUANTS/blob/develop/docs/output.md), [pyQUEST](https://github.com/cancerit/pyQUEST)
* Example `meta.csv`, `meta_consequence.tsv` and `config.ini` files: `/lustre/scratch124/humgen/projects_v2/5utr_sge/users/vb9/sge_analysis/files_to_run_metapipeline/logs`

--- 
### STEP ONE: Load the sge-metapipeline

#### Requirements:
* Access to the `HGI/common/sge-metapipeline` module

#### Running the script:
Run: 
```bash
module load HGI/common/sge-metapipeline/v0.1.1
```

Once loaded, there are two basic commands, `fetch` and `execute`, each with a `--help` option:
```bash
sge-metapipeline fetch --help
sge-metapipeline execute --help
```

#### Notes:
* `fetch` retrieves sample information from the MLWH and creates the initial Sample Manifest
* `execute` runs the `irods_to_lustre`, `QUANTS` and `MAVEQC` sub-commands

--- 
### STEP TWO: Fetch sample manifest and run iRODS to lustre

#### Requirements:
* Study ID and run ID for the samples of interest
* `config.ini` file with warehouse credentials

#### Running the script:
First, load the iRODS module and initiate a session:
```bash
module load ISG/IRODS/1.0
iinit
```
*`iinit` will prompt you to enter your Sanger password.*

Fetch sample information from the MLWH (example):
```bash
sge-metapipeline fetch --study-id 7885 --run-id 49656 --warehouse-creds config.ini
```
Modify the resulting `manifest_fetched.tsv` if needed, to include only the desired samples.

Then run iRODS to lustre (example):
```bash
module load badger/samtools/1.20
sge-metapipeline execute irods_to_lustre --manifest-file manifest_fetched.tsv --output irods_to_lustre --pipeline-config config.ini
```

#### Output:
`manifest_fetched.irods_to_lustre.tsv`, along with merged cram/fastq files on lustre, to be used as input for QUANTS.

#### Notes:
* You may need to load `samtools` (`badger/samtools/1.20`) before running `irods_to_lustre`

--- 
### STEP THREE: Run QUANTS

#### Requirements:
* `meta.csv` file (see example file location above)
* `manifest_fetched.irods_to_lustre.tsv`, populated with the path to the `meta.csv` file
* pipeline `config.ini` file

#### Running the script:
Run (example):
```bash
sge-metapipeline execute QUANTS --manifest-file irods_to_lustre/sg5_manifest_fetched.irods_to_lustre.tsv --output sg5_SLC2A1_QUANTS --pipeline-config sge_metapipeline_vb9.config.ini
```

#### Output:
See expected output structure: [QUANTS output.md](https://github.com/cancerit/QUANTS/blob/develop/docs/output.md), [pyQUEST](https://github.com/cancerit/pyQUEST)

--- 
### STEP FOUR: Run MAVEQC

### Background:
QC is performed using an adapted version of [MAVE-QC](https://github.com/wtsi-hgi/MAVEQC), with the adapted script [run_maveqc_VB.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/run_maveqc_VB.R).

#### Requirements:
* `HGI/common/sge-metapipeline` module
* `sg5_manifest_fetched.irods_to_lustre.quants.tsv` manifest generated by QUANTS
* pipeline `config.ini` file
* [run_maveqc_VB.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/run_maveqc_VB.R)

#### Running the script:
First, load the sge-metapipeline module (if not already loaded):
```bash
module load HGI/common/sge-metapipeline/v0.1.1
```

Then run MAVEQC (example):
```bash
sge-metapipeline execute MAVEQC --manifest-file sg5_SLC2A1_QUANTS/sg5_manifest_fetched.irods_to_lustre.quants.tsv --output sg5_SLC2A1_MAVEQC --pipeline-config sge_metapipeline_vb9.config.ini
```
This will generate a folder called **'input'**. Once this has been generated, the job can be killed. Delete the *contents* of the output folder (but not the output folder itself).

Then run the adapted script:
```bash
Rscript Code/run_maveqc_VB.R <path/to/input_dir> <path/to/output_dir>
```
*where `<path/to/input_dir>` is the folder one level above the 'input' folder generated above, and `<path/to/output_dir>` is the desired output location.*

The script takes exactly two positional arguments (no flags) — `<input_dir>` and `<output_dir>` — and expects:
* `<input_dir>/input/sample_sheet.tsv`, containing (tab-delimited, header required) the columns: `sample_name`, `replicate`, `condition`, `ref_time_point`, `library_independent_count`, `library_dependent_count`, `valiant_meta`, `vep_anno`, `adapt5`, `adapt3`, `per_r1_adaptor`, `per_r2_adaptor`, `library_name`, `library_type` — all generated automatically as part of the `sge-metapipeline execute MAVEQC` step above
* exactly one `ref_time_point` value (matching one of the `condition` values), used as the reference for the DESeq2 comparisons
* R packages: `configr`, `vroom`, `data.table`, `Ckmeans.1d.dp`, `gplots`, `ggplot2`, `plotly`, `ggcorrplot`, `corrplot`, `see`, `ggbeeswarm`, `reactable`, `htmltools`, `sparkline`, `dendextend`, `reshape2`, `gtools`, `DESeq2`, `DEGreport`, `apeglm` (the script will attempt to auto-install any missing packages)

#### Output:
The script creates three subfolders in `<output_dir>`:

* **`plasmid_qc/`** and **`screen_qc/`** — sample-level QC plots and stats tables (read counts, missing variants, mapping rates, etc., against the thresholds in the auto-generated `config.yaml`)
* **`experiment_qc/`** — the DESeq2-based experiment QC, including:
  * sample correlation and PCA plots
  * bee swarm and position plots of log2FoldChange (the plots you assess in STEP FIVE to choose a timepoint/reference comparison)
  * a positional loess diagnostic plot
  * `library_deseq2_results_<comparison>.tsv` and `all_deseq2_results_<comparison>.tsv` per condition-vs-reference comparison (the latter are the files used going into STEP FIVE/SIX)
  * `library_normalized_counts.tsv` and `all_normalized_counts.tsv`

A `config.yaml` is also written to `<output_dir>`, recording the QC thresholds used for that run (e.g. Gini coefficient, minimum total reads, missing variant %, mapping %, LFC cutoffs). These are currently hardcoded defaults inside the script itself (regenerated fresh each run) rather than read from an existing file — to change a threshold, edit the default value in the script (search for `maveqc_config <- list(...)`) and re-run.

#### Notes:
* `run_maveqc_VB.R` is an adapted version of the original [MAVE-QC](https://github.com/wtsi-hgi/MAVEQC) script.

--- 
### STEP FIVE: Select timepoint and reference comparison

#### Requirements:
* `experiment_qc` output folder generated by `run_maveqc_VB.R`

#### Running the script:
No script to run here — assess the bee swarm and position plots in the `experiment_qc` output folder for each condition, and select the desired timepoint and reference comparison (e.g. `D16vD4`, `D16vPlasmid`).

#### Output:
A chosen `all_deseq2_results_condition_<timepoint>_vs_<reference>.tsv` file (or files, if combining guides) to take forward to STEP SIX.

--- 
### STEP SIX: Combine or extract guide-level results

#### Background:
* If there are **2 guides for 1 gene**, run [merge_and_combine.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/merge_and_combine.R) to combine both guides into a single positional/adjusted result.
* If there is **only 1 guide**, run [extract_single_guide.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/extract_single_guide.R) to extract and filter the result table.

#### Requirements:
* R package `data.table` (required); `ggplot2` (optional, only needed if `--plot_prefix` is used)
* the `all_deseq2_results_condition_*.tsv` file(s) selected in STEP FIVE

#### Option A — Two guides: `merge_and_combine.R`

Merges two guide/library result TSVs and computes a combined log2FoldChange (LFC) via inverse-variance weighting.

| Flag | Required? | Description |
| :--- | :--- | :--- |
| `-f1 <file1.tsv>` | Yes | Guide/library 1 results TSV |
| `-f2 <file2.tsv>` | Yes | Guide/library 2 results TSV |
| `-o <out.tsv>` | No | Output TSV (default: `combined_output.tsv`) |
| `--mode <adjusted\|positional>` | No | Result type to combine (default: `adjusted`) |
| `--effect <raw\|shrunk>` | No | Effect size type (default: `raw`) |
| `--drop_prefix_fields <int>` | No | Number of underscore-delimited fields to drop from `oligo_name` when building the merge key (default: `2`) |
| `--pos_col1 <name>` / `--pos_col2 <name>` | No | Position column name in file1/file2 (default: `position`) |
| `--pam_regions1 "start-end,..."` / `--pam_regions2 "start-end,..."` | No | PAM ranges (genomic coordinates) to flag per library, e.g. `"42958767-42958789"` |
| `--pam_col1 <name>` / `--pam_col2 <name>` | No | Alternative to `--pam_regions*`: an existing column in file1/file2 flagging PAM overlap |
| `--exclude_pam <TRUE\|FALSE>` | No | If `TRUE`, sets weights of PAM-overlapping variants to 0 when combining (default: `TRUE`) |
| `--pcut <double>` | No | FDR cutoff used to assign `combined_status` (default: `0.05`) |
| `--plot_prefix <prefix>` | No | If set, writes diagnostic plots (requires `ggplot2`): `<prefix>_lfc_scatter.png`, `<prefix>_combined_SE_hist.png`, `<prefix>_pos_discordant.tsv` |

Run (example):
```bash
Rscript Code/merge_and_combine.R -f1 SLC2A1/D4_ref/sg5/output/experiment_qc/all_deseq2_results_condition_Day16_vs_Day4.tsv -f2 SLC2A1/D4_ref/sg6/output/experiment_qc/all_deseq2_results_condition_Day16_vs_Day4.tsv -o SLC2A1/D4_ref/merged_positional/SLC2A1_merged_positional_D16vD4.tsv --mode positional --effect raw --pam_regions1 "42958767-42958789" --pam_regions2 "42958622-42958644" --plot_prefix sg5_sg6_D16vD4_positional_diagnostics
```

##### Output:
* The specified `-o` TSV, containing both guides' sequences, oligo names, positions, LFCs, SEs and statuses (suffixed `_1`/`_2`), a reconciled `position`/`position_source`, and the combined `combined_LFC`, `combined_SE`, `combined_Z`, `combined_p`, `combined_FDR` and `combined_status` columns.
* If `--plot_prefix` is set: an LFC scatter plot (guide 1 vs guide 2), a combined SE histogram, and a table of variants with discordant positions between the two guides.

#### Option B — One guide: `extract_single_guide.R`

Extracts a compact, self-describing results table from a single guide TSV.

| Flag | Required? | Description |
| :--- | :--- | :--- |
| `-f <input.tsv>` | Yes | Input guide results TSV |
| `-o <output.tsv>` | No | Output TSV (default: `single_output.tsv`) |
| `--mode <adjusted\|positional>` | No | Result type to extract (default: `positional`) |
| `--effect <raw\|shrunk>` | No | Effect size type (default: `raw`) |
| `--drop_prefix_fields <int>` | No | Number of underscore-delimited fields to drop from `oligo_name` when building the variant key (default: `2`) |
| `--pos_col <name>` | No | Position column name (default: `position`) |
| `--pam_regions "start-end,..."` | No | PAM ranges (genomic coordinates) to flag, e.g. `"33543125-33543147"` |
| `--pam_col <name>` | No | Alternative to `--pam_regions`: an existing column flagging PAM overlap |
| `--exclude_pam <TRUE\|FALSE>` | No | If `TRUE`, PAM-overlapping variants are flagged and written to a separate audit file (they are retained, not removed, in the main output) (default: `TRUE`) |
| `--removed_out <removed.tsv>` | No | Path for the PAM-overlapping audit file (default: `<output_basename>_pam_flagged.tsv`) |

Run (example):
```bash
Rscript Code/extract_single_guide.R -f SON/Plasmid_ref/sg10/output/experiment_qc/all_deseq2_results_condition_Day16_vs_Plasmid.tsv -o SON/Plasmid_ref/extracted_positional/SON_extracted_positional_D16vPlasmid.tsv --mode positional --effect raw --pam_regions "33543125-33543147"
```

##### Output:
* The specified `-o` TSV, containing `variant_key`, `oligo_name`, `consequence`, `sequence`, position, `pam_flag`, and the LFC/SE/p-value/FDR/stat columns corresponding to the chosen `--mode`/`--effect`.
* If any variants are flagged as PAM-overlapping, an additional audit TSV of those variants (retained in the main output, but also listed separately).

#### Notes:
* Both scripts exclude variants with `consequence` of `"Others"` or `"Backup_gRNA"` before output.
* Both scripts support `--help` for the full built-in usage message.
* For `--mode positional`, both scripts use `pos_total_se_raw`/`pos_total_se_shrunk` as the SE column, falling back to `lfcSE_raw`/`lfcSE_shrunk` (with a warning) if those aren't present in the input (e.g. outputs generated before `pos_total_se_*` was added). For `--mode adjusted`, `lfcSE_raw`/`lfcSE_shrunk` is always used.

--- 
### STEP SEVEN: Fit a Gaussian Mixture Model — [gmm.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/gmm.R)

#### Background:
Fits a 2-component Gaussian Mixture Model (via `mclust`) to the combined or extracted LFC values from STEP SIX, to classify variants as **depleted** or **no impact** (enriched variants are left unchanged).

#### Requirements:
* R packages: `data.table`, `mclust`, `ggplot2`, `ragg`
* the combined (`merge_and_combine.R`) or extracted (`extract_single_guide.R`) output TSV from STEP SIX
* [gmm.R](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/Code/gmm.R)

| Flag | Required? | Description |
| :--- | :--- | :--- |
| `-f <input.tsv>` | Yes | Combined or single-guide results TSV from STEP SIX |
| `-o <output.tsv>` | No | Output TSV (default: `<input_basename>_with_GMM.tsv`, same directory as input) |
| `--lfc_col <name>` | No | LFC column to fit on (auto-detected if omitted: `combined_LFC`, or one of the `adj_log2FoldChange_*`/`pos_adj_log2FoldChange_*` columns) |
| `--status_col <name>` | No | Status column to use (auto-detected from `--lfc_col` if omitted) |
| `--seed <int>` | No | RNG seed (default: `1`) |
| `--plot_prefix <path>` | No | Diagnostic plot prefix (default: next to output, `<out_basename>_GMM`) |

Run (example):
```bash
Rscript Code/gmm.R -f SLC2A1/D4_ref/merged_positional/SLC2A1_merged_positional_D16vD4.tsv -o SLC2A1/D4_ref/merged_positional/SLC2A1_merged_positional_D16vD4_with_GMM.tsv
```

#### Output:
* The specified `-o` TSV: all original columns, plus `GMM_fit_set`, `GMM_cluster`, `GMM_label`, `GMM_prob_cluster1`/`GMM_prob_cluster2`, `GMM_sign_override`, `GMM_status`, and `GMM_status_changed`.
* A diagnostic plot (`<prefix>_plot.png`) showing variants coloured by fitted GMM cluster, with cluster mean lines.

#### Notes:
* The model is fit only on variants with a status of `depleted` or `no impact`; `enriched` variants are excluded from fitting and left unchanged.
* The number of components (`G`) is fixed at 2 and is not user-configurable.
* Any variant whose fitted component would be labelled "depleted" but which has a positive LFC is automatically reassigned to "no impact" (flagged via `GMM_sign_override`), since a positive LFC cannot represent depletion.

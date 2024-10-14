# 5-UTR-SGE-Pilot

## HDR Oligo Design 

Step ONE: design HDR oligo libraries using [Mutator.VB.py](https://github.com/vb9Sanger/5-UTR/blob/main/Mutator.VB.py)

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



  

  

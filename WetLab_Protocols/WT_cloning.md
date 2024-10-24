# WT Vector Cloning 

## Homology Arm Insert

### PCR amplify HA region (twist region + homology arms)

| Reagent | Volume |
| :--- | :--- |
| WT Hap1 gDNA	| 100ng |
| HA Primer F (10uM)	| 1.5ul |
| HA Primer R (10uM)	| 1.5ul |
| Kapa Hifi MM |	25ul |
| H2O	| to 50ul |

1.	Run PCR protocol under Rebecca’s name 
2.	Run on gel
3.	Excise bands 
4.	[Gel purify](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/gel_purification.md) (elute in 12ul EB buffer) 
5.	Nanodrop 

### Digest with NotI and SbfI:

| Reagent | Volume |
| :--- | :--- |
| Purified insert |	11ul |
| SbfI	| 1ul |
|NotI |	1ul |
|Cutsmart Buffer |	2ul |
|H2O |	to 20ul |

1.	Keep at 37C for at least 1hr, 80C for 20min



## WT Backbone:

| Reagent | Volume |
| :--- | :--- |
| Uncut pmin plasmid |	4000ng |
| SbfI	| 1ul |
| NotI |	1ul |
| Cutsmart Buffer |	2ul |
| H2O |	to 20ul |

1.	Keep at 37C for at least 1hr, 80C for 20min
2.	Run on a 1% gel 
3.	Excise band 
4.	[Gel purify](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/gel_purification.md) (Elute in 12ul EB buffer)
5.	Nanodrop


## Backbone and insert ligation: 

* First: Calculate insert concentrations needed for 3:1 and 1:1 ratios of insert to vector on NEBiocalculator 
* Assumes 1ul vector 

| Reagent | Volume |
| :--- | :--- |
|DNA (insert)	| * |
|Purified / digested Backbone |	1ul |
| Ligase |	1ul |
| Ligase buffer |	2ul | 
| H2O	| to 20ul |

*NEBiocalculator output 

1.	Incubate at room temp for 10min 
2.	Heat inactivate for 10min at 65C
3.	Place on ice 


### Transformation: 

1. Thaw a tube of NEB 10-beta E.coli cells on ice
2. Add 2ul plasmid DNA to 50ul of cells
3. Flick 4-5x to mix
4. Place mix on ice for 30min
5. Heat shock at 42C for exactly 30sec
6. Place on ice for 5min
7. Pipette 950ul of RT outgrowth medium into mix
8. Place at 37C for 60min and shake at 250rpm
9. Warm plates to 37C
10. Mix cells thoroughly by flicking / inverting
11. Spread 1 plate with 40ul of mix, and another with 400ul evenly
12. Incubate overnight at 37C

### Miniprep - (QIAGEN - Cat. No. 28106)
* Create 1:1000 LB:AMP (e.g. 400ul AMP in 400ml LB broth)
* Add 3ml of LB + AMP to tubes with clear caps
* Select a single colony with a P200
* Dump whole tip into tube
* At end of day: put tubes in shaker overnight at 37C

* Next day: remove tubes from shaker
* Swirl, then add 2ml to a 2ml tube 
* Miniprep according to manufarturer’s protocol

  
### Diagnostic Digest 
* Prior to sequencing, perform a diagnostic digest on 2ul of samples 

| Reagent | Volume |
| :--- | :--- |
| DNA	| 2ul |
| NspI	| 1ul |
| Cutsmart buffer |	2ul |
| H2O	| to 20ul |
* 37C for 1hr

# Notes:
* NspI selected based off of expexted band patterns on simulated Snapgene gel
* Looked for RE that would produce identifiable band patterns for all samples 

### Sequence 
* Sanger sequence to validate constructs
* Use minBBseq and HDR WT primers F/R (four tubes per sample)
* Check traces carefully to check for introduced variants (accounting for Hap1 variants)





# SGE HDR Plasmid DNA library NovaSeq Sequencing Prep

## STEP ONE: Maxiprep (QIAGEN - Cat. No. 12963)

### Notes:
* Maxiprep plasmid libraries according to manufacturer’s protocols
* Create two sequential elutions in 200ul EB buffer

## STEP TWO: Test ILL1 qPCR

For each 50ul Reaction:

| Reagent | Volume |
| :--- | :--- |
|HDR Library (Maxiprep)|	25ng|
|F ILL1 Primer (10uM)|	1.5ul|
|R ILL1 Primer (10uM)|	1.5ul|
|Kapa Hifi MM	|25ul|
|20x Evagreen |	2.5ul|
|100x ROX|	0.5ul|
|H2O| 	To 50ul |

### Notes: 
* Amplicon size 250-350bp
* Make up enough master mix for 5 reactions (to test 4 different temperatures)
* Perform 4 reactions across 4 temperatures between 57C and 63C


### PCR Conditions: 

| Temperature       | Duration           |
| ------------- |:-------------:|
|95C	|3min|
|98C	|20s|
|57-63C	|15s|
|72C |	30s|
|Cycles |	30|
|72C |30s|
|4C	|Forever| 

### Notes: 
* Analyse MC plots for each targeton to select optimal temperature and cycle number 

## Step THREE: True ILL1 qPCR
 	
For each 50ul Reaction:

| Reagent | Volume |
| :--- | :--- |
|HDR Library (Maxiprep)|	25ng|
|F ILL1 Primer (10uM)|	1.5ul|
|R ILL1 Primer (10uM)|	1.5ul|
|Kapa Hifi MM	|25ul|
|20x Evagreen |	2.5ul|
|100x ROX|	0.5ul|
|H2O| 	To 50ul |

### Notes: 
* Make up enough master mix for 6 reactions (4x 50ul reactions + another 50uL reaction with Evagreen to trace the amplification with qPCR machine).


### PCR Conditions:

| Temperature       | Duration           |
| ------------- |:-------------:|
|95C|	3min|
|98C	|20s|
|__C	|15s|
|72C |	30s|
|Cycles |	___|
|72C |	30s|
|4C|	Forever |

### Notes: 
* Analyse MC plots for each targeton to select optimal temperature and cycle number 


## STEP FOUR: PCR purification 

1. Merge 4x reactions into one Eppendorf tube.
2. [PCR purify](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/PCR_purification.md) the reaction through Qiagen PCR column (not minelute).
  * Elute DNA in 30ul of EB buffer for each column
3. Nanodrop.
4. Qubit. dsDNA HS kit.
5. Dilute samples to ~ 10ng/ul.
6. Nanodrop.
7. Qubit. 


## STEP FIVE: ILL2-PCR

### Background:
* Addition of Illumina adapters

For each 100ul Reaction:

| Reagent | Volume |
| :--- | :--- |
|HDR Library (ILL1 purified)|	50ng|
|F and R ILL2 Primers (10uM)|	6ul |
|Kapa Hifi MM	|50ul|
|H2O| 	To 100ul | 


### PCR Conditions:

| Temperature       | Duration           |
| ------------- |:-------------:|
|95C	|3min|
|98C	|20s|
|59C	|15s|
|72C |	30s|
|Cycles |	7|
|72C| 	30s|
|4C	|Forever |


## STEP SIX: Ampure bead purification and pooling
1. [Bead purify](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/bead_purification.md) ILL2 pcr product on magnetic rack
2. Qubit.
3. Normalise samples to 20ul of 10ng/ul in EB Buffer
4. Use 10uL of the diluted library for 1% gel electrophoresis. 
  * Should see a single band at ~400bp
  * Each library should have similar intensity if your Qubit was done correctly.
5.	If the libraries look good on the gel, pool 5uL of each diluted library into one tube. (10ng/uL pooled library)
6.	Submit for sequencing or proceed to KAPA library quantification (50,000X and 500,000X dilution on the 10ng/uL pooled library)

# NGS Library preparation (gDNA)
## Version 4, 260624, Malin Andersson (reviewed by Sebastian Gerety)

## Background:
This protocol is for pooled library preparation for Novaseq 6000 using gDNA extracted following an SGE screen. 
In the following protocol, “replicate” refers to each replicate at each timepoint in a targeton (e.g. 3 replicates per time point, with 4 timepoints = 12 replicates in a targeton). 

At least one of the two ILL0 primers should be purposely designed outside the HDR homology arm. Thus, ILL0 should not amplify the HDR plasmid, if any is present in the genomic prep.

It is possible to adapt this for semi-high-throughput using the Qiagen vacuum system and 96w PCR purification kits. If you are working with 48+ samples, this will be much faster than batches of spin columns. However, the yield might be lower as not all elution can be recovered from the plate (generally ~75% of elution can be recovered). The cost of the plates also does not make them suitable for lower sample numbers.

Previous optimisation of this protocol found that for various genomic loci, KAPA HiFi gave better amplification than NEB Q5. The maximum amount of gDNA one can add per reaction is dependent on the primer pair. Most primer pairs can tolerate 1μg to 2μg gDNA per 100μL reaction, but this will vary with product length.


## Step 1: ILL0 optimisation qPCR

### Materials:
* gDNA
* Primer F: (short primer, 18-30bp) [target-specific sequence]
* Primer R: (short primer, 18-30bp) [target-specific sequence]
* KAPA HiFi Hot Start Ready Mix (Roche, KK2602)
* EvaGreen 20X (Biotium, 31000BT)
* ROX 50uM (Thermo Fisher, 13273189)
* Nuclease-free water

### qPCR:
Use qPCR to find out the primer annealing temperature and optimal PCR cycle number. Overcycling can cause non-specific amplification, and skew representation of variants.

*Note: If you are using your SGE samples, be mindful of the amount of gDNA available, as you’ll need 4μg for the actual PCR later!*

1. Prepare 7 reactions/primer pair according to the table below
2. Run the program below for 35 cycles
3. Using the multi-component plot, determine the optimal annealing temperature and cycle number

|        | For each 20μL reaction           |For 7 reactions|
| ------------- |:-------------:|:-------------:|
|KAPA HiFi Ready Mix |10μL|70μL|
|10uM Primer F|0.6μL|4.2μL|
|10uM Primer R|0.6μL|4.2μL|
|gDNA|200ng|1.4μg|
|EvaGreen 20X|1μL|7μL|
|Rox 100X|0.2μL|1.4μL|
|Add water to|20μL|140μL|

#### PCR Conditions: 

| Temperature       | Duration           |
| ------------- |:-------------:|
|95C| 3min|
|98C| 20sec|
|grad*|15sec|
|72C| 1min 30sec**|
|Back to (2)|35 cycles|
|72C| 1min 30sec**|
|4C| forever|
 
*Recommended range: 56-66C as a starting point, in 2C steps (6 reactions/sample)

**If the amplicon is >1500bp, increase the extension time to 2 minutes.

### Determining annealing temperature and cycle number
Use the multicomponent plot to determine the optimal temperature and cycle number. In the image below, each green trace represents an annealing temperature for a specific primer pair. The trace furthest to the left has the highest fluorescence at the start of the plateau (“best DNA yield”). To determine the cycle number, we want to pick it as low as possible in the linear range (see reference image further down) or at the end of the exponential (aka log-linear) range. The estimate is 20-25 cycles for amplification.


### Gel electrophoresis: (Optional)
Run 1% TAE gel (120V, 45 min) to ensure a specific DNA band is present at the correct size.


## Step 2: ILL0
(Amplicon size 1000-2000 bp)

### Materials:
* gDNA
* Primer F: (short primer, 18-30bp) [target-specific sequence]
* Primer R: (short primer, 18-30bp) [target-specific sequence]
* KAPA HiFi Hot Start Ready Mix (Roche, KK2602)
* Nuclease-free water

### PCR: 
1. For each replicate, prepare the reaction mix as outlined in the table below.
2. Split the 200μL into 4x 50μL reactions in a 96wp PCR plate.
3. Run the PCR program as outlined below.
4. Once the program is completed, store at -20C or proceed to purification steps.

| Reagent     | 4x 50μL reactions          |
| ------------- |:-------------:|
|KAPA HiFi Ready Mix|100μL|
|10uM Primer F|6μL|
|10uM Primer R|6μL|
|gDNA|4μg|
|Add water to|200μL|


#### ILL0 PCR conditions: 
| Temperature       | Duration           |
| ------------- |:-------------:|
|95C| 3min|
|98C| 20sec|
|Pre-optimised temp*|15sec|
|72C| 1min 30sec**|
|Back to (2)|X* cycles|
|72C| 1m 30 sec**|
|4C| forever|

*Optimised in step 1

**If the amplicon is >1500bp, increase the extension time to 2 minutes.


### Purifications: 

#### Materials:
* ExoI + buffer (NEB, M0293L)
* NaOAc 3M (Sodium acetate) (Invitrogen, AM9740)
* Qiaquick PCR Purification Kit (Qiagen, 28106, UBW order code: BCPU0005)
* Qiagen Minelute PCR Purification Kit (Qiagen, 28006, UBW order code: BCMI0004)
* Qiaquick 96 well PCR Purification Kit (vacuum system) (Qiagen, 28181)

#### Method: 
1. Pool the 4 reactions together.
2. Purify PCR products using Qiagen PCR Purification kit (standard column or 96 well plate kit). Add 3M NaOAC to adjust the pH. Usually, 10% of the reaction volume is good enough (ie, 20μL of NaOAc per 200μL PCR reaction)
  * For column purification, use 2 columns per replicate. Follow the Qiagen protocol for PCR purification. Elute each column in 40μL EB buffer, then pool the two (total volume 80μL).
  * For 96 well plate purification, elute in 80μL EB buffer. Follow the Qiagen protocol for the Qiaquick 96 well kit.
3. Nanodrop. Expected yield: 50-200ng/μL. 260/280: 1.70-1.90. 260/230: >2.00
4. Gel electrophoresis, 1% agarose, TAE, 120 V, 45 minutes. You should see one distinct band.


## Step 3: ILL1 optimisation qPCR

### Materials
* Primers (IDT or Sigma, PAGE purified)
* Reverse/i7 index adapter: GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT[Target-specific sequence]
* Forward/i5 index adapter: ACACTCTTTCCCTACACGACGCTCTTCCGATCT[Target-specific sequence]
* KAPA HiFi Hot Start Ready Mix (Roche, KK2602)
* EvaGreen 20X (Biotium, 31000BT)
* ROX 50uM (Thermo Fisher, 13273189)
* Nuclease-free water
* Purified ILL0 product

### qPCR: 
1. Prepare 7 reactions/primer pair according to the table below
2. Run the program below for 25 cycles
3. Using the multi-component plot, determine the optimal annealing temperature and cycle number as performed above for ILL0.

|        | For each 20μL reaction           |For 7 reactions|
| ------------- |:-------------:|:-------------:|
|KAPA HiFi Ready Mix |10μL|70μL|
|10uM Primer F|0.6μL|4.2μL|
|10uM Primer R|0.6μL|4.2μL|
|purified ILL0 DNA |10ng|70ng|
|EvaGreen 20X|1μL|7μL|
|Rox 100X|0.2μL|1.4μL|
|Add water to|20μL|140μL|

#### PCR Conditions: 

| Temperature       | Duration           |
| ------------- |:-------------:|
|95C| 3min|
|98C| 20sec|
|grad*|15sec|
|72C| 30sec|
|Back to (2)|25 cycles|
|72C| 30sec|
|4C| forever|

*Recommended range: 56-66C as a starting point, in 2C steps (6 reactions/sample).*

#### Notes: 
* The estimate for ILL1 cycle number is 6-12.


## Step 4: ILL1
(Amplicon size 250-350bp)

### Materials:
* Primers (IDT or Sigma, PAGE purified)
* Reverse/i7 index adapter: GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT[Target-specific sequence]
* Forward/i5 index adapter:
* ACACTCTTTCCCTACACGACGCTCTTCCGATCT[Target-specific sequence]
* KAPA HiFi Hot Start Ready Mix (Roche, KK2602)
* Nuclease-free water
* Purified ILL0 product

### PCR:
1. For each replicate, prepare the reaction mix as outlined in the table below.
2. Split the 200μL into 4x 50μL reactions in a 96wp PCR plate.
3. Run the PCR program as outlined below.
4. Once the program is completed, store at -20C or proceed to purification steps.

| Reagent     | 4x 50μL reactions          |
| ------------- |:-------------:|
|KAPA HiFi Ready Mix|100μL|
|10uM Primer F|6μL|
|10uM Primer R|6μL|
|purified ILL0 DNA|100ng|
|Add water to|200μL|

#### ILL1 PCR conditions: 
| Temperature       | Duration           |
| ------------- |:-------------:|
|95C| 3min|
|98C| 20sec|
|Pre-optimised temp*|15sec|
|72C| 30sec|
|Back to (2)|X* cycles|
|72C| 30sec|
|4C| forever|

*Optimised in step 3


### Purifications: 

#### Materials:
* ExoI + buffer (NEB, M0293L)
* NaOAc 3M (Sodium acetate) (Invitrogen, AM9740)
* Qiaquick PCR Purification Kit (Qiagen, 28106, UBW order code: BCPU0005)
* Qiagen Minelute PCR Purification Kit (Qiagen, 28006, UBW order code: BCMI0004)
* Qiaquick 96 well PCR Purification Kit (vacuum system) (Qiagen, 28181)
* Qubit dsDNA High Sense kit (Invitrogen, Q32854)

#### Method:
1. Pool the 4 reactions together.
2. Purify PCR products using Qiagen PCR Purification kit (standard column or 96 well plate kit). Add 3M NaOAC to adjust the pH. Usually, 10% of the reaction volume is good enough (ie, 20μL of NaOAc per 200μL PCR reaction)
  * For column purification, use 2 columns per replicate. Follow the Qiagen protocol for PCR purification. Elute each column in 40μL EB buffer.
  * For 96 well plate purification, elute in 80μL EB buffer. Follow the Qiagen protocol for the Qiaquick 96 well kit.
3. Digest the purified DNA with ExoI. Digest ssDNA (e.g. oligos) at 37C for 20 min, then inactivate ExoI by heating to 80C for 20 min. For 80μL elutions, make a 60μL reaction (6μL ExoI buffer, 2μL ExoI, 52μL purified DNA).
4. Purify ExoI digested DNA using Qiagen kit (Minelute column kit or 96 well plate kit). Add 3M NaOAC to adjust the pH. Usually, 10% of the reaction volume is good enough (ie, 6μL of NaOAc per 60μL ExoI digested DNA).
  * For column purification, use 1 MinElute column per sample. Elute each column in 15μL EB buffer.
  * For plate purification, elute in 80μL EB buffer. Transfer 60μL from the collection tubes to a 96w PCR plate for storage.
5. Measure the DNA concentration using Qubit High Sense dsDNA kit. Expected yield: 1-3μg with the plate kit. Might be higher if using MinElute columns.


## Step 5: ILL2

### Materials:
* Primers: NEBNext® Multiplex Oligos for Illumina® (96 Unique Dual Index Primer Pairs) (NEB, E6440S)
* KAPA HiFi Hot Start Ready Mix (Roche, KK2602)
* Nuclease-free water
* Purified ILL1 product

*Note: Each replicate included in the pool must have a unique index! Do not use the same well multiple times in the same pool.*

### PCR:
1. Prepare a 50μL reaction for each sample in a 96 well PCR plate.
2. Run the PCR with the below conditions.
3. Purify the DNA using Ampure beads directly in the 96 well plate using the protocol below.

| Reagent     | 1x 50μL reaction|
| ------------- |:-------------:|
|KAPA HiFi Ready Mix|25μL|
|NEB UDI primers|3μL|
|ILL1 DNA|25ng|
|Add water to|50μL|

#### ILL2 PCR conditions: 
| Temperature       | Duration           |
| ------------- |:-------------:|
|95C| 3min|
|98C| 20sec|
|59C|15sec|
|72C| 30sec|
|Back to (2)|7 cycles|
|72C| 30sec|
|4C| forever|


## Step 6: Bead purification

### Materials:
* Ampure XP beads (Beckman Coulter, A63881)
* 96 well magnet
* Ethanol, absolute

### Method: 
1. Wipe the bench and pipettes thoroughly with RNAse away.
2. Warm the Ampure XP bead to room temp (put at RT for ~30min)
3. Vortex Ampure XP bead to resuspend.
4. Add 45μl (0.9X) resuspended beads to the 50μL PCR reaction. Mix well by pipetting up and down at least 10 times. Be careful to expel all the liquid out of the tip during the last mix.
5. Incubate samples on the benchtop for 15 minutes at room temperature. Prepare enough fresh 80% ethanol for the wash later.
6. Place the plate on a magnetic stand to separate the beads from the supernatant.
7. After 5 minutes, carefully remove 20μL of the supernatant from the top. The removed supernatant should look clear. Seal the plate and spin the plate at 300rcf for 20 sec. Most of the bead should be at the bottom. Remove the seal and place on the magnetic rack.
8. After 5 minutes, remove the remaining supernatant. Seal and quick-spin the plate. Place the plate on the magnetic rack again. Remove all remaining solution in the wells using a long reach p10 pipette tip (you do not have to wait for 5 min after putting the plate on the magnetic rack this time - the beads should stick on the wall instantly).
9. Add 200μl of 80% freshly prepared ethanol to the plate while in the magnetic stand. Incubate at room temperature for 30 seconds, and then carefully remove and discard the supernatant. Be careful not to disturb the beads that contain DNA targets.
10. Repeat Step 9 once for a total of two washes. 
11. Seal the plate and spin the plate at 300rcf for 20 sec. Remove the seal and place the plate on the magnetic rack. Remove all traces of ethanol with a long reach p10 pipette tip. Spin the plate again if necessary.
12. Air dry the beads for up to 5 minutes while the plate is on the magnetic stand. Start counting the 5 min after the first plate spinning step from Step 11.
*Caution: Do not over-dry the beads. This may result in lower recovery of DNA. Elute the samples when the beads are still dark brown and glossy looking. When the beads turn lighter brown and start to crack they are too dry.*
13. Remove the tube/plate from the magnetic stand. Elute the DNA target from the beads by adding 33 μl of Qiagen EB buffer. Keep a separate EB buffer reagent that you use only for post-bead NGS prep, to avoid RNA contamination.
14. Mix well by pipetting up and down 10 times. Incubate for 10 minutes at room temperature.
15. Spin the plate at 300rcf for 1 min (to remove air bubbles at the bottom of the well).
16. Place the tube/plate on the magnetic stand. After 5 minutes (or when the solution is clear, it only takes 1 min usually), transfer 30 μl to a new plate (ensure that you do not take out any beads, if possible) and store at –20°C or proceed to Qubit.


## Step 7: QC

### Materials:
* Purified ILL2 product
* Qubit dsDNA High Sense kit (Invitrogen, Q32854)
* DNA ladder, diluted in water to 50ng/μL
* Qubit

### Method: 
1. Prepare a master mix of Buffer+dye. You need (number of samples) + 5 reactions (i.e. if you have 10 samples, you need to prepare 15 reactions master mix). Each reaction consists of 200μL of Qubit HS buffer and 1μL of Qubit HS dye.
2. Label the qubit assay tube. Always label the cap only. Do not label anything at the side of the tubes. Pipette 190μL of the master mix to Tube 1 and Tube 2 (these are for the standard provided by the kit). Pipette 199μL of master mix to Tube 3 (this is for the DNA ladder which has been diluted to 50ng/μL, I used this as additional control). Pipette 199μL of master mix to the other tubes for your samples.
3. Add 10μL of each standard to Tube 1 and Tube 2. Add 1μL of diluted DNA ladder to Tube 3. Add 1μL of your sample to other tubes.
4. Vortex the tubes and quick spin. Incubate for at least 10min at RT (preferably in dark. For example, just leave it in the centrifuge) before assayed by the Qubit reader. 

### Normalisation:
After finding out the concentrations, for each library, prepare 20μL of 10ng/μL library (i.e. 200ng in 20μL). Use Qiagen EB buffer to do the dilution.

### Gel electrophoresis: 
Use 10μL of the diluted library for 2% gel electrophoresis. 120V, 45min, TAE.

*Note: You should see a single band at ~400bp. Each library should have similar intensity if your Qubit was done correctly.*


## Step 8: Pooling
If the libraries look good on the gel, pool 5μL of each diluted library into one tube (10ng/μL library = 50ng).

### Submit the sample for sequencing:
* If you are submitting a pool using the SGE Novaseq study pipeline, label the pool with the barcode provided by DNAP and drop it off as per their instructions.
* If you are submitting the library for walk-up (or if you have been asked to quantify it for another reason), quantify the pooled library using KAPA Quant qPCR kit (Roche, for Illumina sequencing - follow the provided protocol). Then follow the instructions for walk-up drop-off.

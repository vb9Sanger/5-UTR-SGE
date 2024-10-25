# TWIST Oligo Amplification 

## Background:

I had two TWIST oligo pools each containing five targets across five genes. 

Their concentrations came at 179ng/ul and 225ng/ul. Tubes were dissolved in 50ul EB Buffer to ensure there was enough sample to play with. (Yielded 3.58ng/ul and 4.50ng/ul, respectively). 


## STEP ONE: Practice qPCR on WT constructs

1.	Perform qPCRs to test TWIST amplication primers (amplify ~300bp target region)

For each 15ul reaction:

| Reagent        | Volume           |
| ------------- |:-------------:|
| ssDNA (WT Construct)	|10ng |
| F TWIST Primer (10uM)|	0.45ul|
|R TWIST Primer (10uM)	|0.45ul|
|Kapa Hifi MM	|7.5ul|
|Evagreen (20X)|	0.75ul|
|ROX (100X)|	0.15ul |
|H2O |	To 15ul| 

### Notes:
* This is done using WT constructs that will not be used
* Ran 4X reactions (15ul each) 
* Remember to include a negative control
* Spin plate at 2400rpm for 1min prior to PCR

### Gradient PCR Conditions: 
(Modified version of “kapa grad ill0 grad”)

| Temperature       | Duration           |
| ------------- |:-------------:|
|95C|	3min|
|98C	|20sec|
|55-63C*|	15sec|
|72C |	30sec|
|Cycles |	35|
|72C |	30s|
|4C	|Forever|  

### Notes:
* I performed four reactions across four temperatures between 55C and 63C. (I selected 55C, 57C, 59C and 61C)
* Hong Kee selected 57C, 59C, 61C, 63C. Retrospectively, I would select the latter.
* When setting up qPCR plates:
    * Avoid excess freeze thaw of ROX
    * Set up on ice
    * Remember two adjacent wells have the same temp, so load every other well (e.g. B2, B4, B6, B8...)
    * Avoid bubbles when dyes are being used 

2. Assess MC plot to confirm amplification 

## STEP TWO: Optimisation qPCRs on oligo pools 

1.	Perform qPCRs to optimize annealing temperature and cycle number for each target

For each 15ul reaction:


| Reagent        | Volume (1X)         |
| ------------- |:-------------:|
|ssDNA (WT Construct)	|10ng|
|F TWIST Primer (10uM)	|0.45ul|
|R TWIST Primer (10uM)	|0.45ul|
|Kapa Hifi MM	|7.5ul|
|Evagreen (20X)	|0.75ul|
|ROX (100X)	|0.15ul |
|H2O |	To 15ul|

### Notes:
* Ran 4X reactions (15ul each)
* Remember to include a negative control
* Spin plate at 2400rpm for 1min

### Gradient PCR Conditions: 
(Modified version of “kapa grad ill0 grad”)

| Temperature       | Duration           |
| ------------- |:-------------:|
|95C|	3min|
|98C	|20sec|
|55-63C*|	15sec|
|72C |	30sec|
|Cycles |	35|
|72C |	30s|
|4C	|Forever|  

### Notes:
* I performed four reactions across four temperatures between 55C and 63C. (I selected 55C, 57C, 59C and 61C)
* Hong Kee selected 57C, 59C, 61C, 63C. Retrospectively, I would select the latter.

2. Assess MC plot to confirm amplification 

### Notes:
* I selected optimal annealing temp by assessing which sample (if any) demonstrated specific product at an early cycle number. If no difference, selected 60C.
* I selected optimal cycle number by looking at the exponential phase of the corresponding multicomponent plot. (I then subtracted two from this to account for increase in DNA for 50ul reactions performed for subsequent true amplification) **Retrospectively unnecessary**
* For two of my targets, my negative control amplified. A gel confirmed likely primer dimers, necessitating gel extraction to purify these samples following true amplification.
* Due to poor amplification of NIPBL, this qPCR was redone with twice the amount of DNA
* Targets where negative control samples amplified were run again at higher temperatures. (65C was selected for these samples in step 3). 


## STEP THREE: True qPCRs to amplify HDR oligos 

To ensure enough amplified product, 15ul reactions were scaled up to 50ul. For each target region, I ran 2X 50ul reactions. 

For each 50ul reaction:


| Reagent        | Volume    (1X)       |
| ------------- |:-------------:|
|ssDNA (TWIST oligo pool)	|4ng|
|F TWIST Primer (10uM)	|1.5ul|
|R TWIST Primer (10uM)	|1.5ul|
|Kapa Hifi MM	|25ul|
|H2O 	|To 50ul| 

### Notes:
* Ran 2X reactions (50ul each)
* Remember to include a negative control
* Spin plate at 2400rpm for 1min

### Gradient PCR Conditions: 
(Modified version of “kapa grad ill0 grad”)

| Temperature       | Duration           |
| ------------- |:-------------:|
|95C|	3min|
|98C	|20sec|
|__C*|	15sec|
|72C |	30sec|
|Cycles |	__*|
|72C |	30s|
|4C	|Forever|  

*Dependent on step 2

### Notes:
* Following true amplification, samples were bioanalysed (to ensure presence of product at expected size) then stored at -80C

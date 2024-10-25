# WT Backbone Linearisation 


## Background:
* Linearised backbone was initially amplified from 50ng of WT Vector. Gibson Assembly of this ligation product demonstrated contamination of uncut vector. Subsequent gradient PCRs testing 0.5ng – 50ng of starting WT Vector were therefore conducted.
* The above starting concentrations were also tested at both 30 and 35 cycles.
* Running all samples on a gel to visualise corresponding amounts of uncut vector yielded the following optimised conditions:
  * 0.25ng-1.0ng of starting WT Vector
  * 30 cycles
* 8 identical 50ul reactions were run in order to generate a comfortable excess of linearised backbone. 

1. conduct PCR to amplify backbone

For each 50ul reaction:

| Reagent        | Volume           |
| ------------- |:-------------:|
|ssDNA (WT Construct)|	0.25ng – 1.0ng|
|F HDR WT Primer (10uM)	|1.5ul|
|R HDR WT Primer (10uM)	|1.5ul|
|Kapa Hifi MM	|25ul|
|H2O| 	To 50ul| 

### Notes:
* Remember to include a negative control 


## PCR Conditions: 

| Temperature       | Duration           |
| ------------- |:-------------:|
|95C	|3min|
|98C	|20s|
|63C	|15s|
|72C |	30s|
|Cycles |	30|
|72C| 	30s|
|4C	|Forever| 

2. Add 1ul of dpnI to each sample (to digest uncut plasmid)
3. Keep all samples at 37C for 1hr and 80C for 20min
4. Run a gel (25ul per well)
5. Excise bands
6. [Gel purify](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/gel_purification.md)
7. Nanodrop

### Notes: 
* Any samples deemed poor quality following purification were tossed
* Remaining samples were merged into one and nanodropped once more 
 
# TWIST Oligo Purification 

## Background:
Whether or not amplified oligos were purified by PCR purification or [gel purification](https://github.com/vb9Sanger/5-UTR-SGE/blob/main/WetLab_Protocols/gel_purification.md) was determined based off optimisation experiments and corresponding presence of primer dimers (or lack thereof). If there was concern of primer dimers, oligos were purified by gel purification. Otherwise, oligos were purified by PCR purification (according to manufacturer’s protocol).

## Notes:
* All samples were eluted in 12ul of EB buffer
* Ideally, should yield between 20-50ng/ul

# NEBuilder Assembly of WT Vector and Amplified Oligos 	

## Background:
Selection of NEBuilder (vs Gibson), starting amount of ligation mix added to each transformation and volume of NEB-5 alpha cells used in each transformation were the result of a series of optimisation experiments. 
	NEBioCalculator was used to calculate starting amounts of WT Vector and insert for all ligations, where a 1:2 ratio of vector to insert was recommended for assembly of 2 fragments.

## Ligations:

For each 20ul reaction:

| Reagent        | Volume           |
| ------------- |:-------------:|
|Linearised WT backbone	|100ng*|
|Purified insert |	*|
|NEBuilder |	10ul|
|H2O| 	To 20ul| 
*Corresponding volume calculated with NEBioCalculator 

## PCR Conditions: 

| Temperature       | Duration           |
| ------------- |:-------------:|
|50C	|15min|
|4C	|Forever |


## Transformations:

1. Thaw a tube of NEB 10-alpha E.coli cells on ice
2. Create the following transformations:
* 8ul of vector + insert ligation in 200ul of NEB-5 alpha cells (2X)
* 2ul of vector + insert ligation in 50ul of NEB-5 alpha cells
* 2ul of vector alone + NEBuilder in 50ul of NEB-5 alpha cells
* 2ul of vector alone without NEBuilder in 50ul of NEB-5 alpha cells
3. Flick 4-5x to mix
4. Place mix on ice for 30min
5. Heat shock at 42C for exactly 30sec
6. Place on ice for 2min
7. Pipette 950ul of SOC into mix
8. Place at 37C for 60min and shake at 650rpm
9. Warm plates to 37C
10. Mix cells thoroughly by flicking / inverting
11. Following incubation, plate out 20ul of all transformations
12. Incubate plates at 37C overnight
13. Merge remaining volumes from 200ul vector + insert transformations
14. Add this to 150ul of fresh LB+Amp and leave to shake (220rpm)
15. Incubate at 37C for 14h (overnight)
16. The next day: create cell pellets and calculate CFUs for all transformations using the following formula:

***count of colonies per plate x  (total volume of transformation)/(volume of transformation plated)=CFUs***
	
 ### Notes:
* For all transformations, the following controls were included:
  * Vector alone + NEBuilder
  * Vector alone without NEBuilder
  * (Insert alone + NEBuilder was also included where there was risk of unwanted product)
* Across all target regions, we aimed for at least 100X coverage (i.e. If there were 1000 oligos in the library, we aimed for at least 100,000 CFUs)




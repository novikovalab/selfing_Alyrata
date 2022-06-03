# selfing_Alyrata
Scripts associated to the paper: 'Selfing Siberian Arabidopsis lyrata is a progenitor of the allopolyploid Arabidopsis kamchatica'. 

PiW: Calculation of nucleotide diversity in 10 Kb windows.


CoGe processing: Working with output of Synmap = synteny analysis of online platform CoGe: https://genomevolution.org/coge/SynMap.pl

  step 1: transforming output to separate files containing start and end of each syntenic block, performed on command line.
	
  step 2: merging starts and ends to create following file format: chrA, startA, endA, chrB, startB, endB. A and B depict two different genome assemblies.
	
  step 3: Visualizing the result in R using the library circlize. 
	
  
 

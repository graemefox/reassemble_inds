# reassemble_inds


Takes the output from Msatallele (the 1 column per locus) file, checks for inconsistent genotypes 
of the same sample/genotype combination and produces various output files for genepop, structure, 
popgenreport, demetics etc. Specific to my samples in that (it deals with inconsistent naming 
which is hard coded to my samples, also the lat/long of my samples are hard coded.)

The renaming section could easily just be commented out without too much trouble. The lat/long is required
for some of the output formats so will need adjusted if attempted with different data.

#!/bin/bash

################################
##### Run the simulations ######
################################



###################################### Launch the simulations :
# First modify the parameters : 
# for the mutations, admixture proportion and population size in the Positive_model.sim
########Average value of the selection coefficient for both population mutation line 12 and 17
########The number and proportion of mutations on sex chromosome and within each population line 97 to 116
######## Admixture proportion, population size and sex-ratio line 131


# Launching SLiM simulations
./SLiM_3.5/slim Positive_model.slim

echo The simulations are done !


###################################### Analysing ancestry along the genome :


python3 Load_ancestry.py

echo done with python

paste breaks.txt ancestry.txt > Output_Ancestry.txt
rm ancestry.txt                                           
rm breaks.txt
rm Slim_output.tree


echo The ancestry is analysed


############## Script to have the ancestry along with the position and the recombination rate




Rscript Averaging_no_overlap.R
rm Output_Ancestry.txt


#!/bin/bash

################################
##### Run the simulations ######
################################


##################################### Generating the input file:
# First modify the values in the R script (Creating_BDMI_script.R) to have the correct scenario (number of BDMI and mean of the fitnes reduction coefficient)



# run the R script to prepare the SLiM script
# first choose the mean value of the distribution of the fitness impact in Creating_BDMI_script.R
Rscript Creating_BDMI_script.R


## Modifying the files to make the slim script :
sed -i 's/*/"/g' Mutations_def
sed -i 's/*/"/g' Genomic_element
sed -i 's/vide/    /g' Mutations_def
sed -i 's/vide/    /g' Genomic_element
sed -i 's/vide/    /g' Fitness
sed -i 's/vide/    /g' add_mutations
cat Debut.txt Mutations_def Genomic_element Creating_populations.txt add_mutations Demographic_model.txt Fitness > Model_BDMI.slim
rm Genomic_element
rm Fitness
rm add_mutations
rm Mutations_def

echo The script is ready

###################################### Launch the simulations :
# First modify the parameters :  % of admixture, population size and sex-ratio in the file Demographic_model.txt
#line 13 and 14




# Launching SLiM simulations
./SLiM_3.5/slim Model_BDMI.slim
#slim Model_BDMI.slim
rm Model_BDMI.slim

mv Slim_output_BDMI.tree Slim_output.tree

#echo The simulations are done !



###################################### Analysing ancestry along the genome :
python3 Load_ancestry_females.py
paste breaks.txt ancestry.txt > Output_Ancestry.txt
rm breaks.txt
rm ancestry.txt
rm Slim_output.tree

echo The ancestry is analysed


############## Script to have the ancestry along with the position and the reco$
Rscript Averaging_no_overlap.R
rm Output_Ancestry.txt

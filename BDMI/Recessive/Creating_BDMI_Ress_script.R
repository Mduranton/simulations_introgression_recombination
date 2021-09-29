#######

### first we have to choose the number of mutations:

Nb_mut=100

###### choose the mean of the distribution of the impact on fitness:

mean = 0.1



### First we have to create the neutral mutation too mark the sex chromosome
Ligne1 <- paste("initializeMutationType(*m1*, 1, *f*, 0);", sep="")
Ligne2 <- paste("m1.mutationStackPolicy = *l*;", sep="")
Ligne3 <- "vide"
m <- rbind(Ligne1, Ligne2, Ligne3)
Definition_of_mutations <- m 


#### now we are defining the name of the mutation within p1 and p2

number <- 2:((Nb_mut*2)+1)
Mutations <- NULL
for(i in 1:(Nb_mut*2)){
  Mutations <- c(Mutations, paste("m",number[i], sep=""))  
}




for(i in 1:length(Mutations)){
  m <- NULL
  Ligne1 <- paste("initializeMutationType(*", Mutations[i], "*, 0.5, *f*, 0);", sep="")
  Ligne2 <- paste(Mutations[i], ".mutationStackPolicy = *l*;", sep="")
  Ligne3 <- "vide"
  m <- rbind(Ligne1, Ligne2, Ligne3)
  Definition_of_mutations <- rbind(Definition_of_mutations, m)
}

Colone1 <- rep("vide", nrow(Definition_of_mutations))
Definition_of_mutations <- cbind(Colone1, Definition_of_mutations)  

write.table (Definition_of_mutations, "Mutations_def", row.names = F, col.names = F, quote=F)





############## Now we have to define the genomic elements with the different mutations we have just created
###### Sex K
deux <- NULL
trois <- "c(1," 

for(t in 1:(length(Mutations)-1)){
  deux <- paste(deux, paste(Mutations[t], ",", sep=""), sep=" ")
  trois <- paste(trois, "1,")
}

deux <- paste(deux, paste(as.character(Mutations[length(Mutations)]), "),", sep=""), sep=" ")
trois <- paste(trois, "1));", sep=" ")

g1 <- paste("vide","initializeGenomicElementType(*g1*, c(m1,", deux, trois, sep=" ")

deux <- NULL
trois <- "c(" 

for(t in 1:(length(Mutations)-1)){
  deux <- paste(deux, paste(Mutations[t], ",", sep=""), sep=" ")
  trois <- paste(trois, "1,")
}

deux <- paste(deux, paste(as.character(Mutations[length(Mutations)]), "),", sep=""), sep=" ")
trois <- paste(trois, "1));", sep=" ")

g2 <- paste("vide","initializeGenomicElementType(*g2*, c(", deux, trois, sep=" ")


genomic_element <- rbind(g1, g2)


write.table (genomic_element, "Genomic_element", row.names = F, col.names = F, quote=F)


############# fitness definition

#### Using an exponential distribution
Selective_coeff <- rexp(5000,1/mean)

Fitness <- NULL


for (i in seq(1,length(Mutations),2)){
  Value <- sample(Selective_coeff, 1)
  Ligne1 <- paste("fitness(", Mutations[i], "){", sep="")
  Ligne2 <- paste("vide", paste("Nb_", Mutations[i], "=asInteger(genome1.countOfMutationsOfType(", Mutations[i], ")) ", "+ asInteger(genome2.countOfMutationsOfType(", Mutations[i], "));" ,sep=""),sep=" ")
  Ligne3 <- paste("vide", paste("Nb_", Mutations[(i+1)], "=asInteger(genome1.countOfMutationsOfType(", Mutations[(i+1)], ")) ", "+ asInteger(genome2.countOfMutationsOfType(", Mutations[(i+1)], "));" ,sep=""),sep=" ")
  #Ligne4 <- paste("vide", paste("Tot = Nb_", Mutations[i], " + Nb_", Mutations[(i+1)], ";", sep=""),sep=" ")
  #Ligne5 <- "vide if (Tot == 4)"
  Ligne4 <- paste("vide", paste("if ((Nb_", Mutations[i], "==2) & (Nb_", Mutations[(i+1)], "==2))", sep=""))
  Ligne6 <- paste("vide", "vide", paste("return (1-", Value, ");", sep=""), sep=" ")
  #Ligne7 <- "vide else if (Tot == 3)"
  #Ligne7 <- paste("vide", paste("else if ((Nb_", Mutations[i], "==1) & (Nb_", Mutations[(i+1)], "==2))", sep=""))
  #Ligne8 <- paste("vide", "vide", paste("return (1 - (1/2*", Value, "));", sep=""), sep=" ")
  #Ligne9 <- paste("vide", paste("else if ((Nb_", Mutations[i], "==2) & (Nb_", Mutations[(i+1)], "==1))", sep=""))
  #Ligne10 <- paste("vide", "vide", paste("return (1 - (1/2*", Value, "));", sep=""), sep=" ")
  #Ligne9 <- "vide else if (Tot == 2)"
  #Ligne11 <- paste("vide", paste("else if ((Nb_", Mutations[i], "==1) & (Nb_", Mutations[(i+1)], "==1))", sep="")) 
  #Ligne12 <- paste("vide", "vide", paste("return (1 - (1/4*", Value, "));", sep=""), sep=" ")
  Ligne13 <- "vide else"
  Ligne14 <- "vide vide return relFitness;"
  Ligne15 <- "}"
  Ligne16 <- "vide"
  Bloc <- rbind(Ligne1, Ligne2, Ligne3, Ligne4, Ligne6,  Ligne13, Ligne14, Ligne15, Ligne16)
  Fitness <- rbind(Fitness, Bloc)
}

write.table (Fitness, "Fitness", row.names = F, col.names = F, quote=F)




########### créer la partie qui permet d'introduire les mutations au hasard
##### une partie pour les mutations qui tombent sur les chromosomes sexuels et une
## pour le reste du génome

mutations_p1 <-as.character(Mutations[seq(2,length(Mutations),2)])
mutations_p2 <-as.character(Mutations[seq(1,length(Mutations),2)])


add_mutations <- NULL


for(i in 1:((length(Mutations)/2)*0.18)){
  num_p1 <- sample(length(mutations_p1),1)
  mut_p1 <- mutations_p1[num_p1]
  mutations_p1 <- mutations_p1[-num_p1]
  num_p2 <- sample(length(mutations_p2),1)
  mut_p2 <- mutations_p2[num_p2]
  mutations_p2 <- mutations_p2[-num_p2]
  p1 <- paste("GenomesP1SexMut.addNewDrawnMutation(", mut_p1, ",rdunif(1,0,22500000));", sep="")
  p2 <- paste("GenomesP2SexMut.addNewDrawnMutation(", mut_p2, ",rdunif(1,0,22500000));", sep="")
  add_mutations <- rbind(add_mutations, p1, p2)  
}


for (i in 1:length(mutations_p1)){
  num_p1 <- sample(length(mutations_p1),1)
  mut_p1 <- mutations_p1[num_p1]
  mutations_p1 <- mutations_p1[-num_p1]
  num_p2 <- sample(length(mutations_p2),1)
  mut_p2 <- mutations_p2[num_p2]
  mutations_p2 <- mutations_p2[-num_p2]
  p1 <- paste("Indp1.genomes.addNewDrawnMutation(", mut_p1, ",rdunif(1,22500002,119000003));", sep="")
  p2 <- paste("Indp2.genomes.addNewDrawnMutation(", mut_p2, ",rdunif(1,22500002,119000003));", sep="")
  add_mutations <- rbind(add_mutations, p1, p2)
}





Colone1 <- rep("vide", nrow(add_mutations))
add_mutations <- cbind(Colone1, add_mutations)

write.table (add_mutations, "add_mutations", row.names = F, col.names = F, quote=F)


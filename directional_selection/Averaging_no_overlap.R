
##### Load the information for all chromosomes

Recomb <- read.table("./Recomb_average.txt")

Ancestry <- read.table("./Output_Ancestry.txt")


options(scipen = 999)

###### modify the recombination file to have the begining and the end of each window

end <- Recomb[,1]
end <- end[-1]
end <- end-1
end <- c(end, (end[length(end)])+100000)
Recomb <- cbind(Recomb[,1], end, Recomb[,2])


#### same thing for the ancestry file to have the starting and ending position of each block
end2 <- Ancestry[,1]
end2 <- end2[-1]
end2 <- end2-1
end2 <- c(end2, 119000003)
Ancestry <- cbind(Ancestry[,1], end2, Ancestry[,2])
Ancestry[1,1] <- 1


Test <- Ancestry[,2]-Ancestry[,1]
Ancestry <- cbind(Ancestry, Test)
Ancestry <- Ancestry[Ancestry[,4]>0,1:3]

#### loop to separate every block that is overlaping two recombination windows
SUP <- NULL
limit <- 1



remooving_lines <- NULL


i=1
while(i!=nrow(Recomb)){
  while(Ancestry[limit,1]<Recomb[i,2] & Ancestry[limit,2]<Recomb[i,2]){limit=limit+1}
  Divided <- Ancestry[limit,]
  if(Divided[1]!=Recomb[i,2]){
    remooving_lines <- c(remooving_lines, limit)
    Newline <- c(Divided[1], Recomb[i,2], Divided[3])
    Newline2 <- c((Recomb[i,2]+1), Divided[2], Divided[3])
    if(Newline2[2]<=Newline2[1]){Newline2<-NULL}
    SUP <- rbind(SUP, Newline, Newline2 )
  }
  count <- 1
  while(SUP[nrow(SUP),2]>Recomb[i,2]+(100000*count)){
    Newline <- c(SUP[nrow(SUP),1], (Recomb[i,2]+(100000*count)) ,SUP[nrow(SUP),3])
    Newline2 <- c((Recomb[i,2]+((100000*count)+1)), SUP[nrow(SUP),2] ,SUP[nrow(SUP),3])
    SUP <- SUP[-nrow(SUP),]
    SUP <- rbind(SUP, Newline, Newline2)
    count <- count + 1
  }
  if(count==1){i=i+1}else{i=i+count}
}


SUP[nrow(SUP),2] <-Recomb[nrow(Recomb),2]



####### Test to verify that there are no window with a negative length
Verif <- SUP[,2]-SUP[,1]
SUP_test <- cbind(SUP, Verif)
SUP_test <- SUP_test[SUP_test[,4]>0,1:3]
if(nrow(SUP)!=nrow(SUP_test)){print("Something is wrong")}


##### add the new block and remoove the old blocks that were overlaping

Ancestry <- Ancestry[-remooving_lines,]

Ancestry <- rbind(SUP, Ancestry)
Ancestry <- Ancestry[order(Ancestry[,1]),]


##### verif 
Test <- is.element(Recomb[,2], Ancestry[,2])
Verif <- cbind(Recomb, Test)
Test <- Verif[Verif[,4]==1,]
if(nrow(Test)!=nrow(Recomb)){
  miss <- Verif[Verif[,4]==0,2]
  miss <- miss[-length(miss)]
  Test <- is.element(miss, Ancestry[,1])
  miss2 <- cbind(miss, Test)
  miss2 <- miss2[miss2[,2]==1,]
  if(length(miss)!=nrow(miss2)){print("Some border are missing")}
}



################
## first measuring the size of each block 


Size <- Ancestry[,2]-Ancestry[,1]
Size <- Size+1
Ancestry <- cbind(Ancestry, Size)

SAFE <- Ancestry

FINAL <- NULL

i=1
for(i in 1:nrow(Recomb)){
  TAB <- Ancestry[Ancestry[,1]>=Recomb[i,1],]
  if(length(TAB)>4){TAB <- TAB[TAB[,2]<=Recomb[i,2],]}
  if(is.array(TAB)){
    Line <- c(Recomb[i,], weighted.mean(TAB[,3], TAB[,4]))
    if(sum(TAB[,4])!=100000){
      #if(sum(TAB[,4])!=99999){print(paste("total length not good on line ", i, ":", sum(TAB[,4])))}
    }  
  }else{
    Line <- c(Recomb[i,], TAB[3])
  }
  FINAL <- rbind(FINAL,Line)
}


write.table(FINAL, "Position_recomb_ancestry.txt", col.names = F, row.names = F, quote=F)








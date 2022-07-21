#
# SELECTION OF INDIVIDUALS FOR HD 700K UNICORN SNP CHIP
# 20th SEPTEMBER 2013 (Updated 2022-07-21!)
# Susan E. Johnston
#
# PEDIGREE METHOD
#
# Using the method outlined in Druet et al 2013 (Heredity) with stepwise selection of individuals
# We also used the ten most prolific sires and dams within the pedigree
#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up working environment            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(reshape)
library(MASS)
library(kinship2)

pedigree <- read.table("data/pedigree200913.txt", header = T)   # Pedigree File

analysisPrefix <- "Unicorn"      # Prefix for saved files. Will be over-written on each iteration!

NoProlificSires <- 10            # How many of the most prolific sires should be included?
NoProlificDams  <- 10            # How many of the most prolific dams should be included?
NoSampReqd <- 20                 # How many more samples are required?

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Prepare the pedigree for analysis     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# This part determines discrete families for Unicorns and selects the largest
# one. If you have a look, you can see that basically all unicorns are related,
# and the rest are just duos and trios that only fit into the pedigree. We
# discarded these here. If you have many small pedigrees and a lot of group
# structure, you might want to modify the script to somehow take into account
# variation between family groups AND pedigrees.

families <- makefamid(pedigree[,1], pedigree[,2], pedigree[,3])
families.stats <- data.frame(table(families))
pedigree <- pedigree[families == which(families.stats$Freq == max(families.stats$Freq)),]

#~~ create a relationship matrix

rshipmat <- kinship(pedigree[,1], pedigree[,2], pedigree[,3])*2

#~~ create a reference dataset

ped.info <- pedigree
names(ped.info)[1] <- "X1"

rm(families, families.stats)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Calculate the average relationship of each individual with rest of population      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

MeanRship <- function(ID){
  ID <- as.character(ID)
  x <- (sum(rshipmat[,ID])-1)/(nrow(rshipmat)-1)
  x
}

ped.info$MeanRship <- sapply(ped.info$X1, MeanRship)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Determine the most prolific dams and sires                                         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Please admire some vintage sorting functions here!!!

sires <- data.frame(table(pedigree$FATHER), stringsAsFactors=F)
dams  <- data.frame(table(pedigree$MOTHER), stringsAsFactors=F)

sires <- as.numeric(as.character(sires[with(sires, order(-Freq)), ][1:NoProlificSires, 1]))
dams  <- as.numeric(as.character(dams [with(dams,  order(-Freq)), ][1:NoProlificDams,  1]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 4. Create a function to estimate "proportion of genome" represented by individuals    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

RShipSelection <- function(rshipmat, samp.ids, mean.rship){
  
  x <- which(dimnames(rshipmat)[[1]] %in% samp.ids)
  
  ped.sub <- rshipmat[x, x]
  
  Am.inv <- tryCatch(ginv(ped.sub), error = function(e) e)   #CREATE GENERALIZED INVERSE MATRIX
  
  if(!inherits(Am.inv, "error")) return(sum(Am.inv %*% mean.rship))
  
  if(inherits(Am.inv, "error")){
    return(0)
  }
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Carry out sequential sampling of individuals                                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ create columns to save the results

ped.info$Order <- NA
ped.info$CumuRship <- NA

#~~ Find the most related individual and score as '1' with cumulative relationship 0

ped.info$Order[which(ped.info$MeanRship == max(ped.info$MeanRship) & is.na(ped.info$Order))] <- 1
ped.info$CumuRship[which(ped.info$Order == 1)] <- 0

#~~ Find the most prolific individuals and score as '1' also

ped.info$Order[which(ped.info$X1 %in% c(sires, dams))] <- 1
ped.info$CumuRship[which(ped.info$X1 %in% c(sires, dams))] <- 0

#~~ set the next order for the iterative procedure

counter <- 2

#~~ Run the loop for the number of samples required.

for(i in 1:NoSampReqd){
  
  message(paste("Running iteration", i, "of", NoSampReqd))
  
  selectedIDs <- ped.info$X1[which(!is.na(ped.info$Order) & ped.info$Order != -9)]   # selected IDs so far
  nonselectedIDs <- data.frame(ID = ped.info$X1[which(is.na(ped.info$Order))])       # non-selected IDs so far
  
  if(nrow(nonselectedIDs) == 0) stop("Reached end of list - no more IDs possible!")  # 'warning' when no more IDs can be sampled before no. required
  nonselectedIDs$Rship <- NA
  
  # Go through as yet unselected IDs to see the change in relationship with the addition of each
  
  for(j in 1:nrow(nonselectedIDs)){             
    
    sampledIDs <- c(selectedIDs, nonselectedIDs$ID[j])
    meanRship  <- ped.info$MeanRship[which(ped.info$X1 %in% sampledIDs)]
    nonselectedIDs$Rship[j] <- RShipSelection(rshipmat, sampledIDs, meanRship)
    rm(sampledIDs, meanRship)
  }
  
  #~~ Which ID(s) increase the relationship captured? (If more than 1, sample one of them)
  
  idsconforming <- nonselectedIDs$ID[which(nonselectedIDs$Rship == max(nonselectedIDs$Rship,na.rm=T))]
  
  if(length(idsconforming) > 1) idsconforming <- sample(idsconforming, 1)
  message(paste0("Adding ID ", idsconforming))
  
  ped.info$Order[which(ped.info$X1 == idsconforming)] <- counter       # assign order
  ped.info$CumuRship[which(ped.info$X1 == idsconforming)] <- max(nonselectedIDs$Rship, na.rm=T)     # record the new relationship value
  
  counter <- counter + 1
  
  rm(selectedIDs, nonselectedIDs, idsconforming, j)
  
}

rm(i, sires, dams, counter)

#~~ save the results after finishing

write.table(ped.info, paste(analysisPrefix, "_Results.txt", sep = ""), row.names = F, sep = "\t", quote = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6. Data checking                                                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ped.info <- read.table(paste(analysisPrefix, "_Results.txt", sep = ""), header = T)

ped.info$Order[which(ped.info$Order == -9)] <- NA

plot(ped.info$Order, ped.info$MeanRship)
plot(ped.info$Order, ped.info$CumuRship)




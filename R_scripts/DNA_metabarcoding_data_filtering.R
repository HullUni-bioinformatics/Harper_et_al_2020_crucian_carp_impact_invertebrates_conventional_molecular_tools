#' ---
#' Title: "DNA metabarcoding data filtering and refinement"
#' Author: "Lynsey Rebecca Harper"
#' Date: "17th August 2018"
#' ---
#' 
#' 
#' Bulk tissue and eDNA samples from  ponds in North Norfolk, East of England, and
#' East Riding of Yorkshire, were screened for freshwater macroinvertebrates.
#' 
#' Here, the raw DNA metabarcoding data are filtered and refined for community analyses.
#' 
#' 
#' ## Prepare working environment
#' 
#' Clear R memory, set working directory and load required packages. Then,
#' load functions for calculating kappa coefficient, plotting model 
#' residuals and testing model fit.
#' 

## Clear memory
rm(list=ls())

## Direct R to folder containing packages on workstation
.libPaths(c("C:\\R\\rlib", .libPaths("rlib")))

## set working directory to the location of the script
# install.packages("rstudioapi") # first time for each computer
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Check working directory
getwd()


## Load required packages
p <- c("ggplot2","ggpubr","munsell","lazyeval","grid","gridExtra","lme4",
       "glmmADMB","coda","MASS","car","scales","AICcmodavg","xtable","gtools",
       "xlsx","reshape2","dplyr","plyr","tidyr","arm","RVAideMemoire","permute",
       "ResourceSelection","bbmle","RColorBrewer", "MuMIn",
       "ggmap","mapproj","geosphere","jpeg","proto","rjson","RgoogleMaps",
       "maps","labeling","ggsn","png","coin","modeltools","mvtnorm","pROC",
       "vegan","multcomp","betapart","adespatial","adegraphics","ade4")
new.packages <- p[!(p %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cran.ma.imperial.ac.uk/",
                                          dependencies=TRUE)
lapply(p, require, character.only = TRUE)


#'
#' To ensure reproducibility, print details about the version of R being
#' used for analysis.
#' 

sessionInfo()


#' ---
#' 
#' ## 1) Raw data processing
#' 
#' Examine the raw data output from metaBEAT.
#' 

## Original BLAST 
DNA.ass.raw <- read.csv("../Data/DNA_metabarcoding_assigned_raw.csv", header=TRUE)
summary(DNA.ass.raw[,c(1:5,length(DNA.ass.raw))])
head(DNA.ass.raw)
names(DNA.ass.raw)
str(DNA.ass.raw)

## Unassigned BLAST
DNA.unass.raw <- read.csv("../Data/DNA_metabarcoding_unassigned_raw.csv", header=TRUE)
summary(DNA.unass.raw[,c(1:5,length(DNA.unass.raw))])
head(DNA.unass.raw)
names(DNA.unass.raw)
str(DNA.unass.raw)

## Make first row header names
names(DNA.ass.raw) <- lapply(DNA.ass.raw[1,], as.character)
names(DNA.unass.raw) <- lapply(DNA.unass.raw[1,], as.character)
DNA.ass.raw <- DNA.ass.raw[-1,]
DNA.unass.raw <- DNA.unass.raw[-1,]

## Reset row names
rownames(DNA.ass.raw) <- NULL
rownames(DNA.unass.raw) <- NULL

## Remove '-nc.blast'
colnames(DNA.ass.raw) <- gsub('-nc.blast', '', colnames(DNA.ass.raw))
colnames(DNA.unass.raw) <- gsub('-nc.blast.blast', '', colnames(DNA.unass.raw))

## Rename first column
colnames(DNA.ass.raw)[1] <- "Assignment"
colnames(DNA.unass.raw)[1] <- "Assignment"

## Remove any taxonomic assignments that aren't metazoa from each dataframe
## using the taxonomy column created during processing with metaBEAT
DNA.ass <- DNA.ass.raw[which(grepl("Metazoa|unassigned", DNA.ass.raw$taxomomy)),]
DNA.unass <- DNA.unass.raw[which(grepl("Metazoa|unassigned", DNA.unass.raw$taxomomy)),]

## Save copies of dataframes for family level analysis
DNA.ass.fam <- DNA.ass
DNA.unass.fam <- DNA.unass

## Remove last column containing taxonomy
DNA.ass <- DNA.ass[-386]
DNA.unass <- DNA.unass[-386]

## Bind data frames
DNA.merged <- rbind(DNA.ass, DNA.unass)

## Merge read counts by taxonomic assignment
## Anything with the same name will be merged and anything unique will
## be retained
DNA.merged[,2:385] <- lapply(DNA.merged[,2:385], function(x) as.numeric(as.character(x)))
DNA.df <- ddply(DNA.merged, .(Assignment), numcolwise(sum))

## Export as .csv file
write.csv(DNA.df, "../Data/DNA_metabarcoding_merged.csv", row.names=FALSE)


#' --- 
#' 
#' ## 2) Refine dataset
#' 
#' Now, we need to further refine the metabarcoding dataset.
#' 
#' 1. Any spurious species must be removed. Use NBN atlas to check species
#'    occurrence records and ensure they match with sampling locations. This
#'    is also a good source for checking current taxonomy.
#' 2. Any Genus or Family assignments containing only one species in 
#'    the UK must be changed to that species.
#' 3. Likewise, for any species assignment which is the only species in 
#'    the UK and also has genus/family reads, read counts from all
#'    assignments must be merged.
#'  
#' A record of these changes will be kept as these changes are made.
#'      

## Inspect new dataframe
summary(DNA.df[,1:6])
names(DNA.df)
str(DNA.df)

#' 
#' First, remove spurious assignments.
#' 
#' - invertebrate_environmental_sample, row 210
#' 

DNA.true <- DNA.df[-c(210),]

## Reset row names of data frame for further indexing
rownames(DNA.true) <- NULL


#'
#' Now, correct species names:
#' 
#' - Chironomidae_sp._BOLD:ACD0662 = Chironomidae, row 155
#' - Chironomidae_sp._BOLD:ACI7830 = Chironomidae, row 156
#' - Chironomidae_sp._PA2_1 = Chironomidae, row 157
#' - Chironomus_sp._BOLD:AAI4299 = *Chironomus*, row 164
#' - Cladopelma_sp._BOLD-2016 = *Cladopelma*, row 168
#' - Colymbetes_sp._BMNH_1425212 = *Colymbetes fuscus*, row 170
#' - Dasybranchus_sp._DH1 = *Dasybranchus*, row 179
#' - Hypoderaeum_sp._Hubei-2014 = *Hypoderaeum*, row 208
#' - Limnesia_(Limnesia)_sp._HP-Hyd013 = *Limnesia*, row 215
#' - Limnesiidae_sp._BOLD:ACJ8659 = Limesiidae, row 217
#' - Macropelopia_sp._G_BA30 = *Macropelopia*, row 219
#' - Macrothrix_sp._HE-364 = *Macrothrix*, row 220
#' - Parachironomus_sp._BOLD:ACB9399 = *Parachironomus*, row 232
#' - Podocopida_sp._BOLD:AAH0903 = *Podocopida*, row 240
#' - Podocopida_sp._BOLD:AAH0910 = *Podocopida*, row 241
#' - Polypedilum_sp._BOLD-2016 = *Polypedilum*, row 246
#' - Procladius_cf._fuscus_BOLD-2016 = *Procladius*, row 249

DNA.true$Assignment <- as.character(DNA.true$Assignment)
DNA.true[155:157, "Assignment"] <- "Chironomidae"
DNA.true[164, "Assignment"] <- "Chironomus"
DNA.true[168, "Assignment"] <- "Cladopelma"
DNA.true[170, "Assignment"] <- "Colymbetes fuscus"
DNA.true[179, "Assignment"] <- "Dasybranchus"
DNA.true[208, "Assignment"] <- "Hypoderaeum"
DNA.true[215, "Assignment"] <- "Limnesia"
DNA.true[217, "Assignment"] <- "Limnesiidae"
DNA.true[219, "Assignment"] <- "Macropelopia"
DNA.true[220, "Assignment"] <- "Macrothrix"
DNA.true[232, "Assignment"] <- "Parachironomus"
DNA.true[240:241, "Assignment"] <- "Podocopida"
DNA.true[246, "Assignment"] <- "Polypedilum"
DNA.true[249, "Assignment"] <- "Procladius"

## Now merge read counts again by taxonomic assignment
## Anything with the same name in the data frame will be merged
DNA.true <- ddply(DNA.true, .(Assignment), numcolwise(sum))

## Now remove the underscore from species names
DNA.true$Assignment <- gsub('_', ' ', DNA.true$Assignment)


#' ---
#' 
#' ## 3) Clean up dataset
#' 
#' Now the data is in a form that can be manipulated easily, it must be 
#' filtered to remove potential contaminants and false positives. 
#' 
#' There are several ways of doing this:
#' 1. Identify highest level of assassin bug DNA contamination across all 
#'    samples
#' 3. Identify the highest level of contamination in positive (non-assassin
#'    bug DNA) and negative (any DNA) controls
#' 4. Identify species-specific thresholds using positive controls, i.e. 
#'    the frequency required to remove a given species from the positive 
#'    control as only assassin bug DNA should be present.
#' 
#' Arguably, species-specific thresholds based on positive controls are 
#' more effective as negative controls have no template DNA for 
#' contaminant DNA to compete with for amplification, thus contaminant
#' DNA amplifies exponentially. However, only 8 positive controls were
#' included on this MiSeq run which may render this approach ineffective.
#'

#########################################################
# OPTION 1: highest level of assassin bug contamination #
#########################################################

## Create copy of dataframe without positive or negative controls for 
## threshold determination
DNA.ass.bug <- DNA.true[,which(grepl("Assignment|-S-|-M-|-L-|-SO-|-UN-|MC0|MC1",
                                     colnames(DNA.true)))]

## Make Assignment column row names
rownames(DNA.ass.bug) <- DNA.ass.bug$Assignment
DNA.ass.bug <- DNA.ass.bug[,-1]

## Calculate total number of reads in samples
DNA.ass.bug <- rbind(DNA.ass.bug, colSums(DNA.ass.bug))

## Create new dataframe containing the frequency of reads in each sample
DNA.ass.bug.freq  <- DNA.ass.bug/c(DNA.ass.bug[273,])
DNA.ass.bug.freq[is.na(DNA.ass.bug.freq)] <- 0

## Add new column to this dataframe containing the maximum frequency of
## DNA from each taxon across samples
DNA.ass.bug.freq$Threshold <- apply(DNA.ass.bug.freq, 1, max)

## Check Threshold column has been created properly
head(DNA.ass.bug.freq[,335:340])

## Combine read frequencies with taxonomic assignment
species <- data.frame(DNA.true$Assignment)
colnames(species) <- "Assignment"
total <- data.frame("Total")
colnames(total) <- "Assignment"
species <- rbind(species, total)
DNA.ass.bug.freq <- cbind(species, DNA.ass.bug.freq)

## Print contamination threshold based on max. level of assassin bug
## DNA contamination
max(DNA.ass.bug.freq[204,341])   # 0.00016%

## In this scenario, any assignments <0.00016% total reads in biological 
## samples would be considered contamination. To examine effect of this
## threshold on data, replace all assignments that are less than or equal 
## to this threshold for a given species with zero.

## Create new dataframe containing the frequency of reads in each sample
## and apply threshold
DNA.ass.bug.test <- DNA.ass.bug.freq
DNA.ass.bug.test[DNA.ass.bug.test <= 0.00016] <- 0

## Now convert back into read counts.
## Remove last row containing frequencies and add the total read 
## counts to convert assignment frequencies back to read counts for 
## all samples.
DNA.ass.bug.test <- DNA.ass.bug.test[-273,-341]
rownames(DNA.ass.bug.test) <- NULL

total.counts <- data.frame(colSums(DNA.ass.bug[-273,]))
total.counts <- data.frame(t(total.counts))
colnames(total.counts) <- gsub('[.]', '-', colnames(total.counts))
total.counts <- cbind(total, total.counts)

## Now convert frequencies back to read counts
DNA.ass.bug.conversion <- smartbind(DNA.ass.bug.test, total.counts)
rownames(DNA.ass.bug.conversion) <- DNA.ass.bug.conversion$Assignment
DNA.ass.bug.conversion <- DNA.ass.bug.conversion[,-1]
DNA.ass.bug.FP <- DNA.ass.bug.conversion*c(DNA.ass.bug.conversion[273,])

## Remove total row, reset row names, recreate Assignment column
DNA.ass.bug.FP <- DNA.ass.bug.FP[-273,]
DNA.ass.bug.FP$Assignment <- rownames(DNA.ass.bug.FP)
DNA.ass.bug.FP <- DNA.ass.bug.FP[,c(340,1:339)]
rownames(DNA.ass.bug.FP) <- NULL

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- DNA.ass.bug[-273,]
temp$Assignment <- DNA.ass.bug.FP$Assignment
temp <- temp[,c(340,1:339)]
temp[,2:340][temp[,2:340] > 0] <- 1
DNA.ass.bug.FP[,2:340][DNA.ass.bug.FP[,2:340] > 0] <- 1
DNA.ass.bug1 <- data.frame(colSums(temp[,2:340]))
DNA.ass.bug2 <- data.frame(colSums(DNA.ass.bug.FP[,2:340]))
DNA.ass.bug.compare <- cbind(DNA.ass.bug1, DNA.ass.bug2)

## Calculate proportion of information lost
DNA.ass.bug.compare$proportion <- DNA.ass.bug.compare[,2]/DNA.ass.bug.compare[,1]*100
DNA.ass.bug.compare[is.na(DNA.ass.bug.compare)] <- 0

## Examine mean and range of proportion of taxa retained
## Exclude WHALL-M-03 as no taxa were detected in this sample
range(DNA.ass.bug.compare[-303,3])
mean(DNA.ass.bug.compare[-303,3])

## This would result in up to 70% taxa being removed.
## On average, 13% taxa are removed from samples.


########################################################
# OPTION 2: highest level of contamination in controls #
########################################################

## Store mock community data, extraction blanks and controls in new 
## dataframes
DNA.neg.controls <- DNA.true %>% select(Assignment, matches("Empty-|ExtBlank|ExtBlack|Negative"))
DNA.pos.controls <- DNA.true %>% select(Assignment, contains("Positive"))

## Positive and negative controls were subset into seperate dataframes
## earlier in script
head(DNA.pos.controls)
head(DNA.neg.controls)

## Check positive controls for highest level of contamination (any DNA
## except assassin bug). Remove column with taxonomic assignment.
DNA.pos <- DNA.pos.controls[,-1]

## Calculate total number of reads in samples
DNA.pos <- rbind(DNA.pos, colSums(DNA.pos))

## Create new dataframe containing the frequency of reads in each sample
DNA.pos.freq <- DNA.pos/c(DNA.pos[273,])

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
DNA.pos.freq$Threshold <- apply(DNA.pos.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(DNA.pos.freq)

## Combine read frequencies with taxonomic assignment
species <- data.frame(DNA.pos.controls$Assignment)
colnames(species) <- "Assignment"
total <- data.frame("Total")
colnames(total) <- "Assignment"
species <- rbind(species, total)
DNA.pos.freq <- cbind(species, DNA.pos.freq)

## Print contamination threshold based on max. level of non-assassin bug
## DNA contamination
## Exclude unassigned from threshold determination due to large number of
## reads, which cannot be distinguished as poorly amplified assassin bug DNA 
## or other taxa
rownames(DNA.pos.freq) <- DNA.pos.freq$Assignment
DNA.pos.freq <- DNA.pos.freq[,-1]
max(DNA.pos.freq[-c(204,269,273),9])   # 0.00037%

## In this scenario, any assignments <0.00037% total reads in biological 
## samples would be considered contamination. To examine effect of this
## threshold on data, replace all assignments that are less than or equal 
## to this threshold for a given species with zero.

## Create new dataframe containing the frequency of reads in each sample
## and apply threshold
DNA.pos.test <- DNA.ass.bug.freq[,-341]
DNA.pos.test[DNA.pos.test <= 0.00037] <- 0

## Now convert back into read counts.
## Remove last row and column containing total frequency and thresholds.
## Add the total read counts to convert assignment frequencies back to 
## read counts for all samples.
DNA.pos.test <- DNA.pos.test[-273,]
rownames(DNA.pos.test) <- NULL

total.counts <- data.frame(colSums(DNA.ass.bug[-273,]))
total.counts <- data.frame(t(total.counts))
colnames(total.counts) <- gsub('[.]', '-', colnames(total.counts))
total.counts <- cbind(total, total.counts)

## Now convert frequencies back to read counts
DNA.pos.conversion <- smartbind(DNA.pos.test, total.counts)
rownames(DNA.pos.conversion) <- DNA.pos.conversion$Assignment
DNA.pos.conversion <- DNA.pos.conversion[,-1]
DNA.pos.FP <- DNA.pos.conversion*c(DNA.pos.conversion[273,])

## Remove total row, reset row names and recreate Assignment column
DNA.pos.FP <- DNA.pos.FP[-273,]
DNA.pos.FP$Assignment <- rownames(DNA.pos.FP)
DNA.pos.FP <- DNA.pos.FP[,c(340,1:339)]
rownames(DNA.pos.FP) <- NULL

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- DNA.ass.bug[-273,]
temp$Assignment <- DNA.pos.FP$Assignment
temp <- temp[,c(340,1:339)]
temp[,2:340][temp[,2:340] > 0] <- 1
DNA.pos.FP[,2:340][DNA.pos.FP[,2:340] > 0] <- 1
DNA.pos1 <- data.frame(colSums(temp[,2:340]))
DNA.pos2 <- data.frame(colSums(DNA.pos.FP[,2:340]))
DNA.pos.compare <- cbind(DNA.pos1, DNA.pos2)

## Calculate proportion of information lost
DNA.pos.compare$proportion <- DNA.pos.compare[,2]/DNA.pos.compare[,1]*100
DNA.pos.compare[is.na(DNA.pos.compare)] <- 0

## Examine mean and range of proportion of taxa retained
## Exclude WHALL-M-03 as no taxa were detected in this sample
range(DNA.pos.compare[-303,3])
mean(DNA.pos.compare[-303,3])

## This would result in up to 70% taxa being removed.
## On average, 24% taxa are removed from samples.


## Check negative controls for highest level of contamination (any DNA). 
## Remove column with taxonomic assignment.
DNA.neg <- DNA.neg.controls[,-1]

## Calculate total number of reads in samples
DNA.neg <- rbind(DNA.neg, colSums(DNA.neg))

## Create new dataframe containing the frequency of reads in each sample
DNA.neg.freq <- DNA.neg/c(DNA.neg[273,])

## Convert NAs to 0s
DNA.neg.freq[is.na(DNA.neg.freq)] <- 0

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
DNA.neg.freq$Threshold <- apply(DNA.neg.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(DNA.neg.freq)

## Combine read frequencies with taxonomic assignment
species <- data.frame(DNA.neg.controls$Assignment)
colnames(species) <- "Assignment"
total <- data.frame("Total")
colnames(total) <- "Assignment"
species <- rbind(species, total)
DNA.neg.freq <- cbind(species, DNA.neg.freq)

## Print contamination threshold based on max. level of DNA contamination
rownames(DNA.neg.freq) <- DNA.neg.freq$Assignment
DNA.neg.freq <- DNA.neg.freq[,-1]
max(DNA.neg.freq[-c(269,273),38])   # 100%

## In this scenario, any assignments <100% total reads in biological 
## samples would be considered contamination. This is too extreme and
## thus cannot be used as a threshold.


#########################################
# OPTION 3: species-specific thresholds #
#########################################

## Check positive controls for highest level of contamination (any DNA
## except assassin bug). Remove column with taxonomic assignment.
DNA.pos <- DNA.pos.controls[,-1]

## Calculate total number of reads in samples
DNA.pos <- rbind(DNA.pos, colSums(DNA.pos))

## Create new dataframe containing the frequency of reads in each sample
DNA.pos.freq <- DNA.pos/c(DNA.pos[273,])

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
DNA.pos.freq$Threshold <- apply(DNA.pos.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(DNA.pos.freq)

## Combine read frequencies with taxonomic assignment
species <- data.frame(DNA.pos.controls$Assignment)
colnames(species) <- "Assignment"
total <- data.frame("Total")
colnames(total) <- "Assignment"
species <- rbind(species, total)
DNA.pos.freq <- cbind(species, DNA.pos.freq)

## Manually change the species threshold value for assassin bug to 0 as 
## this will be a true contaminant in the dataset and no false positive 
## threshold required.
DNA.pos.freq$Threshold[204] <- 0

## Check this manual edit has occurred
head(DNA.pos.freq[204,10])

## In this scenario, any assignments less than threshold in biological 
## samples would be considered contamination. To examine effect of this
## threshold on data, replace all assignments that are less than or equal 
## to this threshold for a given species with zero.

## Apply thresholds
DNA.SS.test <- DNA.ass.bug.freq[,-341]
DNA.SS.test[DNA.SS.test <= DNA.pos.freq$Threshold] <- 0

## Now convert back into read counts.
## Remove last row and column containing total frequency and thresholds.
## Add the total read counts to convert assignment frequencies back to 
## read counts for all samples.
DNA.SS.test <- DNA.SS.test[-273,]
rownames(DNA.SS.test) <- NULL

total.counts <- data.frame(colSums(DNA.ass.bug[-273,]))
total.counts <- data.frame(t(total.counts))
colnames(total.counts) <- gsub('[.]', '-', colnames(total.counts))
total.counts <- cbind(total, total.counts)

## Now convert frequencies back to read counts
DNA.SS.conversion <- smartbind(DNA.SS.test, total.counts)
rownames(DNA.SS.conversion) <- DNA.SS.conversion$Assignment
DNA.SS.conversion <- DNA.SS.conversion[,-1]
DNA.SS.FP <- DNA.SS.conversion*c(DNA.SS.conversion[273,])

## Remove total row, reset row names and recreate Assignment column
DNA.SS.FP <- DNA.SS.FP[-273,]
DNA.SS.FP$Assignment <- rownames(DNA.SS.FP)
DNA.SS.FP <- DNA.SS.FP[,c(340,1:339)]
rownames(DNA.SS.FP) <- NULL

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- DNA.ass.bug[-273,]
temp$Assignment <- DNA.SS.FP$Assignment
temp <- temp[,c(340,1:339)]
temp[,2:340][temp[,2:340] > 0] <- 1
DNA.SS.FP[,2:340][DNA.SS.FP[,2:340] > 0] <- 1
DNA.SS1 <- data.frame(colSums(temp[,2:340]))
DNA.SS2 <- data.frame(colSums(DNA.SS.FP[,2:340]))
DNA.SS.compare <- cbind(DNA.SS1, DNA.SS2)

## Calculate proportion of information lost
DNA.SS.compare$proportion <- DNA.SS.compare[,2]/DNA.SS.compare[,1]*100
DNA.SS.compare[is.na(DNA.SS.compare)] <- 0

## Examine mean and range of proportion of taxa retained
## Exclude WHALL-M-03 as no taxa were detected in this sample
range(DNA.SS.compare[-303,3])
mean(DNA.SS.compare[-303,3])

## This would result in up to 25% taxa being removed.
## On average, 2% taxa are removed from samples.


###########
# SUMMARY # 
###########

## Tidy dataframes
DNA.ass.bug.compare$Sample <- rownames(DNA.ass.bug.compare)
DNA.ass.bug.compare <- DNA.ass.bug.compare[,c(4,1:3)]
rownames(DNA.ass.bug.compare) <- NULL
colnames(DNA.ass.bug.compare)[2:4] <- c("sp_richness_NT",
                                        "sp_richness_TA",
                                        "prop_TA")
DNA.ass.bug.compare$type <- "assassin bug"

DNA.pos.compare$Sample <- rownames(DNA.pos.compare)
DNA.pos.compare <- DNA.pos.compare[,c(4,1:3)]
rownames(DNA.pos.compare) <- NULL
colnames(DNA.pos.compare)[2:4] <- c("sp_richness_NT",
                                    "sp_richness_TA",
                                    "prop_TA")
DNA.pos.compare$type <- "positive"

DNA.SS.compare$Sample <- rownames(DNA.SS.compare)
DNA.SS.compare <- DNA.SS.compare[,c(4,1:3)]
rownames(DNA.SS.compare) <- NULL
colnames(DNA.SS.compare)[2:4] <- c("sp_richness_NT",
                                   "sp_richness_TA",
                                   "prop_TA")
DNA.SS.compare$type <- "species-specific"

## Combine dataframes
DNA.threshold.df <- rbind(DNA.ass.bug.compare, DNA.pos.compare, DNA.SS.compare)

## Plot number of species retained after thresholds applied
p1 <- ggplot(dat=DNA.threshold.df, aes(x=Sample, y=sp_richness_TA, fill=type))
p1 <- p1 + geom_bar(stat="identity", position=position_dodge())
p1 <- p1 + facet_grid(type ~ .)
p1

## Plot proportion of species retained after thresholds applied
p1 <- ggplot(dat=DNA.threshold.df, aes(x=Sample, y=prop_TA, fill=type))
p1 <- p1 + geom_bar(stat="identity", position=position_dodge())
p1 <- p1 + scale_y_continuous(limits=c(0,100))
p1 <- p1 + facet_grid(type ~ .)
p1

## Going forward the assassin bug threshold will be used as this is
## stringent but also retains the majority of biological information.
## Recreate dataframe with threshold applied.
DNA.ass.bug.FP <- DNA.ass.bug.conversion*c(DNA.ass.bug.conversion[273,])
DNA.ass.bug.FP <- DNA.ass.bug.FP[-273,]
DNA.ass.bug.FP$Assignment <- rownames(DNA.ass.bug.FP)
DNA.ass.bug.FP <- DNA.ass.bug.FP[,c(340,1:339)]
rownames(DNA.ass.bug.FP) <- NULL


#' ---
#'
#' ## 4) Spurious assigments
#'
#' Remove positive control and any non-invertebrate assignments. Also remove 
#' invertebrate assignments higher than family level as these are too coarse to 
#' compare with microscopy. Assignments to be removed are as follows:
#' 
#' - Arthropoda, row 17
#' - Aves, row 22
#' - Bdelloidea, row 23
#' - Coleoptera, row 48
#' - Decapoda, row 70
#' - Hemiptera, row 115
#' - *Homo sapiens*, row 122
#' - *Hydra* (cnidaria), row 125
#' - *Hydra oligactis* (cnidaria), row 126
#' - *Hydra viridissima* (cnidaria), row 127
#' - Insecta, row 143
#' - *Lissotriton vulgaris* (smooth newt), row 160
#' - Metazoa, row 167
#' - *Platymeris biguttata*, row 204
#' - Podocopida, rows 206
#' - Rallidae (waterfowl), row 227
#' - Rhabditida, row 228
#' - *Scardinius erythrophthalmus* (Rudd), row 234
#' - Trombidiformes, row 267
#' - unassigned, row 269
#' 
#' 

DNA.refine <- DNA.ass.bug.FP[-c(17,22:23,48,70,115,122,125:127,143,160,167,
                                204,206,227:228,234,267,269),]

## Reset row names of data frame for further indexing
rownames(DNA.refine) <- NULL


#' ---
#' 
#' ## 5) Dataframes for downstream analyses
#' 
#' Now, create dataframe for analysing variation in PCR/sequencing replicates
#' for each bulk tissue sample.
#' 

## Transpose data frame
DNA.pcr.rep <- setNames(data.frame(t(DNA.refine[,-1])), DNA.refine[,1])

## Make row names a column in data frame
DNA.pcr.rep$Sample <- rownames(DNA.pcr.rep)

## Create new column containing PCR/sequencing replicate number
DNA.pcr.rep$Replicate <- c(1,2,3)

## Move new columns to start of dataframe
DNA.pcr.rep <- DNA.pcr.rep[,c(253:254,1:252)]

## Reset row names
rownames(DNA.pcr.rep) <- NULL

## Store Sample and Replicate columns in new dataframe as metadata
## Remove numbers after sample ID
DNA.pcr.rep.metadata <- DNA.pcr.rep[,1:2]
DNA.pcr.rep.metadata$Sample <- gsub('-01|-02|-03', '', DNA.pcr.rep.metadata$Sample)

## Save data in new dataframe for PCR/sequencing replicates to be pooled
DNA.to.pool <- DNA.pcr.rep

## Remove Replicate column from replicates dataframe
## Make Sample column the row names and remove from dataframe
rownames(DNA.pcr.rep) <- DNA.pcr.rep$Sample
DNA.pcr.rep <- DNA.pcr.rep[,-c(1:2)]

## Write dataframes to .csv files to be imported for analysis in a different script
write.csv(DNA.pcr.rep, "../Data/DNA_PCR_replicates.csv")
write.csv(DNA.pcr.rep.metadata, "../Data/DNA_PCR_replicates_metadata.csv", row.names=FALSE)


#'
#' Now, make dataframe for analysing variation in size categories for each
#' bulk tissue sample, and a dataframe for analysing variation in mock communities.
#' Pool the PCR/sequencing replicates according to size category/mock community.
#' 

## Make dataframes for size categories and mock communities
DNA.size.cat <- DNA.to.pool[which(!grepl("MC0|MC1", DNA.to.pool$Sample)),]
mock.comm <- DNA.to.pool[which(grepl("MC0|MC1", DNA.to.pool$Sample)),]


## Mock community dataframe
## Split mock community dataframe into two dataframes: one for annealing
## temperature gradient experiment, and one for different scenarios 
## experiment
mock.comm.grad <- mock.comm[which(grepl("-4|-5", mock.comm$Sample)),]
mock.comm <- mock.comm[which(!grepl("-4|-5", mock.comm$Sample)),]

## For annealing temperature gradient dataframe: 
## - create new columns containing community and annealing temperature 
## - move new columns to start of dataframe after Sample column
## - remove replicate numbers after temperature in new column
mock.comm.grad$Community <- rep(c("MC01","MC08"), each=18)
mock.comm.grad$Temperature <- gsub("^.*?-","", mock.comm.grad$Sample)
mock.comm.grad <- mock.comm.grad[,c(1,255:256,2:254)]
mock.comm.grad$Temperature <- gsub('-01|-02|-03', '', mock.comm.grad$Temperature)

## Store Sample, Community, Temperature and Replicate columns in new dataframe as 
## metadata
mock.comm.grad.metadata <- mock.comm.grad[,1:4]

## Make Sample column the row names, and remove metadata columns from temperature 
## gradient dataframe
rownames(mock.comm.grad) <- mock.comm.grad$Sample
mock.comm.grad <- mock.comm.grad[,-c(1:4)]

## Write dataframes to .csv files to be imported for analysis in a different script
write.csv(mock.comm.grad, "../Data/DNA_mock_community_gradient.csv")
write.csv(mock.comm.grad.metadata, "../Data/DNA_mock_community_gradient_metadata.csv",
          row.names=FALSE)


## For different community scenarios dataframe: 
## - store Sample and Replicate columns in new dataframe as metadata
## - remove replicate number after sample ID
## - make Sample column the row names of dataframe with assignments
## - remove the Sample and Replicate columns
mock.comm.metadata <- mock.comm[,1:2]
mock.comm.metadata$Sample <- gsub('-01|-02|-03', '', mock.comm.metadata$Sample)
rownames(mock.comm) <- mock.comm$Sample
mock.comm <- mock.comm[,-c(1:2)]

## Write dataframes to .csv files to be imported for analysis in a different script
write.csv(mock.comm, "../Data/DNA_mock_community.csv")
write.csv(mock.comm.metadata, "../Data/DNA_mock_community_metadata.csv", row.names=FALSE)


## Size category dataframe
## Remove PCR/sequencing replicate numbers after sample ID
DNA.size.cat$Sample <- gsub('-01|-02|-03', '', DNA.size.cat$Sample)

## Remove Replicate column
DNA.size.cat <- DNA.size.cat[,-2]

## Merge read counts by size category/mock community
## Anything with the same name in the data frame will be merged
## First, have to make sure that data to be summed is read as numeric
DNA.size.cat[,2:253] <- lapply(DNA.size.cat[,2:253], function(x) as.numeric(as.character(x)))
DNA.size.cat <- ddply(DNA.size.cat, .(Sample), numcolwise(sum))

## Create new column containing size category and move to start of dataframe
## after Sample colum
DNA.size.cat$Category <- gsub("^.*?-","", DNA.size.cat$Sample)
DNA.size.cat <- DNA.size.cat[,c(1,254,2:253)]

## Store Sample and Replicate columns in new dataframe as metadata
## Remove size catogory after sample ID
DNA.size.cat.metadata <- DNA.size.cat[,1:2]
DNA.size.cat.metadata$Sample <- gsub('-L|-M|-S|-SO|-UN', '', 
                                     DNA.size.cat.metadata$Sample)

## Save data in new dataframe for size categories to be pooled
DNA.pooled <- DNA.size.cat

## Make Sample column the row names 
## Remove Sample and Category columns from size categories dataframe
rownames(DNA.size.cat) <- DNA.size.cat$Sample
DNA.size.cat <- DNA.size.cat[,-c(1:2)]

## Write dataframes to .csv files to be imported for analysis in a different script
write.csv(DNA.size.cat, "../Data/DNA_size_categories.csv")
write.csv(DNA.size.cat.metadata, "../Data/DNA_size_categories_metadata.csv", row.names=FALSE)


#'
#' Finally, to assess how presence of crucian carp affects macroinvertebrate diversity,
#' we need to make dataframe where the communities in S/M/L bulk tissue samples are
#' pooled for each pond.
#' 

## First, subset dataframe for small, medium and large size categories
DNA.pooled <- DNA.pooled[which(!grepl("-SO|-UN", DNA.pooled$Sample)),]

## Remove column containing size category
DNA.pooled <- DNA.pooled[,-2]

## Remove numbers after sample ID
DNA.pooled$Sample <- gsub('-L|-M|-S', '', DNA.pooled$Sample)

## Merge read counts by sample
## Anything with the same name in the data frame will be merged
## First, have to make sure that data to be summed is read as numeric
DNA.pooled[,2:253] <- lapply(DNA.pooled[,2:253], function(x) as.numeric(as.character(x)))
DNA.pooled <- ddply(DNA.pooled, .(Sample), numcolwise(sum))

## Export as .csv file
write.csv(DNA.pooled, "../Data/DNA_metabarcoding_pooled.csv", row.names=FALSE)

## The metadata for the pooled DNA samples is the same as that which will be used for 
## netting and is contained in a manually prepared .csv file.


#' ---
#' 
#' 6) Dataframe for family level analysis
#' 
#' Now make dataframe containing family level assignments for S/M/L size categories
#' pooled by pond. This will be used to analyse the impact of crucian carp on 
#' invertebrate diversity at family level and compared to results obtained at species
#' level.
#' 

## Split taxonomy column into columns according to taxonomic rank
DNA.ass.fam <- DNA.ass.fam %>% separate(taxomomy, into = paste("id", 1:18, sep = ""))
DNA.unass.fam <- DNA.unass.fam %>% separate(taxomomy, into = paste("id", 1:18, sep = ""))

## Delete new columns except one containing family assignment
DNA.ass.fam <- DNA.ass.fam %>% select(Assignment, matches("-|id13"))
DNA.unass.fam <- DNA.unass.fam %>% select(Assignment, matches("-|id13"))

## Rename column containing family assignments and move to start of dataframe after
## column containing original assignment
names(DNA.ass.fam)[names(DNA.ass.fam) == "id13"] <- "Family"
names(DNA.unass.fam)[names(DNA.unass.fam) == "id13"] <- "Family"
DNA.ass.fam <- DNA.ass.fam[,c(1,386,2:385)]
DNA.unass.fam <- DNA.unass.fam[,c(1,386,2:385)]

## NB: some assignments higher than family level are still present but these have
## been given NA values. We need to retain these until the false positive threshold 
## has been applied

## Bind data frames
DNA.fam <- rbind(DNA.ass.fam, DNA.unass.fam)

## Replace any NA values with the original taxonomic assignment
DNA.fam$Family[is.na(DNA.fam$Family)] <- as.character(DNA.fam$Assignment[is.na(DNA.fam$Family)])

## Merge read counts by original taxonomic assignment
## Anything with the same name will be merged and anything unique will
## be retained
DNA.fam[,3:386] <- lapply(DNA.fam[,3:386], function(x) as.numeric(as.character(x)))
DNA.fam <- ddply(DNA.fam, .(Assignment, Family), numcolwise(sum))

## Remove spurious assignments
DNA.fam <- DNA.fam[-c(210),]

## Reset row names of data frame for further indexing
rownames(DNA.fam) <- NULL

#'
#' Now, correct names of original assignments as done earlier:
#' 
#' - Chironomidae_sp._BOLD:ACD0662 = Chironomidae, row 155
#' - Chironomidae_sp._BOLD:ACI7830 = Chironomidae, row 156
#' - Chironomidae_sp._PA2_1 = Chironomidae, row 157
#' - Chironomus_sp._BOLD:AAI4299 = *Chironomus*, row 164
#' - Cladopelma_sp._BOLD-2016 = *Cladopelma*, row 168
#' - Colymbetes_sp._BMNH_1425212 = *Colymbetes fuscus*, row 170
#' - Dasybranchus_sp._DH1 = *Dasybranchus*, row 179
#' - Hypoderaeum_sp._Hubei-2014 = *Hypoderaeum*, row 208
#' - Limnesia_(Limnesia)_sp._HP-Hyd013 = *Limnesia*, row 215
#' - Limnesiidae_sp._BOLD:ACJ8659 = Limesiidae, row 217
#' - Macropelopia_sp._G_BA30 = *Macropelopia*, row 219
#' - Macrothrix_sp._HE-364 = *Macrothrix*, row 220
#' - Parachironomus_sp._BOLD:ACB9399 = *Parachironomus*, row 232
#' - Podocopida_sp._BOLD:AAH0903 = *Podocopida*, row 240
#' - Podocopida_sp._BOLD:AAH0910 = *Podocopida*, row 241
#' - Polypedilum_sp._BOLD-2016 = *Polypedilum*, row 246
#' - Procladius_cf._fuscus_BOLD-2016 = *Procladius*, row 249

DNA.fam$Assignment <- as.character(DNA.fam$Assignment)
DNA.fam[155:157, "Assignment"] <- "Chironomidae"
DNA.fam[164, "Assignment"] <- "Chironomus"
DNA.fam[168, "Assignment"] <- "Cladopelma"
DNA.fam[170, "Assignment"] <- "Colymbetes fuscus"
DNA.fam[179, "Assignment"] <- "Dasybranchus"
DNA.fam[208, "Assignment"] <- "Hypoderaeum"
DNA.fam[215, "Assignment"] <- "Limnesia"
DNA.fam[217, "Assignment"] <- "Limnesiidae"
DNA.fam[219, "Assignment"] <- "Macropelopia"
DNA.fam[220, "Assignment"] <- "Macrothrix"
DNA.fam[232, "Assignment"] <- "Parachironomus"
DNA.fam[240:241, "Assignment"] <- "Podocopida"
DNA.fam[246, "Assignment"] <- "Polypedilum"
DNA.fam[249, "Assignment"] <- "Procladius"

## Also correct 'unknown' family assignment for Podocopida
DNA.fam[240:241, "Family"] <- "Podocopida"
DNA.fam[257:259, "Family"] <- "Philodinidae"

## Now merge read counts again by taxonomic assignment
## Anything with the same name in the data frame will be merged
DNA.fam <- ddply(DNA.fam, .(Assignment, Family), numcolwise(sum))

## Apply false positive threshold (0.00016%) to data
## Create dataframe containing only read counts for biological samples
fam.counts <- DNA.fam[,which(grepl("-S-|-M-|-L-|-SO-|-UN-|MC0|MC1",
                                   colnames(DNA.fam)))]

## Calculate total number of reads per sample
fam.freq <- rbind(fam.counts, colSums(fam.counts))

## Create new dataframe containing the frequency of reads in each sample
fam.freq <- fam.freq/c(fam.freq[273,])
fam.freq[is.na(fam.freq)] <- 0

## Apply false positive threshold
fam.freq[fam.freq <= 0.00016] <- 0

## Now convert back into read counts.
## Remove last row containing frequencies and add the total read 
## counts to convert assignment frequencies back to read counts for 
## all samples.
fam.freq <- fam.freq[-273,]
rownames(fam.freq) <- NULL

total.counts <- data.frame(colSums(fam.counts))
total.counts <- data.frame(t(total.counts))
colnames(total.counts) <- gsub('[.]', '-', colnames(total.counts))

## Now convert frequencies back to read counts
fam.conversion <- smartbind(fam.freq, total.counts)
fam.FP <- fam.conversion*c(fam.conversion[273,])

## Remove total row, reset row names, recreate Assignment column
fam.FP <- fam.FP[-273,]
fam.FP <- cbind(fam.FP, DNA.fam[,1:2])
fam.FP <- fam.FP[,c(340:341,1:339)]
rownames(fam.FP) <- NULL


#'
#' Remove positive control and any non-invertebrate assignments. Also remove 
#' invertebrate assignments higher than family level as these are too coarse 
#' to compare with microscopy. Assignments to be removed are as follows:
#' 
#' - Arthropoda, row 17
#' - Aves, row 22
#' - Bdelloidea, row 23
#' - Chromadorea, row 41
#' - Coleoptera, row 48
#' - Decapoda, row 70
#' - Diptera, row 74
#' - Hemiptera, row 115
#' - *Homo sapiens*, row 122
#' - *Hydra* (cnidaria), row 125
#' - *Hydra oligactis* (cnidaria), row 126
#' - *Hydra viridissima* (cnidaria), row 127
#' - Insecta, row 143
#' - *Lissotriton vulgaris* (smooth newt), row 160
#' - Metazoa, row 167
#' - *Platymeris bigutta*, row 204
#' - Podocopida, rows 206
#' - Rallidae (waterfowl), row 227
#' - Rhabditida, row 228
#' - *Scardinius erythrophthalmus* (Rudd), row 234
#' - Spongillida, row 241
#' - Trombidiformes, row 267
#' - unassigned, row 269
#' 
#' 

DNA.fam.refine <- fam.FP[-c(17,22:23,41,48,70,74,115,122,125:127,143,
                            160,167,204,206,227:228,234,241,267,269),]

## Reset row names of data frame for further indexing
rownames(DNA.fam.refine) <- NULL

## Now the column containing original taxonomic assignments can be removed, and the
## family assignments merged
DNA.fam.refine <- DNA.fam.refine[,-1]
DNA.fam.refine <- ddply(DNA.fam.refine, .(Family), numcolwise(sum))

## Now create dataframe that only contains S/M/L size categories to pooled for each
## pond
DNA.fam.pool <- DNA.fam.refine[,which(grepl("Family|-S-|-M-|-L-", colnames(DNA.fam.refine)))]

## Transpose data frame
DNA.fam.pool <- setNames(data.frame(t(DNA.fam.pool[,-1])), DNA.fam.pool[,1])

## Make row names first column in data frame
DNA.fam.pool$Sample <- rownames(DNA.fam.pool)
DNA.fam.pool <- DNA.fam.pool[,c(78,1:77)]

## Remove size category and replicate number from sample ID
DNA.fam.pool$Sample <- gsub("-S|-M|-L|-01|-02|-03", "", DNA.fam.pool$Sample)

## Pool data for each pond
DNA.fam.pool <- ddply(DNA.fam.pool, .(Sample), numcolwise(sum))

## Write dataframe to .csv file to be used for analysis in different script
write.csv(DNA.fam.pool, "../Data/DNA_metabarcoding_site_by_family.csv", row.names=FALSE)

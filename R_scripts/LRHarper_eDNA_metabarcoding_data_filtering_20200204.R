#' ---
#' Title: "eDNA metabarcoding data filtering and refinement"
#' Author: "Lynsey Rebecca Harper"
#' Date: "4th February 2020"
#' ---
#' 
#' 
#' Bulk tissue and eDNA samples from  ponds in North Norfolk, East of England, and
#' East Riding of Yorkshire, were screened for freshwater macroinvertebrates.
#' 
#' Here, the raw eDNA metabarcoding data are filtered and refined for community analyses.
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
eDNA.ass.raw <- read.csv("../Data/eDNA_metabarcoding_assigned_raw.csv", header=TRUE)
summary(eDNA.ass.raw[,c(1:5,length(eDNA.ass.raw))])
head(eDNA.ass.raw)
names(eDNA.ass.raw)
str(eDNA.ass.raw)

## Unassigned BLAST
eDNA.unass.raw <- read.csv("../Data/eDNA_metabarcoding_unassigned_raw.csv", header=TRUE)
summary(eDNA.unass.raw[,c(1:5,length(eDNA.unass.raw))])
head(eDNA.unass.raw)
names(eDNA.unass.raw)
str(eDNA.unass.raw)

## Make first row header names
names(eDNA.ass.raw) <- lapply(eDNA.ass.raw[1,], as.character)
names(eDNA.unass.raw) <- lapply(eDNA.unass.raw[1,], as.character)
eDNA.ass.raw <- eDNA.ass.raw[-1,]
eDNA.unass.raw <- eDNA.unass.raw[-1,]

## Reset row names
rownames(eDNA.ass.raw) <- NULL
rownames(eDNA.unass.raw) <- NULL

## Remove '-nc.blast'
colnames(eDNA.ass.raw) <- gsub('-nc.blast', '', colnames(eDNA.ass.raw))
colnames(eDNA.unass.raw) <- gsub('-nc.blast.blast', '', colnames(eDNA.unass.raw))

## Rename first column
colnames(eDNA.ass.raw)[1] <- "Assignment"
colnames(eDNA.unass.raw)[1] <- "Assignment"

## Rename WOFA1 as DAN, THW as THW16, and WOFA3 as WHALL for consistency 
## with DNA metabarcoding
colnames(eDNA.ass.raw) <- gsub("WOFA1","DAN", colnames(eDNA.ass.raw))
colnames(eDNA.ass.raw) <- gsub("WOFA3","WHALL", colnames(eDNA.ass.raw))
colnames(eDNA.ass.raw) <- gsub("THW","THW16", colnames(eDNA.ass.raw))
colnames(eDNA.unass.raw) <- gsub("WOFA1","DAN", colnames(eDNA.unass.raw))
colnames(eDNA.unass.raw) <- gsub("WOFA3","WHALL", colnames(eDNA.unass.raw))
colnames(eDNA.unass.raw) <- gsub("THW","THW16", colnames(eDNA.unass.raw))

## Create new dataframe for calculating the proportional read counts 
## for each taxon downstream.
## Duplicate raw read data
raw1 <- eDNA.ass.raw
raw2 <- eDNA.unass.raw

## Remove unassigned from reference database BLAST results as unassigned
## reads were extracted for BLAST against entire NCBI nucleotide database.
## Also remove last column containing taxonomy.
raw1 <- raw1[-85,-341]
raw2 <- raw2[,-341]

## Bind data frames
raw.merged <- rbind(raw1, raw2)

## Merge read counts by taxonomic assignment
## Anything with the same name will be merged and anything unique will
## be retained
raw.merged[,2:340] <- lapply(raw.merged[,2:340], function(x) as.numeric(as.character(x)))
raw.merged <- ddply(raw.merged, .(Assignment), numcolwise(sum))

## Make Assignment column row names
rownames(raw.merged) <- raw.merged$Assignment
raw.merged <- raw.merged[,-1]

## Calculate total number of reads in samples/controls
raw.merged <- rbind(raw.merged, colSums(raw.merged))

## Make new dataframe containing sample ID and total number of reads 
## per sample
raw.total <- raw.merged[515,]
raw.total$Assignment <- "Total"
raw.total <- raw.total[,c(340,1:339)]
rownames(raw.total) <- NULL

## Remove any taxonomic assignments that aren't metazoa from each dataframe
## using the taxonomy column created during processing with metaBEAT
eDNA.ass <- eDNA.ass.raw[which(grepl("Metazoa", eDNA.ass.raw$taxomomy)),]
eDNA.unass <- eDNA.unass.raw[which(grepl("Metazoa", eDNA.unass.raw$taxomomy)),]

## Save copies of dataframes for family level analysis
eDNA.ass.fam <- eDNA.ass
eDNA.unass.fam <- eDNA.unass

## Remove last column containing taxonomy
eDNA.ass <- eDNA.ass[,-341]
eDNA.unass <- eDNA.unass[,-341]

## Bind data frames
eDNA.merged <- rbind(eDNA.ass, eDNA.unass)

## Merge read counts by taxonomic assignment
## Anything with the same name will be merged and anything unique will
## be retained
eDNA.merged[,2:340] <- lapply(eDNA.merged[,2:340], function(x) as.numeric(as.character(x)))
eDNA.df <- ddply(eDNA.merged, .(Assignment), numcolwise(sum))

## Export as .csv file
write.csv(eDNA.df, "../Data/eDNA_metabarcoding_merged.csv", row.names=FALSE)


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
summary(eDNA.df[,1:6])
names(eDNA.df)
str(eDNA.df)

#' 
#' First, remove spurious assignments.
#' 
#' - invertebrate_environmental_sample, row 191
#' - uncultured_Philodina, row 293
#' 

eDNA.true <- eDNA.df[-c(191,293),]

## Reset row names of data frame for further indexing
rownames(eDNA.true) <- NULL


#'
#' Now, correct species names:
#' 
#' - Alona_sp._1_NA = *Alona*, row 4
#' - Acroperus_sp._2_NA = *Acroperus*, row 86
#' - Anystidae_sp._BOLD:AAM7961 = Anystidae, row 91
#' - Canthocamptidae_sp._BOLD:ACJ8158 = Canthocamptidae, row 104
#' - Chaetonotus_sp._1_TK-2012 = *Chaetonotus*, row 117
#' - Chaitophorus_sp._3_LMH-2015 = *Chaitophorus*, row 118
#' - Chironomidae_sp._BOLD:ACD0662 = Chironomidae, row 121
#' - Chironomidae_sp._BOLD:ACI7830 = Chironomidae, row 122
#' - Chironomidae_sp._PA2_1 = Chironomidae, row 123
#' - Chironomus_sp._BOLD-2016 = *Chironomus*, row 129
#' - Chironomus_sp._BOLD:AAI4299 = *Chironomus*, row 130
#' - Chydorus_sp._PS-2013 = *Chydorus*, row 134
#' - Cricotopus_sp._22ES = *Cricotopus*, row 141
#' - Cricotopus_sp._8ES = *Cricotopus*, row 142
#' - Daphnia_sp._BOLD:ACW5340 = *Daphnia*, row 153
#' - Dasybranchus_sp._DH1 = *Dasybranchus*, row 154
#' - Heterolepidoderma_sp._2_TK-2012 = *Heterolepidoderma*, row 181
#' - Hypsibius_cf._dujardini_DS-2016 = *Hypsibius*, row 189
#' - Krenopelopia_sp._BOLD:AAI2213 = *Krenopelopia*, row 195
#' - Macrothrix_sp._HE-364 = *Macrothrix*, row 206
#' - Maxillopoda_sp._BOLD:ACW5478 = Maxillopoda, row 207
#' - Maxillopoda_sp._BOLD:ACW5664 = Maxillopoda, row 208
#' - Pionidae_sp._BOLD:ACE2606 = Pionidae, row 229
#' - Podocopida_sp._BOLD:AAG1450 = *Podocopida*, row 236
#' - Podocopida_sp._BOLD:AAH0893 = *Podocopida*, row 237
#' - Podocopida_sp._BOLD:AAH0903 = *Podocopida*, row 238
#' - Podocopida_sp._BOLD:AAH0908 = *Podocopida*, row 239
#' - Podocopida_sp._BOLD:AAH0910 = *Podocopida*, row 240
#' - Polyarthra_dolichoptera_complex_sp._UO-2013 = *Polyarthra dolichoptera*, row 242
#' - Polyarthra_sp._EM-2017 = *Polyarthra*, row 243
#' - Polyarthra_sp._UO-2013 = *Polyarthra*, row 244
#' - Polyarthra_sp._WM-2017a = *Polyarthra*, row 245
#' - Polypedilum_sp._BOLD-2016 = *Polypedilum*, row 249
#' - Procladius_cf._fuscus_BOLD-2016 = *Procladius*, row 254
#' - Psocoptera_sp._BOLD:AAN8452 = *Psocoptera*, row 258
#' - Smittia_sp._8ES = *Smittia*, row 272
#' - Synchaeta_cf._tremula/oblonga_UO-2012 = *Synchaeta*, row 277
#' - Tarsonemidae_sp._BOLD:ABV3248 = Tarsonemidae, row 282
#' - Trombidiformes_sp._BOLD:ACI3657 = Trombidiformes, row 289
#' - Tydeidae_sp._BOLD:ACI5267 = Tydeidae, row 290

eDNA.true$Assignment <- as.character(eDNA.true$Assignment)
eDNA.true[4, "Assignment"] <- "Alona"
eDNA.true[86, "Assignment"] <- "Acroporus"
eDNA.true[91, "Assignment"] <- "Anystidae"
eDNA.true[104, "Assignment"] <- "Canthocamptidae"
eDNA.true[117, "Assignment"] <- "Chaetonotus"
eDNA.true[118, "Assignment"] <- "Chaitophorus"
eDNA.true[121:123, "Assignment"] <- "Chironomidae"
eDNA.true[129:130, "Assignment"] <- "Chironomus"
eDNA.true[134, "Assignment"] <- "Chydorus"
eDNA.true[141:142, "Assignment"] <- "Cricotopus"
eDNA.true[153, "Assignment"] <- "Daphnia"
eDNA.true[154, "Assignment"] <- "Dasybranchus"
eDNA.true[181, "Assignment"] <- "Heterolepidoderma"
eDNA.true[189, "Assignment"] <- "Hypsibius"
eDNA.true[195, "Assignment"] <- "Krenopelopia"
eDNA.true[206, "Assignment"] <- "Macrothrix"
eDNA.true[207:208, "Assignment"] <- "Maxillopoda"
eDNA.true[229, "Assignment"] <- "Pionidae"
eDNA.true[236:240, "Assignment"] <- "Podocopida"
eDNA.true[242, "Assignment"] <- "Polyarthra dolichoptera"
eDNA.true[243:245, "Assignment"] <- "Polyarthra"
eDNA.true[249, "Assignment"] <- "Polypedilum"
eDNA.true[254, "Assignment"] <- "Procladius"
eDNA.true[258, "Assignment"] <- "Psocoptera"
eDNA.true[272, "Assignment"] <- "Smittia"
eDNA.true[277, "Assignment"] <- "Synchaeta"
eDNA.true[282, "Assignment"] <- "Tarsonemidae"
eDNA.true[289, "Assignment"] <- "Trombidiformes"
eDNA.true[290, "Assignment"] <- "Tydeidae"

## Now merge read counts again by taxonomic assignment
## Anything with the same name in the data frame will be merged
eDNA.true <- ddply(eDNA.true, .(Assignment), numcolwise(sum))

## Now remove the underscore from species names
eDNA.true$Assignment <- gsub("_", " ", eDNA.true$Assignment)


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
#' 2. Identify the highest level of contamination in positive (non-assassin
#'    bug DNA) and negative (any DNA) controls
#' 3. Identify species-specific thresholds using positive controls, i.e. 
#'    the frequency required to remove a given species from the positive 
#'    control as only assassin bug DNA should be present.
#' 
#' Arguably, species-specific thresholds based on positive controls are 
#' more effective as negative controls have no template DNA for 
#' contaminant DNA to compete with for amplification, thus contaminant
#' DNA amplifies exponentially. However, only 11 positive controls were
#' included on this MiSeq run which may render this approach ineffective.
#'

#########################################################
# OPTION 1: highest level of assassin bug contamination #
#########################################################

## Bind total read counts to refined dataset
eDNA.ass.bug <- rbind(eDNA.true, raw.total)

## Create copy of dataframe without positive or negative controls for 
## threshold determination
eDNA.ass.bug <- eDNA.ass.bug[,which(!grepl("Blank|Empty|Negative|Positive",
                                           colnames(eDNA.ass.bug)))]

## Make Assignment column row names
rownames(eDNA.ass.bug) <- eDNA.ass.bug$Assignment
eDNA.ass.bug <- eDNA.ass.bug[,-1]

## Create new dataframe containing the frequency of reads in each sample
eDNA.ass.bug.freq <- eDNA.ass.bug/c(eDNA.ass.bug[271,])
eDNA.ass.bug.freq[is.na(eDNA.ass.bug.freq)] <- 0

## Add new column to this dataframe containing the maximum frequency of
## DNA from each taxon across samples
eDNA.ass.bug.freq$Threshold <- apply(eDNA.ass.bug.freq, 1, max)

## Check Threshold column has been created properly
head(eDNA.ass.bug.freq[,230:231])

## Combine read frequencies with taxonomic assignment
eDNA.ass.bug.freq <- tibble:::rownames_to_column(eDNA.ass.bug.freq, "Assignment")

## Print contamination threshold based on max. level of assassin bug
## DNA contamination
max(eDNA.ass.bug.freq[202,232])   # 0%

## In this scenario, 0% of reads in biological samples would be considered 
## contamination, thus 0% of taxa would be removed and all biological 
## information would be retained.


########################################################
# OPTION 2: highest level of contamination in controls #
########################################################

## Store mock community data, extraction blanks and controls in new 
## dataframes
eDNA.neg.controls <- eDNA.true %>% select(Assignment, matches("Empty|Blank|Negative"))
eDNA.pos.controls <- eDNA.true %>% select(Assignment, contains("Positive"))

## Duplicate positive control dataframe
eDNA.pos <- eDNA.pos.controls

## Make Assignment column row names
rownames(eDNA.pos) <- eDNA.pos$Assignment
eDNA.pos <- eDNA.pos[,-1]

## Extract positive controls from raw totals dataframe
pos.total <- raw.total[,which(grepl("Positive", colnames(raw.total)))]

## Bind total read counts to positive controls
eDNA.pos <- rbind(eDNA.pos, pos.total)

## Create new dataframe containing the frequency of reads in each sample
eDNA.pos.freq <- eDNA.pos/c(eDNA.pos[271,])
eDNA.pos.freq[is.na(eDNA.pos.freq)] <- 0

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
eDNA.pos.freq$Threshold <- apply(eDNA.pos.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(eDNA.pos.freq)

## Combine read frequencies with taxonomic assignment
eDNA.pos.freq <- tibble:::rownames_to_column(eDNA.pos.freq, "Assignment")

## Print contamination threshold based on max. level of non-assassin bug
## DNA contamination
rownames(eDNA.pos.freq) <- eDNA.pos.freq$Assignment
eDNA.pos.freq <- eDNA.pos.freq[,-1]
max(eDNA.pos.freq[-c(202,271),12])   # 0.0007530717

## In this scenario, any assignments <0.07530717% total reads in biological 
## samples would be considered contamination. To examine effect of this
## threshold on data, replace all assignments that are less than or equal 
## to this threshold for a given species with zero.

## Create new dataframe containing the frequency of reads in each sample
## and apply threshold
eDNA.pos.test <- eDNA.ass.bug.freq[,-232]
eDNA.pos.test[eDNA.pos.test <= 0.0007530717] <- 0

## Now convert back into read counts.
## Remove last row and column containing total frequency and thresholds.
## Add the total read counts to convert assignment frequencies back to 
## read counts for all samples.
eDNA.pos.test <- eDNA.pos.test[-271,]
rownames(eDNA.pos.test) <- NULL
sample.total <- raw.total[,which(!grepl("Blank|Empty|Negative|Positive",
                                        colnames(raw.total)))]
eDNA.pos.conversion <- rbind(eDNA.pos.test, sample.total)

## Now convert frequencies back to read counts
rownames(eDNA.pos.conversion) <- eDNA.pos.conversion$Assignment
eDNA.pos.conversion <- eDNA.pos.conversion[,-1]
eDNA.pos.FP <- eDNA.pos.conversion*c(eDNA.pos.conversion[271,])

## Remove total row, reset row names and recreate Assignment column
eDNA.pos.FP <- eDNA.pos.FP[-271,]
eDNA.pos.FP <- tibble:::rownames_to_column(eDNA.pos.FP, "Assignment")

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- eDNA.ass.bug[-271,]
temp$Assignment <- eDNA.pos.FP$Assignment
temp <- temp[,c(231,1:230)]
temp[,2:231][temp[,2:231] > 0] <- 1
eDNA.pos.FP[,2:231][eDNA.pos.FP[,2:231] > 0] <- 1
eDNA.pos1 <- data.frame(colSums(temp[,2:231]))
eDNA.pos2 <- data.frame(colSums(eDNA.pos.FP[,2:231]))
eDNA.pos.compare <- cbind(eDNA.pos1, eDNA.pos2)

## Calculate proportion of information lost
eDNA.pos.compare$proportion <- eDNA.pos.compare[,2]/eDNA.pos.compare[,1]*100
eDNA.pos.compare[is.na(eDNA.pos.compare)] <- 0

## Examine mean and range of proportion of taxa retained
## Exclude LDUN2-04, SKEY1-1-02, THW16-03 and THW16-05, WADD17-1-08 as no taxa 
## were detected in these samples
range(eDNA.pos.compare[-c(19,62,123,125,148),3])
mean(eDNA.pos.compare[-c(19,62,123,125,148),3])

## This would result in up to 78% taxa being removed.
## On average, 35% taxa are removed from samples.


## Check negative controls for highest level of contamination (any DNA). 
## Duplicate negative controls dataframe
eDNA.neg <- eDNA.neg.controls

## Make Assignment column row names
rownames(eDNA.neg) <- eDNA.neg$Assignment
eDNA.neg <- eDNA.neg[,-1]

## Extract positive controls from raw totals dataframe
neg.total <- raw.total[,which(grepl("Empty|Blank|Negative", 
                                    colnames(raw.total)))]

## Bind total read counts to negative controls
eDNA.neg <- rbind(eDNA.neg, neg.total)

## Create new dataframe containing the frequency of reads in each sample
eDNA.neg.freq <- eDNA.neg/c(eDNA.neg[271,])

## Convert NAs to 0s
eDNA.neg.freq[is.na(eDNA.neg.freq)] <- 0

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
eDNA.neg.freq$Threshold <- apply(eDNA.neg.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(eDNA.neg.freq)

## Combine read frequencies with taxonomic assignment
eDNA.neg.freq <- tibble:::rownames_to_column(eDNA.neg.freq, "Assignment")

## Print contamination threshold based on max. level of DNA contamination
rownames(eDNA.neg.freq) <- eDNA.neg.freq$Assignment
eDNA.neg.freq <- eDNA.neg.freq[,-1]
max(eDNA.neg.freq[-271,99])   # 0.9951373

## In this scenario, any assignments <99.51% total reads in biological 
## samples would be considered contamination. This is too extreme and
## thus cannot be used as a threshold. High levels of contamination 
## have occurred in empty tag combinations included on the sequencing
## run, indicative of tag jumping. Many of the full process blanks for
## the eDNA work are also contaminated. If thresholds were based on PCR 
## negative controls alone, the contamination threshold would be 0%, 


#########################################
# OPTION 3: species-specific thresholds #
#########################################

## Create new dataframe containing the frequency of reads in each sample
eDNA.pos.freq <- eDNA.pos/c(eDNA.pos[271,])

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
eDNA.pos.freq$Threshold <- apply(eDNA.pos.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(eDNA.pos.freq)

## Combine read frequencies with taxonomic assignment
eDNA.pos.freq <- tibble:::rownames_to_column(eDNA.pos.freq, "Assignment")

## Manually change the species threshold value for assassin bug to 0 as 
## this will be a true contaminant in the dataset and no false positive 
## threshold required.
eDNA.pos.freq$Threshold[202] <- 0

## Check this manual edit has occurred
head(eDNA.pos.freq[202,13])

## In this scenario, any assignments less than threshold in biological 
## samples would be considered contamination. To examine effect of this
## threshold on data, replace all assignments that are less than or equal 
## to this threshold for a given species with zero.

## Apply thresholds
eDNA.SS.test <- eDNA.ass.bug.freq[,-232]
eDNA.SS.test[eDNA.SS.test <= eDNA.pos.freq$Threshold] <- 0

## Now convert back into read counts.
## Remove last row and column containing total frequency and thresholds.
## Add the total read counts to convert assignment frequencies back to 
## read counts for all samples.
eDNA.SS.test <- eDNA.SS.test[-271,]
rownames(eDNA.SS.test) <- NULL
eDNA.SS.conversion <- rbind(eDNA.SS.test, sample.total)

## Now convert frequencies back to read counts
rownames(eDNA.SS.conversion) <- eDNA.SS.conversion$Assignment
eDNA.SS.conversion <- eDNA.SS.conversion[,-1]
eDNA.SS.FP <- eDNA.SS.conversion*c(eDNA.SS.conversion[271,])

## Remove total row, reset row names and recreate Assignment column
eDNA.SS.FP <- eDNA.SS.FP[-271,]
eDNA.SS.FP <- tibble:::rownames_to_column(eDNA.SS.FP, "Assignment")

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- eDNA.ass.bug[-271,]
temp$Assignment <- eDNA.SS.FP$Assignment
temp <- temp[,c(231,1:230)]
temp[,2:231][temp[,2:231] > 0] <- 1
eDNA.SS.FP[,2:231][eDNA.SS.FP[,2:231] > 0] <- 1
eDNA.SS1 <- data.frame(colSums(temp[,2:231]))
eDNA.SS2 <- data.frame(colSums(eDNA.SS.FP[,2:231]))
eDNA.SS.compare <- cbind(eDNA.SS1, eDNA.SS2)

## Calculate proportion of information lost
eDNA.SS.compare$proportion <- eDNA.SS.compare[,2]/eDNA.SS.compare[,1]*100
eDNA.SS.compare[is.na(eDNA.SS.compare)] <- 0

## Examine mean and range of proportion of taxa retained
## Exclude LDUN2-04, SKEY1-1-02, THW-03 and THW-05, WADD17-1-08 as no taxa 
## were detected in these samples
range(eDNA.SS.compare[-c(19,62,123,125,148),3])
mean(eDNA.SS.compare[-c(19,62,123,125,148),3])

## This would result in up to 9% taxa being removed.
## On average, 1% taxa are removed from samples.


###########
# SUMMARY # 
###########

## Tidy dataframes
eDNA.pos.compare$Sample <- rownames(eDNA.pos.compare)
eDNA.pos.compare <- eDNA.pos.compare[,c(4,1:3)]
rownames(eDNA.pos.compare) <- NULL
colnames(eDNA.pos.compare)[2:4] <- c("sp_richness_NT",
                                     "sp_richness_TA",
                                     "prop_TA")
eDNA.pos.compare$type <- "positive"

eDNA.SS.compare$Sample <- rownames(eDNA.SS.compare)
eDNA.SS.compare <- eDNA.SS.compare[,c(4,1:3)]
rownames(eDNA.SS.compare) <- NULL
colnames(eDNA.SS.compare)[2:4] <- c("sp_richness_NT",
                                    "sp_richness_TA",
                                    "prop_TA")
eDNA.SS.compare$type <- "species-specific"

## Combine dataframes
eDNA.threshold.df <- rbind(eDNA.pos.compare, eDNA.SS.compare)

## Plot number of species retained after thresholds applied
p1 <- ggplot(dat=eDNA.threshold.df, aes(x=Sample, y=sp_richness_TA, fill=type))
p1 <- p1 + geom_bar(stat="identity", position=position_dodge())
p1 <- p1 + labs(x="Sample", 
                y="Detections remaining after threshold application")
p1 <- p1 + theme_bw()
p1 <- p1 + theme(panel.background = element_rect(fill = 'white'),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.y = element_text(colour="black"),
                 legend.position = "none",
                 text = element_text(size=20))
p1 <- p1 + facet_grid(type ~ .)
p1

## Plot proportion of species retained after thresholds applied
p1 <- ggplot(dat=eDNA.threshold.df, aes(x=Sample, y=prop_TA, fill=type))
p1 <- p1 + geom_bar(stat="identity", position=position_dodge())
p1 <- p1 + scale_y_continuous(limits=c(0,100))
p1 <- p1 + labs(x="Sample", 
                y="Detections remaining after threshold application (%)")
p1 <- p1 + theme_bw()
p1 <- p1 + theme(panel.background = element_rect(fill = 'white'),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.y = element_text(colour="black"),
                 legend.position = "none",
                 text = element_text(size=20))
p1 <- p1 + facet_grid(type ~ .)
p1

## Going forward the species-specific threshold will be used. The assassin
## bug threshold indicated we did not have any contamination in the dataset,
## but applying the species-specific threshold as well will retain biological
## information and add an extra level of confidence to the data.
## Recreate dataframe with threshold applied.
eDNA.SS.FP <- eDNA.SS.conversion*c(eDNA.SS.conversion[271,])
eDNA.SS.FP <- eDNA.SS.FP[-271,]
eDNA.SS.FP <- tibble:::rownames_to_column(eDNA.SS.FP, "Assignment")


#' ---
#'
#' ## 4) Spurious assigments
#'
#' Remove positive control and any non-invertebrate assignments. Also remove 
#' invertebrate assignments higher than family level as these are too coarse 
#' to compare with microscopy. Assignments to be removed are as follows:
#' 
#' - Arthropoda, row 14
#' - *Arvicola amphibius* (water vole), row 15
#' - Bdelloidea, row 20
#' - *Carassius carassius* (crucian carp), row 31
#' - Chromadorea, row 54
#' - Columbidae, row 60
#' - Cyclopoida, row 74
#' - *Gallinula chloropus* (common moorhen), row 110
#' - *Gasterosteus aculeatus* (three-spined stickleback), row 111
#' - Hominidae, row 127
#' - *Homo sapiens* (human), row 128
#' - Hymenoptera, row 133
#' - Insecta, row 138
#' - *Lissotriton helveticus* (palmate newt), row 156
#' - *Lissotriton vulgaris* (smooth newt), row 157
#' - Maxillopoda, row 163
#' - Metazoa, row 165
#' - *Mimeoma maculata* (beetle not found in UK), row 170
#' - Muridae, row 171
#' - *Pica pica* (magpie), row 195
#' - Platyhelminthes, row 200
#' - *Platymeris biguttata*, row 202
#' - Podocopida, row 207
#' - Psocoptera, row 225
#' - *Pungitius* (ninespine stickleback), row 227
#' - *Scardinius erythrophthalmus* (Rudd), row 239
#' - Trombidiformes, row 265
#' 

eDNA.refine <- eDNA.SS.FP[-c(14:15,20,31,54,60,74,110:111,127:128,133,
                             138,156:157,163,165,170:171,195,200,202,207,
                             225,227,239,265),]

## Reset row names of data frame for further indexing
rownames(eDNA.refine) <- NULL

## Create separate dataframes for eDNA samples processed with standard workflow,
## and eDNA samples which had 12 PCR/sequencing replicates performed for them
eDNA.samples <- eDNA.refine[,which(!grepl("SKEY1|WADD17",colnames(eDNA.refine)))]
eDNA.exp <- eDNA.refine[,which(grepl("Assignment|SKEY1|WADD17",colnames(eDNA.refine)))]


#' ---
#' 
#' ## 5) Dataframes for downstream analyses
#' 
#' Make dataframe for analysing variation in PCR/sequencing replicates.
#' 

## Transpose data frame
eDNA.pcr.rep <- setNames(data.frame(t(eDNA.exp[,-1])), eDNA.exp[,1])

## Make row names a column in data frame
eDNA.pcr.rep$Sample <- rownames(eDNA.pcr.rep)

## Create new column containing biological replicate number
eDNA.pcr.rep$pcr_rep <- c(1:12)

## Move new columns to start of dataframe
eDNA.pcr.rep <- eDNA.pcr.rep[,c(244:245,1:243)]

## Reset row names
rownames(eDNA.pcr.rep) <- NULL

## Store Sample and Replicate columns in new dataframe as metadata
## Remove numbers after sample ID
eDNA.pcr.rep.metadata <- eDNA.pcr.rep[,1:2]
eDNA.pcr.rep.metadata$Sample <- gsub("-01|-02|-03|-04|-05|-06|-07|-08|-09|-10|-11|-12", "", 
                                     eDNA.pcr.rep.metadata$Sample)

## Save data in new dataframe for PCR/sequencing replicates to be pooled
eDNA.to.pool <- eDNA.pcr.rep

## Remove Replicate column from replicates dataframe
## Make Sample column the row names and remove from dataframe
rownames(eDNA.pcr.rep) <- eDNA.pcr.rep$Sample
eDNA.pcr.rep <- eDNA.pcr.rep[,-c(1:2)]

## Write dataframes to .csv files to be imported for analysis in a different script
write.csv(eDNA.pcr.rep, "../Data/eDNA_PCR_replicates.csv")
write.csv(eDNA.pcr.rep.metadata, "../Data/eDNA_PCR_replicates_metadata.csv", row.names=FALSE)


#'
#' Now, make dataframe for analysing variation in biological replicates.
#' Pool the PCR/sequencing replicates according to biological replicate
#' and combine with other biological replicates.
#' 

## Remove PCR/sequencing replicate numbers after sample ID
eDNA.to.pool$Sample <- gsub("-01|-02|-03|-04|-05|-06|-07|-08|-09|-10|-11|-12", "",
                            eDNA.to.pool$Sample)

## Remove Replicate column
eDNA.to.pool <- eDNA.to.pool[,-2]

## Merge read counts by biological replicate
## Anything with the same name in the data frame will be merged
## First, have to make sure that data to be summed is read as numeric
eDNA.to.pool[,2:244] <- lapply(eDNA.to.pool[,2:244], function(x) as.numeric(as.character(x)))
eDNA.to.pool <- ddply(eDNA.to.pool, .(Sample), numcolwise(sum))

## Create new column containing biological replicate number
eDNA.to.pool$bio_rep <- c(1:5)

## Move to start of dataframe after Sample column
eDNA.to.pool <- eDNA.to.pool[,c(1,245,2:244)]

## Make dataframe for biological replicates
eDNA.bio.rep <- setNames(data.frame(t(eDNA.samples[,-1])), eDNA.samples[,1])

## Make row names a column in data frame
eDNA.bio.rep$Sample <- rownames(eDNA.bio.rep)

## Create new column containing biological replicate number
eDNA.bio.rep$bio_rep <- c(1:5)

## Move new columns to start of dataframe
eDNA.bio.rep <- eDNA.bio.rep[,c(244:245,1:243)]

## Reset row names
rownames(eDNA.bio.rep) <- NULL

## Combine data for all biological replicates into single dataframe
eDNA.bio.rep <- rbind(eDNA.bio.rep, eDNA.to.pool)

## Make bioligical replicate numbers consistent
eDNA.bio.rep$Sample <- gsub("-1","-01", eDNA.bio.rep$Sample)
eDNA.bio.rep$Sample <- gsub("-2","-02", eDNA.bio.rep$Sample)
eDNA.bio.rep$Sample <- gsub("-3","-03", eDNA.bio.rep$Sample)
eDNA.bio.rep$Sample <- gsub("-4","-04", eDNA.bio.rep$Sample)
eDNA.bio.rep$Sample <- gsub("-5","-05", eDNA.bio.rep$Sample)

## Sort biological replicates by sample ID
eDNA.bio.rep <- eDNA.bio.rep[order(eDNA.bio.rep$Sample),] 

## Store Sample and Replicate columns in new dataframe as metadata
## Remove numbers after sample ID
eDNA.bio.rep.metadata <- eDNA.bio.rep[,1:2]
eDNA.bio.rep.metadata$Sample <- gsub("-01|-02|-03|-04|-05", "", 
                                     eDNA.bio.rep.metadata$Sample)

## Save data in new dataframe for biological replicates to be pooled
eDNA.pooled <- eDNA.bio.rep

## Remove Replicate column from replicates dataframe
## Make Sample column the row names and remove from dataframe
rownames(eDNA.bio.rep) <- eDNA.bio.rep$Sample
eDNA.bio.rep <- eDNA.bio.rep[,-c(1:2)]

## Write dataframes to .csv files to be imported for analysis in a different script
write.csv(eDNA.bio.rep, "../Data/eDNA_biological_replicates.csv")
write.csv(eDNA.bio.rep.metadata, "../Data/eDNA_biological_replicates_metadata.csv", row.names=FALSE)


#'
#' Now, make dataframe for analysing communities in pooled eDNA samples 
#' from ponds with and without crucian carp.
#' 

## Remove column containing biological replicate number
eDNA.pooled <- eDNA.pooled[,-2]

## Remove numbers after sample ID
eDNA.pooled$Sample <- gsub("-01|-02|-03|-04|-05", "", eDNA.pooled$Sample)

## Merge read counts by sample
## Anything with the same name in the data frame will be merged
## First, have to make sure that data to be summed is read as numeric
eDNA.pooled[,2:244] <- lapply(eDNA.pooled[,2:244], function(x) as.numeric(as.character(x)))
eDNA.pooled <- ddply(eDNA.pooled, .(Sample), numcolwise(sum))

## For comparison with the DNA metabarcoding data, we must remove
## eDNA samples from pooled dataframe that do not have corresponding 
## bulk tissue samples: GUES1, MYST, OTOM, POFA4, PYES2, WOFA2
eDNA.pooled <- eDNA.pooled[which(!grepl("GUES1|MYST|OTOM|POFA4|PYES2|WOFA2", 
                                        eDNA.pooled$Sample)),]
rownames(eDNA.pooled) <- NULL

## Export as .csv file
write.csv(eDNA.pooled, "../Data/eDNA_metabarcoding_pooled.csv", row.names=FALSE)

## The metadata for the pooled eDNA samples is the same as that which will be used for 
## netting and is contained in a manually prepared .csv file.


#' ---
#' 
#' 6) Dataframe for family level analysis
#' 
#' Now make dataframe containing family level assignments for biological replicates
#' pooled by pond. This will be used to analyse the impact of crucian carp on 
#' invertebrate diversity at family level and compared to results obtained at species
#' level.
#' 

## Split taxonomy column into columns according to taxonomic rank
eDNA.ass.fam <- eDNA.ass.fam %>% separate(taxomomy, into = paste("id", 1:18, sep = ""))
eDNA.unass.fam <- eDNA.unass.fam %>% separate(taxomomy, into = paste("id", 1:18, sep = ""))

## Delete new columns except one containing family assignment
eDNA.ass.fam <- eDNA.ass.fam %>% select(Assignment, matches("-|id13"))
eDNA.unass.fam <- eDNA.unass.fam %>% select(Assignment, matches("-|id13"))

## Rename column containing family assignments and move to start of dataframe after
## column containing original assignment
names(eDNA.ass.fam)[names(eDNA.ass.fam) == "id13"] <- "Family"
names(eDNA.unass.fam)[names(eDNA.unass.fam) == "id13"] <- "Family"
eDNA.ass.fam <- eDNA.ass.fam[,c(1,341,2:340)]
eDNA.unass.fam <- eDNA.unass.fam[,c(1,341,2:340)]

## NB: some assignments higher than family level are still present but these have
## been given NA values. We need to retain these until the false positive threshold 
## has been applied

## Bind data frames
eDNA.fam <- rbind(eDNA.ass.fam, eDNA.unass.fam)

## Replace any NA values with the original taxonomic assignment
eDNA.fam$Family[is.na(eDNA.fam$Family)] <- as.character(eDNA.fam$Assignment[is.na(eDNA.fam$Family)])

## Merge read counts by original taxonomic assignment
## Anything with the same name will be merged and anything unique will
## be retained
eDNA.fam[,3:341] <- lapply(eDNA.fam[,3:341], function(x) as.numeric(as.character(x)))
eDNA.fam <- ddply(eDNA.fam, .(Assignment, Family), numcolwise(sum))

## Remove spurious assignments
eDNA.fam <- eDNA.fam[-c(191,293),]

## Reset row names of data frame for further indexing
rownames(eDNA.fam) <- NULL

#' 
#' Now, correct names of original assignments as done earlier:
#' 
#' - Alona_sp._1_NA = *Alona*, row 4
#' - Acroperus_sp._2_NA = *Acroperus*, row 86
#' - Anystidae_sp._BOLD:AAM7961 = Anystidae, row 91
#' - Canthocamptidae_sp._BOLD:ACJ8158 = Canthocamptidae, row 104
#' - Chaetonotus_sp._1_TK-2012 = *Chaetonotus*, row 117
#' - Chaitophorus_sp._3_LMH-2015 = *Chaitophorus*, row 118
#' - Chironomidae_sp._BOLD:ACD0662 = Chironomidae, row 121
#' - Chironomidae_sp._BOLD:ACI7830 = Chironomidae, row 122
#' - Chironomidae_sp._PA2_1 = Chironomidae, row 123
#' - Chironomus_sp._BOLD-2016 = *Chironomus*, row 129
#' - Chironomus_sp._BOLD:AAI4299 = *Chironomus*, row 130
#' - Chydorus_sp._PS-2013 = *Chydorus*, row 134
#' - Cricotopus_sp._22ES = *Cricotopus*, row 141
#' - Cricotopus_sp._8ES = *Cricotopus*, row 142
#' - Daphnia_sp._BOLD:ACW5340 = *Daphnia*, row 153
#' - Dasybranchus_sp._DH1 = *Dasybranchus*, row 154
#' - Heterolepidoderma_sp._2_TK-2012 = *Heterolepidoderma*, row 181
#' - Hypsibius_cf._dujardini_DS-2016 = *Hypsibius*, row 189
#' - Krenopelopia_sp._BOLD:AAI2213 = *Krenopelopia*, row 195
#' - Macrothrix_sp._HE-364 = *Macrothrix*, row 206
#' - Maxillopoda_sp._BOLD:ACW5478 = Maxillopoda, row 207
#' - Maxillopoda_sp._BOLD:ACW5664 = Maxillopoda, row 208
#' - Pionidae_sp._BOLD:ACE2606 = Pionidae, row 229
#' - Podocopida_sp._BOLD:AAG1450 = *Podocopida*, row 236
#' - Podocopida_sp._BOLD:AAH0893 = *Podocopida*, row 237
#' - Podocopida_sp._BOLD:AAH0903 = *Podocopida*, row 238
#' - Podocopida_sp._BOLD:AAH0908 = *Podocopida*, row 239
#' - Podocopida_sp._BOLD:AAH0910 = *Podocopida*, row 240
#' - Polyarthra_dolichoptera_complex_sp._UO-2013 = *Polyarthra dolichoptera*, row 242
#' - Polyarthra_sp._EM-2017 = *Polyarthra*, row 243
#' - Polyarthra_sp._UO-2013 = *Polyarthra*, row 244
#' - Polyarthra_sp._WM-2017a = *Polyarthra*, row 245
#' - Polypedilum_sp._BOLD-2016 = *Polypedilum*, row 249
#' - Procladius_cf._fuscus_BOLD-2016 = *Procladius*, row 254
#' - Psocoptera_sp._BOLD:AAN8452 = *Psocoptera*, row 258
#' - Smittia_sp._8ES = *Smittia*, row 272
#' - Synchaeta_cf._tremula/oblonga_UO-2012 = *Synchaeta*, row 277
#' - Tarsonemidae_sp._BOLD:ABV3248 = Tarsonemidae, row 282
#' - Trombidiformes_sp._BOLD:ACI3657 = Trombidiformes, row 289
#' - Tydeidae_sp._BOLD:ACI5267 = Tydeidae, row 290

eDNA.fam$Assignment <- as.character(eDNA.fam$Assignment)
eDNA.fam[4, "Assignment"] <- "Alona"
eDNA.fam[86, "Assignment"] <- "Acroporus"
eDNA.fam[91, "Assignment"] <- "Anystidae"
eDNA.fam[104, "Assignment"] <- "Canthocamptidae"
eDNA.fam[117, "Assignment"] <- "Chaetonotus"
eDNA.fam[118, "Assignment"] <- "Chaitophorus"
eDNA.fam[121:123, "Assignment"] <- "Chironomidae"
eDNA.fam[129:130, "Assignment"] <- "Chironomus"
eDNA.fam[134, "Assignment"] <- "Chydorus"
eDNA.fam[141:142, "Assignment"] <- "Cricotopus"
eDNA.fam[153, "Assignment"] <- "Daphnia"
eDNA.fam[154, "Assignment"] <- "Dasybranchus"
eDNA.fam[181, "Assignment"] <- "Heterolepidoderma"
eDNA.fam[189, "Assignment"] <- "Hypsibius"
eDNA.fam[195, "Assignment"] <- "Krenopelopia"
eDNA.fam[206, "Assignment"] <- "Macrothrix"
eDNA.fam[207:208, "Assignment"] <- "Maxillopoda"
eDNA.fam[229, "Assignment"] <- "Pionidae"
eDNA.fam[236:240, "Assignment"] <- "Podocopida"
eDNA.fam[242, "Assignment"] <- "Polyarthra dolichoptera"
eDNA.fam[243:245, "Assignment"] <- "Polyarthra"
eDNA.fam[249, "Assignment"] <- "Polypedilum"
eDNA.fam[254, "Assignment"] <- "Procladius"
eDNA.fam[258, "Assignment"] <- "Psocoptera"
eDNA.fam[272, "Assignment"] <- "Smittia"
eDNA.fam[277, "Assignment"] <- "Synchaeta"
eDNA.fam[282, "Assignment"] <- "Tarsonemidae"
eDNA.fam[289, "Assignment"] <- "Trombidiformes"
eDNA.fam[290, "Assignment"] <- "Tydeidae"

## Also correct 'unknown' family assignment for Podocopida
eDNA.fam[207:208, "Family"] <- "Maxillopoda"
eDNA.fam[236:240, "Family"] <- "Podocopida"
eDNA.fam[258, "Family"] <- "Psocoptera"
eDNA.fam[263:268, "Family"] <- "Philodinidae"
eDNA.fam[289, "Family"] <- "Trombidiformes"

## Now merge read counts again by taxonomic assignment
## Anything with the same name in the data frame will be merged
eDNA.fam <- ddply(eDNA.fam, .(Assignment, Family), numcolwise(sum))

## Apply species-specific false positive thresholds to data
## Create dataframe containing only read counts for biological samples
fam.counts <- eDNA.fam[,which(!grepl("Assignment|Family|Blank|Empty|Negative|Positive",
                                     colnames(eDNA.fam)))]

## Calculate total number of reads per sample
fam.freq <- rbind(fam.counts, sample.total[,-1])

## Create new dataframe containing the frequency of reads in each sample
fam.freq <- fam.freq/c(fam.freq[271,])
fam.freq[is.na(fam.freq)] <- 0

## Apply species-specific false positive thresholds
fam.freq[fam.freq <= eDNA.pos.freq$Threshold] <- 0

## Now convert back into read counts.
## Remove last row containing frequencies and add the total read 
## counts to convert assignment frequencies back to read counts for 
## all samples.
fam.freq <- fam.freq[-271,]
rownames(fam.freq) <- NULL
fam.conversion <- rbind(fam.freq, sample.total[,-1])
fam.FP <- fam.conversion*c(fam.conversion[271,])

## Remove total row, reset row names, recreate Assignment column
fam.FP <- fam.FP[-271,]
fam.FP <- cbind(eDNA.fam[,1:2], fam.FP)


#' 
#' Remove any non-invertebrate assignments. Also remove invertebrate assignments 
#' higher than family level as these are too coarse to compare with microscopy. 
#' Assignments to be removed are as follows:
#' 
#' - Arthropoda, row 14
#' - *Arvicola amphibius* (water vole), row 15
#' - Bdelloidea, row 20
#' - *Carassius carassius* (crucian carp), row 31
#' - Chromadorea, row 54
#' - Columbidae, row 60
#' - Cyclopoida, row 74
#' - Diptera, row 90
#' - *Gallinula chloropus* (common moorhen), row 110
#' - *Gasterosteus aculeatus* (three-spined stickleback), row 111
#' - Hominidae, row 127
#' - *Homo sapiens* (human), row 128
#' - Hymenoptera, row 133
#' - Insecta, row 138
#' - *Lissotriton helveticus* (palmate newt), row 156
#' - *Lissotriton vulgaris* (smooth newt), row 157
#' - Maxillopoda, row 163
#' - Metazoa, row 165
#' - *Mimeoma maculata* (beetle not found in UK), row 170
#' - Muridae, row 171
#' - *Pica pica* (magpie), row 195
#' - Platyhelminthes, row 200
#' - *Platymeris biguttata*, row 202
#' - Ploima, row 205
#' - Podocopida, row 207
#' - Psocoptera, row 225
#' - *Pungitius* (ninespine stickleback), row 227
#' - *Scardinius erythrophthalmus* (Rudd), row 239
#' - Trombidiformes, row 265
#' 

eDNA.fam.refine <- fam.FP[-c(14:15,20,31,54,60,74,90,110:111,127:128,133,
                             138,156:157,163,165,170:171,195,200,202,205,
                             207,225,227,239,265),]

## Reset row names of data frame for further indexing
rownames(eDNA.fam.refine) <- NULL

## Now the column containing original taxonomic assignments can be removed, and the
## family assignments merged
eDNA.fam.refine <- eDNA.fam.refine[,-1]
eDNA.fam.refine <- ddply(eDNA.fam.refine, .(Family), numcolwise(sum))

## Now create dataframe of pooled eDNA data for each pond
## Transpose data frame
eDNA.fam.pool <- setNames(data.frame(t(eDNA.fam.refine[,-1])), eDNA.fam.refine[,1])

## Make row names first column in data frame
eDNA.fam.pool$Sample <- rownames(eDNA.fam.pool)
eDNA.fam.pool <- eDNA.fam.pool[,c(98,1:97)]

## Remove biological and PCR replicate number from sample ID
eDNA.fam.pool$Sample <- gsub("-01|-02|-03|-04|-05|-06|-07|-08|-09|-10|-11|-12", "", 
                             eDNA.fam.pool$Sample)
eDNA.fam.pool$Sample <- gsub("-1|-2|-3|-4|-5", "", eDNA.fam.pool$Sample)

## Pool data for each pond
eDNA.fam.pool <- ddply(eDNA.fam.pool, .(Sample), numcolwise(sum))

## For comparison with the DNA metabarcoding data, we must remove
## eDNA samples from pooled dataframe that do not have corresponding 
## bulk tissue samples: GUES1, MYST, OTOM, POFA4, PYES2, WOFA2
eDNA.fam.pool <- eDNA.fam.pool[which(!grepl("GUES1|MYST|OTOM|POFA4|PYES2|WOFA2", 
                                            eDNA.fam.pool$Sample)),]
rownames(eDNA.fam.pool) <- NULL

## Write dataframe to .csv file to be used for analysis in different script
write.csv(eDNA.fam.pool, "../Data/eDNA_metabarcoding_site_by_family.csv", row.names=FALSE)

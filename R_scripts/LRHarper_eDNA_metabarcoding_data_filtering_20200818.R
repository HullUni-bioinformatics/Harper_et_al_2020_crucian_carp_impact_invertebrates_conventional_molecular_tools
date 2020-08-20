#' ---
#' Title: "eDNA metabarcoding data filtering and refinement"
#' Author: "Lynsey Rebecca Harper"
#' Date: "18th August 2020"
#' ---
#' 
#' 
#' Bulk tissue and eDNA samples from  ponds in North Norfolk, East of 
#' England, and East Riding of Yorkshire, were screened for freshwater 
#' macroinvertebrates.
#' 
#' Here, the raw eDNA metabarcoding data are filtered and refined for 
#' community analyses.
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
p <- c("plyr","tidyverse","reshape2","ggpubr","taxize")
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
## for each taxon downstream. Duplicate raw read data:
raw1 <- eDNA.ass.raw
raw2 <- eDNA.unass.raw

## Remove unassigned from reference database BLAST results as unassigned
## reads were extracted for BLAST against entire NCBI nucleotide database.
## Also remove last column containing taxonomy.
raw1 <- raw1[-85,-341]
raw2 <- raw2[,-341]

## Bind data frames
raw.merged <- rbind(raw1, raw2)

## Merge read counts by taxonomic assignment. Anything with the same name 
## will be merged and anything unique will be retained.
raw.merged[,2:340] <- lapply(raw.merged[,2:340], function(x) as.numeric(as.character(x)))
raw.merged <- ddply(raw.merged, .(Assignment), numcolwise(sum))

## Make Assignment column row names
raw.merged <- column_to_rownames(raw.merged, "Assignment")

## Calculate total number of reads in samples/controls
raw.merged <- rbind(raw.merged, colSums(raw.merged))

## Make new dataframe containing sample ID and total number of reads 
## per sample
raw.total <- raw.merged[515,]
raw.total$Assignment <- "Total"
raw.total <- raw.total[,c(340,1:339)]
rownames(raw.total) <- NULL

## Remove any taxonomic assignments that aren't metazoa from each 
## dataframe using the taxonomy column created during processing with 
## metaBEAT
eDNA.ass <- eDNA.ass.raw[which(grepl("Metazoa", eDNA.ass.raw$taxomomy)),]
eDNA.unass <- eDNA.unass.raw[which(grepl("Metazoa", eDNA.unass.raw$taxomomy)),]

## Remove last column containing taxonomy
eDNA.ass <- eDNA.ass[,-341]
eDNA.unass <- eDNA.unass[,-341]

## Bind data frames
eDNA.merged <- rbind(eDNA.ass, eDNA.unass)

## Merge read counts by taxonomic assignment. Anything with the same name 
## will be merged and anything unique will be retained.
eDNA.merged[,2:340] <- lapply(eDNA.merged[,2:340], function(x) as.numeric(as.character(x)))
eDNA.df <- ddply(eDNA.merged, .(Assignment), numcolwise(sum))

## Export as .csv file
write.csv(eDNA.df, 
          "../Data/eDNA_metabarcoding_merged_20200818.csv", 
          row.names=FALSE)



#' --- 
#' 
#' ## 2) Refine dataset
#' 
#' Now, we need to further refine the metabarcoding dataset.
#' 
#' 1. Any spurious species must be removed. Use NBN atlas to check species
#'    occurrence records and ensure they match with sampling locations. 
#'    This is also a good source for checking current taxonomy.
#' 2. Any genus or family assignments containing only one species in 
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

## First, remove spurious assignments:
## 
## - invertebrate_environmental_sample, row 191
## - uncultured_Philodina, row 293
## 

eDNA.true <- eDNA.df[-c(191,293),]

## Reset row names of data frame for further indexing
rownames(eDNA.true) <- NULL

## Now, correct species names:
## 
## - Alona_sp._1_NA = Alona, row 4
## - Abacarus_lolii = Abacarus, row 85
## - Acroperus_sp._2_NA = Acroperus, row 86
## - Anystidae_sp._BOLD:AAM7961 = Anystidae, row 91
## - Canthocamptidae_sp._BOLD:ACJ8158 = Canthocamptidae, row 104
## - Chaetonotus_aemilianus = Chaetonotus, row 113
## - Chaetonotus_antrumus = Chaetonotus, row 114
## - Chaetonotus_daphnes = Chaetonotus, row 115
## - Chaetonotus_similis = Chaetonotus, row 116
## - Chaetonotus_sp._1_TK-2012 = Chaetonotus, row 117
## - Chaitophorus_sp._3_LMH-2015 = Chaitophorus, row 118
## - Chironomidae_sp._BOLD:ACD0662 = Chironomidae, row 121
## - Chironomidae_sp._BOLD:ACI7830 = Chironomidae, row 122
## - Chironomidae_sp._PA2_1 = Chironomidae, row 123
## - Chironomus_sp._BOLD-2016 = Chironomus, row 129
## - Chironomus_sp._BOLD:AAI4299 = Chironomus, row 130
## - Chydorus_brevilabris = Chydorus, row 133
## - Chydorus_sp._PS-2013 = Chydorus, row 134
## - Cricotopus_sp._22ES = Cricotopus, row 141
## - Cricotopus_sp._8ES = Cricotopus, row 142
## - Daphnia_sp._BOLD:ACW5340 = Daphnia, row 153
## - Dasybranchus_sp._DH1 = Dasybranchus, row 154
## - Enochrus_ater = Enochrus, row 166
## - Haplothrips_tenuipennis = Haplothrips, row 177
## - Heterolepidoderma_macrops = Heterolepidoderma, row 180
## - Heterolepidoderma_sp._2_TK-2012 = Heterolepidoderma, row 181
## - Hypsibius_cf._dujardini_DS-2016 = Hypsibius, row 189
## - Krenopelopia_sp._BOLD:AAI2213 = Krenopelopia, row 195
## - Macrothrix_sp._HE-364 = Macrothrix, row 206
## - Maxillopoda_sp._BOLD:ACW5478 = Maxillopoda, row 207
## - Maxillopoda_sp._BOLD:ACW5664 = Maxillopoda, row 208
## - Mimeoma_maculata = Scarabaeidae, row 213
## - Pionidae_sp._BOLD:ACE2606 = Pionidae, row 229
## - Podocopida_sp._BOLD:AAG1450 = Podocopida, row 236
## - Podocopida_sp._BOLD:AAH0893 = Podocopida, row 237
## - Podocopida_sp._BOLD:AAH0903 = Podocopida, row 238
## - Podocopida_sp._BOLD:AAH0908 = Podocopida, row 239
## - Podocopida_sp._BOLD:AAH0910 = Podocopida, row 240
## - Polyarthra_dolichoptera_complex_sp._UO-2013 = Polyarthra dolichoptera, row 242
## - Polyarthra_sp._EM-2017 = Polyarthra, row 243
## - Polyarthra_sp._UO-2013 = Polyarthra, row 244
## - Polyarthra_sp._WM-2017a = Polyarthra, row 245
## - Polypedilum_sp._BOLD-2016 = Polypedilum, row 249
## - Procladius_cf._fuscus_BOLD-2016 = Procladius, row 254
## - Psocoptera_sp._BOLD:AAN8452 = Psocoptera, row 258
## - Smittia_sp._8ES = Smittia, row 272
## - Stenostomum_sthenum = Stenostomum, row 274
## - Stylochaeta_scirtetica = Stylochaeta, row 275
## - Synchaeta_cf._tremula/oblonga_UO-2012 = Synchaeta, row 277
## - Tarsonemidae_sp._BOLD:ABV3248 = Tarsonemidae, row 282
## - Trombidiformes_sp._BOLD:ACI3657 = Trombidiformes, row 289
## - Tydeidae_sp._BOLD:ACI5267 = Tydeidae, row 290

eDNA.true$Assignment <- as.character(eDNA.true$Assignment)
eDNA.true[4, "Assignment"] <- "Alona"
eDNA.true[85, "Assignment"] <- "Abacarus"
eDNA.true[86, "Assignment"] <- "Acroperus"
eDNA.true[91, "Assignment"] <- "Anystidae"
eDNA.true[104, "Assignment"] <- "Canthocamptidae"
eDNA.true[113:117, "Assignment"] <- "Chaetonotus"
eDNA.true[118, "Assignment"] <- "Chaitophorus"
eDNA.true[121:123, "Assignment"] <- "Chironomidae"
eDNA.true[129:130, "Assignment"] <- "Chironomus"
eDNA.true[133:134, "Assignment"] <- "Chydorus"
eDNA.true[141:142, "Assignment"] <- "Cricotopus"
eDNA.true[153, "Assignment"] <- "Daphnia"
eDNA.true[154, "Assignment"] <- "Dasybranchus"
eDNA.true[166, "Assignment"] <- "Enochrus"
eDNA.true[177, "Assignment"] <- "Haplothrips"
eDNA.true[180:181, "Assignment"] <- "Heterolepidoderma"
eDNA.true[189, "Assignment"] <- "Hypsibius"
eDNA.true[195, "Assignment"] <- "Krenopelopia"
eDNA.true[206, "Assignment"] <- "Macrothrix"
eDNA.true[207:208, "Assignment"] <- "Maxillopoda"
eDNA.true[213, "Assignment"] <- "Scarabaeidae"
eDNA.true[229, "Assignment"] <- "Pionidae"
eDNA.true[236:240, "Assignment"] <- "Podocopida"
eDNA.true[242, "Assignment"] <- "Polyarthra dolichoptera"
eDNA.true[243:245, "Assignment"] <- "Polyarthra"
eDNA.true[249, "Assignment"] <- "Polypedilum"
eDNA.true[254, "Assignment"] <- "Procladius"
eDNA.true[258, "Assignment"] <- "Psocoptera"
eDNA.true[272, "Assignment"] <- "Smittia"
eDNA.true[274, "Assignment"] <- "Stenostomum"
eDNA.true[275, "Assignment"] <- "Stylochaeta"
eDNA.true[277, "Assignment"] <- "Synchaeta"
eDNA.true[282, "Assignment"] <- "Tarsonemidae"
eDNA.true[289, "Assignment"] <- "Trombidiformes"
eDNA.true[290, "Assignment"] <- "Tydeidae"

## Now merge read counts again by taxonomic assignment. Anything with 
## the same name in the data frame will be merged.
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
#'    samples.
#' 2. Identify the highest level of contamination in PCR positive 
#'    (non-assassin bug DNA) and negative (any DNA) controls.
#' 3. Identify species-specific thresholds using positive controls, i.e. 
#'    the frequency required to remove a given species from the positive 
#'    control as only assassin bug DNA should be present.
#' 
#' Arguably, species-specific thresholds based on PCR positive controls 
#' are more effective as PCR negative controls have no template DNA for 
#' contaminant DNA to compete with for amplification, thus contaminant
#' DNA amplifies exponentially. However, only 11 PCR positive controls were
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
eDNA.ass.bug <- column_to_rownames(eDNA.ass.bug, "Assignment")

## Create new dataframe containing the frequency of reads in each sample
eDNA.ass.bug.freq <- eDNA.ass.bug/c(eDNA.ass.bug[265,])
eDNA.ass.bug.freq[is.na(eDNA.ass.bug.freq)] <- 0

## Add new column to this dataframe containing the maximum frequency of
## DNA from each taxon across samples
eDNA.ass.bug.freq$Threshold <- apply(eDNA.ass.bug.freq, 1, max)

## Check Threshold column has been created properly
head(eDNA.ass.bug.freq[,230:231])

## Combine read frequencies with taxonomic assignment
eDNA.ass.bug.freq <- rownames_to_column(eDNA.ass.bug.freq, "Assignment")

## Print contamination threshold based on max. level of assassin bug
## DNA contamination
max(eDNA.ass.bug.freq[195,232])   # 0%

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
eDNA.pos <- column_to_rownames(eDNA.pos, "Assignment")

## Extract positive controls from raw totals dataframe
pos.total <- raw.total[,which(grepl("Positive", colnames(raw.total)))]

## Bind total read counts to positive controls
eDNA.pos <- rbind(eDNA.pos, pos.total)

## Create new dataframe containing the frequency of reads in each sample
eDNA.pos.freq <- eDNA.pos/c(eDNA.pos[265,])
eDNA.pos.freq[is.na(eDNA.pos.freq)] <- 0

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
eDNA.pos.freq$Threshold <- apply(eDNA.pos.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(eDNA.pos.freq)

## Combine read frequencies with taxonomic assignment
eDNA.pos.freq <- rownames_to_column(eDNA.pos.freq, "Assignment")

## Print contamination threshold based on max. level of non-assassin bug
## DNA contamination
eDNA.pos.freq <- column_to_rownames(eDNA.pos.freq, "Assignment")
max(eDNA.pos.freq[-c(195,265),12])   # 0.0007530717

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
eDNA.pos.test <- eDNA.pos.test[-265,]
sample.total <- raw.total[,which(!grepl("Blank|Empty|Negative|Positive",
                                        colnames(raw.total)))]
eDNA.pos.conversion <- rbind(eDNA.pos.test, sample.total)

## Now convert frequencies back to read counts
eDNA.pos.conversion <- column_to_rownames(eDNA.pos.conversion, "Assignment")
eDNA.pos.FP <- eDNA.pos.conversion*c(eDNA.pos.conversion[265,])

## Remove total row, reset row names and recreate Assignment column
eDNA.pos.FP <- eDNA.pos.FP[-265,]
eDNA.pos.FP <- rownames_to_column(eDNA.pos.FP, "Assignment")

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- eDNA.ass.bug[-265,]
temp$Assignment <- eDNA.pos.FP$Assignment
temp <- temp[,c(231,1:230)]
temp[,2:231][temp[,2:231] > 0] <- 1
eDNA.pos.FP[,2:231][eDNA.pos.FP[,2:231] > 0] <- 1
eDNA.pos1 <- data.frame(colSums(temp[,2:231]))
eDNA.pos2 <- data.frame(colSums(eDNA.pos.FP[,2:231]))
eDNA.pos.compare <- cbind(eDNA.pos1, eDNA.pos2)

## Calculate proportion of information retained
eDNA.pos.compare$proportion <- eDNA.pos.compare[,2]/eDNA.pos.compare[,1]*100
eDNA.pos.compare[is.na(eDNA.pos.compare)] <- 0

## Examine mean and range of proportion of taxa retained
## Exclude LDUN2-04, SKEY1-1-02, THW16-03 and THW16-05, and WADD17-1-08 
## as no taxa were detected in these samples
eDNA.pos.compare <- rownames_to_column(eDNA.pos.compare, "sample_ID")
range(eDNA.pos.compare[-c(19,62,123,125,148),4])
mean(eDNA.pos.compare[-c(19,62,123,125,148),4])

## This would result in up to 78% taxa being removed.
## On average, 34% taxa are removed from samples.


## Check negative controls for highest level of contamination (any DNA). 
## Duplicate negative controls dataframe
eDNA.neg <- eDNA.neg.controls

## Make Assignment column row names
eDNA.neg <- column_to_rownames(eDNA.neg, "Assignment")

## Extract positive controls from raw totals dataframe
neg.total <- raw.total[,which(grepl("Empty|Blank|Negative", 
                                    colnames(raw.total)))]

## Bind total read counts to negative controls
eDNA.neg <- rbind(eDNA.neg, neg.total)

## Create new dataframe containing the frequency of reads in each sample
eDNA.neg.freq <- eDNA.neg/c(eDNA.neg[265,])

## Convert NAs to 0s
eDNA.neg.freq[is.na(eDNA.neg.freq)] <- 0

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
eDNA.neg.freq$Threshold <- apply(eDNA.neg.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(eDNA.neg.freq)

## Combine read frequencies with taxonomic assignment
eDNA.neg.freq <- rownames_to_column(eDNA.neg.freq, "Assignment")

## Print contamination threshold based on max. level of DNA contamination
eDNA.neg.freq <- column_to_rownames(eDNA.neg.freq, "Assignment")
max(eDNA.neg.freq[-c(195,265),99])   # 0.4107143

## In this scenario, any assignments <41.07% total reads in biological 
## samples would be considered contamination. To examine effect of this
## threshold on data, replace all assignments that are less than or equal 
## to this threshold for a given species with zero.

## Create new dataframe containing the frequency of reads in each sample
## and apply threshold
eDNA.neg.test <- eDNA.ass.bug.freq[,-232]
eDNA.neg.test[eDNA.neg.test <= 0.4107143] <- 0

## Now convert back into read counts.
## Remove last row and column containing total frequency and thresholds.
## Add the total read counts to convert assignment frequencies back to 
## read counts for all samples.
eDNA.neg.test <- eDNA.neg.test[-265,]
eDNA.neg.conversion <- rbind(eDNA.neg.test, sample.total)
eDNA.neg.conversion <- column_to_rownames(eDNA.neg.conversion, "Assignment")
eDNA.neg.FP <- eDNA.neg.conversion*c(eDNA.neg.conversion[265,])

## Remove total row, reset row names and recreate Assignment column
eDNA.neg.FP <- eDNA.neg.FP[-265,]
eDNA.neg.FP <- rownames_to_column(eDNA.neg.FP, "Assignment")

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- eDNA.ass.bug[-265,]
temp$Assignment <- eDNA.neg.FP$Assignment
temp <- temp[,c(231,1:230)]
temp[,2:231][temp[,2:231] > 0] <- 1
eDNA.neg.FP[,2:231][eDNA.neg.FP[,2:231] > 0] <- 1
eDNA.neg1 <- data.frame(colSums(temp[,2:231]))
eDNA.neg2 <- data.frame(colSums(eDNA.neg.FP[,2:231]))
eDNA.neg.compare <- cbind(eDNA.neg1, eDNA.neg2)

## Calculate proportion of information retained
eDNA.neg.compare$proportion <- eDNA.neg.compare[,2]/eDNA.neg.compare[,1]*100
eDNA.neg.compare[is.na(eDNA.neg.compare)] <- 0

## Examine mean and range of proportion of taxa retained
## Exclude LDUN2-04, SKEY1-1-02, THW16-03 and THW16-05, and WADD17-1-08 
## as no taxa were detected in these samples
eDNA.neg.compare <- rownames_to_column(eDNA.neg.compare, "sample_ID")
range(eDNA.neg.compare[-c(19,62,123,125,148),4])
mean(eDNA.neg.compare[-c(19,62,123,125,148),4])

## This would result in up to 100% taxa being removed.
## On average, 99% taxa are removed from samples.


#########################################
# OPTION 3: species-specific thresholds #
#########################################

## Create new dataframe containing the frequency of reads in each sample
eDNA.pos.freq <- eDNA.pos/c(eDNA.pos[265,])

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
eDNA.pos.freq$Threshold <- apply(eDNA.pos.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(eDNA.pos.freq)

## Combine read frequencies with taxonomic assignment
eDNA.pos.freq <- rownames_to_column(eDNA.pos.freq, "Assignment")

## Manually change the species threshold value for assassin bug to 0 as 
## this will be a true contaminant in the dataset and no false positive 
## threshold required.
eDNA.pos.freq$Threshold[195] <- 0

## Check this manual edit has occurred
head(eDNA.pos.freq[195,13])

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
eDNA.SS.test <- eDNA.SS.test[-265,]
eDNA.SS.conversion <- rbind(eDNA.SS.test, sample.total)
eDNA.SS.conversion <- column_to_rownames(eDNA.SS.conversion, "Assignment")
eDNA.SS.FP <- eDNA.SS.conversion*c(eDNA.SS.conversion[265,])

## Remove total row, reset row names and recreate Assignment column
eDNA.SS.FP <- eDNA.SS.FP[-265,]
eDNA.SS.FP <- rownames_to_column(eDNA.SS.FP, "Assignment")

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- eDNA.ass.bug[-265,]
temp$Assignment <- eDNA.SS.FP$Assignment
temp <- temp[,c(231,1:230)]
temp[,2:231][temp[,2:231] > 0] <- 1
eDNA.SS.FP[,2:231][eDNA.SS.FP[,2:231] > 0] <- 1
eDNA.SS1 <- data.frame(colSums(temp[,2:231]))
eDNA.SS2 <- data.frame(colSums(eDNA.SS.FP[,2:231]))
eDNA.SS.compare <- cbind(eDNA.SS1, eDNA.SS2)

## Calculate proportion of information retained
eDNA.SS.compare$proportion <- eDNA.SS.compare[,2]/eDNA.SS.compare[,1]*100
eDNA.SS.compare[is.na(eDNA.SS.compare)] <- 0

## Examine mean and range of proportion of taxa retained
## Exclude LDUN2-04, SKEY1-1-02, THW-03 and THW-05, WADD17-1-08 as no taxa 
## were detected in these samples
eDNA.SS.compare <- rownames_to_column(eDNA.SS.compare, "sample_ID")
range(eDNA.SS.compare[-c(19,62,123,125,148),4])
mean(eDNA.SS.compare[-c(19,62,123,125,148),4])

## This would result in up to 9% taxa being removed.
## On average, 1% taxa are removed from samples.


###########
# SUMMARY # 
###########

## Tidy dataframes
colnames(eDNA.pos.compare)[2:4] <- c("spp_richness_NT",
                                     "spp_richness_TA",
                                     "prop_TA")
eDNA.pos.compare$type <- "positive"

colnames(eDNA.neg.compare)[2:4] <- c("spp_richness_NT",
                                     "spp_richness_TA",
                                     "prop_TA")
eDNA.neg.compare$type <- "negative"

colnames(eDNA.SS.compare)[2:4] <- c("spp_richness_NT",
                                    "spp_richness_TA",
                                    "prop_TA")
eDNA.SS.compare$type <- "species-specific"

## Combine dataframes
eDNA.threshold.df <- rbind(eDNA.pos.compare, 
                           eDNA.neg.compare,
                           eDNA.SS.compare)

## Plot number of species retained after thresholds applied
p1a <- ggplot(eDNA.threshold.df, aes(x=sample_ID, 
                                     y=spp_richness_TA, 
                                     fill=type)) +
        geom_bar(stat="identity", position=position_dodge()) + 
        labs(x="Sample", 
             y="Detections remaining after threshold application") + 
        theme_bw() + 
        theme(panel.background = element_rect(fill = 'white'),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.y = element_text(colour="black"),
              legend.position = "none",
              text = element_text(size=20)) + 
        facet_grid(type ~ .)
p1a

## Plot proportion of species retained after thresholds applied
p1b <- ggplot(eDNA.threshold.df, aes(x=sample_ID,
                                     y=prop_TA, 
                                     fill=type)) + 
        geom_bar(stat="identity", position=position_dodge()) + 
        scale_y_continuous(limits=c(0,100)) + 
        labs(x="Sample", 
             y="Detections remaining after threshold application (%)") + 
        theme_bw() + 
        theme(panel.background = element_rect(fill = 'white'),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.y = element_text(colour="black"),
              legend.position = "none",
              text = element_text(size=20)) + 
        facet_grid(type ~ .)
p1b

## Going forward the species-specific threshold will be used. The assassin
## bug threshold indicated we did not have any contamination in the dataset,
## but applying the species-specific threshold as well will retain biological
## information and add an extra level of confidence to the data.
## Recreate dataframe with threshold applied.
eDNA.SS.FP <- eDNA.SS.conversion*c(eDNA.SS.conversion[265,])
eDNA.SS.FP <- eDNA.SS.FP[-265,]
eDNA.SS.FP <- rownames_to_column(eDNA.SS.FP, "Assignment")



#' ---
#'
#' ## 4) Spurious assigments
#'
#' Remove PCR positive control and any non-invertebrate assignments. Also 
#' remove invertebrate assignments higher than family level as these are 
#' too coarse to compare with microscopy. Assignments to be removed are 
#' as follows:
#' 
#' - Arthropoda, row 14
#' - *Arvicola amphibius* (water vole), row 15
#' - Bdelloidea, row 20
#' - *Carassius carassius* (crucian carp), row 31
#' - Chromadorea, row 50
#' - Columbidae, row 55
#' - *Cristatella mucedo* (bryozoan), row 63
#' - Cyclopoida, row 69
#' - Diptera, row 85
#' - *Gallinula chloropus* (common moorhen), row 105
#' - *Gasterosteus aculeatus* (three-spined stickleback), row 106
#' - Hominidae, row 121
#' - *Homo sapiens* (human), row 122
#' - *Hydra* (cnidaria), row 124
#' - *Hydra oligactis* (cnidaria), row 125
#' - *Hydra viridissima* (cnidaria), row 126
#' - Hymenoptera, row 127
#' - Insecta, row 132
#' - *Lissotriton helveticus* (palmate newt), row 150
#' - *Lissotriton vulgaris* (smooth newt), row 151
#' - Maxillopoda, row 157
#' - Metazoa, row 159
#' - Muridae, row 164
#' - *Pica pica* (magpie), row 188
#' - Platyhelminthes, row 193
#' - *Platymeris biguttata* (PCR positive control), row 195
#' - Ploima, row 198
#' - *Plumatella fungosa* (bryozoan), row 199
#' - Podocopida, row 200
#' - Psocoptera, row 218
#' - *Pungitius* (ninespine stickleback), row 220
#' - *Scardinius erythrophthalmus* (Rudd), row 233
#' - Trombidiformes, row 259
#' 

eDNA.refine <- eDNA.SS.FP[-c(14:15,20,31,50,55,63,69,85,105:106,121:122,
                             124:127,132,150:151,157,159,164,188,193,
                             195,198:200,218,220,233,259),]

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
eDNA.pcr.rep <- rownames_to_column(eDNA.pcr.rep, "Sample")

## Create new column containing biological replicate number
eDNA.pcr.rep$pcr_rep <- c(1:12)

## Move new columns to start of dataframe
eDNA.pcr.rep <- eDNA.pcr.rep[,c(1,233,2:232)]

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
write.csv(eDNA.pcr.rep, 
          "../Data/eDNA_PCR_replicates_20200818.csv")
write.csv(eDNA.pcr.rep.metadata, 
          "../Data/eDNA_PCR_replicates_metadata_20200818.csv", 
          row.names=FALSE)


## Now, make dataframe for analysing variation in biological replicates.
## Pool the PCR/sequencing replicates according to biological replicate
## and combine with other biological replicates.

## Remove PCR/sequencing replicate numbers after sample ID
eDNA.to.pool$Sample <- gsub("-01|-02|-03|-04|-05|-06|-07|-08|-09|-10|-11|-12", "",
                            eDNA.to.pool$Sample)

## Remove Replicate column
eDNA.to.pool <- eDNA.to.pool[,-2]

## Merge read counts by biological replicate
## Anything with the same name in the data frame will be merged
## First, have to make sure that data to be summed is read as numeric
eDNA.to.pool[,2:232] <- lapply(eDNA.to.pool[,2:232], function(x) as.numeric(as.character(x)))
eDNA.to.pool <- ddply(eDNA.to.pool, .(Sample), numcolwise(sum))

## Create new column containing biological replicate number
eDNA.to.pool$bio_rep <- c(1:5)

## Move to start of dataframe after Sample column
eDNA.to.pool <- eDNA.to.pool[,c(1,233,2:232)]

## Make dataframe for biological replicates
eDNA.bio.rep <- setNames(data.frame(t(eDNA.samples[,-1])), eDNA.samples[,1])

## Make row names a column in data frame
eDNA.bio.rep <- rownames_to_column(eDNA.bio.rep, "Sample")

## Create new column containing biological replicate number
eDNA.bio.rep$bio_rep <- c(1:5)

## Move new columns to start of dataframe
eDNA.bio.rep <- eDNA.bio.rep[,c(1,233,2:232)]

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
write.csv(eDNA.bio.rep, 
          "../Data/eDNA_biological_replicates_20200818.csv")
write.csv(eDNA.bio.rep.metadata, 
          "../Data/eDNA_biological_replicates_metadata_20200818.csv", 
          row.names=FALSE)


## Now, make dataframe for analysing communities in pooled eDNA samples 
## from ponds with and without crucian carp.

## Remove column containing biological replicate number
eDNA.pooled <- eDNA.pooled[,-2]

## Remove numbers after sample ID
eDNA.pooled$Sample <- gsub("-01|-02|-03|-04|-05", "", eDNA.pooled$Sample)

## Merge read counts by sample
## Anything with the same name in the data frame will be merged
## First, have to make sure that data to be summed is read as numeric
eDNA.pooled[,2:232] <- lapply(eDNA.pooled[,2:232], function(x) as.numeric(as.character(x)))
eDNA.pooled <- ddply(eDNA.pooled, .(Sample), numcolwise(sum))

## For comparison with the DNA metabarcoding data, we must remove
## eDNA samples from pooled dataframe that do not have corresponding 
## bulk tissue samples: GUES1, MYST, OTOM, POFA4, PYES2, WOFA2
eDNA.pooled <- eDNA.pooled[which(!grepl("GUES1|MYST|OTOM|POFA4|PYES2|WOFA2", 
                                        eDNA.pooled$Sample)),]
rownames(eDNA.pooled) <- NULL

## Export as .csv file
write.csv(eDNA.pooled, 
          "../Data/eDNA_metabarcoding_pooled_20200818.csv", 
          row.names=FALSE)

## The metadata for the pooled eDNA samples is the same as that which will be used for 
## netting and is contained in a manually prepared .csv file.



#' ---
#' 
#' 6) Dataframe for family level analysis
#' 
#' Now make dataframe containing family level assignments for biological 
#' replicates pooled by pond. This will be used to analyse the impact of 
#' crucian carp on invertebrate diversity at family level and compared to 
#' results obtained at species level.
#' 

## Store taxa in a vector
taxa <- colnames(eDNA.pooled[,-1])

## Provide NCBI API key
## Sys.setenv(ENTREZ_KEY = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

## Using taxize package, fetch Family information for each taxa
families <- tax_name(sci=taxa, get="family", db="ncbi")

## Correct family for Potamopyrgus antipodarum
families[178, "family"] <- "Hydrobiidae"

## Add missing family information for certain taxa
families[194:199, "family"] <- "Philodinidae"

## Now merge taxonomy with refined species data (DNA.pooled)
eDNA.fam <- setNames(data.frame(t(eDNA.pooled[,-1])), eDNA.pooled[,1])
eDNA.fam <- rownames_to_column(eDNA.fam, "query")
eDNA.fam <- merge(families, eDNA.fam, by="query")

## Remove query and db columns
eDNA.fam <- eDNA.fam[,-c(1:2)]

## Merge read counts by family. Anything with the same name will be 
## merged and anything unique will be retained.
eDNA.fam[2:19] <- lapply(eDNA.fam[,2:19], function(x) as.numeric(as.character(x)))
eDNA.fam <- ddply(eDNA.fam, .(family), numcolwise(sum))

## Transpose data frame
eDNA.fam.pooled <- setNames(data.frame(t(eDNA.fam[,-1])), eDNA.fam[,1])
eDNA.fam.pooled <- rownames_to_column(eDNA.fam.pooled, "Sample")

## Write dataframe to .csv file to be used for analysis in different script
write.csv(eDNA.fam.pooled, 
          "../Data/eDNA_metabarcoding_site_by_family_20200818.csv",
          row.names=FALSE)


#' ---
#' Title: "DNA metabarcoding data filtering and refinement"
#' Author: "Lynsey Rebecca Harper"
#' Date: "18th August 2020"
#' ---
#' 
#' 
#' Bulk tissue and eDNA samples from  ponds in North Norfolk, East of 
#' England, and East Riding of Yorkshire, were screened for freshwater 
#' macroinvertebrates.
#' 
#' Here, the raw DNA metabarcoding data are filtered and refined for 
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

## Create new dataframe for calculating the proportional read counts 
## for each taxon downstream. Duplicate raw read data:
raw1 <- DNA.ass.raw
raw2 <- DNA.unass.raw

## Remove unassigned reads from reference database BLAST results as 
## unassigned reads were extracted for BLAST against entire NCBI nucleotide 
## database. Also remove last column containing taxonomy.
raw1 <- raw1[-141,-386]
raw2 <- raw2[,-386]

## Bind data frames
raw.merged <- rbind(raw1, raw2)

## Merge read counts by taxonomic assignment. Anything with the same name 
## will be merged and anything unique will be retained.
raw.merged[,2:385] <- lapply(raw.merged[,2:385], function(x) as.numeric(as.character(x)))
raw.merged <- ddply(raw.merged, .(Assignment), numcolwise(sum))

## Make Assignment column row names
raw.merged <- column_to_rownames(raw.merged, "Assignment")

## Calculate total number of reads in samples/controls
raw.merged <- rbind(raw.merged, colSums(raw.merged))

## Make new dataframe containing sample ID and total number of reads 
## per sample
raw.total <- raw.merged[336,]
raw.total$Assignment <- "Total"
raw.total <- raw.total[,c(385,1:384)]
rownames(raw.total) <- NULL

## Remove any taxonomic assignments that aren't metazoa from each 
## dataframe using the taxonomy column created during processing with 
## metaBEAT
DNA.ass <- DNA.ass.raw[which(grepl("Metazoa", DNA.ass.raw$taxomomy)),]
DNA.unass <- DNA.unass.raw[which(grepl("Metazoa", DNA.unass.raw$taxomomy)),]

## Remove last column containing taxonomy
DNA.ass <- DNA.ass[,-386]
DNA.unass <- DNA.unass[,-386]

## Bind data frames
DNA.merged <- rbind(DNA.ass, DNA.unass)

## Merge read counts by taxonomic assignment. Anything with the same name 
## will be merged and anything unique will be retained.
DNA.merged[,2:385] <- lapply(DNA.merged[,2:385], function(x) as.numeric(as.character(x)))
DNA.df <- ddply(DNA.merged, .(Assignment), numcolwise(sum))

## Export as .csv file
write.csv(DNA.df, 
          "../Data/DNA_metabarcoding_merged_20200818.csv", 
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
summary(DNA.df[,1:6])
names(DNA.df)
str(DNA.df)

## First, remove spurious assignments:
##
## - invertebrate_environmental_sample, row 209
## 

DNA.true <- DNA.df[-209,]

## Reset row names of data frame for further indexing
rownames(DNA.true) <- NULL

## Now, correct species names:
## 
## - Callicorixa_audeni = Callicorixa, row 17
## - Kurzia_media = Kurzia, row 82
## - Bichromomyia_flaviscutellata = Psychodidae, row 149
## - Chironomidae_sp._BOLD:ACD0662 = Chironomidae, row 154
## - Chironomidae_sp._BOLD:ACI7830 = Chironomidae, row 155
## - Chironomidae_sp._PA2_1 = Chironomidae, row 156
## - Chironomus_muratensis = Chironomus, row 159
## - Chironomus_sp._BOLD:AAI4299 = Chironomus, row 163
## - Cladopelma_sp._BOLD-2016 = Cladopelma, row 167
## - Colymbetes_sp._BMNH_1425212 = Colymbetes fuscus, row 169
## - Dasybranchus_sp._DH1 = Dasybranchus, row 178
## - Enochrus_ater = Enochrus, row 185
## - Erythromma_humerale = Erythromma, row 188
## - Glossiphonia_concolor = Glossiphonia, row 193
## - Haplothrips_tenuipennis = Haplothrips, row 197
## - Helobdella_modesta = Helobdella, row 198
## - Hypoderaeum_sp._Hubei-2014 = Hypoderaeum, row 207
## - Limnesia_(Limnesia)_sp._HP-Hyd013 = Limnesia, row 214
## - Limnesiidae_sp._BOLD:ACJ8659 = Limesiidae, row 216
## - Macropelopia_sp._G_BA30 = Macropelopia, row 218
## - Macrothrix_sp._HE-364 = Macrothrix, row 219
## - Parachironomus_major = Parachironomus, row 229
## - Parachironomus_sp._BOLD:ACB9399 = Parachironomus, row 231
## - Podocopida_sp._BOLD:AAH0903 = Podocopida, row 239
## - Podocopida_sp._BOLD:AAH0910 = Podocopida, row 240
## - Polypedilum_sp._BOLD-2016 = Polypedilum, row 245
## - Procladius_cf._fuscus_BOLD-2016 = Procladius, row 248
## - Stenostomum_sthenum = Stenostomum, row 262
## - Stylochaeta_scirtetica = Stylochaeta, row 264
## - Thermonectus = Dytiscidae, row 276

DNA.true$Assignment <- as.character(DNA.true$Assignment)
DNA.true[17, "Assignment"] <- "Callicorixa"
DNA.true[82, "Assignment"] <- "Kurzia"
DNA.true[149, "Assignment"] <- "Psychodidae"
DNA.true[154:156, "Assignment"] <- "Chironomidae"
DNA.true[c(159,163), "Assignment"] <- "Chironomus"
DNA.true[167, "Assignment"] <- "Cladopelma"
DNA.true[169, "Assignment"] <- "Colymbetes fuscus"
DNA.true[178, "Assignment"] <- "Dasybranchus"
DNA.true[185, "Assignment"] <- "Enochrus"
DNA.true[188, "Assignment"] <- "Erythromma"
DNA.true[193, "Assignment"] <- "Glossiphonia"
DNA.true[197, "Assignment"] <- "Haplothrips"
DNA.true[198, "Assignment"] <- "Helobdella"
DNA.true[207, "Assignment"] <- "Hypoderaeum"
DNA.true[214, "Assignment"] <- "Limnesia"
DNA.true[216, "Assignment"] <- "Limnesiidae"
DNA.true[218, "Assignment"] <- "Macropelopia"
DNA.true[219, "Assignment"] <- "Macrothrix"
DNA.true[c(229,231), "Assignment"] <- "Parachironomus"
DNA.true[239:240, "Assignment"] <- "Podocopida"
DNA.true[245, "Assignment"] <- "Polypedilum"
DNA.true[248, "Assignment"] <- "Procladius"
DNA.true[262, "Assignment"] <- "Stenostomum"
DNA.true[264, "Assignment"] <- "Stylochaeta"
DNA.true[276, "Assignment"] <- "Dytiscidae"

## Now merge read counts again by taxonomic assignment. Anything with 
## the same name in the data frame will be merged.
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
#' DNA amplifies exponentially. However, only 8 PCR positive controls were
#' included on this MiSeq run which may render this approach ineffective.
#'

#########################################################
# OPTION 1: highest level of assassin bug contamination #
#########################################################

## Bind total read counts to refined dataset
DNA.ass.bug <- rbind(DNA.true, raw.total)

## Remove positive and negative controls for threshold determination
DNA.ass.bug <- DNA.ass.bug[,which(grepl("Assignment|-S-|-M-|-L-|-SO-|-UN-|MC0|MC1",
                                        colnames(DNA.ass.bug)))]

## Make Assignment column row names
DNA.ass.bug <- column_to_rownames(DNA.ass.bug, "Assignment")

## Create new dataframe containing the frequency of reads in each sample
DNA.ass.bug.freq  <- DNA.ass.bug/c(DNA.ass.bug[267,])
DNA.ass.bug.freq[is.na(DNA.ass.bug.freq)] <- 0

## Add new column to this dataframe containing the maximum frequency of
## DNA from each taxon across samples
DNA.ass.bug.freq$Threshold <- apply(DNA.ass.bug.freq, 1, max)

## Check Threshold column has been created properly
head(DNA.ass.bug.freq[,335:340])

## Combine read frequencies with taxonomic assignment
DNA.ass.bug.freq <- rownames_to_column(DNA.ass.bug.freq, "Assignment")

## Print contamination threshold based on max. level of assassin bug
## DNA contamination
max(DNA.ass.bug.freq[199,341])   # 0.0002361414

## In this scenario, any assignments <0.02361414% total reads in biological 
## samples would be considered contamination. To examine effect of this
## threshold on data, replace all assignments that are less than or equal 
## to this threshold for a given species with zero.

## Create new dataframe containing the frequency of reads in each sample
## and apply threshold
DNA.ass.bug.test <- DNA.ass.bug.freq
DNA.ass.bug.test[DNA.ass.bug.test <= 0.0002361414] <- 0

## Now convert back into read counts. Remove last row containing 
## frequencies and add the total read counts to convert assignment 
## frequencies back to read counts for all samples.
DNA.ass.bug.test <- DNA.ass.bug.test[-267,-341]
rownames(DNA.ass.bug.test) <- NULL
sample.total <- raw.total[,which(grepl("Assignment|-S-|-M-|-L-|-SO-|-UN-|MC0|MC1",
                                       colnames(raw.total)))]
DNA.ass.bug.conversion <- rbind(DNA.ass.bug.test, sample.total)

## Now convert frequencies back to read counts
DNA.ass.bug.conversion <- column_to_rownames(DNA.ass.bug.conversion, "Assignment")
DNA.ass.bug.FP <- DNA.ass.bug.conversion*c(DNA.ass.bug.conversion[267,])

## Remove total row, reset row names, and recreate Assignment column
DNA.ass.bug.FP <- DNA.ass.bug.FP[-267,]
DNA.ass.bug.FP <- rownames_to_column(DNA.ass.bug.FP, "Assignment")

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- DNA.ass.bug[-267,]
temp$Assignment <- DNA.ass.bug.FP$Assignment
temp <- temp[,c(340,1:339)]
temp[,2:340][temp[,2:340] > 0] <- 1
DNA.ass.bug.FP[,2:340][DNA.ass.bug.FP[,2:340] > 0] <- 1
DNA.ass.bug1 <- data.frame(colSums(temp[,2:340]))
DNA.ass.bug2 <- data.frame(colSums(DNA.ass.bug.FP[,2:340]))
DNA.ass.bug.compare <- cbind(DNA.ass.bug1, DNA.ass.bug2)

## Calculate proportion of information retained
DNA.ass.bug.compare$proportion <- DNA.ass.bug.compare[,2]/DNA.ass.bug.compare[,1]*100
DNA.ass.bug.compare[is.na(DNA.ass.bug.compare)] <- 0

## Examine mean and range of proportion of taxa retained
## Exclude WHALL-M-03 as no taxa were detected in this sample
DNA.ass.bug.compare <- rownames_to_column(DNA.ass.bug.compare, "sample_ID")
range(DNA.ass.bug.compare[-303,4])
mean(DNA.ass.bug.compare[-303,4])

## This would result in up to 62.5% taxa being removed.
## On average, 15% taxa are removed from samples.


########################################################
# OPTION 2: highest level of contamination in controls #
########################################################

## Store controls in new dataframes
DNA.neg.controls <- DNA.true %>% select(Assignment, matches("Empty-|ExtBlank|ExtBlack|Negative"))
DNA.pos.controls <- DNA.true %>% select(Assignment, contains("Positive"))

## Duplicate positive control dataframe
DNA.pos <- DNA.pos.controls

## Make Assignment column row names
DNA.pos <- column_to_rownames(DNA.pos, "Assignment")

## Extract positive controls from raw totals dataframe
pos.total <- raw.total[,which(grepl("Positive", colnames(raw.total)))]

## Bind total read counts to positive controls
DNA.pos <- rbind(DNA.pos, pos.total)

## Create new dataframe containing the frequency of reads in each sample
DNA.pos.freq <- DNA.pos/c(DNA.pos[267,])

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
DNA.pos.freq$Threshold <- apply(DNA.pos.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(DNA.pos.freq)

## Combine read frequencies with taxonomic assignment and make note of 
## which row contains assassin bug (Platymeris biguttata)
DNA.pos.freq <- rownames_to_column(DNA.pos.freq, "Assignment")

## Print contamination threshold based on max. level of non-assassin bug
## DNA contamination
DNA.pos.freq <- column_to_rownames(DNA.pos.freq, "Assignment")
max(DNA.pos.freq[-c(199,267),9])   # 0.0003758551

## In this scenario, any assignments <0.03758551% total reads in biological 
## samples would be considered contamination. To examine effect of this
## threshold on data, replace all assignments that are less than or equal 
## to this threshold for a given species with zero.

## Create new dataframe containing the frequency of reads in each sample
## and apply threshold
DNA.pos.test <- DNA.ass.bug.freq[,-341]
DNA.pos.test[DNA.pos.test <= 0.0003758551] <- 0

## Now convert frequencies back into read counts.
## Remove last row and column containing total frequency and thresholds.
## Add the total read counts to convert assignment frequencies back to 
## read counts for all samples.
DNA.pos.test <- DNA.pos.test[-267,]
DNA.pos.conversion <- rbind(DNA.pos.test, sample.total)
DNA.pos.conversion <- column_to_rownames(DNA.pos.conversion, "Assignment")
DNA.pos.FP <- DNA.pos.conversion*c(DNA.pos.conversion[267,])

## Remove total row, reset row names and recreate Assignment column
DNA.pos.FP <- DNA.pos.FP[-267,]
DNA.pos.FP <- rownames_to_column(DNA.pos.FP, "Assignment")

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- DNA.ass.bug[-267,]
temp$Assignment <- DNA.pos.FP$Assignment
temp <- temp[,c(340,1:339)]
temp[,2:340][temp[,2:340] > 0] <- 1
DNA.pos.FP[,2:340][DNA.pos.FP[,2:340] > 0] <- 1
DNA.pos1 <- data.frame(colSums(temp[,2:340]))
DNA.pos2 <- data.frame(colSums(DNA.pos.FP[,2:340]))
DNA.pos.compare <- cbind(DNA.pos1, DNA.pos2)

## Calculate proportion of information retained
DNA.pos.compare$proportion <- DNA.pos.compare[,2]/DNA.pos.compare[,1]*100
DNA.pos.compare[is.na(DNA.pos.compare)] <- 0

## Examine mean and range of proportion of taxa retained
## Exclude WHALL-M-03 as no taxa were detected in this sample
DNA.pos.compare <- rownames_to_column(DNA.pos.compare, "sample_ID")
range(DNA.pos.compare[-303,4])
mean(DNA.pos.compare[-303,4])

## This would result in up to 84% taxa being removed.
## On average, 21% taxa are removed from samples.


## Check negative controls for highest level of contamination (any DNA). 
## Duplicate negative controls dataframe
DNA.neg <- DNA.neg.controls

## Make Assignment column row names
DNA.neg <- column_to_rownames(DNA.neg, "Assignment")

## Extract positive controls from raw totals dataframe
neg.total <- raw.total[,which(grepl("Empty-|ExtBlank|ExtBlack|Negative", 
                                    colnames(raw.total)))]

## Bind total read counts to negative controls
DNA.neg <- rbind(DNA.neg, neg.total)

## Create new dataframe containing the frequency of reads in each sample
DNA.neg.freq <- DNA.neg/c(DNA.neg[267,])

## Convert NAs to 0s
DNA.neg.freq[is.na(DNA.neg.freq)] <- 0

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
DNA.neg.freq$Threshold <- apply(DNA.neg.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(DNA.neg.freq)

## Combine read frequencies with taxonomic assignment
DNA.neg.freq <- rownames_to_column(DNA.neg.freq, "Assignment")

## Print contamination threshold based on max. level of DNA contamination
DNA.neg.freq <- column_to_rownames(DNA.neg.freq, "Assignment")
max(DNA.neg.freq[-c(199,267),38])   # 100%

## In this scenario, any assignments <100% total reads in biological 
## samples would be considered contamination. This is too extreme and
## thus cannot be used as a threshold.


#########################################
# OPTION 3: species-specific thresholds #
#########################################

## Check positive controls for highest level of contamination belonging to
## each taxonomic assignment. Create new dataframe containing the frequency 
## of reads in each positive control.
DNA.pos.freq <- DNA.pos/c(DNA.pos[267,])

## Add new column to this dataframe containing the maximum frequency for
## DNA from each assignment across all controls
DNA.pos.freq$Threshold <- apply(DNA.pos.freq, 1, max)

## Check maxmimum frequency column has been created properly
head(DNA.pos.freq)

## Combine read frequencies with taxonomic assignment and make note of 
## which row contains assassin bug (Platymeris biguttata)
DNA.pos.freq <- rownames_to_column(DNA.pos.freq, "Assignment")

## Manually change the species threshold value for assassin bug to 0 as 
## this will be a true contaminant in the dataset and no false positive 
## threshold required.
DNA.pos.freq$Threshold[199] <- 0

## Check this manual edit has occurred
head(DNA.pos.freq[199,10])

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
DNA.SS.test <- DNA.SS.test[-267,]
rownames(DNA.SS.test) <- NULL
DNA.SS.conversion <- rbind(DNA.SS.test, sample.total)

## Now convert frequencies back to read counts
DNA.SS.conversion <- column_to_rownames(DNA.SS.conversion, "Assignment")
DNA.SS.FP <- DNA.SS.conversion*c(DNA.SS.conversion[267,])

## Remove total row, reset row names and recreate Assignment column
DNA.SS.FP <- DNA.SS.FP[-267,]
DNA.SS.FP <- rownames_to_column(DNA.SS.FP, "Assignment")

## Compare taxon richness of samples with and without threshold
## Convert sequence read counts to presence-absence
temp <- DNA.ass.bug[-267,]
temp$Assignment <- DNA.SS.FP$Assignment
temp <- temp[,c(340,1:339)]
temp[,2:340][temp[,2:340] > 0] <- 1
DNA.SS.FP[,2:340][DNA.SS.FP[,2:340] > 0] <- 1
DNA.SS1 <- data.frame(colSums(temp[,2:340]))
DNA.SS2 <- data.frame(colSums(DNA.SS.FP[,2:340]))
DNA.SS.compare <- cbind(DNA.SS1, DNA.SS2)

## Calculate proportion of information retained
DNA.SS.compare$proportion <- DNA.SS.compare[,2]/DNA.SS.compare[,1]*100
DNA.SS.compare[is.na(DNA.SS.compare)] <- 0

## Examine mean and range of proportion of taxa retained
## Exclude WHALL-M-03 as no taxa were detected in this sample
DNA.SS.compare <- rownames_to_column(DNA.SS.compare, "sample_ID")
range(DNA.SS.compare[-303,4])
mean(DNA.SS.compare[-303,4])

## This would result in up to 15% taxa being removed.
## On average, 1% taxa are removed from samples.


###########
# SUMMARY # 
###########

## Tidy dataframes
colnames(DNA.ass.bug.compare)[2:4] <- c("spp_richness_NT",
                                        "spp_richness_TA",
                                        "prop_TA")
DNA.ass.bug.compare$type <- "assassin bug"

colnames(DNA.pos.compare)[2:4] <- c("spp_richness_NT",
                                    "spp_richness_TA",
                                    "prop_TA")
DNA.pos.compare$type <- "positive"

colnames(DNA.SS.compare)[2:4] <- c("spp_richness_NT",
                                   "spp_richness_TA",
                                   "prop_TA")
DNA.SS.compare$type <- "species-specific"

## Combine dataframes
DNA.threshold.df <- rbind(DNA.ass.bug.compare, 
                          DNA.pos.compare, 
                          DNA.SS.compare)

## Plot number of species retained after thresholds applied
p1a <- ggplot(dat=DNA.threshold.df, aes(x=sample_ID, 
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
p1b <- ggplot(dat=DNA.threshold.df, aes(x=sample_ID, 
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

## Going forward the assassin bug threshold will be used as this is
## stringent but also retains the majority of biological information.
## Recreate dataframe with threshold applied.
DNA.ass.bug.FP <- DNA.ass.bug.conversion*c(DNA.ass.bug.conversion[267,])
DNA.ass.bug.FP <- DNA.ass.bug.FP[-267,]
DNA.ass.bug.FP <- rownames_to_column(DNA.ass.bug.FP, "Assignment")



#' ---
#'
#' ## 4) Spurious assigments
#'
#' Remove PCR positive control and any non-invertebrate assignments. Also 
#' remove invertebrate assignments higher than family level as these are 
#' too coarse to compare with microscopy. Assignments to be removed are 
#' as follows:
#' 
#' - Arthropoda, row 17
#' - Aves, row 22
#' - Bdelloidea, row 23
#' - Chromadorea, row 39
#' - Coleoptera, row 46
#' - Decapoda, row 68
#' - Diptera, row 72
#' - Hemiptera, row 111
#' - *Homo sapiens*, row 118
#' - *Hydra* (cnidaria), row 121
#' - *Hydra oligactis* (cnidaria), row 122
#' - *Hydra viridissima* (cnidaria), row 123
#' - Insecta, row 139
#' - *Lissotriton vulgaris* (smooth newt), row 156
#' - Metazoa, row 163
#' - *Platymeris biguttata* (PCR positive control), row 199
#' - *Plumatella fungosa* (bryozoan), row 200
#' - Podocopida, row 201
#' - Rallidae (waterfowl), row 223
#' - Rhabditida, row 224
#' - *Scardinius erythrophthalmus* (Rudd), row 230
#' - Spongillida, row 237
#' - Trombidiformes, row 262
#' 

DNA.refine <- DNA.ass.bug.FP[-c(17,22:23,39,46,68,72,111,118,121:123,139,
                                156,163,199:201,223:224,230,237,262),]

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
DNA.pcr.rep <- rownames_to_column(DNA.pcr.rep, "Sample")

## Create new column containing PCR/sequencing replicate number
DNA.pcr.rep$Replicate <- c(1,2,3)

## Move new columns to start of dataframe
DNA.pcr.rep <- DNA.pcr.rep[,c(1,245,2:244)]

## Store Sample and Replicate columns in new dataframe as metadata
## Remove numbers after sample ID
DNA.pcr.rep.metadata <- DNA.pcr.rep[,1:2]
DNA.pcr.rep.metadata$Sample <- gsub("-01|-02|-03", "", DNA.pcr.rep.metadata$Sample)

## Save data in new dataframe for PCR/sequencing replicates to be pooled
DNA.to.pool <- DNA.pcr.rep

## Make Sample column the row names, and remove metadata columns from 
## dataframe
rownames(DNA.pcr.rep) <- DNA.pcr.rep$Sample
DNA.pcr.rep <- DNA.pcr.rep[,-c(1:2)]

## Write dataframes to .csv files to be imported for analysis in a different script
write.csv(DNA.pcr.rep, 
          "../Data/DNA_PCR_replicates_20200818.csv")
write.csv(DNA.pcr.rep.metadata, 
          "../Data/DNA_PCR_replicates_metadata_20200818.csv", 
          row.names=FALSE)


## Now, make dataframe for analysing variation in size categories for each
## bulk tissue sample, and a dataframe for analysing variation in mock 
## communities. Pool the PCR/sequencing replicates according to size 
## category/mock community.

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
mock.comm.grad <- mock.comm.grad[,c(1,246:247,2:245)]
mock.comm.grad$Temperature <- gsub("-01|-02|-03", "", mock.comm.grad$Temperature)

## Store Sample, Community, Temperature and Replicate columns in new dataframe as 
## metadata
mock.comm.grad.metadata <- mock.comm.grad[,1:4]

## Make Sample column the row names, and remove metadata columns from temperature 
## gradient dataframe
rownames(mock.comm.grad) <- mock.comm.grad$Sample
mock.comm.grad <- mock.comm.grad[,-c(1:4)]

## Write dataframes to .csv files to be imported for analysis in a different script
write.csv(mock.comm.grad, 
          "../Data/DNA_mock_community_gradient_20200818.csv")
write.csv(mock.comm.grad.metadata, 
          "../Data/DNA_mock_community_gradient_metadata_20200818.csv",
          row.names=FALSE)


## For different community scenarios dataframe: 
## - store Sample and Replicate columns in new dataframe as metadata
## - remove replicate number after sample ID
## - make Sample column the row names of dataframe with assignments
## - remove the Sample and Replicate columns
mock.comm.metadata <- mock.comm[,1:2]
mock.comm.metadata$Sample <- gsub("-01|-02|-03", "", mock.comm.metadata$Sample)
rownames(mock.comm) <- mock.comm$Sample
mock.comm <- mock.comm[,-c(1:2)]

## Write dataframes to .csv files to be imported for analysis in a different script
write.csv(mock.comm, 
          "../Data/DNA_mock_community_20200818.csv")
write.csv(mock.comm.metadata, 
          "../Data/DNA_mock_community_metadata_20200818.csv", 
          row.names=FALSE)


## Size category dataframe
## Remove PCR/sequencing replicate numbers after sample ID
DNA.size.cat$Sample <- gsub("-01|-02|-03", "", DNA.size.cat$Sample)

## Remove Replicate column
DNA.size.cat <- DNA.size.cat[,-2]

## Merge read counts by size category/mock community
## Anything with the same name in the data frame will be merged
## First, have to make sure that data to be summed is read as numeric
DNA.size.cat[,2:244] <- lapply(DNA.size.cat[,2:244], function(x) as.numeric(as.character(x)))
DNA.size.cat <- ddply(DNA.size.cat, .(Sample), numcolwise(sum))

## Create new column containing size category and move to start of dataframe
## after Sample colum
DNA.size.cat$Category <- gsub("^.*?-","", DNA.size.cat$Sample)
DNA.size.cat <- DNA.size.cat[,c(1,245,2:244)]

## Store Sample and Replicate columns in new dataframe as metadata
## Remove size catogory after sample ID
DNA.size.cat.metadata <- DNA.size.cat[,1:2]
DNA.size.cat.metadata$Sample <- gsub("-L|-M|-S|-SO|-UN", "", 
                                     DNA.size.cat.metadata$Sample)

## Save data in new dataframe for size categories to be pooled
DNA.pooled <- DNA.size.cat

## Make Sample column the row names 
## Remove Sample and Category columns from size categories dataframe
rownames(DNA.size.cat) <- DNA.size.cat$Sample
DNA.size.cat <- DNA.size.cat[,-c(1:2)]

## Write dataframes to .csv files to be imported for analysis in a different script
write.csv(DNA.size.cat, 
          "../Data/DNA_size_categories_20200818.csv")
write.csv(DNA.size.cat.metadata, 
          "../Data/DNA_size_categories_metadata_20200818.csv", 
          row.names=FALSE)


## Finally, to assess how presence of crucian carp affects macroinvertebrate 
## diversity, we need to make dataframe where the communities in S/M/L bulk 
## tissue samples are pooled for each pond. 

## First, subset dataframe for small, medium and large size categories
DNA.pooled <- DNA.pooled[which(!grepl("-SO|-UN", DNA.pooled$Sample)),]

## Remove column containing size category
DNA.pooled <- DNA.pooled[,-2]

## Remove numbers after sample ID
DNA.pooled$Sample <- gsub("-L|-M|-S", "", DNA.pooled$Sample)

## Merge read counts by sample
## Anything with the same name in the data frame will be merged
## First, have to make sure that data to be summed is read as numeric
DNA.pooled[,2:244] <- lapply(DNA.pooled[,2:244], function(x) as.numeric(as.character(x)))
DNA.pooled <- ddply(DNA.pooled, .(Sample), numcolwise(sum))

## Export as .csv file
write.csv(DNA.pooled, 
          "../Data/DNA_metabarcoding_pooled_20200818.csv", 
          row.names=FALSE)

## The metadata for the pooled DNA samples is the same as that which will 
## be used for netting and is contained in a manually prepared .csv file.



#' ---
#' 
#' 6) Dataframe for family level analysis
#' 
#' Now make dataframe containing family level assignments for S/M/L size 
#' categories pooled by pond. This will be used to analyse the impact of 
#' crucian carp on invertebrate diversity at family level and compared to 
#' results obtained at species level.
#' 

## Store taxa in a vector
taxa <- colnames(DNA.pooled[,-1])

## Provide NCBI API key
## Sys.setenv(ENTREZ_KEY = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

## Using taxize package, fetch Family information for each taxa
families <- tax_name(query=taxa, get="family", db="ncbi")

## Correct family for Potamopyrgus antipodarum
families[189, "family"] <- "Hydrobiidae"

## Add family information for rotifers
families[207:209, "family"] <- "Philodinidae"

## Now merge taxonomy with refined species data (DNA.pooled)
DNA.fam <- setNames(data.frame(t(DNA.pooled[,-1])), DNA.pooled[,1])
DNA.fam <- rownames_to_column(DNA.fam, "query")
DNA.fam <- merge(families, DNA.fam, by="query")

## Remove query and db columns
DNA.fam <- DNA.fam[,-c(1:2)]

## Merge read counts by family. Anything with the same name will be 
## merged and anything unique will be retained.
DNA.fam[2:19] <- lapply(DNA.fam[,2:19], function(x) as.numeric(as.character(x)))
DNA.fam <- ddply(DNA.fam, .(family), numcolwise(sum))

## Transpose data frame
DNA.fam.pooled <- setNames(data.frame(t(DNA.fam[,-1])), DNA.fam[,1])
DNA.fam.pooled <- rownames_to_column(DNA.fam.pooled, "Sample")

## Write dataframe to .csv file to be used for analysis in different script
write.csv(DNA.fam.pooled, 
          "../Data/DNA_metabarcoding_site_by_family_20200818.csv", 
          row.names=FALSE)


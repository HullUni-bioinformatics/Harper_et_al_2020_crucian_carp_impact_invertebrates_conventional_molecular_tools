#' ---
#' Title: "Assessing the impact of the threatened crucian carp (*Carassius carassius*) on pond invertebrate diversity - a comparison of conventional and molecular tools"
#' Author: "Lynsey Rebecca Harper"
#' Date: "18th August 2020"
#' ---
#' 
#' 
#' Samples from 20 ponds in North Norfolk, East of England, and 3 ponds 
#' in Selby, East Riding of Yorkshire, were screened for freshwater 
#' invertebrates.
#' 
#' Pond-netting and eDNA metabarcoding were performed in 14 ponds with 
#' crucian carp (*Carassius carassius*) and 10 ponds without crucian 
#' carp. DNA metabarcoding was performed on netted samples from 18 of 
#' 24 ponds. Concordance between methods for invertebrate assessment, 
#' and the effect of crucian carp on alpha and beta diversity will be 
#' examined.
#' 
#' Indices of alpha diversity, such as Fishers Alpha and Shannon diversity, 
#' cannot be used as they take abundance into account. We do not know 
#' whether sequencing data accurately reflects specimen abundance, thus we 
#' treat all data as presence-absence.
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
p <- c("plyr","tidyverse","ggpubr","reshape2","taxize","gtools",
       "vegan","multcomp","RVAideMemoire","arm","Matrix","lme4",
       "car","betapart","ade4","adespatial","adegraphics")
new.packages <- p[!(p %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, 
                                          repos = "https://cran.ma.imperial.ac.uk/",
                                          dependencies=TRUE)
lapply(p, require, character.only = TRUE)

## Load custom functions
f <- c("CheckResidsFunction.R", "OverdispersalFunction.R", 
       "CheckConvergenceFunction.R","HighstatLibV6.R","MyLibrary.r")
lapply(f, source)

#'
#' To ensure reproducibility, print details about the version of R being
#' used for analysis.
#' 

sessionInfo()



###################
# IMPORT DATASETS #
###################

#------------------------#
# ENVIRONMENTAL METADATA #
#------------------------#

## Import data and check structure 
metadata <- read.csv("../Data/Site_metadata.csv", header=TRUE)
summary(metadata)
head(metadata)
names(metadata)
str(metadata)


#------------------------------#
# SWEEP-NETTING AND MICROSCOPY #
#------------------------------#

## Import species-level data and check structure
net.df <- read.csv("../Data/Netting_site_by_species.csv", 
                   row.names=1, 
                   header=TRUE)
summary(net.df)
head(net.df)
names(net.df)
str(net.df)

## Remove '_' from species names
colnames(net.df) <- gsub("_", " ", colnames(net.df))

## Create copy of dataframe without macroinvertebrates that have no 
## reference sequences for method comparison as metabarcoding sequences 
## belonging to these species would not have been taxonomically assigne
net.to.remove <- c("Callicorixa praeusta","Corixa affinis", 
                   "Corixa dentipes","Corixa panzeri",
                   "Hesperocorixa castanea","Hesperocorixa moesta",
                   "Micronecta scholtzi","Plea minutissima",
                   "Sigara distincta", "Sigara dorsalis",
                   "Sigara lateralis","Sigara scotti",
                   "Pisidium milium","Pisidium pseudosphaerium",
                   "Pisidium pulchellum","Sphaerium corneum",
                   "Valvata macrostoma","Erythromma viridulum",
                   "Orthetrum cancellatum")
net.noref <- net.df[,!(colnames(net.df) %in% net.to.remove)]

## Create copy of net.noref to be used to make family-level dataframe
net.fam.noref <- net.noref

## Create vector of taxa to fetch family information for from NCBI
net.to.fetch <- c(colnames(net.fam.noref))

## Remove genus and family-level assignments
net.df <- net.df[,which(grepl(" ", colnames(net.df)))]
net.noref <- net.noref[,which(grepl(" ", colnames(net.noref)))]

## Remove empty taxonomic assignments as these are problematic for vegan
net.df <- net.df[!sapply(net.df, function(x) all(x == 0))]
net.noref <- net.noref[!sapply(net.noref, function(x) all(x == 0))]

## Import family-level data and check structure
net.fam <- read.csv("../Data/Netting_site_by_family.csv", 
                    row.names=1, 
                    header=TRUE)
summary(net.fam)
head(net.fam)
names(net.fam)
str(net.fam)

## Create copy of dataframe without macroinvertebrates that have no 
## reference sequences for method comparison as metabarcoding sequences 
## belonging to these species would not have been taxonomically assigned

## Provide NCBI API key
## Sys.setenv(ENTREZ_KEY = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

## Using taxize package, fetch family information for each taxon
net.families <- tax_name(sci=net.to.fetch, get="family", db="ncbi")

## Correct family for Potamopyrgus antipodarum to match metaBEAT
## output
net.families[79, "family"] <- "Hydrobiidae"

## Add family information where missing
net.families[62, "family"] <- "Dytiscidae"
net.families[65, "family"] <- "Lymnaeidae"
net.families[73, "family"] <- "Sphaeriidae"

## Remove rows containing taxa higher than family-level
net.families <- net.families[complete.cases(net.families), ]

## Now merge taxonomy with net.fam.noref
net.fam.noref <- data.frame(t(net.fam.noref))
net.fam.noref <- rownames_to_column(net.fam.noref, "query")
net.fam.noref <- merge(net.families, net.fam.noref, by="query")

## Remove query and db columns
net.fam.noref <- net.fam.noref[,-c(1:2)]

## Merge read counts by family. Anything with the same name will be 
## merged and anything unique will be retained.
net.fam.noref[2:19] <- lapply(net.fam.noref[,2:19], function(x) as.numeric(as.character(x)))
net.fam.noref <- ddply(net.fam.noref, .(family), numcolwise(sum))

## Transpose data frame
net.fam.noref <- setNames(data.frame(t(net.fam.noref[,-1])), net.fam.noref[,1])

## Remove empty taxonomic assignments as these are problematic for vegan
net.fam <- net.fam[!sapply(net.fam, function(x) all(x == 0))]
net.fam.noref <- net.fam.noref[!sapply(net.fam.noref, function(x) all(x == 0))]


#-------------------#
# DNA METABARCODING #
#-------------------#

## Import species-level data and check structure
DNA.df <- read.csv("../Data/DNA_metabarcoding_pooled_20200818.csv", 
                   row.names=1, 
                   header=TRUE)
summary(DNA.df)
head(DNA.df)
names(DNA.df)
str(DNA.df)

## Remove '.' from sample names
colnames(DNA.df) <- gsub("[.]", " ", colnames(DNA.df))

## Create copy of dataframe without meiofauna for method comparison 
## as these would not have been captured by netting due to net mesh 
## size of 1 mm
DNA.to.remove <- c("Ascomorpha ecaudis","Aulodrilus pluriseta",
                   "Australapatemon","Canthocamptus",
                   "Cypridopsis vidua","Dactylobiotus parthenogeneticus",
                   "Eucyclops","Eucyclops serrulatus",
                   "Keratella","Keratella cochlearis",
                   "Keratella quadrata","Kurzia",
                   "Lecane bulla","Lecane closterocerca","Macrothrix",
                   "Mytilina mucronata","Mytilina ventralis",
                   "Naididae","Nais","Nais christinae","Nais communis",
                   "Nais elinguis","Nais variabilis",
                   "Philodina citrina","Philodina megalotrocha",
                   "Polymerurus rhomboides","Potamothrix bavaricus",
                   "Potamothrix hammoniensis","Rotaria",
                   "Rotaria magnacalcarata","Rotaria rotatoria",
                   "Stenostomum","Stylaria lacustris",
                   "Stylochaeta","Synchaeta",
                   "Synchaeta pectinata","Testudinella caeca",
                   "Trichocerca similis","Trichotria tetractis")
DNA.nomicro <- DNA.df[,!(colnames(DNA.df) %in% DNA.to.remove)]

## Create copy of DNA.nomicro to be used to make family-level dataframe
DNA.fam.nomicro <- DNA.nomicro

## Create vector of taxa to fetch family information for from NCBI
DNA.to.fetch <- c(colnames(DNA.fam.nomicro))

## Before progressing with analysis using vegan, we have to alter the 
## dataframe to remove genus and family-level assignments
DNA.df <- DNA.df[,which(grepl(" ", colnames(DNA.df)))]
DNA.nomicro <- DNA.nomicro[,which(grepl(" ", colnames(DNA.nomicro)))]

## Remove empty taxonomic assignments as these are problematic for vegan
DNA.df <- DNA.df[!sapply(DNA.df, function(x) all(x == 0))]
DNA.nomicro <- DNA.nomicro[!sapply(DNA.nomicro, function(x) all(x == 0))]

## Import family-level data and check structure
DNA.fam <- read.csv("../Data/DNA_metabarcoding_site_by_family_20200818.csv", 
                    row.names=1, 
                    header=TRUE)
summary(DNA.fam)
head(DNA.fam)
names(DNA.fam)
str(DNA.fam)

## Using taxize package, fetch family information for each taxon
DNA.families <- tax_name(sci=DNA.to.fetch, get="family", db="ncbi")

## Correct family for Potamopyrgus antipodarum to match metaBEAT
## output
DNA.families[163, "family"] <- "Hydrobiidae"

## Now merge taxonomy with DNA.fam.nomicro
DNA.fam.nomicro <- data.frame(t(DNA.fam.nomicro))
DNA.fam.nomicro <- rownames_to_column(DNA.fam.nomicro, "query")
DNA.fam.nomicro <- merge(DNA.families, DNA.fam.nomicro, by="query")

## Remove query and db columns
DNA.fam.nomicro <- DNA.fam.nomicro[,-c(1:2)]

## Merge read counts by family. Anything with the same name will be 
## merged and anything unique will be retained.
DNA.fam.nomicro[2:19] <- lapply(DNA.fam.nomicro[,2:19], function(x) as.numeric(as.character(x)))
DNA.fam.nomicro <- ddply(DNA.fam.nomicro, .(family), numcolwise(sum))

## Transpose data frame
DNA.fam.nomicro <- setNames(data.frame(t(DNA.fam.nomicro[,-1])), DNA.fam.nomicro[,1])

## Remove empty taxonomic assignments as these are problematic for vegan
DNA.fam <- DNA.fam[!sapply(DNA.fam, function(x) all(x == 0))]
DNA.fam.nomicro <- DNA.fam.nomicro[!sapply(DNA.fam.nomicro, function(x) all(x == 0))]


#--------------------#
# eDNA METABARCODING #
#--------------------#

## Import species-level data and check structure
eDNA.df <- read.csv("../Data/eDNA_metabarcoding_pooled_20200818.csv", 
                    row.names=1, 
                    header=TRUE)
summary(eDNA.df)
head(eDNA.df)
names(eDNA.df)
str(eDNA.df)

## Remove '.' from sample names
colnames(eDNA.df) <- gsub("[.]", " ", colnames(eDNA.df))

## Create copy of dataframe for method comparison without meiofauna
## as these would not have been captured by netting due to net mesh 
## size of 1 mm
eDNA.to.remove <- c("Abacarus", "Acanthocyclops robustus",
                    "Acroperus","Adineta","Adineta vaga",
                    "Alona","Anystidae","Ascomorpha ecaudis",
                    "Asplanchna sieboldi","Brachionus",
                    "Brachionus calyciflorus","Brachionus quadridentatus",      
                    "Brachionus urceolaris","Canthocamptidae",
                    "Canthocamptus staphylinus","Chaetonotus",
                    "Chydorus","Cypridopsis vidua",
                    "Dactylobiotus parthenogeneticus",
                    "Dermatophagoides pteronyssinus",
                    "Dissotrocha macrostyla","Eubosmina coregoni",
                    "Eucyclops","Eucyclops macrurus",
                    "Eucyclops serrulatus","Eucypris virens",
                    "Eudiaptomus gracilis","Floscularia melicerta",
                    "Fridericia perrieri","Graptoleberis testudinaria",
                    "Heterolepidoderma","Hypsibius","Ilyocryptus agilis",
                    "Keratella","Keratella cochlearis",
                    "Keratella quadrata","Lecane","Lecane closterocerca",
                    "Lepidochaetus zelinkai","Lepidodermella squamata",
                    "Limnocythere inopinata","Macrothrix",
                    "Mytilina mucronata","Mytilina ventralis",
                    "Naididae","Nais","Nais christinae","Nais communis",
                    "Nais elinguis","Nais variabilis",
                    "Philodina citrina","Philodina megalotrocha",
                    "Platyias quadricornis","Plectus aquatilis",
                    "Pleuroxus denticulatus","Polyarthra",
                    "Polyarthra dolichoptera","Polymerurus rhomboides",
                    "Proales daphnicola","Proales fallaciosa",
                    "Rotaria","Rotaria macrura",
                    "Rotaria magnacalcarata","Rotaria neptunia",
                    "Rotaria neptunoida","Rotaria rotatoria",
                    "Scapholeberis mucronata","Stenostomum",
                    "Stylaria lacustris","Stylochaeta",
                    "Synchaeta","Synchaeta pectinata",
                    "Tarsonemidae","Testudinella patina",
                    "Thrips","Thrips tabaci","Trichocerca similis",
                    "Tydeidae","Tyrophagus putrescentiae")
eDNA.nomicro <- eDNA.df[,!(colnames(eDNA.df) %in% eDNA.to.remove)]

## Create copy of eDNA.nomicro to be used to make family-level dataframe
eDNA.fam.nomicro <- eDNA.nomicro

## Create vector of taxa to fetch family information for from NCBI
eDNA.to.fetch <- c(colnames(eDNA.fam.nomicro))

## Before progressing with analysis using vegan, we have to alter the 
## dataframe to contain only species-level assignments
eDNA.df <- eDNA.df[,which(grepl(" ", colnames(eDNA.df)))]
eDNA.nomicro <- eDNA.nomicro[,which(grepl(" ", colnames(eDNA.nomicro)))]

## Remove empty taxonomic assignments as these are problematic for vegan
eDNA.df <- eDNA.df[!sapply(eDNA.df, function(x) all(x == 0))]
eDNA.nomicro <- eDNA.nomicro[!sapply(eDNA.nomicro, function(x) all(x == 0))]

## Import family-level data and check structure
eDNA.fam <- read.csv("../Data/eDNA_metabarcoding_site_by_family_20200818.csv", 
                     row.names=1, 
                     header=TRUE)
summary(eDNA.fam)
head(eDNA.fam)
names(eDNA.fam)
str(eDNA.fam)

## Using taxize package, fetch family information for each taxon
eDNA.families <- tax_name(sci=eDNA.to.fetch, get="family", db="ncbi")

## Correct family for Potamopyrgus antipodarum to match metaBEAT
## output
eDNA.families[120, "family"] <- "Hydrobiidae"

## Now merge taxonomy with DNA.fam.nomicro
eDNA.fam.nomicro <- data.frame(t(eDNA.fam.nomicro))
eDNA.fam.nomicro <- rownames_to_column(eDNA.fam.nomicro, "query")
eDNA.fam.nomicro <- merge(eDNA.families, eDNA.fam.nomicro, by="query")

## Remove query and db columns
eDNA.fam.nomicro <- eDNA.fam.nomicro[,-c(1:2)]

## Merge read counts by family. Anything with the same name will be 
## merged and anything unique will be retained.
eDNA.fam.nomicro[2:19] <- lapply(eDNA.fam.nomicro[,2:19], function(x) as.numeric(as.character(x)))
eDNA.fam.nomicro <- ddply(eDNA.fam.nomicro, .(family), numcolwise(sum))

## Transpose data frame
eDNA.fam.nomicro <- setNames(data.frame(t(eDNA.fam.nomicro[,-1])), eDNA.fam.nomicro[,1])

## Remove empty taxonomic assignments as these are problematic for vegan
eDNA.fam <- eDNA.fam[!sapply(eDNA.fam, function(x) all(x == 0))]
eDNA.fam.nomicro <- eDNA.fam.nomicro[!sapply(eDNA.fam.nomicro, function(x) all(x == 0))]


#-------------------------------------#
# COMBINED DATASETS (NO TAXA REMOVED) #
#-------------------------------------#

## First, replicate all species-level dataframes
netting <- net.df
DNA.meta <- DNA.df
eDNA.meta <- eDNA.df

## Add columns containing pond ID, whether crucian carp are present in 
## pond, and which tool was used to generate data 
netting$Sample <- rownames(netting)
netting$Crucian <- metadata$Crucian
netting$Method <- c("Netting")
netting <- netting[,c(92:94,1:91)]
rownames(netting) <- NULL

DNA.meta$Sample <- rownames(DNA.meta)
DNA.meta$Crucian <- metadata$Crucian
DNA.meta$Method <- c("DNA")
DNA.meta <- DNA.meta[,c(132:134,1:131)]
rownames(DNA.meta) <- NULL

eDNA.meta$Sample <- rownames(eDNA.meta)
eDNA.meta$Crucian <- metadata$Crucian
eDNA.meta$Method <- c("eDNA")
eDNA.meta <- eDNA.meta[,c(146:148,1:145)]
rownames(eDNA.meta) <- NULL

## Combine data frames
comb.df <- smartbind(netting, DNA.meta, eDNA.meta)

## Replace NAs with 0s
comb.df[is.na(comb.df)] <- 0

## Take the first three columns and store in new dataframe as metadata
method.metadata <- comb.df[,1:3]
rownames(method.metadata) <- NULL

## Add methodology to pond ID in dataframe of assignments
comb.df <- within(comb.df, id <- paste(Sample, Method, sep="_"))
comb.df <- comb.df[,c(253,1:252)]

## Remove metadata columns and rename ID column
comb.df <- comb.df[,-c(2:4)]
colnames(comb.df)[1] <- "Sample"

## Save dataframe in current form for comparing monitoring tools later
## Remove empty taxonomic assignments
method.df <- comb.df
rownames(method.df) <- method.df$Sample
method.df <- method.df[,-1]
method.df <- method.df[!sapply(method.df, function(x) all(x == 0))]

## Now we can pool samples in comb.df
## Remove methodology from sample name
comb.df$Sample <- gsub("_Netting|_DNA|_eDNA", "", comb.df$Sample)
comb.df[,2:250] <- lapply(comb.df[,2:250], function(x) as.numeric(as.character(x)))
comb.df <- ddply(comb.df, .(Sample), numcolwise(sum))

## Make Sample column row names
comb.df <- column_to_rownames(comb.df, "Sample")

## Remove empty taxonomic assignments
comb.df <- comb.df[!sapply(comb.df, function(x) all(x == 0))]

## Make copy of dataframe for analysing individual invertebrate orders 
## in ponds with and without crucian carp later
comb.order.spp <- comb.df

## Now, replicate all family-level dataframes
netting.fam <- net.fam
DNA.meta.fam <- DNA.fam
eDNA.meta.fam <- eDNA.fam

## Add columns containing pond ID, whether crucian carp are present in pond,
## and which tool was used to generate data
netting.fam$Sample <- rownames(netting.fam)
netting.fam$Crucian <- metadata$Crucian
netting.fam$Method <- c("Netting")
netting.fam <- netting.fam[,c(39:41,1:38)]
rownames(netting.fam) <- NULL

DNA.meta.fam$Sample <- rownames(DNA.meta.fam)
DNA.meta.fam$Crucian <- metadata$Crucian
DNA.meta.fam$Method <- c("DNA")
DNA.meta.fam <- DNA.meta.fam[,c(56:58,1:55)]
rownames(DNA.meta.fam) <- NULL

eDNA.meta.fam$Sample <- rownames(eDNA.meta.fam)
eDNA.meta.fam$Crucian <- metadata$Crucian
eDNA.meta.fam$Method <- c("eDNA")
eDNA.meta.fam <- eDNA.meta.fam[,c(91:93,1:90)]
rownames(eDNA.meta.fam) <- NULL

## Combine data frames
comb.fam.df <- smartbind(netting.fam, DNA.meta.fam, eDNA.meta.fam)

## Replace NAs with 0s
comb.fam.df[is.na(comb.fam.df)] <- 0

## Take the first three columns and store in new dataframe as metadata
method.fam.metadata <- comb.fam.df[,1:3]
rownames(method.fam.metadata) <- NULL

## Add methodology to pond ID in dataframe of assignments
comb.fam.df <- within(comb.fam.df, id <- paste(Sample, Method, sep="_"))
comb.fam.df <- comb.fam.df[,c(113,1:112)]

## Remove metadata columns and rename ID column
comb.fam.df <- comb.fam.df[,-c(2:4)]
colnames(comb.fam.df)[1] <- "Sample"

## Save dataframe in current form for comparing monitoring tools later
## Remove empty taxonomic assignments
method.fam.df <- comb.fam.df
rownames(method.fam.df) <- method.fam.df$Sample
method.fam.df <- method.fam.df[,-1]
method.fam.df <- method.fam.df[!sapply(method.fam.df, function(x) all(x == 0))]

## Now we can pool samples in comb.fam.df
## Remove methodology from sample name
comb.fam.df$Sample <- gsub("_Netting|_DNA|_eDNA", "", comb.fam.df$Sample)
comb.fam.df[,2:110] <- lapply(comb.fam.df[,2:110], function(x) as.numeric(as.character(x)))
comb.fam.df <- ddply(comb.fam.df, .(Sample), numcolwise(sum))

## Make Sample column row names
comb.fam.df <- column_to_rownames(comb.fam.df, "Sample")

## Remove empty taxonomic assignments
comb.fam.df <- comb.fam.df[!sapply(comb.fam.df, function(x) all(x == 0))]

## Make copy of dataframe for analysing individual invertebrate orders 
## in ponds with and without crucian carp later
comb.order.fam <- comb.fam.df


#------------------------------------------------------------------#
# COMBINED DATASETS (MEIOFAUNA/TAXA MISSING REF SEQUENCES REMOVED) #
#------------------------------------------------------------------#

## First, replicate all species-level dataframes
netting.noref <- net.noref
DNA.meta.nomicro <- DNA.nomicro
eDNA.meta.nomicro <- eDNA.nomicro

## Add columns containing pond ID, whether crucian carp are present in pond,
## and which tool was used to generate data
netting.noref$Sample <- rownames(netting.noref)
netting.noref$Crucian <- metadata$Crucian
netting.noref$Method <- c("Netting")
netting.noref <- netting.noref[,c(73:75,1:72)]
rownames(netting.noref) <- NULL

DNA.meta.nomicro$Sample <- rownames(DNA.meta.nomicro)
DNA.meta.nomicro$Crucian <- metadata$Crucian
DNA.meta.nomicro$Method <- c("DNA")
DNA.meta.nomicro <- DNA.meta.nomicro[,c(119:121,1:118)]
rownames(DNA.meta.nomicro) <- NULL

eDNA.meta.nomicro$Sample <- rownames(eDNA.meta.nomicro)
eDNA.meta.nomicro$Crucian <- metadata$Crucian
eDNA.meta.nomicro$Method <- c("eDNA")
eDNA.meta.nomicro <- eDNA.meta.nomicro[,c(94:96,1:93)]
rownames(eDNA.meta.nomicro) <- NULL

## Combine data frames
comb.nomicro <- smartbind(netting.noref, DNA.meta.nomicro, eDNA.meta.nomicro)

## Replace NAs with 0s
comb.nomicro[is.na(comb.nomicro)] <- 0

## Take the first three columns and store in new dataframe as metadata
method.nomicro.metadata <- comb.nomicro[,1:3]
rownames(method.nomicro.metadata) <- NULL

## Add methodology to pond ID in dataframe of assignments
comb.nomicro <- within(comb.nomicro, id <- paste(Sample, Method, sep="_"))
comb.nomicro <- comb.nomicro[,c(180,1:179)]

## Remove metadata columns and rename ID column
comb.nomicro <- comb.nomicro[,-c(2:4)]
colnames(comb.nomicro)[1] <- "Sample"

## Save dataframe in current form for comparing monitoring tools later
## Remove empty taxonomic assignments
method.nomicro <- comb.nomicro
rownames(method.nomicro) <- method.nomicro$Sample
method.nomicro <- method.nomicro[,-1]
method.nomicro <- method.nomicro[!sapply(method.nomicro, function(x) all(x == 0))]

## Now we can pool samples in comb.nomicro
## Remove methodology from sample name
comb.nomicro$Sample <- gsub("_Netting|_DNA|_eDNA", "", comb.nomicro$Sample)
comb.nomicro[,2:177] <- lapply(comb.nomicro[,2:177], function(x) as.numeric(as.character(x)))
comb.nomicro <- ddply(comb.nomicro, .(Sample), numcolwise(sum))

## Make Sample column row names
comb.nomicro <- column_to_rownames(comb.nomicro, "Sample")

## Remove empty taxonomic assignments
comb.nomicro <- comb.nomicro[!sapply(comb.nomicro, function(x) all(x == 0))]

## Make copy of dataframe for analysing individual invertebrate orders in 
## ponds with and without crucian carp later
comb.nomicro.order.spp <- comb.nomicro

## Now, replicate all family-level dataframes
netting.fam.noref <- net.fam.noref
DNA.meta.fam.nomicro <- DNA.fam.nomicro
eDNA.meta.fam.nomicro <- eDNA.fam.nomicro

## Add columns containing pond ID, whether crucian carp are present in pond,
## and which tool was used to generate data
netting.fam.noref$Sample <- rownames(netting.fam.noref)
netting.fam.noref$Crucian <- metadata$Crucian
netting.fam.noref$Method <- c("Netting")
netting.fam.noref <- netting.fam.noref[,c(37:39,1:36)]
rownames(netting.fam.noref) <- NULL

DNA.meta.fam.nomicro$Sample <- rownames(DNA.meta.fam.nomicro)
DNA.meta.fam.nomicro$Crucian <- metadata$Crucian
DNA.meta.fam.nomicro$Method <- c("DNA")
DNA.meta.fam.nomicro <- DNA.meta.fam.nomicro[,c(47:49,1:46)]
rownames(DNA.meta.fam.nomicro) <- NULL

eDNA.meta.fam.nomicro$Sample <- rownames(eDNA.meta.fam.nomicro)
eDNA.meta.fam.nomicro$Crucian <- metadata$Crucian
eDNA.meta.fam.nomicro$Method <- c("eDNA")
eDNA.meta.fam.nomicro <- eDNA.meta.fam.nomicro[,c(62:64,1:61)]
rownames(eDNA.meta.fam.nomicro) <- NULL

## Combine data frames
comb.fam.nomicro <- smartbind(netting.fam.noref, DNA.meta.fam.nomicro, eDNA.meta.fam.nomicro)

## Replace NAs with 0s
comb.fam.nomicro[is.na(comb.fam.nomicro)] <- 0

## Take the first three columns and store in new dataframe as metadata
method.fam.nomicro.metadata <- comb.fam.nomicro[,1:3]
rownames(method.fam.nomicro.metadata) <- NULL

## Add methodology to pond ID in dataframe of assignments
comb.fam.nomicro <- within(comb.fam.nomicro, id <- paste(Sample, Method, sep="_"))
comb.fam.nomicro <- comb.fam.nomicro[,c(81,1:80)]

## Remove metadata columns and rename ID column
comb.fam.nomicro <- comb.fam.nomicro[,-c(2:4)]
colnames(comb.fam.nomicro)[1] <- "Sample"

## Save dataframe in current form for comparing monitoring tools later
method.fam.nomicro <- comb.fam.nomicro
rownames(method.fam.nomicro) <- method.fam.nomicro$Sample
method.fam.nomicro <- method.fam.nomicro[,-1]
method.fam.nomicro <- method.fam.nomicro[!sapply(method.fam.nomicro, function(x) all(x == 0))]

## Now we can pool samples in comb.fam.df
## Remove methodology from sample name
comb.fam.nomicro$Sample <- gsub("_Netting|_DNA|_eDNA", "", comb.fam.nomicro$Sample)
comb.fam.nomicro[,2:78] <- lapply(comb.fam.nomicro[,2:78], function(x) as.numeric(as.character(x)))
comb.fam.nomicro <- ddply(comb.fam.nomicro, .(Sample), numcolwise(sum))

## Make Sample column row names
comb.fam.nomicro <- column_to_rownames(comb.fam.nomicro, "Sample")

## Remove empty taxonomic assignments
comb.fam.nomicro <- comb.fam.nomicro[!sapply(comb.fam.nomicro, function(x) all(x == 0))]

## Make copy of dataframe for analysing individual invertebrate orders 
## in ponds with and without crucian carp later
comb.nomicro.order.fam <- comb.fam.nomicro



#######################################
# METHOD COMPARISON (NO TAXA REMOVED) #
#######################################

#--------------------------#
# VISUALISE METHOD OVERLAP #
#--------------------------#

## Tranform data from counts (specimens or reads) to presence-absence by
## rows (ponds)
net.spp.pa <- decostand(net.df, method = "pa", MARGIN = 1)
head(net.spp.pa[,1:3], n = 3)

net.fam.pa <- decostand(net.fam, method = "pa", MARGIN = 1)
head(net.fam.pa[,1:3], n = 3)

DNA.spp.pa <- decostand(DNA.df, method = "pa", MARGIN = 1)
head(DNA.spp.pa[,1:3], n = 3)

DNA.fam.pa <- decostand(DNA.fam, method = "pa", MARGIN = 1)
head(DNA.fam.pa[,1:3], n = 3)

eDNA.spp.pa <- decostand(eDNA.df, method = "pa", MARGIN = 1)
head(eDNA.spp.pa[,1:3], n = 3)

eDNA.fam.pa <- decostand(eDNA.fam, method = "pa", MARGIN = 1)
head(eDNA.fam.pa[,1:3], n = 3)


## METHOD OVERLAP
## Produce venn diagrams of overlap in taxonomic assignments each method

## Code to install bioconductor and limma on first use:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install("limma")
library("limma")

## Species-level
## Create dataframe for each method containing species name and total 
## number of ponds it was found in
net.total <- data.frame(colnames(net.spp.pa), colSums(net.spp.pa))
colnames(net.total) <- c("Species", "Netting")

DNA.total <- data.frame(colnames(DNA.spp.pa), colSums(DNA.spp.pa))
colnames(DNA.total) <- c("Species", "DNA metabarcoding")

eDNA.total <- data.frame(colnames(eDNA.spp.pa), colSums(eDNA.spp.pa))
colnames(eDNA.total) <- c("Species", "eDNA metabarcoding")

## Make list of dataframes and bind them together
## Empty rows will be filled with NA values using rbind.fill() function
spp.totals <- list(net.total, DNA.total, eDNA.total)
spp.df <- do.call(rbind.fill, spp.totals)

## Replace NA values with 0 and merge rows by species
spp.df[is.na(spp.df)] <- 0
spp.df <- ddply(spp.df, .(Species), numcolwise(sum))

## Write as .csv file and add column containing taxonomic order for 
## species assignments
#write.csv(spp.df, "../Data/venn_species_20200818.csv", 
#          row.names=FALSE)

## Reimport file
spp.df <- read.csv("../Data/venn_species_edited_20200818.csv", 
                   header=TRUE)

## Venn diagram of total species detections across methods
## Take species detections by method and make new presence-absence dataframe
venn.df <- spp.df[,3:5]
venn.df[venn.df > 0] <- 1

## Make venn diagram
spp <- vennCounts(venn.df)
svg(filename="./venn_diagrams_inc-all-taxa/total_species_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

## Venn diagram of species detections by taxonomic group across methods
## Take species detections by method and make new presence-absence dataframe
groups.df <- spp.df[,-1]
groups.df[groups.df > 0] <- 1

## Subset new dataframe by major invertebrate groups and make individual
## venn diagrams
Annelida <- subset(groups.df, Group == "Annelida")
Annelida <- Annelida[,-1]
spp <- vennCounts(Annelida)
svg(filename="./venn_diagrams_inc-all-taxa/Annelida_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Arachnida <- subset(groups.df, Group == "Arachnida")
Arachnida <- Arachnida[,-1]
spp <- vennCounts(Arachnida)
svg(filename="./venn_diagrams_inc-all-taxa/Arachnida_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Coleoptera <- subset(groups.df, Group == "Coleoptera")
Coleoptera <- Coleoptera[,-1]
spp <- vennCounts(Coleoptera)
svg(filename="./venn_diagrams_inc-all-taxa/Coleoptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Collembola <- subset(groups.df, Group == "Collembola")
Collembola <- Collembola[,-1]
spp <- vennCounts(Collembola)
svg(filename="./venn_diagrams_inc-all-taxa/Collembola_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Crustacea <- subset(groups.df, Group == "Crustacea")
Crustacea <- Crustacea[,-1]
spp <- vennCounts(Crustacea)
svg(filename="./venn_diagrams_inc-all-taxa/Crustacea_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Diptera <- subset(groups.df, Group == "Diptera")
Diptera <- Diptera[,-1]
spp <- vennCounts(Diptera)
svg(filename="./venn_diagrams_inc-all-taxa/Diptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Ephemeroptera <- subset(groups.df, Group == "Ephemeroptera")
Ephemeroptera <- Ephemeroptera[,-1]
spp <- vennCounts(Ephemeroptera)
svg(filename="./venn_diagrams_inc-all-taxa/Ephemeroptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Gastrotricha <- subset(groups.df, Group == "Gastrotricha")
Gastrotricha <- Gastrotricha[,-1]
spp <- vennCounts(Gastrotricha)
svg(filename="./venn_diagrams_inc-all-taxa/Gastrotricha_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Hemiptera <- subset(groups.df, Group == "Hemiptera")
Hemiptera <- Hemiptera[,-1]
spp <- vennCounts(Hemiptera)
svg(filename="./venn_diagrams_inc-all-taxa/Hemiptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Hirudinea <- subset(groups.df, Group == "Hirudinea")
Hirudinea <- Hirudinea[,-1]
spp <- vennCounts(Hirudinea)
svg(filename="./venn_diagrams_inc-all-taxa/Hirudinea_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Hymenoptera <- subset(groups.df, Group == "Hymenoptera")
Hymenoptera <- Hymenoptera[,-1]
spp <- vennCounts(Hymenoptera)
svg(filename="./venn_diagrams_inc-all-taxa/Hymenoptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Lepidoptera <- subset(groups.df, Group == "Lepidoptera")
Lepidoptera <- Lepidoptera[,-1]
spp <- vennCounts(Lepidoptera)
svg(filename="./venn_diagrams_inc-all-taxa/Lepidoptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Megaloptera <- subset(groups.df, Group == "Megaloptera")
Megaloptera <- Megaloptera[,-1]
spp <- vennCounts(Megaloptera)
svg(filename="./venn_diagrams_inc-all-taxa/Megaloptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Mollusca <- subset(groups.df, Group == "Mollusca")
Mollusca <- Mollusca[,-1]
spp <- vennCounts(Mollusca)
svg(filename="./venn_diagrams_inc-all-taxa/Mollusca_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Nematoda <- subset(groups.df, Group == "Nematoda")
Nematoda <- Nematoda[,-1]
spp <- vennCounts(Nematoda)
svg(filename="./venn_diagrams_inc-all-taxa/Nematoda_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Odonata <- subset(groups.df, Group == "Odonata")
Odonata <- Odonata[,-1]
spp <- vennCounts(Odonata)
svg(filename="./venn_diagrams_inc-all-taxa/Odonata_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Platyhelminthes <- subset(groups.df, Group == "Platyhelminthes")
Platyhelminthes <- Platyhelminthes[,-1]
spp <- vennCounts(Platyhelminthes)
svg(filename="./venn_diagrams_inc-all-taxa/Platyhelminthes_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Psocoptera <- subset(groups.df, Group == "Psocoptera")
Psocoptera <- Psocoptera[,-1]
spp <- vennCounts(Psocoptera)
svg(filename="./venn_diagrams_inc-all-taxa/Psocoptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Rotifera <- subset(groups.df, Group == "Rotifera")
Rotifera <- Rotifera[,-1]
spp <- vennCounts(Rotifera)
svg(filename="./venn_diagrams_inc-all-taxa/Rotifera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Tardigrada <- subset(groups.df, Group == "Tardigrada")
Tardigrada <- Tardigrada[,-1]
spp <- vennCounts(Tardigrada)
svg(filename="./venn_diagrams_inc-all-taxa/Tardigrada_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Trichoptera <- subset(groups.df, Group == "Trichoptera")
Trichoptera <- Trichoptera[,-1]
spp <- vennCounts(Trichoptera)
svg(filename="./venn_diagrams_inc-all-taxa/Trichoptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()


## Family-level
## Create dataframe for each method containing family name and total number
## of ponds it was found in
net.fam.total <- data.frame(colnames(net.fam.pa), colSums(net.fam.pa))
colnames(net.fam.total) <- c("Family", "Netting")

DNA.fam.total <- data.frame(colnames(DNA.fam.pa), colSums(DNA.fam.pa))
colnames(DNA.fam.total) <- c("Family", "DNA metabarcoding")

eDNA.fam.total <- data.frame(colnames(eDNA.fam.pa), colSums(eDNA.fam.pa))
colnames(eDNA.fam.total) <- c("Family", "eDNA metabarcoding")

## Make list of dataframes and bind them together
## Empty rows will be filled with NA values using rbind.fill() function
fam.totals <- list(net.fam.total, DNA.fam.total, eDNA.fam.total)
fam.df <- do.call(rbind.fill, fam.totals)

## Replace NA values with 0 and merge rows by species
fam.df[is.na(fam.df)] <- 0
fam.df <- ddply(fam.df, .(Family), numcolwise(sum))

## Write as .csv file and add column containing taxonomic order for species
## assignments
#write.csv(fam.df, "../Data/venn_family_20200818.csv", 
#          row.names=FALSE)

## Reimport file
fam.df <- read.csv("../Data/venn_family_edited_20200818.csv", 
                   header=TRUE)

## Venn diagram of total family detections across methods
## Take species detections by method and make new presence-absence dataframe
venn.fam.df <- fam.df[,3:5]
venn.fam.df[venn.fam.df > 0] <- 1

## Make venn diagram
fam <- vennCounts(venn.fam.df)
svg(filename="./venn_diagrams_inc-all-taxa/total_family_method.svg", width = 13, height = 10)
vennDiagram(fam, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

## Venn diagram of species detections by taxonomic group across methods
## Take species detections by method and make new presence-absence dataframe
groups.fam.df <- fam.df[,-1]
groups.fam.df[groups.fam.df > 0] <- 1

## Subset new dataframe by major invertebrate groups and make individual
## venn diagrams
Annelida <- subset(groups.fam.df, Group == "Annelida")
Annelida <- Annelida[,-1]
spp <- vennCounts(Annelida)
svg(filename="./venn_diagrams_inc-all-taxa/Annelida_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Arachnida <- subset(groups.fam.df, Group == "Arachnida")
Arachnida <- Arachnida[,-1]
spp <- vennCounts(Arachnida)
svg(filename="./venn_diagrams_inc-all-taxa/Arachnida_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Coleoptera <- subset(groups.fam.df, Group == "Coleoptera")
Coleoptera <- Coleoptera[,-1]
spp <- vennCounts(Coleoptera)
svg(filename="./venn_diagrams_inc-all-taxa/Coleoptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Collembola <- subset(groups.fam.df, Group == "Collembola")
Collembola <- Collembola[,-1]
spp <- vennCounts(Collembola)
svg(filename="./venn_diagrams_inc-all-taxa/Collembola_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Crustacea <- subset(groups.fam.df, Group == "Crustacea")
Crustacea <- Crustacea[,-1]
spp <- vennCounts(Crustacea)
svg(filename="./venn_diagrams_inc-all-taxa/Crustacea_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Diptera <- subset(groups.fam.df, Group == "Diptera")
Diptera <- Diptera[,-1]
spp <- vennCounts(Diptera)
svg(filename="./venn_diagrams_inc-all-taxa/Diptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Ephemeroptera <- subset(groups.fam.df, Group == "Ephemeroptera")
Ephemeroptera <- Ephemeroptera[,-1]
spp <- vennCounts(Ephemeroptera)
svg(filename="./venn_diagrams_inc-all-taxa/Ephemeroptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Gastrotricha <- subset(groups.fam.df, Group == "Gastrotricha")
Gastrotricha <- Gastrotricha[,-1]
spp <- vennCounts(Gastrotricha)
svg(filename="./venn_diagrams_inc-all-taxa/Gastrotricha_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Hemiptera <- subset(groups.fam.df, Group == "Hemiptera")
Hemiptera <- Hemiptera[,-1]
spp <- vennCounts(Hemiptera)
svg(filename="./venn_diagrams_inc-all-taxa/Hemiptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Hirudinea <- subset(groups.fam.df, Group == "Hirudinea")
Hirudinea <- Hirudinea[,-1]
spp <- vennCounts(Hirudinea)
svg(filename="./venn_diagrams_inc-all-taxa/Hirudinea_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Hymenoptera <- subset(groups.fam.df, Group == "Hymenoptera")
Hymenoptera <- Hymenoptera[,-1]
spp <- vennCounts(Hymenoptera)
svg(filename="./venn_diagrams_inc-all-taxa/Hymenoptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Lepidoptera <- subset(groups.fam.df, Group == "Lepidoptera")
Lepidoptera <- Lepidoptera[,-1]
spp <- vennCounts(Lepidoptera)
svg(filename="./venn_diagrams_inc-all-taxa/Lepidoptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Megaloptera <- subset(groups.fam.df, Group == "Megaloptera")
Megaloptera <- Megaloptera[,-1]
spp <- vennCounts(Megaloptera)
svg(filename="./venn_diagrams_inc-all-taxa/Megaloptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Mollusca <- subset(groups.fam.df, Group == "Mollusca")
Mollusca <- Mollusca[,-1]
spp <- vennCounts(Mollusca)
svg(filename="./venn_diagrams_inc-all-taxa/Mollusca_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Nematoda <- subset(groups.fam.df, Group == "Nematoda")
Nematoda <- Nematoda[,-1]
spp <- vennCounts(Nematoda)
svg(filename="./venn_diagrams_inc-all-taxa/Nematoda_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Odonata <- subset(groups.fam.df, Group == "Odonata")
Odonata <- Odonata[,-1]
spp <- vennCounts(Odonata)
svg(filename="./venn_diagrams_inc-all-taxa/Odonata_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Platyhelminthes <- subset(groups.fam.df, Group == "Platyhelminthes")
Platyhelminthes <- Platyhelminthes[,-1]
spp <- vennCounts(Platyhelminthes)
svg(filename="./venn_diagrams_inc-all-taxa/Platyhelminthes_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Psocoptera <- subset(groups.fam.df, Group == "Psocoptera")
Psocoptera <- Psocoptera[,-1]
spp <- vennCounts(Psocoptera)
svg(filename="./venn_diagrams_inc-all-taxa/Psocoptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Rotifera <- subset(groups.fam.df, Group == "Rotifera")
Rotifera <- Rotifera[,-1]
spp <- vennCounts(Rotifera)
svg(filename="./venn_diagrams_inc-all-taxa/Rotifera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Tardigrada <- subset(groups.fam.df, Group == "Tardigrada")
Tardigrada <- Tardigrada[,-1]
spp <- vennCounts(Tardigrada)
svg(filename="./venn_diagrams_inc-all-taxa/Tardigrada_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Thysanoptera <- subset(groups.fam.df, Group == "Thysanoptera")
Thysanoptera <- Thysanoptera[,-1]
spp <- vennCounts(Thysanoptera)
svg(filename="./venn_diagrams_inc-all-taxa/Thysanoptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Trichoptera <- subset(groups.fam.df, Group == "Trichoptera")
Trichoptera <- Trichoptera[,-1]
spp <- vennCounts(Trichoptera)
svg(filename="./venn_diagrams_inc-all-taxa/Trichoptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()


#------------------------#
# SPECIES-LEVEL ANALYSIS #
#------------------------#

## ALPHA DIVERSITY
## Tranform data from counts (specimens or reads) to presence-absence by
## rows (ponds)
method.pa <- decostand(method.df, method = "pa", MARGIN = 1)
head(method.pa[,1:3], n = 3)

## Calculate basic species richness for each pond using each method
method.richness <- specnumber(method.pa)
method.richness

## Statistically compare and plot alpha diversity
## Create data frame
method.alpha <- data.frame(method.richness)

## Add sample metadata
method.alpha <- cbind(method.metadata[,1:3], method.alpha)
str(method.alpha)
method.alpha <- method.alpha %>% mutate_if(is.character, as.factor)

## Reset row names of data frame for further indexing
rownames(method.alpha) <- NULL

## Statistically compare whether method influences alpha diversity
method.regression <- glm.nb(method.richness ~ Method, data=method.alpha)
summary(method.regression)
anova(method.regression, test = "Chi")
drop1(method.regression, test = "Chi")

## Test difference between groups
summary(glht(glm.nb(method.richness ~ Method, data=method.alpha), 
             linfct=mcp(Method="Tukey")))

## Significant difference between means of all three methods

## Check model meets GLM assumptions
## Test for overdispersion
55.033/51
1-pchisq(55.033, df=51)  # not overdispersed

## Plot the fitted data against the observed data
plot(method.alpha$method.richness ~ fitted(method.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(method.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(method.regression)
shapiro.test(sresid) # P = 0.1063

## Residuals are normally distributed therefore model reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ method.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ method.alpha$Method, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(method.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(method.regression)
summary(influence)
CD <- cooks.distance(method.regression)
plot(CD ~ sresid)

## Plot species richness
m1 <- ggplot(method.alpha, aes(x=Method, y=method.richness)) + 
  geom_jitter(aes(colour=Method), 
              cex=2, pch=16, width=0.2, 
              show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  scale_y_continuous(limits=c(0,60), breaks=seq(0,60,10)) + 
  annotate("text", 
           x = c("DNA","eDNA","Netting"), y = 58, 
           label = c("a","a","b"), cex=10) + 
  labs(title="(a) Standard methods",
       subtitle="(i) Species-level",
       x="", y=expression(alpha~Diversity)) + 
  scale_colour_manual(values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = "white"),
        axis.line.x = element_line(colour = "black", size=0.5, linetype="solid"),
        axis.line.y = element_line(colour = "black", size=0.5, linetype="solid"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0),
        plot.subtitle = element_text(face="bold", hjust=0),
        legend.position = "none",
        text = element_text(size=20))
m1

## Species richness appears to be highest with eDNA metabarcoding, followed by
## DNA metabarcoding, then netting.


## BETA DIVERSITY
## Now look at beta diversity in ponds with different methods.
## The beta.pair() function computes three dissimilarity metrics and uses 
## the argument index.family to calculate beta diversity according to 
## Sorensen or Jaccard index of total dissimilarity. The function returns 
## three matrices containing the pairwise between-site values of each 
## component of beta diversity (turnover, nestedness, total beta). The 
## dissimilarity matrices yielded by beta.pair are objects of class dist 
## and can be used for further analyses, e.g. NMDS, cluster analysis.
## (NB: Jaccard index is used to calculate beta diversity from 
## presence-absence data. Abundance data would use Bray-curtis 
## dissimilarity)

## Pairwise between-site values of each component of beta diversity
method.dist <- beta.pair(method.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
method.bd.turn <- betadisper(method.dist$beta.jtu, method.metadata$Method)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(method.metadata, method.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(method.bd.turn$distances, method.metadata$Method, mean)

## Compute variance per group
tapply(method.bd.turn$distances, method.metadata$Method, var)

## Ordination plot of turnover partition
plot(method.bd.turn)

## Boxplot of turnover partition
boxplot(method.bd.turn, xlab="Method", xaxt="n", bty="n")
axis(side=1, at=c(1,2,3), labels=c("DNA","eDNA","Netting"))

## Plots indicate that there is some difference in turnover partition
## of beta diversity between methods. 
## Statistically check whether turnover is different between methods 
## using standard parametric anova or permutation tests.
anova(method.bd.turn)     # No significant difference between methods
permutest(method.bd.turn) # No significant difference between methods

## Analyse pairwise differences between groups (methods) using 
## parametric Tukey's HSD test.
TukeyHSD(method.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
method.comm.turn <- metaMDS(method.dist$beta.jtu, 
                            dist="jaccard", 
                            k=2,
                            maxit=999,
                            trymax=1000,
                            wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
method.comm.turn$stress
stressplot(method.comm.turn)

## plot site scores as text
ordiplot(method.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
method.NMDS1 <- method.comm.turn$points[,1]
method.NMDS2 <- method.comm.turn$points[,2]
method.turn.NMDS <- data.frame(NMDS1=method.NMDS1, 
                               NMDS2=method.NMDS2,
                               Method = method.metadata$Method)

## Check data
head(method.turn.NMDS)

## Plot data frame
m2a <- ggplot(method.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Method)) + 
  geom_point() + 
  stat_ellipse() + 
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.6,0.6)) + 
  labs(title="(a) Turnover", 
       subtitle="(i) Standard methods (species-level)",
       x="",y="NMDS2") + 
  scale_colour_manual(name="Method",
                      values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position = "bottom",
        legend.key=element_blank(),
        legend.key.size = unit(2, 'lines'))
m2a

## Statistically check difference in spatial turnover of communities
method.turn.anosim <- anosim(method.dist$beta.jtu, method.metadata$Method)

## Inspect results
method.turn.anosim
summary(method.turn.anosim)
plot(method.turn.anosim)

## There appears to be a significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
method.turn.adonis <- adonis(method.dist$beta.jtu ~ Method, method.metadata)

## Inspect results
## no summary() or plot() methods included
method.turn.adonis

## Again result is significant. There is a substantial difference in 
## species replacement (i.e. turnover) between monitoring tools.
## Therefore, species produced by one method are substituted by species 
## with another method.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
method.bd.nest <- betadisper(method.dist$beta.jne, method.metadata$Method)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
method.nest <- with(method.metadata, method.bd.nest)
method.nest

## Compute mean distance to centroid per group
tapply(method.bd.nest$distances, method.metadata$Method, mean)

## Compute variance per group
tapply(method.bd.nest$distances, method.metadata$Method, var)

## Ordination plot of turnover partition
plot(method.bd.nest)

## Boxplot of turnover partition
boxplot(method.bd.nest, xlab="Method", xaxt="n", bty="n")
axis(side=1, at=c(1,2,3), labels=c("DNA","eDNA","Netting"))

## Plots indicates that there is no difference in nestedness partition
## of beta diversity between methods 
## Statistically check whether nestedness is different between methods
## using standard parametric anova or permutation tests.
anova(method.bd.nest)     # No significant difference between methods
permutest(method.bd.nest) # No significant difference between methods

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(method.bd.nest)  # No significant difference between methods

## Ordination of beta diversity partitioned by nestedness:
method.comm.nest <- metaMDS(method.dist$beta.jne, 
                            dist="jaccard", 
                            k=2,
                            maxit=999,
                            trymax=1000,
                            wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
method.comm.nest$stress
stressplot(method.comm.nest)

## plot site scores as text
ordiplot(method.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
method.NMDS1 <- method.comm.nest$points[,1]
method.NMDS2 <- method.comm.nest$points[,2]
method.nest.NMDS <- data.frame(NMDS1=method.NMDS1, 
                               NMDS2=method.NMDS2,
                               Method = method.metadata$Method)

## Check data
head(method.nest.NMDS)

## Plot data frame
m2b <- ggplot(method.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Method)) + 
  geom_point() + 
  stat_ellipse() + 
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.6,0.6)) + 
  labs(title="(b) Nestedness", 
       subtitle="", 
       x="", y="") + 
  scale_colour_manual(name="Method",
                      values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
m2b

## Statistically check difference in nestedness of communities
method.nest.anosim <- anosim(method.dist$beta.jne, method.metadata$Method)

## Inspect results
method.nest.anosim
summary(method.nest.anosim)
plot(method.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
method.nest.adonis <- adonis(method.dist$beta.jne ~ Method, method.metadata)

## Inspect results
## no summary() or plot() methods included
method.nest.adonis

## Result is not significant. There is no substantial difference in species 
## loss or gain (i.e. nestedness) between methods.


## 3. TOTAL BETA DIVERSITY
method.bd.total <- betadisper(method.dist$beta.jac, method.metadata$Method)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(method.metadata, method.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(method.bd.total$distances, method.metadata$Method, mean)

## Compute variance per group
tapply(method.bd.total$distances, method.metadata$Method, var)

## Ordination plot of total beta diversity
plot(method.bd.total)

## Boxplot of total beta diversity
boxplot(method.bd.total, xlab="Method", xaxt="n", bty="n")
axis(side=1, at=c(1,2,3), labels=c("DNA","eDNA","Netting"))

## Plots indicates that there is no difference in total beta diversity 
## between methods.
## Statistically check whether total beta is different between methods
## using standard parametric anova or permutation tests.
anova(method.bd.total)     # No significant difference between methods
permutest(method.bd.total) # No significant difference between methods

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(method.bd.total)  # No significant difference between methods

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
method.comm.total <- metaMDS(method.dist$beta.jac, 
                             dist="jaccard", 
                             k=2,
                             maxit=999,
                             trymax=1000,
                             wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
method.comm.total$stress
stressplot(method.comm.total)

## plot site scores as text
ordiplot(method.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
method.NMDS1 <- method.comm.total$points[,1]
method.NMDS2 <- method.comm.total$points[,2]
method.total.NMDS <- data.frame(NMDS1=method.NMDS1,
                                NMDS2=method.NMDS2,
                                Method=method.metadata$Method)

## Check data
head(method.total.NMDS)

## Plot data frame
m2c <- ggplot(method.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Method)) +
  geom_point() + 
  stat_ellipse() + 
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.6,0.6)) + 
  labs(title=expression(bold("(c) Total"~beta~"Diversity")),
       subtitle="", 
       x="", y="") + 
  scale_colour_manual(name="Method",
                      values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
m2c

## Statistically check difference in total beta diversity of ponds
method.total.anosim <- anosim(method.dist$beta.jac, method.metadata$Method)

## Inspect results
method.total.anosim
summary(method.total.anosim)
plot(method.total.anosim)

## There appears to be a significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
method.total.adonis <- adonis(method.dist$beta.jac ~ Method, method.metadata)

## Inspect results
## no summary() or plot() methods included
method.total.adonis

## Again result is significant. There is substantial variation in overall 
## species composition of samples with different monitoring tools. Are 
## methods more or less complementary at family level?


#-----------------------#
# FAMILY-LEVEL ANALYSIS #
#-----------------------#

## ALPHA DIVERSITY
## Tranform data from counts (specimens or reads) to presence-absence by 
## rows (ponds)
method.fam.pa <- decostand(method.fam.df, method = "pa", MARGIN = 1)
head(method.fam.pa[,1:3], n = 3)

## Basic richness in each pond with each method
method.fam.richness <- specnumber(method.fam.pa)
method.fam.richness

## Statistically compare and plot alpha diversity
## Create data frame
method.fam.alpha <- data.frame(method.fam.richness)

## Add sample metadata
method.fam.alpha <- cbind(method.fam.metadata[,1:3], method.fam.alpha)
str(method.fam.alpha)
method.fam.alpha <- method.fam.alpha %>% mutate_if(is.character, as.factor)

## Reset row names of data frame for further indexing
rownames(method.fam.alpha) <- NULL

## Statistically compare whether method influences alpha diversity
method.fam.regression <- glm.nb(method.fam.richness ~ Method, data=method.fam.alpha)
summary(method.fam.regression)
anova(method.fam.regression, test = "Chi")
drop1(method.fam.regression, test = "Chi")

## Test difference between groups
summary(glht(glm.nb(method.fam.richness ~ Method, data=method.fam.alpha), 
             linfct=mcp(Method="Tukey")))

## Significant difference between means of all three methods

## Check model meets GLM assumptions
## Test for overdispersion
54.633/51
1-pchisq(54.633, df=51)  # not overdispersed

## Plot the fitted data against the observed data
plot(method.fam.alpha$method.fam.richness ~ fitted(method.fam.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(method.fam.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(method.fam.regression)
shapiro.test(sresid) # 0.9082

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ method.fam.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ method.fam.alpha$Method, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(method.fam.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(method.fam.regression)
summary(influence)
CD <- cooks.distance(method.fam.regression)
plot(CD ~ sresid)

## Plot family richness
m3 <- ggplot(method.fam.alpha, aes(x=Method, y=method.fam.richness)) + 
  geom_jitter(aes(colour=Method), cex=2, pch=16, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  scale_y_continuous(limits=c(0,60), breaks=seq(0,60,10)) + 
  annotate("text", 
           x = c("DNA","eDNA","Netting"), y = 58, 
           label = c("a","b","c"), cex=10) + 
  labs(subtitle="(ii) Family-level",
       x="Method", y=expression(alpha~Diversity)) + 
  scale_colour_manual(values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0),
        legend.position = "none",
        text = element_text(size=20))
m3

## Family richness appears to be highest with eDNA metabarcoding, followed by
## DNA metabarcoding, then netting.


## BETA DIVERSITY
## Now look at beta diversity in ponds with different methods
method.fam.dist <- beta.pair(method.fam.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
method.fam.bd.turn <- betadisper(method.fam.dist$beta.jtu, method.fam.metadata$Method)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(method.fam.metadata, method.fam.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(method.fam.bd.turn$distances, method.metadata$Method, mean)

## Compute variance per group
tapply(method.fam.bd.turn$distances, method.metadata$Method, var)

## Ordination plot of turnover partition
plot(method.fam.bd.turn)

## Boxplot of turnover partition
boxplot(method.fam.bd.turn, xlab="Method", xaxt="n", bty="n")
axis(side=1, at=c(1,2,3), labels=c("DNA","eDNA","Netting"))

## Plots indicate that there is some difference in turnover partition
## of beta diversity between methods. 
## Statistically check whether turnover is different between methods 
## using standard parametric anova or permutation tests.
anova(method.fam.bd.turn)     # No significant difference between methods
permutest(method.fam.bd.turn) # No significant difference between methods

## Analyse pairwise differences between groups (methods) using 
## parametric Tukey's HSD test.
TukeyHSD(method.fam.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
method.fam.comm.turn <- metaMDS(method.fam.dist$beta.jtu,
                                dist="jaccard", 
                                k=2,
                                maxit=999,
                                trymax=1000,
                                wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
method.fam.comm.turn$stress
stressplot(method.fam.comm.turn)

## plot site scores as text
ordiplot(method.fam.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
method.NMDS1 <- method.fam.comm.turn$points[,1]
method.NMDS2 <- method.fam.comm.turn$points[,2]
method.fam.turn.NMDS <- data.frame(NMDS1=method.NMDS1, 
                                   NMDS2=method.NMDS2,
                                   Method = method.fam.metadata$Method)

## Check data
head(method.fam.turn.NMDS)

## Plot data frame
m4a <- ggplot(method.fam.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Method)) + 
  geom_point() + 
  stat_ellipse() + 
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.6,0.6)) + 
  labs(subtitle="(ii) Standard methods (family-level)",
       x="",y="NMDS2") + 
  scale_colour_manual(name="Method",
                      values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position = "none",
        legend.key=element_blank(),
        legend.key.size = unit(2, 'lines'))
m4a

## Statistically check difference in spatial turnover of communities
method.fam.turn.anosim <- anosim(method.fam.dist$beta.jtu, method.fam.metadata$Method)

## Inspect results
method.fam.turn.anosim
summary(method.fam.turn.anosim)
plot(method.fam.turn.anosim)

## There appears to be a significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
method.fam.turn.adonis <- adonis(method.fam.dist$beta.jtu ~ Method, method.fam.metadata)

## Inspect results
## no summary() or plot() methods included
method.fam.turn.adonis

## Again result is significant. There is a substantial difference in 
## family replacement (i.e. turnover) between monitoring tools.
## Therefore, families produced by one method are substituted by families 
## with another method.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
method.fam.bd.nest <- betadisper(method.fam.dist$beta.jne, method.fam.metadata$Method)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
method.nest <- with(method.fam.metadata, method.fam.bd.nest)
method.nest

## Compute mean distance to centroid per group
tapply(method.fam.bd.nest$distances, method.metadata$Method, mean)

## Compute variance per group
tapply(method.fam.bd.nest$distances, method.metadata$Method, var)

## Ordination plot of turnover partition
plot(method.fam.bd.nest)

## Boxplot of turnover partition
boxplot(method.fam.bd.nest, xlab="Method", xaxt="n", bty="n")
axis(side=1, at=c(1,2,3), labels=c("DNA","eDNA","Netting"))

## Plots indicates that there is no difference in nestedness partition
## of beta diversity between methods 
## Statistically check whether nestedness is different between methods
## using standard parametric anova or permutation tests.
anova(method.fam.bd.nest)     # No significant difference between methods
permutest(method.fam.bd.nest) # No significant difference between methods

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(method.fam.bd.nest)  # No significant difference between methods

## Ordination of beta diversity partitioned by nestedness:
method.fam.comm.nest <- metaMDS(method.dist$beta.jne,
                                dist="jaccard", 
                                k=2,
                                maxit=999,
                                trymax=1000,
                                wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
method.fam.comm.nest$stress
stressplot(method.fam.comm.nest)

## plot site scores as text
ordiplot(method.fam.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
method.NMDS1 <- method.fam.comm.nest$points[,1]
method.NMDS2 <- method.fam.comm.nest$points[,2]
method.fam.nest.NMDS <- data.frame(NMDS1=method.NMDS1, 
                                   NMDS2=method.NMDS2,
                                   Method = method.fam.metadata$Method)

## Check data
head(method.fam.nest.NMDS)

## Plot data frame
m4b <- ggplot(method.fam.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Method)) + 
  geom_point() + 
  stat_ellipse() + 
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.6,0.6)) + 
  labs(subtitle="", x="", y="") + 
  scale_colour_manual(name="Method",
                      values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
m4b

## Statistically check difference in nestedness of communities
method.fam.nest.anosim <- anosim(method.fam.dist$beta.jne, method.fam.metadata$Method)

## Inspect results
method.fam.nest.anosim
summary(method.fam.nest.anosim)
plot(method.fam.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
method.fam.nest.adonis <- adonis(method.fam.dist$beta.jne ~ Method, method.fam.metadata)

## Inspect results
## no summary() or plot() methods included
method.fam.nest.adonis

## Result is not significant. There is no substantial difference in family 
## loss or gain (i.e. nestedness) between methods.


## 3. TOTAL BETA DIVERSITY
method.fam.bd.total <- betadisper(method.fam.dist$beta.jac, method.fam.metadata$Method)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(method.fam.metadata, method.fam.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(method.fam.bd.total$distances, method.metadata$Method, mean)

## Compute variance per group
tapply(method.fam.bd.total$distances, method.metadata$Method, var)

## Ordination plot of total beta diversity
plot(method.fam.bd.total)

## Boxplot of total beta diversity
boxplot(method.fam.bd.total, xlab="Method", xaxt="n", bty="n")
axis(side=1, at=c(1,2,3), labels=c("DNA","eDNA","Netting"))

## Plots indicates that there is no difference in total beta diversity 
## between methods.
## Statistically check whether total beta is different between methods
## using standard parametric anova or permutation tests.
anova(method.fam.bd.total)     # No significant difference between methods
permutest(method.fam.bd.total) # No significant difference between methods

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(method.fam.bd.total)  # Significant difference between methods

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
method.fam.comm.total <- metaMDS(method.fam.dist$beta.jac,
                                 dist="jaccard", 
                                 k=2,
                                 maxit=999,
                                 trymax=1000,
                                 wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
method.fam.comm.total$stress
stressplot(method.fam.comm.total)

## plot site scores as text
ordiplot(method.fam.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
method.NMDS1 <- method.fam.comm.total$points[,1]
method.NMDS2 <- method.fam.comm.total$points[,2]
method.fam.total.NMDS <- data.frame(NMDS1=method.NMDS1,
                                    NMDS2=method.NMDS2,
                                    Method=method.fam.metadata$Method)

## Check data
head(method.fam.total.NMDS)

## Plot data frame
m4c <- ggplot(method.fam.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Method)) + 
  geom_point() + 
  stat_ellipse() + 
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.6,0.6)) + 
  labs(subtitle="", x="", y="") + 
  scale_colour_manual(name="Method",
                      values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
m4c

## Statistically check difference in total beta diversity of ponds
method.fam.total.anosim <- anosim(method.fam.dist$beta.jac, method.fam.metadata$Method)

## Inspect results
method.fam.total.anosim
summary(method.fam.total.anosim)
plot(method.fam.total.anosim)

## There appears to be a significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
method.fam.total.adonis <- adonis(method.fam.dist$beta.jac ~ Method, method.fam.metadata)

## Inspect results
## no summary() or plot() methods included
method.fam.total.adonis

## Again result is significant. There is substantial variation in overall 
## family composition of samples with different monitoring tools. Netting
## and DNA metabarcoding are more similar than eDNA metabarcoding at either
## taxonomic rank.



####################################################################
# METHOD COMPARISON (MEIOFAUNA/TAXA MISSING REF SEQUENCES REMOVED) #
####################################################################

#--------------------------#
# VISUALISE METHOD OVERLAP #
#--------------------------#

## Tranform data from read counts to presence-absence by
## rows (ponds)
net.noref.spp.pa <- decostand(net.noref, method = "pa", MARGIN = 1)
head(net.noref.spp.pa[,1:3], n = 3)

net.noref.fam.pa <- decostand(net.fam.noref, method = "pa", MARGIN = 1)
head(net.noref.fam.pa[,1:3], n = 3)

DNA.nomicro.spp.pa <- decostand(DNA.nomicro, method = "pa", MARGIN = 1)
head(DNA.nomicro.spp.pa[,1:3], n = 3)

DNA.nomicro.fam.pa <- decostand(DNA.fam.nomicro, method = "pa", MARGIN = 1)
head(DNA.nomicro.fam.pa[,1:3], n = 3)

eDNA.nomicro.spp.pa <- decostand(eDNA.nomicro, method = "pa", MARGIN = 1)
head(eDNA.nomicro.spp.pa[,1:3], n = 3)

eDNA.nomicro.fam.pa <- decostand(eDNA.fam.nomicro, method = "pa", MARGIN = 1)
head(eDNA.nomicro.fam.pa[,1:3], n = 3)


## Produce venn diagrams of overlap in taxonomic assignments each method
## Species-level
## Create dataframe for each method containing species name and total number
## of ponds it was found in
net.noref.total <- data.frame(colnames(net.noref.spp.pa), colSums(net.noref.spp.pa))
colnames(net.noref.total) <- c("Species", "Netting")

DNA.nomicro.total <- data.frame(colnames(DNA.nomicro.spp.pa), colSums(DNA.nomicro.spp.pa))
colnames(DNA.nomicro.total) <- c("Species", "DNA metabarcoding")

eDNA.nomicro.total <- data.frame(colnames(eDNA.nomicro.spp.pa), colSums(eDNA.nomicro.spp.pa))
colnames(eDNA.nomicro.total) <- c("Species", "eDNA metabarcoding")

## Make list of dataframes and bind them together
## Empty rows will be filled with NA values using rbind.fill() function
nomicro.spp.totals <- list(net.noref.total, DNA.nomicro.total, eDNA.nomicro.total)
nomicro.spp.df <- do.call(rbind.fill, nomicro.spp.totals)

## Replace NA values with 0 and merge rows by species
nomicro.spp.df[is.na(nomicro.spp.df)] <- 0
nomicro.spp.df <- ddply(nomicro.spp.df, .(Species), numcolwise(sum))

## Write as .csv file and add column containing taxonomic order for species
## assignments
#write.csv(nomicro.spp.df, "../Data/venn_species_nomicro_20200818.csv", 
#          row.names=FALSE)

## Reimport file
nomicro.spp.df <- read.csv("../Data/venn_species_nomicro_edited_20200818.csv", 
                           header=TRUE)

## Venn diagram of total species detections across methods
## Take species detections by method and make new presence-absence dataframe
venn.df <- nomicro.spp.df[,3:5]
venn.df[venn.df > 0] <- 1

## Make venn diagram
spp <- vennCounts(venn.df)
svg(filename="./venn_diagrams_exc-some-taxa/total_species_nomicro_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

## Venn diagram of species detections by taxonomic group across methods
## Take species detections by method and make new presence-absence dataframe
nomicro.groups.df <- nomicro.spp.df[,-1]
nomicro.groups.df[nomicro.groups.df > 0] <- 1

## Subset new dataframe by major invertebrate groups and make individual
## venn diagrams
Annelida <- subset(nomicro.groups.df, Group == "Annelida")
Annelida <- Annelida[,-1]
spp <- vennCounts(Annelida)
svg(filename="./venn_diagrams_exc-some-taxa/Annelida_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Arachnida <- subset(nomicro.groups.df, Group == "Arachnida")
Arachnida <- Arachnida[,-1]
spp <- vennCounts(Arachnida)
svg(filename="./venn_diagrams_exc-some-taxa/Arachnida_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Coleoptera <- subset(nomicro.groups.df, Group == "Coleoptera")
Coleoptera <- Coleoptera[,-1]
spp <- vennCounts(Coleoptera)
svg(filename="./venn_diagrams_exc-some-taxa/Coleoptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Collembola <- subset(nomicro.groups.df, Group == "Collembola")
Collembola <- Collembola[,-1]
spp <- vennCounts(Collembola)
svg(filename="./venn_diagrams_exc-some-taxa/Collembola_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Crustacea <- subset(nomicro.groups.df, Group == "Crustacea")
Crustacea <- Crustacea[,-1]
spp <- vennCounts(Crustacea)
svg(filename="./venn_diagrams_exc-some-taxa/Crustacea_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Diptera <- subset(nomicro.groups.df, Group == "Diptera")
Diptera <- Diptera[,-1]
spp <- vennCounts(Diptera)
svg(filename="./venn_diagrams_exc-some-taxa/Diptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Ephemeroptera <- subset(nomicro.groups.df, Group == "Ephemeroptera")
Ephemeroptera <- Ephemeroptera[,-1]
spp <- vennCounts(Ephemeroptera)
svg(filename="./venn_diagrams_exc-some-taxa/Ephemeroptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Hemiptera <- subset(nomicro.groups.df, Group == "Hemiptera")
Hemiptera <- Hemiptera[,-1]
spp <- vennCounts(Hemiptera)
svg(filename="./venn_diagrams_exc-some-taxa/Hemiptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Hirudinea <- subset(nomicro.groups.df, Group == "Hirudinea")
Hirudinea <- Hirudinea[,-1]
spp <- vennCounts(Hirudinea)
svg(filename="./venn_diagrams_exc-some-taxa/Hirudinea_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Hymenoptera <- subset(nomicro.groups.df, Group == "Hymenoptera")
Hymenoptera <- Hymenoptera[,-1]
spp <- vennCounts(Hymenoptera)
svg(filename="./venn_diagrams_exc-some-taxa/Hymenoptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Lepidoptera <- subset(nomicro.groups.df, Group == "Lepidoptera")
Lepidoptera <- Lepidoptera[,-1]
spp <- vennCounts(Lepidoptera)
svg(filename="./venn_diagrams_exc-some-taxa/Lepidoptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Megaloptera <- subset(nomicro.groups.df, Group == "Megaloptera")
Megaloptera <- Megaloptera[,-1]
spp <- vennCounts(Megaloptera)
svg(filename="./venn_diagrams_exc-some-taxa/Megaloptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Mollusca <- subset(nomicro.groups.df, Group == "Mollusca")
Mollusca <- Mollusca[,-1]
spp <- vennCounts(Mollusca)
svg(filename="./venn_diagrams_exc-some-taxa/Mollusca_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Odonata <- subset(nomicro.groups.df, Group == "Odonata")
Odonata <- Odonata[,-1]
spp <- vennCounts(Odonata)
svg(filename="./venn_diagrams_exc-some-taxa/Odonata_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Platyhelminthes <- subset(nomicro.groups.df, Group == "Platyhelminthes")
Platyhelminthes <- Platyhelminthes[,-1]
spp <- vennCounts(Psocoptera)
svg(filename="./venn_diagrams_exc-some-taxa/Plathelminthes_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Psocoptera <- subset(nomicro.groups.df, Group == "Psocoptera")
Psocoptera <- Psocoptera[,-1]
spp <- vennCounts(Psocoptera)
svg(filename="./venn_diagrams_exc-some-taxa/Psocoptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Rotifera <- subset(nomicro.groups.df, Group == "Rotifera")
Rotifera <- Rotifera[,-1]
spp <- vennCounts(Rotifera)
svg(filename="./venn_diagrams_exc-some-taxa/Rotifera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Trichoptera <- subset(nomicro.groups.df, Group == "Trichoptera")
Trichoptera <- Trichoptera[,-1]
spp <- vennCounts(Trichoptera)
svg(filename="./venn_diagrams_exc-some-taxa/Trichoptera_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()


## Family-level
## Create dataframe for each method containing family name and total number
## of ponds it was found in
net.noref.fam.total <- data.frame(colnames(net.noref.fam.pa), colSums(net.noref.fam.pa))
colnames(net.noref.fam.total) <- c("Family", "Netting")

DNA.nomicro.fam.total <- data.frame(colnames(DNA.nomicro.fam.pa), colSums(DNA.nomicro.fam.pa))
colnames(DNA.nomicro.fam.total) <- c("Family", "DNA metabarcoding")

eDNA.nomicro.fam.total <- data.frame(colnames(eDNA.nomicro.fam.pa), colSums(eDNA.nomicro.fam.pa))
colnames(eDNA.nomicro.fam.total) <- c("Family", "eDNA metabarcoding")

## Make list of dataframes and bind them together
## Empty rows will be filled with NA values using rbind.fill() function
nomicro.fam.totals <- list(net.noref.fam.total, DNA.nomicro.fam.total, eDNA.nomicro.fam.total)
nomicro.fam.df <- do.call(rbind.fill, nomicro.fam.totals)

## Replace NA values with 0 and merge rows by species
nomicro.fam.df[is.na(nomicro.fam.df)] <- 0
nomicro.fam.df <- ddply(nomicro.fam.df, .(Family), numcolwise(sum))

## Write as .csv file and add column containing taxonomic order for species
## assignments
#write.csv(nomicro.fam.df, "../Data/venn_family_nomicro_20200818.csv", 
#          row.names=FALSE)

## Reimport file
nomicro.fam.df <- read.csv("../Data/venn_family_nomicro_edited_20200818.csv", 
                           header=TRUE)

## Venn diagram of total family detections across methods
## Take species detections by method and make new presence-absence dataframe
venn.fam.df <- nomicro.fam.df[,3:5]
venn.fam.df[venn.fam.df > 0] <- 1

## Make venn diagram
fam <- vennCounts(venn.fam.df)
svg(filename="./venn_diagrams_exc-some-taxa/total_family_nomicro_method.svg", width = 13, height = 10)
vennDiagram(fam, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

## Venn diagram of species detections by taxonomic group across methods
## Take species detections by method and make new presence-absence dataframe
nomicro.groups.fam.df <- nomicro.fam.df[,-1]
nomicro.groups.fam.df[nomicro.groups.fam.df > 0] <- 1

## Subset new dataframe by major invertebrate groups and make individual
## venn diagrams
Annelida <- subset(nomicro.groups.fam.df, Group == "Annelida")
Annelida <- Annelida[,-1]
spp <- vennCounts(Annelida)
svg(filename="./venn_diagrams_exc-some-taxa/Annelida_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Arachnida <- subset(nomicro.groups.fam.df, Group == "Arachnida")
Arachnida <- Arachnida[,-1]
spp <- vennCounts(Arachnida)
svg(filename="./venn_diagrams_exc-some-taxa/Arachnida_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Coleoptera <- subset(nomicro.groups.fam.df, Group == "Coleoptera")
Coleoptera <- Coleoptera[,-1]
spp <- vennCounts(Coleoptera)
svg(filename="./venn_diagrams_exc-some-taxa/Coleoptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Collembola <- subset(nomicro.groups.fam.df, Group == "Collembola")
Collembola <- Collembola[,-1]
spp <- vennCounts(Collembola)
svg(filename="./venn_diagrams_exc-some-taxa/Collembola_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Crustacea <- subset(nomicro.groups.fam.df, Group == "Crustacea")
Crustacea <- Crustacea[,-1]
spp <- vennCounts(Crustacea)
svg(filename="./venn_diagrams_exc-some-taxa/Crustacea_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Diptera <- subset(nomicro.groups.fam.df, Group == "Diptera")
Diptera <- Diptera[,-1]
spp <- vennCounts(Diptera)
svg(filename="./venn_diagrams_exc-some-taxa/Diptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Ephemeroptera <- subset(nomicro.groups.fam.df, Group == "Ephemeroptera")
Ephemeroptera <- Ephemeroptera[,-1]
spp <- vennCounts(Ephemeroptera)
svg(filename="./venn_diagrams_exc-some-taxa/Ephemeroptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Hemiptera <- subset(nomicro.groups.fam.df, Group == "Hemiptera")
Hemiptera <- Hemiptera[,-1]
spp <- vennCounts(Hemiptera)
svg(filename="./venn_diagrams_exc-some-taxa/Hemiptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Hirudinea <- subset(nomicro.groups.fam.df, Group == "Hirudinea")
Hirudinea <- Hirudinea[,-1]
spp <- vennCounts(Hirudinea)
svg(filename="./venn_diagrams_exc-some-taxa/Hirudinea_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Hymenoptera <- subset(nomicro.groups.fam.df, Group == "Hymenoptera")
Hymenoptera <- Hymenoptera[,-1]
spp <- vennCounts(Hymenoptera)
svg(filename="./venn_diagrams_exc-some-taxa/Hymenoptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Lepidoptera <- subset(nomicro.groups.fam.df, Group == "Lepidoptera")
Lepidoptera <- Lepidoptera[,-1]
spp <- vennCounts(Lepidoptera)
svg(filename="./venn_diagrams_exc-some-taxa/Lepidoptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Megaloptera <- subset(nomicro.groups.fam.df, Group == "Megaloptera")
Megaloptera <- Megaloptera[,-1]
spp <- vennCounts(Megaloptera)
svg(filename="./venn_diagrams_exc-some-taxa/Megaloptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Mollusca <- subset(nomicro.groups.fam.df, Group == "Mollusca")
Mollusca <- Mollusca[,-1]
spp <- vennCounts(Mollusca)
svg(filename="./venn_diagrams_exc-some-taxa/Mollusca_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Odonata <- subset(nomicro.groups.fam.df, Group == "Odonata")
Odonata <- Odonata[,-1]
spp <- vennCounts(Odonata)
svg(filename="./venn_diagrams_exc-some-taxa/Odonata_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Platyhelminthes <- subset(nomicro.groups.fam.df, Group == "Platyhelminthes")
Platyhelminthes <- Platyhelminthes[,-1]
spp <- vennCounts(Platyhelminthes)
svg(filename="./venn_diagrams_exc-some-taxa/Platyhelminthes_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Psocoptera <- subset(nomicro.groups.fam.df, Group == "Psocoptera")
Psocoptera <- Psocoptera[,-1]
spp <- vennCounts(Psocoptera)
svg(filename="./venn_diagrams_exc-some-taxa/Psocoptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Rotifera <- subset(nomicro.groups.fam.df, Group == "Rotifera")
Rotifera <- Rotifera[,-1]
spp <- vennCounts(Rotifera)
svg(filename="./venn_diagrams_exc-some-taxa/Rotifera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Thysanoptera <- subset(nomicro.groups.fam.df, Group == "Thysanoptera")
Thysanoptera <- Thysanoptera[,-1]
spp <- vennCounts(Thysanoptera)
svg(filename="./venn_diagrams_exc-some-taxa/Thysanoptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()

Trichoptera <- subset(nomicro.groups.fam.df, Group == "Trichoptera")
Trichoptera <- Trichoptera[,-1]
spp <- vennCounts(Trichoptera)
svg(filename="./venn_diagrams_exc-some-taxa/Trichoptera_family_method.svg", width = 13, height = 10)
vennDiagram(spp, names = c("Netting", "DNA metabarcoding", "eDNA metabarcoding"), 
            cex = 2, lwd = 2, counts.col = "black",
            circle.col = c("limegreen","purple","orange"))
dev.off()


#------------------------#
# SPECIES-LEVEL ANALYSIS #
#------------------------#

## ALPHA DIVERSITY
## Tranform data from counts (specimens or reads) to presence-absence by
## rows (ponds)
nomicro.pa <- decostand(method.nomicro, method = "pa", MARGIN = 1)
head(nomicro.pa[,1:3], n = 3)

## Basic richness
nomicro.richness <- specnumber(nomicro.pa)
nomicro.richness

## Statistically compare and plot alpha diversity
## Create data frame
nomicro.alpha <- data.frame(nomicro.richness)

## Add sample metadata
nomicro.alpha <- cbind(method.metadata[,1:3], nomicro.alpha)
str(nomicro.alpha)
nomicro.alpha <- nomicro.alpha %>% mutate_if(is.character, as.factor)

## Reset row names of data frame for further indexing
rownames(nomicro.alpha) <- NULL

## Statistically compare whether method influences alpha diversity
nomicro.regression <- glm.nb(nomicro.richness ~ Method, data=nomicro.alpha)
summary(nomicro.regression)
anova(nomicro.regression, test = "Chi")
drop1(nomicro.regression, test = "Chi")

## Test difference between groups
summary(glht(glm.nb(nomicro.richness ~ Method, data=nomicro.alpha), 
             linfct=mcp(Method="Tukey")))

## Significant difference between means of:
## eDNA-DNA
## Netting-DNA

## Check model meets GLM assumptions
## Test for overdispersion
57.708/51
1-pchisq(57.708, df=51)  # not overdispersed

## Plot the fitted data against the observed data
plot(nomicro.alpha$nomicro.richness ~ fitted(nomicro.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(nomicro.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(nomicro.regression)
shapiro.test(sresid) # P = 0.9735

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ method.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ method.alpha$Method, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(method.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(method.regression)
summary(influence)
CD <- cooks.distance(method.regression)
plot(CD ~ sresid)

## Plot species richness
m5 <- ggplot(nomicro.alpha, aes(x=Method, y=nomicro.richness)) + 
  geom_jitter(aes(colour=Method), cex=2, pch=16, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  scale_y_continuous(limits=c(0,60), breaks=seq(0,60,10)) + 
  annotate("text", 
           x = c("DNA","eDNA","Netting"), y = 58, 
           label = c("a","b","b"), cex=10) + 
  labs(title="(b) Unbiased methods",
       subtitle="",
       x="", y="") + 
  scale_colour_manual(values=c("purple","orange","limegreen")) +
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0),
        plot.subtitle = element_text(face="bold", hjust=0),
        legend.position = "none",
        text = element_text(size=20))
m5

## Species richness appears to be highest with DNA metabarcoding, followed by
## netting and eDNA metabarcoding.


## BETA DIVERSITY
## Pairwise between-site values of each component of beta diversity
nomicro.dist <- beta.pair(nomicro.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
nomicro.bd.turn <- betadisper(nomicro.dist$beta.jtu, method.nomicro.metadata$Method)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(method.metadata, nomicro.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(nomicro.bd.turn$distances, method.metadata$Method, mean)

## Compute variance per group
tapply(nomicro.bd.turn$distances, method.metadata$Method, var)

## Ordination plot of turnover partition
plot(nomicro.bd.turn)

## Boxplot of turnover partition
boxplot(nomicro.bd.turn, xlab="Method", xaxt="n", bty="n")
axis(side=1, at=c(1,2,3), labels=c("DNA","eDNA","Netting"))

## Plots indicate that there is some difference in turnover partition
## of beta diversity between methods. 
## Statistically check whether turnover is different between methods 
## using standard parametric anova or permutation tests.
anova(nomicro.bd.turn)     # No significant difference between methods
permutest(nomicro.bd.turn) # No significant difference between methods

## Analyse pairwise differences between groups (methods) using 
## parametric Tukey's HSD test.
TukeyHSD(nomicro.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
nomicro.comm.turn <- metaMDS(nomicro.dist$beta.jtu, 
                             dist="jaccard", 
                             k=2,
                             maxit=999,
                             trymax=1000,
                             wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
nomicro.comm.turn$stress
stressplot(nomicro.comm.turn)

## plot site scores as text
ordiplot(nomicro.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
nomicro.NMDS1 <- nomicro.comm.turn$points[,1]
nomicro.NMDS2 <- nomicro.comm.turn$points[,2]
nomicro.turn.NMDS <- data.frame(NMDS1=nomicro.NMDS1, 
                                NMDS2=nomicro.NMDS2,
                                Method = method.nomicro.metadata$Method)

## Check data
head(nomicro.turn.NMDS)

## Plot data frame
m6a <- ggplot(nomicro.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Method)) + 
  geom_point() + 
  stat_ellipse() + 
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.6,0.6)) + 
  labs(subtitle="(iii) Unbiased methods (species-level)",
       x="",y="NMDS2") + 
  scale_colour_manual(name="Method",
                      values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black"),
        text = element_text(size=20),
        legend.position = "bottom",
        legend.key=element_blank(),
        legend.key.size = unit(2, 'lines'))
m6a

## Statistically check difference in spatial turnover of communities
nomicro.turn.anosim <- anosim(nomicro.dist$beta.jtu, method.nomicro.metadata$Method)

## Inspect results
nomicro.turn.anosim
summary(nomicro.turn.anosim)
plot(nomicro.turn.anosim)

## There appears to be a significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
nomicro.turn.adonis <- adonis(nomicro.dist$beta.jtu ~ Method, method.nomicro.metadata)

## Inspect results
## no summary() or plot() methods included
nomicro.turn.adonis

## Again result is significant. There is a substantial difference in 
## species replacement (i.e. turnover) between monitoring tools.
## Therefore, species produced by one method are substituted by species 
## with another method.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
nomicro.bd.nest <- betadisper(nomicro.dist$beta.jne, method.nomicro.metadata$Method)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
nomicro.nest <- with(method.nomicro.metadata, nomicro.bd.nest)
nomicro.nest

## Compute mean distance to centroid per group
tapply(nomicro.bd.nest$distances, method.nomicro.metadata$Method, mean)

## Compute variance per group
tapply(nomicro.bd.nest$distances, method.nomicro.metadata$Method, var)

## Ordination plot of turnover partition
plot(nomicro.bd.nest)

## Boxplot of turnover partition
boxplot(nomicro.bd.nest, xlab="Method", xaxt="n", bty="n")
axis(side=1, at=c(1,2,3), labels=c("DNA","eDNA","Netting"))

## Plots indicates that there is no difference in nestedness partition
## of beta diversity between methods 
## Statistically check whether nestedness is different between methods
## using standard parametric anova or permutation tests.
anova(nomicro.bd.nest)     # No significant difference between methods
permutest(nomicro.bd.nest) # No significant difference between methods

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(nomicro.bd.nest)  # No significant difference between methods

## Ordination of beta diversity partitioned by nestedness:
nomicro.comm.nest <- metaMDS(nomicro.dist$beta.jne, 
                             dist="jaccard", 
                             k=2,
                             maxit=999,
                             trymax=1000,
                             wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
nomicro.comm.nest$stress
stressplot(nomicro.comm.nest)

## plot site scores as text
ordiplot(nomicro.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
nomicro.NMDS1 <- nomicro.comm.nest$points[,1]
nomicro.NMDS2 <- nomicro.comm.nest$points[,2]
nomicro.nest.NMDS <- data.frame(NMDS1=nomicro.NMDS1, 
                                NMDS2=nomicro.NMDS2,
                                Method = method.nomicro.metadata$Method)

## Check data
head(nomicro.nest.NMDS)

## Plot data frame
m6b <- ggplot(nomicro.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Method)) + 
  geom_point() + 
  stat_ellipse() + 
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.6,0.6)) + 
  labs(subtitle="", x="", y="") + 
  scale_colour_manual(name="Method",
                      values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black"),
        text = element_text(size=20),
        legend.position="none")
m6b

## Statistically check difference in nestedness of communities
nomicro.nest.anosim <- anosim(nomicro.dist$beta.jne, method.nomicro.metadata$Method)

## Inspect results
nomicro.nest.anosim
summary(nomicro.nest.anosim)
plot(nomicro.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
nomicro.nest.adonis <- adonis(nomicro.dist$beta.jne ~ Method, method.nomicro.metadata)

## Inspect results
## no summary() or plot() methods included
nomicro.nest.adonis

## Result is not significant. There is no substantial difference in species 
## loss or gain (i.e. nestedness) between methods.


## 3. TOTAL BETA DIVERSITY
nomicro.bd.total <- betadisper(nomicro.dist$beta.jac, method.nomicro.metadata$Method)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(method.nomicro.metadata, nomicro.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(nomicro.bd.total$distances, method.nomicro.metadata$Method, mean)

## Compute variance per group
tapply(nomicro.bd.total$distances, method.nomicro.metadata$Method, var)

## Ordination plot of total beta diversity
plot(nomicro.bd.total)

## Boxplot of total beta diversity
boxplot(nomicro.bd.total, xlab="Method", xaxt="n", bty="n")
axis(side=1, at=c(1,2,3), labels=c("DNA","eDNA","Netting"))

## Plots indicates that there is no difference in total beta diversity 
## between methods.
## Statistically check whether total beta is different between methods
## using standard parametric anova or permutation tests.
anova(nomicro.bd.total)     # No significant difference between methods
permutest(nomicro.bd.total) # No significant difference between methods

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(nomicro.bd.total)  # No significant difference between methods

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
nomicro.comm.total <- metaMDS(nomicro.dist$beta.jac, 
                              dist="jaccard", 
                              k=2,
                              maxit=999,
                              trymax=1000,
                              wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
nomicro.comm.total$stress
stressplot(nomicro.comm.total)

## plot site scores as text
ordiplot(nomicro.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
nomicro.NMDS1 <- nomicro.comm.total$points[,1]
nomicro.NMDS2 <- nomicro.comm.total$points[,2]
nomicro.total.NMDS <- data.frame(NMDS1=nomicro.NMDS1,
                                 NMDS2=nomicro.NMDS2,
                                 Method=method.nomicro.metadata$Method)

## Check data
head(nomicro.total.NMDS)

## Plot data frame
m6c <- ggplot(nomicro.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Method)) + 
  geom_point() + 
  stat_ellipse() + 
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.6,0.6)) + 
  labs(subtitle="", x="", y="") + 
  scale_colour_manual(name="Method",
                      values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black"),
        text = element_text(size=20),
        legend.position="none")
m6c

## Statistically check difference in total beta diversity of ponds
nomicro.total.anosim <- anosim(nomicro.dist$beta.jac, method.nomicro.metadata$Method)

## Inspect results
nomicro.total.anosim
summary(nomicro.total.anosim)
plot(nomicro.total.anosim)

## There appears to be a significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
nomicro.total.adonis <- adonis(nomicro.dist$beta.jac ~ Method, method.nomicro.metadata)

## Inspect results
## no summary() or plot() methods included
nomicro.total.adonis

## Again result is significant. There is substantial variation in overall 
## species composition of samples with different monitoring tools. Are methods
## more or less complementary at family level?


#-----------------------#
# FAMILY-LEVEL ANALYSIS #
#-----------------------#

## ALPHA DIVERSITY
## Tranform data from counts (specimens or reads) to presence-absence by 
## rows (ponds)
nomicro.fam.pa <- decostand(method.fam.nomicro, method = "pa", MARGIN = 1)
head(nomicro.fam.pa[,1:3], n = 3)

## Basic richness
nomicro.fam.richness <- specnumber(nomicro.fam.pa)
nomicro.fam.richness

## Statistically compare and plot alpha diversity
## Create data frame
nomicro.fam.alpha <- data.frame(nomicro.fam.richness)

## Add sample metadata
nomicro.fam.alpha <- cbind(method.fam.nomicro.metadata[,1:3], nomicro.fam.alpha)
str(nomicro.fam.alpha)
nomicro.fam.alpha <- nomicro.fam.alpha %>% mutate_if(is.character, as.factor)

## Reset row names of data frame for further indexing
rownames(nomicro.fam.alpha) <- NULL

## Statistically compare whether method influences alpha diversity
nomicro.fam.regression <- glm.nb(nomicro.fam.richness ~ Method, data=nomicro.fam.alpha)
summary(nomicro.fam.regression)
anova(nomicro.fam.regression, test = "Chi")
drop1(nomicro.fam.regression, test = "Chi")

## Test difference between groups
summary(glht(glm.nb(nomicro.fam.richness ~ Method, data=nomicro.fam.alpha), 
             linfct=mcp(Method="Tukey")))

## No significant difference between means of any method.

## Check model meets GLM assumptions
## Test for overdispersion
57.840/51
1-pchisq(57.840, df=51)  # not overdispersed

## Plot the fitted data against the observed data
plot(nomicro.fam.alpha$nomicro.fam.richness ~ fitted(nomicro.fam.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(nomicro.fam.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(nomicro.fam.regression)
shapiro.test(sresid) # 0.745

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ nomicro.fam.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ nomicro.fam.alpha$Method, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(nomicro.fam.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(nomicro.fam.regression)
summary(influence)
CD <- cooks.distance(nomicro.fam.regression)
plot(CD ~ sresid)

## Plot species richness
m7 <- ggplot(nomicro.fam.alpha, aes(x=Method, y=nomicro.fam.richness)) + 
  geom_jitter(aes(colour=Method), cex=2, pch=16, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  scale_y_continuous(limits=c(0,60), breaks=seq(0,60,10)) + 
  annotate("text", 
           x = c("DNA","eDNA","Netting"), y = 58, 
           label = c("a","a","a"), cex=10) + 
  labs(subtitle="",
       x="Method", y="") + 
  scale_colour_manual(values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0),
        legend.position = "none",
        text = element_text(size=20))
m7

## Family richness appears to be comparable across methods.


## BETA DIVERSITY
nomicro.fam.dist <- beta.pair(nomicro.fam.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
nomicro.fam.bd.turn <- betadisper(nomicro.fam.dist$beta.jtu, method.fam.nomicro.metadata$Method)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(method.fam.nomicro.metadata, nomicro.fam.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(nomicro.fam.bd.turn$distances, method.fam.nomicro.metadata$Method, mean)

## Compute variance per group
tapply(nomicro.fam.bd.turn$distances, method.fam.nomicro.metadata$Method, var)

## Ordination plot of turnover partition
plot(nomicro.fam.bd.turn)

## Boxplot of turnover partition
boxplot(nomicro.fam.bd.turn, xlab="Method", xaxt="n", bty="n")
axis(side=1, at=c(1,2,3), labels=c("DNA","eDNA","Netting"))

## Plots indicate that there is some difference in turnover partition
## of beta diversity between methods. 
## Statistically check whether turnover is different between methods 
## using standard parametric anova or permutation tests.
anova(nomicro.fam.bd.turn)     # Significant difference between methods
permutest(nomicro.fam.bd.turn) # Significant difference between methods

## Analyse pairwise differences between groups (methods) using 
## parametric Tukey's HSD test.
TukeyHSD(nomicro.fam.bd.turn)  

## No significant difference between methods

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
nomicro.fam.comm.turn <- metaMDS(nomicro.fam.dist$beta.jtu, 
                                 dist="jaccard", 
                                 k=2,
                                 maxit=999,
                                 trymax=1000,
                                 wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
nomicro.fam.comm.turn$stress
stressplot(nomicro.fam.comm.turn)

## plot site scores as text
ordiplot(nomicro.fam.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
nomicro.NMDS1 <- nomicro.fam.comm.turn$points[,1]
nomicro.NMDS2 <- nomicro.fam.comm.turn$points[,2]
nomicro.fam.turn.NMDS <- data.frame(NMDS1=nomicro.NMDS1, 
                                    NMDS2=nomicro.NMDS2,
                                    Method = method.fam.nomicro.metadata$Method)

## Check data
head(nomicro.fam.turn.NMDS)

## Plot data frame
m8a <- ggplot(nomicro.fam.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Method)) + 
  geom_point() + 
  stat_ellipse() + 
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.6,0.6)) + 
  labs(subtitle="(iv) Unbiased methods (family-level)", 
       x="NMDS1",y="NMDS2") + 
  scale_colour_manual(name="Method",
                      values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position = "none",
        legend.key=element_blank(),
        legend.key.size = unit(2, 'lines'))
m8a

## Statistically check difference in spatial turnover of communities
method.fam.turn.anosim <- anosim(method.fam.dist$beta.jtu, method.fam.nomicro.metadata$Method)

## Inspect results
method.fam.turn.anosim
summary(method.fam.turn.anosim)
plot(method.fam.turn.anosim)

## There appears to be a significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
method.fam.turn.adonis <- adonis(method.fam.dist$beta.jtu ~ Method, method.fam.nomicro.metadata)

## Inspect results
## no summary() or plot() methods included
method.fam.turn.adonis

## Again result is significant. There is a substantial difference in 
## family replacement (i.e. turnover) between monitoring tools.
## Therefore, families produced by one method are substituted by families
## with another method.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
nomicro.fam.bd.nest <- betadisper(nomicro.fam.dist$beta.jne, method.fam.nomicro.metadata$Method)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
nomicro.nest <- with(method.fam.nomicro.metadata, nomicro.fam.bd.nest)
nomicro.nest

## Compute mean distance to centroid per group
tapply(nomicro.fam.bd.nest$distances, method.fam.nomicro.metadata$Method, mean)

## Compute variance per group
tapply(nomicro.fam.bd.nest$distances, method.fam.nomicro.metadata$Method, var)

## Ordination plot of turnover partition
plot(nomicro.fam.bd.nest)

## Boxplot of turnover partition
boxplot(nomicro.fam.bd.nest, xlab="Method", xaxt="n", bty="n")
axis(side=1, at=c(1,2,3), labels=c("DNA","eDNA","Netting"))

## Plots indicates that there is no difference in nestedness partition
## of beta diversity between methods 
## Statistically check whether nestedness is different between methods
## using standard parametric anova or permutation tests.
anova(nomicro.fam.bd.nest)     # No significant difference between methods
permutest(nomicro.fam.bd.nest) # No significant difference between methods

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(nomicro.fam.bd.nest)  # No significant difference between methods

## Ordination of beta diversity partitioned by nestedness:
nomicro.fam.comm.nest <- metaMDS(nomicro.fam.dist$beta.jne, 
                                 dist="jaccard", 
                                 k=2,
                                 maxit=999,
                                 trymax=1000,
                                 wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
nomicro.fam.comm.nest$stress
stressplot(nomicro.fam.comm.nest)

## plot site scores as text
ordiplot(nomicro.fam.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
nomicro.NMDS1 <- nomicro.fam.comm.nest$points[,1]
nomicro.NMDS2 <- nomicro.fam.comm.nest$points[,2]
nomicro.fam.nest.NMDS <- data.frame(NMDS1=nomicro.NMDS1, 
                                    NMDS2=nomicro.NMDS2,
                                    Method = method.fam.nomicro.metadata$Method)

## Check data
head(nomicro.fam.nest.NMDS)

## Plot data frame
m8b <- ggplot(nomicro.fam.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Method)) + 
  geom_point() + 
  stat_ellipse() + 
  coord_cartesian(xlim=c(-0.8,0.8), ylim=c(-0.6,0.6)) + 
  labs(subtitle="", x="NMDS1", y="") + 
  scale_colour_manual(name="Method",
                      values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
m8b

## Statistically check difference in nestedness of communities
nomicro.fam.nest.anosim <- anosim(nomicro.fam.dist$beta.jne, method.fam.nomicro.metadata$Method)

## Inspect results
nomicro.fam.nest.anosim
summary(nomicro.fam.nest.anosim)
plot(nomicro.fam.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
nomicro.fam.nest.adonis <- adonis(nomicro.fam.dist$beta.jne ~ Method, method.fam.nomicro.metadata)

## Inspect results
## no summary() or plot() methods included
nomicro.fam.nest.adonis

## Result is not significant. There is no substantial difference in family 
## loss or gain (i.e. nestedness) between methods.


## 3. TOTAL BETA DIVERSITY
nomicro.fam.bd.total <- betadisper(nomicro.fam.dist$beta.jac, method.fam.nomicro.metadata$Method)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(method.fam.nomicro.metadata, nomicro.fam.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(nomicro.fam.bd.total$distances, method.fam.nomicro.metadata$Method, mean)

## Compute variance per group
tapply(nomicro.fam.bd.total$distances, method.fam.nomicro.metadata$Method, var)

## Ordination plot of total beta diversity
plot(nomicro.fam.bd.total)

## Boxplot of total beta diversity
boxplot(nomicro.fam.bd.total, xlab="Method", xaxt="n", bty="n")
axis(side=1, at=c(1,2,3), labels=c("DNA","eDNA","Netting"))

## Plots indicates that there is no difference in total beta diversity 
## between methods.
## Statistically check whether total beta is different between methods
## using standard parametric anova or permutation tests.
anova(nomicro.fam.bd.total)     # Significant difference between methods
permutest(nomicro.fam.bd.total) # Significant difference between methods

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(nomicro.fam.bd.total)  

## No significant difference between any methods

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
nomicro.fam.comm.total <- metaMDS(nomicro.fam.dist$beta.jac,
                                  dist="jaccard", 
                                  k=2,
                                  maxit=999,
                                  trymax=1000,
                                  wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
nomicro.fam.comm.total$stress
stressplot(nomicro.fam.comm.total)

## plot site scores as text
ordiplot(nomicro.fam.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
nomicro.NMDS1 <- nomicro.fam.comm.total$points[,1]
nomicro.NMDS2 <- nomicro.fam.comm.total$points[,2]
nomicro.fam.total.NMDS <- data.frame(NMDS1=nomicro.NMDS1,
                                     NMDS2=nomicro.NMDS2,
                                     Method=method.fam.nomicro.metadata$Method)

## Check data
head(nomicro.fam.total.NMDS)

## Plot data frame
m8c <- ggplot(nomicro.fam.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Method)) + 
  geom_point() + 
  stat_ellipse() + 
  coord_cartesian(xlim=c(-0.6,0.6), ylim=c(-0.6,0.6)) + 
  labs(subtitle="", x="NMDS1", y="") + 
  scale_colour_manual(name="Method",
                      values=c("purple","orange","limegreen")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
m8c

## Statistically check difference in total beta diversity of ponds
method.fam.total.anosim <- anosim(method.fam.dist$beta.jac, method.fam.nomicro.metadata$Method)

## Inspect results
method.fam.total.anosim
summary(method.fam.total.anosim)
plot(method.fam.total.anosim)

## There appears to be a significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
method.fam.total.adonis <- adonis(method.fam.dist$beta.jac ~ Method, method.fam.nomicro.metadata)

## Inspect results
## no summary() or plot() methods included
method.fam.total.adonis

## Again result is significant. There is substantial variation in overall 
## family composition of samples with different monitoring tools. Netting
## and DNA metabarcoding are more similar than eDNA metabarcoding with 
## either method.



#####################################################################
# INFLUENCE OF CRUCIAN CARP ON POND INVERTEBRATES (NO TAXA REMOVED) #
#####################################################################

################################
# SWEEP-NETTING AND MICROSCOPY #
################################

#------------------------#
# SPECIES-LEVEL ANALYSIS #
#------------------------#

## BASIC SUMMARIES
## Total number of individuals per pond:
sum.of.rows <- apply(net.df, 1, sum)
sort(sum.of.rows, decreasing = TRUE) 
sum(sum.of.rows)

## Total number of individuals per species across all sites:
sum.of.columns <- apply(net.df, 2, sum)

## Actual number of individuals:
sort(sum.of.columns, decreasing = TRUE)
sum(sum.of.columns)

## Proportional number of individuals:
sort(round(sum.of.columns/sum(sum.of.columns)*100, 2), decreasing = TRUE)

## Number of ponds where each species occurs:
net.spec.pres <- apply(net.df > 0, 2, sum) 
sort(net.spec.pres, decreasing = TRUE)


## ALPHA DIVERSITY
## Species richness
net.sample.richness <- specnumber(net.spp.pa)
net.sample.richness

## Statistically compare and plot alpha diversity
## Create data frame
net.alpha <- data.frame(net.sample.richness)

## add metadata from external file
net.alpha <- cbind(metadata[,1:2], net.alpha)

## Reset row names of data frame for further indexing
rownames(net.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
net.alpha.regression <- glm.nb(net.sample.richness ~ Crucian, data=net.alpha)
summary(net.alpha.regression)
anova(net.alpha.regression, test = "Chi")
drop1(net.alpha.regression, test = "Chi")
summary(glht(glm.nb(net.sample.richness ~ Crucian, data=net.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
18.253/16
1-pchisq(18.253, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(net.alpha$net.sample.richness ~ fitted(net.alpha.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(net.alpha.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(net.alpha.regression)
shapiro.test(sresid) # P = 0.5825

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ net.alpha.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ net.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(net.alpha.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(net.alpha.regression)
summary(influence)
CD <- cooks.distance(net.alpha.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot species richness
p1 <- ggplot(net.alpha, aes(x=Crucian, y=net.sample.richness)) + 
  geom_jitter(aes(colour=Crucian), width=0.2, cex=2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,30)) +
  labs(title="(a) Species-level", 
       subtitle="(i) Netting and microscopy",
       x="", y=expression(alpha~diversity)) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) +
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p1

## Species richness does not appear to be higher in ponds without crucian
## carp, indicating this species may have a negligible impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot species accumulation curve across ponds
plot(specaccum(net.spp.pa), 
     xlab = "Number of ponds", 
     ylab = "Number of species",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## BETA DIVERSITY
net.multi <- beta.multi(net.spp.pa, index.family="jaccard")
print(net.multi)

## The majority of total beta diversity arises from species turnover
## rather than nestedness.
## Pairwise between-site values of each component of beta diversity
net.dist <- beta.pair(net.spp.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
net.bd.turn <- betadisper(net.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(metadata, net.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(net.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(net.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(net.bd.turn)

## Boxplot of turnover partition
boxplot(net.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is some difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(net.bd.turn)     # No significant difference between ponds
permutest(net.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(net.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
net.comm.turn <- metaMDS(net.dist$beta.jtu,
                         dist="jaccard", 
                         k=2,
                         maxit=999,
                         trymax=1000,
                         wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
net.comm.turn$stress
stressplot(net.comm.turn)

## plot site scores as text
ordiplot(net.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
net.NMDS1 <- net.comm.turn$points[,1]
net.NMDS2 <- net.comm.turn$points[,2]
net.turn.NMDS <- data.frame(NMDS1=net.NMDS1, 
                            NMDS2=net.NMDS2,
                            Crucian = metadata$Crucian)

## Check data
head(net.turn.NMDS)

## Plot data frame
p2a <- ggplot(net.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(x="", y="NMDS2",
       title="(a) Turnover",
       subtitle=" (i) Netting and microscopy") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position = "bottom",
        legend.key=element_blank(),
        legend.key.size = unit(2, 'lines'))
p2a

## Statistically check difference in spatial turnover of communities
net.turn.anosim <- anosim(net.dist$beta.jtu, metadata$Crucian)

## Inspect results
net.turn.anosim
summary(net.turn.anosim)
plot(net.turn.anosim)

## There appears to be a borderline significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
net.turn.adonis <- adonis(net.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
net.turn.adonis

## Now, result is significant. Netting indicates there is substantial 
## difference in species replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## species in ponds with no fish are substituted by species in ponds
## with crucian carp.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
net.bd.nest <- betadisper(net.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.nest <- with(metadata, net.bd.nest)
mod.nest

## Compute mean distance to centroid per group
tapply(net.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(net.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of nestedness partition
plot(net.bd.nest)

## Boxplot of nestedness partition
boxplot(net.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is slight difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(net.bd.nest)     # Significant difference between ponds
permutest(net.bd.nest) # Significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(net.bd.nest)  # Significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
net.comm.nest <- metaMDS(net.dist$beta.jne, 
                         dist="jaccard", 
                         k=2,
                         maxit=999,
                         trymax=1000,
                         wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
net.comm.nest$stress
stressplot(net.comm.nest)

## plot site scores as text
ordiplot(net.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
net.NMDS1 <- net.comm.nest$points[,1]
net.NMDS2 <- net.comm.nest$points[,2]
net.nest.NMDS <- data.frame(NMDS1=net.NMDS1, 
                            NMDS2=net.NMDS2,
                            Crucian = metadata$Crucian)

## Check data
head(net.nest.NMDS)

## Plot data frame
p2b <- ggplot(net.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(title="(b) Nestedness", subtitle="", x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p2b

## Statistically check difference in nestedness of communities
net.nest.anosim <- anosim(net.dist$beta.jne, metadata$Crucian)

## Inspect results
net.nest.anosim
summary(net.nest.anosim)
plot(net.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
net.nest.adonis <- adonis(net.dist$beta.jne ~ Crucian, metadata)
net.nest.adonis

## Again result is not significant. Netting indicates there is no 
## substantial difference in species loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## 3. TOTAL BETA DIVERSITY
net.bd.total <- betadisper(net.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(metadata, net.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(net.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(net.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(net.bd.total)

## Boxplot of total beta diversity
boxplot(net.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is little difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(net.bd.total)     # No significant difference between ponds
permutest(net.bd.total) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(net.bd.total)  # No significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
net.comm.total <- metaMDS(net.dist$beta.jac, 
                          dist="jaccard", 
                          k=2,
                          maxit=999,
                          trymax=1000,
                          wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
net.comm.total$stress
stressplot(net.comm.total)

## plot site scores as text
ordiplot(net.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
net.NMDS1 <- net.comm.total$points[,1]
net.NMDS2 <- net.comm.total$points[,2]
net.total.NMDS <- data.frame(NMDS1=net.NMDS1,
                             NMDS2=net.NMDS2,
                             Crucian = metadata$Crucian)

## Check data
head(net.total.NMDS)

## Plot data frame
p2c <- ggplot(net.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(title=expression(bold("(c) Total"~beta~"Diversity")), 
       subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p2c

## Statistically check difference in total beta diversity of ponds
net.total.anosim <- anosim(net.dist$beta.jac, metadata$Crucian)

## Inspect results
net.total.anosim
summary(net.total.anosim)
plot(net.total.anosim)

## There appears to be borderline significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
net.total.adonis <- adonis(net.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
net.total.adonis

## Result is significant. Netting indicates there is substantial variation in 
## overall species composition of ponds with crucian carp and ponds without 
## crucian carp


#-----------------------#
# FAMILY-LEVEL ANALYSIS #
#-----------------------#

## BASIC SUMMARIES
## Total number of individuals per pond:
sum.of.rows <- apply(net.fam, 1, sum)
sort(sum.of.rows, decreasing = TRUE) 
sum(sum.of.rows)

## Total number of individuals per family across all sites:
sum.of.columns <- apply(net.fam, 2, sum)

## Actual number of individuals:
sort(sum.of.columns, decreasing = TRUE)
sum(sum.of.columns)

## Proportional number of individuals:
sort(round(sum.of.columns/sum(sum.of.columns)*100, 2), decreasing = TRUE)

## Number of ponds where each species occurs:
net.fam.pres <- apply(net.fam > 0, 2, sum) 
sort(net.fam.pres, decreasing = TRUE)


## ALPHA DIVERSITY
## Family richness
net.fam.richness <- specnumber(net.fam.pa)
net.fam.richness

## Statistically compare and plot alpha diversity
## Create data frame
net.fam.alpha <- data.frame(net.fam.richness)

## add metadata from external file
net.fam.alpha <- cbind(metadata[,1:2], net.fam.alpha)

## Reset row names of data frame for further indexing
rownames(net.fam.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
net.fam.regression <- glm.nb(net.fam.richness ~ Crucian, data=net.fam.alpha)
summary(net.fam.regression)
anova(net.fam.regression, test = "Chi")
drop1(net.fam.regression, test = "Chi")
summary(glht(glm(net.fam.richness ~ Crucian, data=net.fam.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
18.666/16
1-pchisq(18.666, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(net.fam.alpha$net.fam.richness ~ fitted(net.fam.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(net.fam.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(net.fam.regression)
shapiro.test(sresid) # P = 0.581

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ net.fam.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ net.fam.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(net.fam.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(net.fam.regression)
summary(influence)
CD <- cooks.distance(net.fam.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot family richness
p3 <- ggplot(net.fam.alpha, aes(x=Crucian, y=net.fam.richness)) + 
  geom_jitter(aes(colour=Crucian), width=0.2, cex=2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,30)) + 
  labs(title="(b) Family-level", 
       subtitle="", 
       x="", y="") + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p3

## Family richness does not appear much higher in ponds without crucian
## carp, indicating this species may have a negligible impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot family accumulation curve across ponds
plot(specaccum(net.fam.pa), 
     xlab = "Number of ponds", 
     ylab = "Number of families",
     ci.type="polygon", ci.col="grey")

## Number of ponds samples was adequate to fully represent invertebrate
## diversity at family-level.


## BETA DIVERSITY
net.fam.multi <- beta.multi(net.fam.pa, index.family="jaccard")
print(net.fam.multi)

## The majority of total beta diversity arises from family turnover
## rather than nestedness.
## Pairwise between-site values of each component of beta diversity
net.fam.dist <- beta.pair(net.fam.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
net.fam.bd.turn <- betadisper(net.fam.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(metadata, net.fam.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(net.fam.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(net.fam.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(net.fam.bd.turn)

## Boxplot of turnover partition
boxplot(net.fam.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is some difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(net.fam.bd.turn)     # No significant difference between ponds
permutest(net.fam.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(net.fam.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
net.fam.comm.turn <- metaMDS(net.fam.dist$beta.jtu, 
                             dist="jaccard", 
                             k=2,
                             maxit=999,
                             trymax=1000,
                             wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
net.fam.comm.turn$stress
stressplot(net.fam.comm.turn)

## plot site scores as text
ordiplot(net.fam.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
net.NMDS1 <- net.fam.comm.turn$points[,1]
net.NMDS2 <- net.fam.comm.turn$points[,2]
net.fam.turn.NMDS <- data.frame(NMDS1=net.NMDS1, 
                                NMDS2=net.NMDS2,
                                Crucian = metadata$Crucian)

## Check data
head(net.fam.turn.NMDS)

## Plot data frame
p4a <- ggplot(net.fam.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(title="(a) Turnover",
       subtitle="(i) Netting and microscopy",
       x="", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) +
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position = "bottom",
        legend.key=element_blank(),
        legend.key.size = unit(2, 'lines'))
p4a

## Statistically check difference in spatial turnover of communities
net.fam.turn.anosim <- anosim(net.fam.dist$beta.jtu, metadata$Crucian)

## Inspect results
net.fam.turn.anosim
summary(net.fam.turn.anosim)
plot(net.fam.turn.anosim)

## There appears to be no significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
net.fam.turn.adonis <- adonis(net.fam.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
net.fam.turn.adonis

## Results are not significant. Netting indicates there is no real 
## difference in family replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## families in ponds with no fish are not substituted by families in ponds
## with crucian carp.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
net.fam.bd.nest <- betadisper(net.fam.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.nest <- with(metadata, net.fam.bd.nest)
mod.nest

## Compute mean distance to centroid per group
tapply(net.fam.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(net.fam.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(net.fam.bd.nest)

## Boxplot of turnover partition
boxplot(net.fam.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is slight difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(net.fam.bd.nest)     # No significant difference between ponds
permutest(net.fam.bd.nest) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(net.fam.bd.nest)  # No significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
net.fam.comm.nest <- metaMDS(net.fam.dist$beta.jne, 
                             dist="jaccard", 
                             k=2,
                             maxit=999,
                             trymax=1000,
                             wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
net.fam.comm.nest$stress
stressplot(net.fam.comm.nest)

## plot site scores as text
ordiplot(net.fam.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
net.NMDS1 <- net.fam.comm.nest$points[,1]
net.NMDS2 <- net.fam.comm.nest$points[,2]
net.fam.nest.NMDS <- data.frame(NMDS1=net.NMDS1, 
                                NMDS2=net.NMDS2,
                                Crucian = metadata$Crucian)

## Check data
head(net.fam.nest.NMDS)

## Plot data frame
p4b <- ggplot(net.fam.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(title="(b) Nestedness", subtitle="", x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p4b

## Statistically check difference in nestedness of communities
net.fam.nest.anosim <- anosim(net.fam.dist$beta.jne, metadata$Crucian)

## Inspect results
net.fam.nest.anosim
summary(net.fam.nest.anosim)
plot(net.fam.nest.anosim)

## There appears to be a significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
net.fam.nest.adonis <- adonis(net.fam.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
net.fam.nest.adonis

## Result is not significant. Netting indicates there is no 
## substantial difference in family loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## 3. TOTAL BETA DIVERSITY
net.fam.bd.total <- betadisper(net.fam.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(metadata, net.fam.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(net.fam.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(net.fam.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(net.fam.bd.total)

## Boxplot of total beta diversity
boxplot(net.fam.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is little difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(net.fam.bd.total)     # No significant difference between ponds
permutest(net.fam.bd.total) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(net.fam.bd.total)  # No significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
net.fam.comm.total <- metaMDS(net.fam.dist$beta.jac, 
                              dist="jaccard", 
                              k=2,
                              maxit=999,
                              trymax=1000,
                              wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
net.fam.comm.total$stress
stressplot(net.fam.comm.total)

## plot site scores as text
ordiplot(net.fam.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
net.NMDS1 <- net.fam.comm.total$points[,1]
net.NMDS2 <- net.fam.comm.total$points[,2]
net.fam.total.NMDS <- data.frame(NMDS1=net.NMDS1,
                                 NMDS2=net.NMDS2,
                                 Crucian = metadata$Crucian)

## Check data
head(net.fam.total.NMDS)

## Plot data frame
p4c <- ggplot(net.fam.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(title=expression(bold("(c) Total"~beta~"Diversity")),
       subtitle="", x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p4c

## Statistically check difference in total beta diversity of ponds
net.fam.total.anosim <- anosim(net.fam.dist$beta.jac, metadata$Crucian)

## Inspect results
net.fam.total.anosim
summary(net.fam.total.anosim)
plot(net.fam.total.anosim)

## There appears to be no significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
net.fam.total.adonis <- adonis(net.fam.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
net.fam.total.adonis

## Result is not significant. Netting indicates there is no variation in 
## overall family composition of ponds with crucian carp and ponds without 
## crucian carp.



#####################
# DNA METABARCODING #
#####################

#------------------------#
# SPECIES-LEVEL ANALYSIS #
#------------------------#

## BASIC SUMMARIES
## Total number of sequences per pond:
sum.of.rows <- apply(DNA.df, 1, sum)
sort(sum.of.rows, decreasing = TRUE) 
sum(sum.of.rows)

## Total number of sequences per species across all sites
sum.of.columns <- apply(DNA.df, 2, sum)

## Actual read counts:
sort(sum.of.columns, decreasing = TRUE)
sum(sum.of.columns)

## Proportional read counts:
sort(round(sum.of.columns/sum(sum.of.columns)*100, 2), decreasing = TRUE)

## Number of ponds where each species occurs
DNA.spec.pres <- apply(DNA.df > 0, 2, sum) 
sort(DNA.spec.pres, decreasing = TRUE)


## ALPHA DIVERSITY
## Basic richness
DNA.spp.richness <- specnumber(DNA.spp.pa)
DNA.spp.richness

## Statistically compare and plot alpha diversity
## Create data frame
DNA.spp.alpha <- data.frame(DNA.spp.richness)

## add metadata from external file
DNA.spp.alpha <- cbind(metadata[,1:2], DNA.spp.alpha)

## Reset row names of data frame for further indexing
rownames(DNA.spp.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
DNA.spp.regression <- glm.nb(DNA.spp.richness ~ Crucian, data=DNA.spp.alpha)
summary(DNA.spp.regression)
anova(DNA.spp.regression, test = "Chi")
drop1(DNA.spp.regression, test = "Chi")
summary(glht(glm(DNA.spp.richness ~ Crucian, data=DNA.spp.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
18.411/16
1-pchisq(18.411, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(DNA.spp.alpha$DNA.spp.richness ~ fitted(DNA.spp.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(DNA.spp.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(DNA.spp.regression)
shapiro.test(sresid) # P = 0.2484

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ DNA.spp.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ DNA.spp.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(DNA.spp.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(DNA.spp.regression)
summary(influence)
CD <- cooks.distance(DNA.spp.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot species richness
p5 <- ggplot(DNA.spp.alpha, aes(x=Crucian, y=DNA.spp.richness)) + 
  geom_jitter(aes(colour=Crucian), cex=2, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  scale_y_continuous(limits=c(0,40)) + 
  labs(subtitle="(ii) DNA metabarcoding",
       x="", y=expression(alpha~diversity)) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p5

## Species richness does not appear higher in ponds without crucian
## carp, indicating this species may have a negligible impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot species accumulation curve across ponds
plot(specaccum(DNA.df), 
     xlab = "Number of ponds", 
     ylab = "Number of species",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## BETA DIVERSITY
DNA.multi <- beta.multi(DNA.spp.pa, index.family="jaccard")
print(DNA.multi)

## The majority of total beta diversity arises from species turnover
## rather than nestedness.
## Pairwise between-site values of each component of beta diversity
DNA.dist <- beta.pair(DNA.spp.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
DNA.bd.turn <- betadisper(DNA.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(metadata, DNA.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(DNA.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(DNA.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(DNA.bd.turn)

## Boxplot of turnover partition
boxplot(DNA.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is no real difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(DNA.bd.turn)     # No significant difference between ponds
permutest(DNA.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(DNA.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
DNA.comm.turn <- metaMDS(DNA.dist$beta.jtu, 
                         dist="jaccard", 
                         k=2,
                         maxit=999,
                         trymax=1000,
                         wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
DNA.comm.turn$stress
stressplot(DNA.comm.turn)

## plot site scores as text
ordiplot(DNA.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
DNA.NMDS1 <- DNA.comm.turn$points[,1]
DNA.NMDS2 <- DNA.comm.turn$points[,2]
DNA.turn.NMDS <- data.frame(NMDS1=DNA.NMDS1, 
                            NMDS2=DNA.NMDS2,
                            Crucian = metadata$Crucian)

## Check data
head(DNA.turn.NMDS)

## Plot data frame
p6a <- ggplot(DNA.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) +
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="(ii) DNA metabarcoding", x="", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, color="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p6a

## Statistically check difference in spatial turnover of communities
DNA.turn.anosim <- anosim(DNA.dist$beta.jtu, metadata$Crucian)

## Inspect results
DNA.turn.anosim
summary(DNA.turn.anosim)
plot(DNA.turn.anosim)

## There appears to be no significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
DNA.turn.adonis <- adonis(DNA.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
DNA.turn.adonis

## Again result is not significant. DNA metabarcoding indicates there is no 
## substantial difference in species replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## species in ponds with no fish are not substituted by species in ponds
## with crucian carp.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
DNA.bd.nest <- betadisper(DNA.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
DNA.nest <- with(metadata, DNA.bd.nest)
DNA.nest

## Compute mean distance to centroid per group
tapply(DNA.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(DNA.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(DNA.bd.nest)

## Boxplot of turnover partition
boxplot(DNA.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is slight difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(DNA.bd.nest)     # Significant difference between ponds
permutest(DNA.bd.nest) # Significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(DNA.bd.nest)  # Significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
DNA.comm.nest <- metaMDS(DNA.dist$beta.jne, 
                         dist="jaccard", 
                         k=2,
                         maxit=999,
                         trymax=1000,
                         wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
DNA.comm.nest$stress
stressplot(DNA.comm.nest)

## plot site scores as text
ordiplot(DNA.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
DNA.NMDS1 <- DNA.comm.nest$points[,1]
DNA.NMDS2 <- DNA.comm.nest$points[,2]
DNA.nest.NMDS <- data.frame(NMDS1=DNA.NMDS1, 
                            NMDS2=DNA.NMDS2,
                            Crucian = metadata$Crucian)

## Check data
head(DNA.nest.NMDS)

## Plot data frame
p6b <- ggplot(DNA.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p6b

## Statistically check difference in nestedness of communities
DNA.nest.anosim <- anosim(DNA.dist$beta.jne, metadata$Crucian)

## Inspect results
DNA.nest.anosim
summary(DNA.nest.anosim)
plot(DNA.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
DNA.nest.adonis <- adonis(DNA.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
DNA.nest.adonis

## Again result is not significant. DNA metabarcoding indicates there is no 
## substantial difference in species loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## 3. TOTAL BETA DIVERSITY
DNA.bd.total <- betadisper(DNA.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(metadata, DNA.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(DNA.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(DNA.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(DNA.bd.total)

## Boxplot of total beta diversity
boxplot(DNA.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is some difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(DNA.bd.total)     # No significant difference between ponds
permutest(DNA.bd.total) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(DNA.bd.total)  # No significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
DNA.comm.total <- metaMDS(DNA.dist$beta.jac, 
                          dist="jaccard", 
                          k=2,
                          maxit=999,
                          trymax=1000,
                          wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
DNA.comm.total$stress
stressplot(DNA.comm.total)

## plot site scores as text
ordiplot(DNA.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
DNA.NMDS1 <- DNA.comm.total$points[,1]
DNA.NMDS2 <- DNA.comm.total$points[,2]
DNA.total.NMDS <- data.frame(NMDS1=DNA.NMDS1,
                             NMDS2=DNA.NMDS2,
                             Crucian = metadata$Crucian)

## Check data
head(DNA.total.NMDS)

## Plot data frame
p6c <- ggplot(DNA.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p6c

## Statistically check difference in total beta diversity of ponds
DNA.total.anosim <- anosim(DNA.dist$beta.jac, metadata$Crucian)

## Inspect results
DNA.total.anosim
summary(DNA.total.anosim)
plot(DNA.total.anosim)

## There appears to be no significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
DNA.total.adonis <- adonis(DNA.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
DNA.total.adonis


#-----------------------#
# FAMILY-LEVEL ANALYSIS #
#-----------------------#

## BASIC SUMMARIES 
## Total number of sequences per pond:
sum.of.rows <- apply(DNA.fam, 1, sum)
sort(sum.of.rows, decreasing = TRUE) 
sum(sum.of.rows)

## Total number of sequences per species across all sites:
sum.of.columns <- apply(DNA.fam, 2, sum)

## Actual read counts:
sort(sum.of.columns, decreasing = TRUE)
sum(sum.of.columns)

## Proportional read counts:
sort(round(sum.of.columns/sum(sum.of.columns)*100, 2), decreasing = TRUE)

## Number of ponds where each species occurs
DNA.fam.pres <- apply(DNA.fam > 0, 2, sum) 
sort(DNA.fam.pres, decreasing = TRUE)


## ALPHA DIVERSITY
## Basic richness
DNA.fam.richness <- specnumber(DNA.fam.pa)
DNA.fam.richness

## Statistically compare and plot alpha diversity
## Create data frame
DNA.fam.alpha <- data.frame(DNA.fam.richness)

## Add sample metadata
DNA.fam.alpha <- cbind(metadata[,1:2], DNA.fam.alpha)

## Reset row names of data frame for further indexing
rownames(DNA.fam.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
DNA.fam.regression <- glm.nb(DNA.fam.richness ~ Crucian, data=DNA.fam.alpha)
summary(DNA.fam.regression)
anova(DNA.fam.regression, test = "Chi")
drop1(DNA.fam.regression, test = "Chi")
summary(glht(glm(DNA.fam.richness ~ Crucian, data=DNA.fam.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
15.825/16
1-pchisq(15.825, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(DNA.fam.alpha$DNA.fam.richness ~ fitted(DNA.fam.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(DNA.fam.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(DNA.fam.regression)
shapiro.test(sresid) # P = 0.5026

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ DNA.fam.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ DNA.fam.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(DNA.fam.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(DNA.fam.regression)
summary(influence)
CD <- cooks.distance(DNA.fam.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot family richness
p7 <- ggplot(DNA.fam.alpha, aes(x=Crucian, y=DNA.fam.richness)) + 
  geom_jitter(aes(colour=Crucian), cex=2, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,40)) + 
  labs(subtitle="", x="", y="") + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p7

## Family richness does not appear higher in ponds without crucian
## carp, indicating this species may have a negligible impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot family accumulation curve across ponds
plot(specaccum(DNA.fam), 
     xlab = "Number of ponds", 
     ylab = "Number of species",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## BETA DIVERSITY
DNA.fam.multi <- beta.multi(DNA.fam.pa, index.family="jaccard")
print(DNA.fam.multi)

## The majority of total beta diversity arises from family turnover
## rather than nestedness.
## Pairwise between-site values of each component of beta diversity
DNA.fam.dist <- beta.pair(DNA.fam.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
DNA.fam.bd.turn <- betadisper(DNA.fam.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(metadata, DNA.fam.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(DNA.fam.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(DNA.fam.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(DNA.fam.bd.turn)

## Boxplot of turnover partition
boxplot(DNA.fam.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is no real difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(DNA.fam.bd.turn)     # No significant difference between ponds
permutest(DNA.fam.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(DNA.fam.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
DNA.fam.comm.turn <- metaMDS(DNA.fam.dist$beta.jtu, 
                             dist="jaccard", 
                             k=2,
                             maxit=999,
                             trymax=1000,
                             wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
DNA.fam.comm.turn$stress
stressplot(DNA.fam.comm.turn)

## plot site scores as text
ordiplot(DNA.fam.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
DNA.NMDS1 <- DNA.fam.comm.turn$points[,1]
DNA.NMDS2 <- DNA.fam.comm.turn$points[,2]
DNA.fam.turn.NMDS <- data.frame(NMDS1=DNA.NMDS1, 
                                NMDS2=DNA.NMDS2,
                                Crucian = metadata$Crucian)

## Check data
head(DNA.fam.turn.NMDS)

## Plot data frame
p8a <- ggplot(DNA.fam.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="(ii) DNA metabarcoding", 
       x="", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p8a

## Statistically check difference in spatial turnover of communities
DNA.fam.turn.anosim <- anosim(DNA.fam.dist$beta.jtu, metadata$Crucian)

## Inspect results
DNA.fam.turn.anosim
summary(DNA.fam.turn.anosim)
plot(DNA.fam.turn.anosim)

## There appears to be no significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
DNA.fam.turn.adonis <- adonis(DNA.fam.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
DNA.fam.turn.adonis

## Again result is not significant. DNA metabarcoding indicates there is no
## substantial difference in family replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## families in ponds with no fish are not substituted by families in ponds
## with crucian carp.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
DNA.fam.bd.nest <- betadisper(DNA.fam.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
DNA.fam.nest <- with(metadata, DNA.fam.bd.nest)
DNA.fam.nest

## Compute mean distance to centroid per group
tapply(DNA.fam.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(DNA.fam.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(DNA.fam.bd.nest)

## Boxplot of turnover partition
boxplot(DNA.fam.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is no difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(DNA.fam.bd.nest)     # No significant difference between ponds
permutest(DNA.fam.bd.nest) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(DNA.fam.bd.nest)  # No significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
DNA.fam.comm.nest <- metaMDS(DNA.fam.dist$beta.jne, 
                             dist="jaccard", 
                             k=2,
                             maxit=999,
                             trymax=1000,
                             wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
DNA.fam.comm.nest$stress
stressplot(DNA.fam.comm.nest)

## plot site scores as text
ordiplot(DNA.fam.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
DNA.NMDS1 <- DNA.fam.comm.nest$points[,1]
DNA.NMDS2 <- DNA.fam.comm.nest$points[,2]
DNA.fam.nest.NMDS <- data.frame(NMDS1=DNA.NMDS1, 
                                NMDS2=DNA.NMDS2,
                                Crucian = metadata$Crucian)

## Check data
head(DNA.fam.nest.NMDS)

## Plot data frame
p8b <- ggplot(DNA.fam.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p8b

## Statistically check difference in nestedness of communities
DNA.fam.nest.anosim <- anosim(DNA.fam.dist$beta.jne, metadata$Crucian)

## Inspect results
DNA.fam.nest.anosim
summary(DNA.fam.nest.anosim)
plot(DNA.fam.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
DNA.fam.nest.adonis <- adonis(DNA.fam.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
DNA.fam.nest.adonis

## Again result is not significant. DNA metabarcoding indicates there is no 
## substantial difference in family loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## 3. TOTAL BETA DIVERSITY
DNA.fam.bd.total <- betadisper(DNA.fam.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(metadata, DNA.fam.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(DNA.fam.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(DNA.fam.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(DNA.fam.bd.total)

## Boxplot of total beta diversity
boxplot(DNA.fam.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is no difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(DNA.fam.bd.total)     # No significant difference between ponds
permutest(DNA.fam.bd.total) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(DNA.fam.bd.total)  # No significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
DNA.fam.comm.total <- metaMDS(DNA.fam.dist$beta.jac, 
                              dist="jaccard", 
                              k=2,
                              maxit=999,
                              trymax=1000,
                              wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
DNA.fam.comm.total$stress
stressplot(DNA.fam.comm.total)

## plot site scores as text
ordiplot(DNA.fam.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
DNA.NMDS1 <- DNA.fam.comm.total$points[,1]
DNA.NMDS2 <- DNA.fam.comm.total$points[,2]
DNA.fam.total.NMDS <- data.frame(NMDS1=DNA.NMDS1,
                                 NMDS2=DNA.NMDS2,
                                 Crucian = metadata$Crucian)

## Check data
head(DNA.fam.total.NMDS)

## Plot data frame
p8c <- ggplot(DNA.fam.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p8c

## Statistically check difference in total beta diversity of ponds
DNA.fam.total.anosim <- anosim(DNA.fam.dist$beta.jac, metadata$Crucian)

## Inspect results
DNA.fam.total.anosim
summary(DNA.fam.total.anosim)
plot(DNA.fam.total.anosim)

## There appears to be no significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
DNA.fam.total.adonis <- adonis(DNA.fam.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
DNA.fam.total.adonis

## Again result is not significant. DNA metabarcoding indicates there is no 
## variation in overall species composition of ponds with crucian carp and 
## ponds without crucian carp.



######################
# eDNA METABARCODING #
######################

#------------------------#
# SPECIES-LEVEL ANALYSIS #
#------------------------#

## BASIC SUMMARIES
## Total number of sequences per pond:
sum.of.rows <- apply(eDNA.df, 1, sum)
sort(sum.of.rows, decreasing = TRUE) 
sum(sum.of.rows)

## Total number of sequences per species across all sites:
sum.of.columns <- apply(eDNA.df, 2, sum)

## Actual read counts:
sort(sum.of.columns, decreasing = TRUE)
sum(sum.of.columns)

## Proportional read counts:
sort(round(sum.of.columns/sum(sum.of.columns)*100, 2), decreasing = TRUE)

## Number of ponds where each species occurs:
eDNA.spec.pres <- apply(eDNA.df > 0, 2, sum) 
sort(eDNA.spec.pres, decreasing = TRUE)


## ALPHA DIVERSITY
## Convert read counts to presence-absence
eDNA.spp.pa <- decostand(eDNA.df, method = "pa", MARGIN = 1)
head(eDNA.spp.pa[,1:3], n = 3)

## Basic richness
eDNA.spp.richness <- specnumber(eDNA.spp.pa)
eDNA.spp.richness

## Statistically compare and plot alpha diversity
## Create data frame
eDNA.spp.alpha <- data.frame(eDNA.spp.richness)

## Add sample metadata
eDNA.spp.alpha <- cbind(metadata[,1:2], eDNA.spp.alpha)

## Reset row names of data frame for further indexing
rownames(eDNA.spp.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
eDNA.spp.regression <- glm.nb(eDNA.spp.richness ~ Crucian, data=eDNA.spp.alpha)
summary(eDNA.spp.regression)
anova(eDNA.spp.regression, test = "Chi")
drop1(eDNA.spp.regression, test = "Chi")
summary(glht(glm.nb(eDNA.spp.richness ~ Crucian, data=eDNA.spp.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
18.724/16
1-pchisq(18.724, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(eDNA.spp.alpha$eDNA.spp.richness ~ fitted(eDNA.spp.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(eDNA.spp.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(eDNA.spp.regression)
shapiro.test(sresid) # P = 0.169

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ eDNA.spp.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ eDNA.spp.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(eDNA.spp.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(eDNA.spp.regression)
summary(influence)
CD <- cooks.distance(eDNA.spp.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot species richness
p9 <- ggplot(eDNA.spp.alpha, aes(x=Crucian, y=eDNA.spp.richness)) + 
  geom_jitter(aes(colour=Crucian), cex=2, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,60), breaks=seq(0,60,10)) + 
  labs(subtitle="(iii) eDNA metabarcoding",
       x="", y=expression(alpha~diversity)) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p9

## Species richness appears to be the same in ponds without crucian
## carp, indicating this species may have little impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot species accumulation curve across ponds
plot(specaccum(eDNA.df), 
     xlab = "Number of ponds", 
     ylab = "Number of species",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## BETA DIVERSITY
eDNA.multi <- beta.multi(eDNA.spp.pa, index.family="jaccard")
print(eDNA.multi)

## The majority of total beta diversity arises from species turnover
## rather than nestedness.
## Pairwise between-site values of each component of beta diversity
eDNA.dist <- beta.pair(eDNA.spp.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
eDNA.bd.turn <- betadisper(eDNA.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(metadata, eDNA.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(eDNA.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(eDNA.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(eDNA.bd.turn)

## Boxplot of turnover partition
boxplot(eDNA.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is substantial difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(eDNA.bd.turn)     # No significant difference between ponds
permutest(eDNA.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(eDNA.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
eDNA.comm.turn <- metaMDS(eDNA.dist$beta.jtu, 
                          dist="jaccard", 
                          k=2,
                          maxit=999,
                          trymax=1000,
                          wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
eDNA.comm.turn$stress
stressplot(eDNA.comm.turn)

## plot site scores as text
ordiplot(eDNA.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
eDNA.NMDS1 <- eDNA.comm.turn$points[,1]
eDNA.NMDS2 <- eDNA.comm.turn$points[,2]
eDNA.turn.NMDS <- data.frame(NMDS1=eDNA.NMDS1,
                             NMDS2=eDNA.NMDS2,
                             Crucian = metadata$Crucian)

## Check data
head(eDNA.turn.NMDS)

## Plot data frame
p10a <- ggplot(eDNA.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="(iii) eDNA metabarcoding", 
       x="", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p10a

## Statistically check difference in spatial turnover of communities
eDNA.turn.anosim <- anosim(eDNA.dist$beta.jtu, metadata$Crucian)

## Inspect results
eDNA.turn.anosim
summary(eDNA.turn.anosim)
plot(eDNA.turn.anosim)

## There appears to be no significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
eDNA.turn.adonis <- adonis(eDNA.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
eDNA.turn.adonis

## Again result is not significant. eDNA metabarcoding indicates there is no 
## substantial difference in species replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## species in ponds with no fish are not substituted by species in ponds
## with crucian carp.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
eDNA.bd.nest <- betadisper(eDNA.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
eDNA.nest <- with(metadata, eDNA.bd.nest)
eDNA.nest

## Compute mean distance to centroid per group
tapply(eDNA.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(eDNA.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(eDNA.bd.nest)

## Boxplot of turnover partition
boxplot(eDNA.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is no difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(eDNA.bd.nest)     # No significant difference between ponds
permutest(eDNA.bd.nest) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(eDNA.bd.nest)  # No significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
eDNA.comm.nest <- metaMDS(eDNA.dist$beta.jne, 
                          dist="jaccard", 
                          k=2,
                          maxit=999,
                          trymax=1000,
                          wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
eDNA.comm.nest$stress
stressplot(eDNA.comm.nest)

## plot site scores as text
ordiplot(eDNA.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
eDNA.NMDS1 <- eDNA.comm.nest$points[,1]
eDNA.NMDS2 <- eDNA.comm.nest$points[,2]
eDNA.nest.NMDS <- data.frame(NMDS1=eDNA.NMDS1, 
                             NMDS2=eDNA.NMDS2,
                             Crucian = metadata$Crucian)

## Check data
head(eDNA.nest.NMDS)

## Plot data frame
p10b <- ggplot(eDNA.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp", 
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p10b

## Statistically check difference in nestedness of communities
eDNA.nest.anosim <- anosim(eDNA.dist$beta.jne, metadata$Crucian)

## Inspect results
eDNA.nest.anosim
summary(eDNA.nest.anosim)
plot(eDNA.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
eDNA.nest.adonis <- adonis(eDNA.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
eDNA.nest.adonis

## Result is not significant. eDNA metabarcoding indicates there is no 
## substantial difference in species loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## 3. TOTAL BETA DIVERSITY
eDNA.bd.total <- betadisper(eDNA.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
eDNA.total <- with(metadata, eDNA.bd.total)
eDNA.total

## Compute mean distance to centroid per group
tapply(eDNA.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(eDNA.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(eDNA.bd.total)

## Boxplot of total beta diversity
boxplot(eDNA.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is substantial difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(eDNA.bd.total)     # Significant difference between ponds
permutest(eDNA.bd.total) # Significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(eDNA.bd.total)  # Significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
eDNA.comm.total <- metaMDS(eDNA.dist$beta.jac, 
                           dist="jaccard", 
                           k=2,
                           maxit=999,
                           trymax=1000,
                           wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
eDNA.comm.total$stress
stressplot(eDNA.comm.total)

## plot site scores as text
ordiplot(eDNA.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
eDNA.NMDS1 <- eDNA.comm.total$points[,1]
eDNA.NMDS2 <- eDNA.comm.total$points[,2]
eDNA.total.NMDS <- data.frame(NMDS1=eDNA.NMDS1,
                              NMDS2=eDNA.NMDS2,
                              Crucian = metadata$Crucian)

## Check data
head(eDNA.total.NMDS)

## Plot data frame
p10c <- ggplot(eDNA.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() +
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p10c

## Statistically check difference in total beta diversity of ponds
eDNA.total.anosim <- anosim(eDNA.dist$beta.jac, metadata$Crucian)

## Inspect results
eDNA.total.anosim
summary(eDNA.total.anosim)
plot(eDNA.total.anosim)

## There appears to be a significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
eDNA.total.adonis <- adonis(eDNA.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
eDNA.total.adonis

## Again result is significant. eDNA metabarcoding indicates there is 
## substantial variation in overall species composition of
## ponds with crucian carp and ponds without crucian carp


#-----------------------#
# FAMILY-LEVEL ANALYSIS #
#-----------------------#

## BASIC SUMMARIES
## Total number of sequences per pond:
sum.of.rows <- apply(eDNA.fam, 1, sum)
sort(sum.of.rows, decreasing = TRUE) 
sum(sum.of.rows)

## Total number of sequences per species across all sites:
sum.of.columns <- apply(eDNA.fam, 2, sum)

## Actual read counts:
sort(sum.of.columns, decreasing = TRUE)
sum(sum.of.columns)

## Proportional read counts:
sort(round(sum.of.columns/sum(sum.of.columns)*100, 2), decreasing = TRUE)

## Number of ponds where each species occurs:
eDNA.fam.pres <- apply(eDNA.fam > 0, 2, sum) 
sort(eDNA.fam.pres, decreasing = TRUE)


## ALPHA DIVERSITY
## Convert read counts to presence-absence
eDNA.fam.pa <- decostand(eDNA.fam, method = "pa", MARGIN = 1)
head(eDNA.fam.pa)

## Basic richness
eDNA.fam.richness <- specnumber(eDNA.fam.pa)
eDNA.fam.richness

## Statistically compare and plot alpha diversity
## Create data frame
eDNA.fam.alpha <- data.frame(eDNA.fam.richness)

## Add sample metadata
eDNA.fam.alpha <- cbind(metadata[,1:2], eDNA.fam.alpha)

## Reset row names of data frame for further indexing
rownames(eDNA.fam.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
eDNA.fam.regression <- glm.nb(eDNA.fam.richness ~ Crucian, data=eDNA.fam.alpha)
summary(eDNA.fam.regression)
anova(eDNA.fam.regression, test = "Chi")
drop1(eDNA.fam.regression, test = "Chi")
summary(glht(glm.nb(eDNA.fam.richness ~ Crucian, data=eDNA.fam.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
18.826/16
1-pchisq(18.826, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(eDNA.fam.alpha$eDNA.fam.richness ~ fitted(eDNA.fam.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(eDNA.fam.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(eDNA.fam.regression)
shapiro.test(sresid) # P = 0.5491

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ eDNA.fam.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ eDNA.fam.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(eDNA.fam.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(eDNA.fam.regression)
summary(influence)
CD <- cooks.distance(eDNA.fam.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot family richness
p11 <- ggplot(eDNA.fam.alpha, aes(x=Crucian, y=eDNA.fam.richness)) + 
  geom_jitter(aes(colour=Crucian), cex=2, width=0.2, show.legend=FALSE) +
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,60), breaks=seq(0,60,10)) + 
  labs(subtitle="",x="", y="") + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) +
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p11

## Family richness does not appear higher in ponds without crucian
## carp, indicating this species may have a negligible impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot family accumulation curve across ponds
plot(specaccum(eDNA.fam), 
     xlab = "Number of ponds", 
     ylab = "Number of species",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## BETA DIVERSITY
eDNA.fam.multi <- beta.multi(eDNA.fam.pa, index.family="jaccard")
print(eDNA.fam.multi)

## The majority of total beta diversity arises from family turnover
## rather than nestedness.
## Pairwise between-site values of each component of beta diversity
eDNA.fam.dist <- beta.pair(eDNA.fam.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
eDNA.fam.bd.turn <- betadisper(eDNA.fam.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
eDNA.turn <- with(metadata, eDNA.fam.bd.turn)
eDNA.turn

## Compute mean distance to centroid per group
tapply(eDNA.fam.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(eDNA.fam.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(eDNA.fam.bd.turn)

## Boxplot of turnover partition
boxplot(eDNA.fam.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is no real difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(eDNA.fam.bd.turn)     # No significant difference between ponds
permutest(eDNA.fam.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(eDNA.fam.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
eDNA.fam.comm.turn <- metaMDS(eDNA.fam.dist$beta.jtu,
                              dist="jaccard", 
                              k=2,
                              maxit=999,
                              trymax=1000,
                              wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
eDNA.fam.comm.turn$stress
stressplot(eDNA.fam.comm.turn)

## plot site scores as text
ordiplot(eDNA.fam.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
eDNA.NMDS1 <- eDNA.fam.comm.turn$points[,1]
eDNA.NMDS2 <- eDNA.fam.comm.turn$points[,2]
eDNA.fam.turn.NMDS <- data.frame(NMDS1=eDNA.NMDS1,
                                 NMDS2=eDNA.NMDS2,
                                 Crucian = metadata$Crucian)

## Check data
head(eDNA.fam.turn.NMDS)

## Plot data frame
p12a <- ggplot(eDNA.fam.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() +
  labs(subtitle="(iii) eDNA metabarcoding", 
       x="", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p12a

## Statistically check difference in spatial turnover of communities
eDNA.fam.turn.anosim <- anosim(eDNA.fam.dist$beta.jtu, metadata$Crucian)

## Inspect results
eDNA.fam.turn.anosim
summary(eDNA.fam.turn.anosim)
plot(eDNA.fam.turn.anosim)

## There appears to be a significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
eDNA.fam.turn.adonis <- adonis(eDNA.fam.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
eDNA.fam.turn.adonis

## Result is not significant. eDNA metabarcoding indicates there is
## no difference in family replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## families in ponds with no fish are not substituted by families in ponds
## with crucian carp.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
eDNA.fam.bd.nest <- betadisper(eDNA.fam.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
eDNA.fam.nest <- with(metadata, eDNA.fam.bd.nest)
eDNA.fam.nest

## Compute mean distance to centroid per group
tapply(eDNA.fam.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(eDNA.fam.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(eDNA.fam.bd.nest)

## Boxplot of turnover partition
boxplot(eDNA.fam.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is slight difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(eDNA.fam.bd.nest)     # No significant difference between ponds
permutest(eDNA.fam.bd.nest) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(eDNA.fam.bd.nest)  # No significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
eDNA.fam.comm.nest <- metaMDS(eDNA.fam.dist$beta.jne, 
                              dist="jaccard", 
                              k=2,
                              maxit=999,
                              trymax=1000,
                              wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
eDNA.fam.comm.nest$stress
stressplot(eDNA.fam.comm.nest)

## plot site scores as text
ordiplot(eDNA.fam.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
eDNA.NMDS1 <- eDNA.fam.comm.nest$points[,1]
eDNA.NMDS2 <- eDNA.fam.comm.nest$points[,2]
eDNA.fam.nest.NMDS <- data.frame(NMDS1=eDNA.NMDS1, 
                                 NMDS2=eDNA.NMDS2,
                                 Crucian = metadata$Crucian)

## Check data
head(eDNA.fam.nest.NMDS)

## Plot data frame
p12b <- ggplot(eDNA.fam.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p12b

## Statistically check difference in nestedness of communities
eDNA.fam.nest.anosim <- anosim(eDNA.fam.dist$beta.jne, metadata$Crucian)

## Inspect results
eDNA.fam.nest.anosim
summary(eDNA.fam.nest.anosim)
plot(eDNA.fam.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
eDNA.fam.nest.adonis <- adonis(eDNA.fam.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
eDNA.fam.nest.adonis

## Again result is not significant. eDNA metabarcoding indicates there is no 
## substantial difference in family loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## 3. TOTAL BETA DIVERSITY
eDNA.fam.bd.total <- betadisper(eDNA.fam.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
eDNA.total <- with(metadata, eDNA.fam.bd.total)
eDNA.total

## Compute mean distance to centroid per group
tapply(eDNA.fam.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(eDNA.fam.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(eDNA.fam.bd.total)

## Boxplot of total beta diversity
boxplot(eDNA.fam.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is large difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(eDNA.fam.bd.total)     # No significant difference between ponds
permutest(eDNA.fam.bd.total) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(eDNA.fam.bd.total)  # No significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
eDNA.fam.comm.total <- metaMDS(eDNA.fam.dist$beta.jac, 
                               dist="jaccard", 
                               k=2,
                               maxit=999,
                               trymax=1000,
                               wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
eDNA.fam.comm.total$stress
stressplot(eDNA.fam.comm.total)

## plot site scores as text
ordiplot(eDNA.fam.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
eDNA.NMDS1 <- eDNA.fam.comm.total$points[,1]
eDNA.NMDS2 <- eDNA.fam.comm.total$points[,2]
eDNA.fam.total.NMDS <- data.frame(NMDS1=eDNA.NMDS1,
                                  NMDS2=eDNA.NMDS2,
                                  Crucian = metadata$Crucian)

## Check data
head(eDNA.fam.total.NMDS)

## Plot data frame
p12c <- ggplot(eDNA.fam.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p12c

## Statistically check difference in total beta diversity of ponds
eDNA.fam.total.anosim <- anosim(eDNA.fam.dist$beta.jac, metadata$Crucian)

## Inspect results
eDNA.fam.total.anosim
summary(eDNA.fam.total.anosim)
plot(eDNA.fam.total.anosim)

## There appears to be a significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
eDNA.fam.total.adonis <- adonis(eDNA.fam.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
eDNA.fam.total.adonis

## Again result is significant. eDNA metabarcoding indicates there 
## is some variation in overall family composition of ponds with 
## crucian carp and ponds without crucian carp.



####################
# METHODS COMBINED #
####################

#------------------------#
# SPECIES-LEVEL ANALYSIS #
#------------------------#

## BASIC SUMMARIES
## Number of ponds where each species occurs
spec.pres <- apply(comb.df > 0, 2, sum) 
sort(spec.pres, decreasing = TRUE)

## Tranform data from counts (specimens or reads) to presence-absence by
## rows (ponds)
comb.pa <- decostand(comb.df, method = "pa", MARGIN = 1)
head(comb.pa[,1:3], n = 3)


## ALPHA DIVERSITY
## Basic richness
comb.richness <- specnumber(comb.pa)
comb.richness

## Statistically compare and plot alpha diversity
## Create data frame
comb.alpha <- data.frame(comb.richness)

## Add sample metadata
comb.alpha <- cbind(metadata[,1:2], comb.alpha)

## Reset row names of data frame for further indexing
rownames(comb.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
comb.regression <- glm.nb(comb.richness ~ Crucian, data=comb.alpha)
summary(comb.regression)
anova(comb.regression, test = "Chi")
drop1(comb.regression, test = "Chi")
summary(glht(glm.nb(comb.richness ~ Crucian, data=comb.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
18.097/16
1-pchisq(18.097, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(comb.alpha$comb.richness ~ fitted(comb.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(comb.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(comb.regression)
shapiro.test(sresid) # P = 0.1786

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ comb.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ comb.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(comb.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(comb.regression)
summary(influence)
CD <- cooks.distance(comb.regression)
plot(CD ~ sresid)

## Plot species richness
p13 <- ggplot(comb.alpha, aes(x=Crucian, y=comb.richness)) + 
  geom_jitter(aes(colour=Crucian), cex=2, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,80), breaks=seq(0,80,10)) + 
  labs(subtitle="(iv) Methods combined",
       x="Crucian carp", y=expression(alpha~diversity)) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p13

## Species richness appears to be the same in ponds with and without crucian
## carp, indicating this species may have negligible impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot species accumulation curve across ponds
plot(specaccum(comb.df), 
     xlab = "Number of ponds", 
     ylab = "Number of species",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## Examine alpha diversity of the major invertebrate groups individually 
## Import data containing order information for invertebrate species and 
## families
spp.order <- read.csv("../Data/Invert_species_order_data_20200818.csv", 
                      header=TRUE)

## Make columns in each dataframe characters rather than factor variables
spp.order$Species <- as.character(spp.order$Species)
spp.order$Group <- as.character(spp.order$Group)

## Transpose dataframe
order.df <- data.frame(t(comb.order.spp))

## Tranform data to presence-absence
order.df[order.df > 0] <- 1

## Make row names a new column in dataframe
order.df <- rownames_to_column(order.df, "Species")

## Add order data to netting data
order.df <- merge(order.df, spp.order, by.x="Species", by.y="Species", all.x=TRUE)
order.df <- order.df[,c(1,20,2:19)]

## Remove species column, and merge data by order
order.df <- order.df[,-1]
order.df <- ddply(order.df, .(Group), numcolwise(sum))

## Make Group column row names and transpose dataframe
order.df <- column_to_rownames(order.df, "Group")
order.dat <- data.frame(t(order.df))

## Make row names column in dataframe
order.dat <- rownames_to_column(order.dat, "Pond")

## Create new column specifying whether pond contains crucian carp
order.dat$Crucian <- factor(ifelse(order.dat$Pond %in% metadata$Site, 
                                   as.character(metadata$Crucian), "NA"))
order.dat <- order.dat[,c(1,23,2:22)]

## Melt dataframe and rename new columns
order.melt <- melt(order.dat, id=c("Pond","Crucian"))
colnames(order.melt)[3:4] <- c("Group","Richness")

## Statistically test for differences in invertebrate species richness
order.alpha.regression <- glm.nb(Richness ~ Group/Crucian, data=order.melt)
summary(order.alpha.regression)
anova(order.alpha.regression, test = "Chi")
drop1(order.alpha.regression, test = "Chi")

## Check model meets GLM assumptions
## Test for overdispersion
312.16/336
1-pchisq(312.16, df=336)  # not overdispersed

## Plot the fitted data against the observed data
plot(order.melt$Richness ~ fitted(order.alpha.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(order.alpha.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(order.alpha.regression)
shapiro.test(sresid) # P = 2.829e-13

## Some deviation from normality as residuals are normally distributed
## therefore model may not be reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ order.alpha.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ order.melt$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity
plot(sresid ~ order.melt$Group, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
vif(order.alpha.regression)

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(order.alpha.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(order.alpha.regression)
summary(influence)
CD <- cooks.distance(order.alpha.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Create additional dataframes to hold text annotations
text1 <- data.frame(label=c("","",".","","","","","","",""),
                    Group=c("Annelida","Arachnida",
                            "Coleoptera","Collembola",
                            "Crustacea","Diptera",
                            "Ephemeroptera","Gastrotricha",
                            "Hemiptera","Hirudinea"))

text2 <- data.frame(label=c("","","","***","","","","","","",""),
                    Group=c("Hymenoptera","Lepidoptera",
                            "Megaloptera","Mollusca",
                            "Nematoda","Odonata",
                            "Platyhelminthes","Psocoptera",
                            "Rotifera","Tardigrada",
                            "Trichoptera"))

## Plot Group richness in crucian carp and non-fish ponds
p14a <- ggplot(order.melt[order.melt$Group %in% c("Annelida","Arachnida",
                                                  "Coleoptera","Collembola",
                                                  "Crustacea","Diptera",
                                                  "Ephemeroptera","Gastrotricha",
                                                  "Hemiptera","Hirudinea"),],
               aes(x=Crucian, y=Richness)) + 
  geom_jitter(aes(colour=Crucian), width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  coord_cartesian(ylim=c(0,15)) + 
  scale_y_continuous(breaks=seq(0,15,1)) + 
  labs(title="(a) Species-level",
       x="", y=expression(alpha~diversity)) + 
  geom_text(data=text1, mapping=aes(x=1.5, y=15, label=label), size=10) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin=unit(c(2,0,2,0), "mm")),
        axis.line.x = element_line(colour="black", size=0.5, linetype="solid"),
        axis.line.y = element_line(colour="black", size=0.5, linetype="solid"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin=unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=20),
        strip.text.x = element_text(size=12),
        legend.position = "none",
        text = element_text(size=20)) + 
  facet_grid(. ~ Group)
p14a

p14b <- ggplot(order.melt[order.melt$Group %in% c("Hymenoptera","Lepidoptera",
                                                  "Megaloptera","Mollusca",
                                                  "Nematoda","Odonata",
                                                  "Platyhelminthes","Psocoptera",
                                                  "Rotifera","Tardigrada",
                                                  "Trichoptera"),],
               aes(x=Crucian, y=Richness)) + 
  geom_jitter(aes(colour=Crucian), width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  coord_cartesian(ylim=c(0,15)) + 
  scale_y_continuous(breaks=seq(0,15,1)) + 
  labs(x="Crucian carp", y=expression(alpha~diversity)) + 
  geom_text(data=text2, mapping=aes(x=1.5, y=15, label=label), size=10) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"), labels=c("Absent","Present")) + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        axis.line.x = element_line(colour="black", size=0.5, linetype="solid"),
        axis.line.y = element_line(colour="black", size=0.5, linetype="solid"),
        axis.title.x = element_text(margin=unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin=unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=20),
        strip.text.x = element_text(size=12),
        legend.position = "none",
        text = element_text(size=20)) + 
  facet_grid(. ~ Group)
p14b


## BETA DIVERSITY
## First, tidy up environmental metadata
## Make new dataframe containing only variables relevant to RDA
env.data <- metadata[,c(1:2,4:10)]

## Make first column containing pond ID the row names of dataframe
env.data <- column_to_rownames(env.data, "Site")

## Code crucian carp presence-absence as numeric (0, 1) for RDA
env.data$Crucian <- gsub("Y", "1", env.data$Crucian)
env.data$Crucian <- gsub("N", "0", env.data$Crucian)
env.data$Crucian <- as.numeric(env.data$Crucian)

## Check environmental data for collinearity
plot(env.data[,1:8], 
     lower.panel = panel.cor, 
     upper.panel = panel.smooth2)

## Percentage of submerged macrophyte cover is strongly collinear (+/- 0.5) 
## with other variables. Percentage of emergent cover and emergent perimeter 
## are collinear. 
## Remove percentages of submerged macrophyte cover and emergent perimeter
## (less informative than emergent cover).
env.data <- env.data[,-c(5,7)]

## Check for collinearity
plot(env.data[,1:6], 
     lower.panel = panel.cor, 
     upper.panel = panel.smooth2)

## Now, examine beta diversity across all ponds
comb.multi <- beta.multi(comb.pa, index.family="jaccard")
print(comb.multi)

## The majority of total beta diversity arises from species turnover
## rather than nestedness.
## Now get pairwise between-site values of each component of beta diversity
comb.dist <- beta.pair(comb.pa, index.family="jaccard")

## Convert each partition of beta diversity to a matrix if you want to 
## write as a csv file
#comb.turn <- as.matrix(dist(net.beta$beta.jtu))
#comb.nest <- as.matrix(dist(net.beta$beta.jne))
#comb.total <- as.matrix(dist(net.beta$beta.jac))


## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)

## First, test whether ponds with and without crucian carp have different
## communities using PERMANOVA and visualise differences with NMDS
comb.bd.turn <- betadisper(comb.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(metadata, comb.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(comb.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(comb.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(comb.bd.turn)

## Boxplot of turnover partition
boxplot(comb.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is substantial difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(comb.bd.turn)     # No significant difference between ponds
permutest(comb.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(comb.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
comb.bd.turn <- metaMDS(comb.dist$beta.jtu, 
                        dist="jaccard", 
                        k=2,
                        maxit=999,
                        trymax=1000,
                        wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
comb.bd.turn$stress
stressplot(comb.bd.turn)

## plot site scores as text
ordiplot(comb.bd.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
comb.NMDS1 <- comb.bd.turn$points[,1]
comb.NMDS2 <- comb.bd.turn$points[,2]
comb.turn.NMDS <- data.frame(NMDS1=comb.NMDS1, 
                             NMDS2=comb.NMDS2,
                             Crucian = metadata$Crucian)

## Check data
head(comb.turn.NMDS)

## Plot data frame
p15a <- ggplot(comb.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() +
  labs(subtitle="(iv) Methods combined", 
       x="NMDS1", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p15a

## Statistically check difference in spatial turnover of communities
comb.turn.anosim <- anosim(comb.dist$beta.jtu, metadata$Crucian)

## Inspect results
comb.turn.anosim
summary(comb.turn.anosim)
plot(comb.turn.anosim)

## There appears to be a significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
comb.turn.adonis <- adonis(comb.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
comb.turn.adonis

## Again result is significant. The combined data indicates there is a
## substantial difference in species replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## species in ponds with no fish are substituted by species in ponds
## with crucian carp.

## Use RDA and variance partitioning to examine the influence of biotic and 
## abiotic variables on turnover partition of beta diversity.
## First, compute principal coordinate decomposition (classical scaling) for 
## turnover distance matrix using the pcoa() function from the ape package. Use 
## lingoes correction to account for negative eigenvalues.
pcoa.jtu <- pcoa(comb.dist$beta.jtu, correction="lingoes", rn=NULL)

## Inspect outputs
pcoa.jtu$correction
pcoa.jtu$note
pcoa.jtu$trace

## The principal component coordinates for each diversity matrix can be written
## to a csv file for inspection
#write.csv(pcoa.jtu$vectors, "../Data/Combined_data_species_turnover_eigenvectors.csv")

## Create groups of variables for RDA
## Biotic:
biotic <- data.frame(env.data$Crucian)
colnames(biotic) <- "Crucian"

## Abiotic:
abiotic <- env.data[,c(2:6)]

## log10 transform abiotic data to remove units
abiotic <- log10(abiotic+1)

## Perform forward model selection for each group of explanatory variables
## Adjusted R-squared will be used as a goodness-of-fit measure. Model selection
## needs to be performed for each group of variables, excluding those containing
## only one variable.
abiotic.rda.null <- rda(pcoa.jtu$vectors ~ 1, abiotic)  # null model
abiotic.rda.all <- rda(pcoa.jtu$vectors ~ ., abiotic)   # model with all variables
step.env <- ordiR2step(abiotic.rda.null, abiotic.rda.all, 
                       direction="forward", perm.max = 1000)

## Summary table and plot:
step.env$anova
plot(step.env)

## Indicates that only pond area should be included from this group of variables
## in model.

## Run RDA with selected biotic and abiotic variables
mod <- varpart(pcoa.jtu$vectors, ~ biotic$Crucian, ~ abiotic$Area)
mod
plot(mod, digits=2, cutoff=0, bg=c("skyblue2","grey60"), 
     Xnames=c("Biotic", "Abiotic"), id.size=1)

## The output indicates that crucian carp presence-absence has a larger
## individual fraction of explained variance. We can also see that there is 
## a fraction of shared variance which indicates that biotic and abiotic 
## variable effects are somewhat dependent on one another.

## Test for significance of the whole model:
mod2 <- rda(pcoa.jtu$vectors ~ biotic$Crucian + abiotic$Area)
summary(mod2)
coef(mod2)       # canonical coefficients
vif.cca(mod2)    # no collinearity present in model
mod2             # 17.98% variance is explained by variables in model
anova(mod2)      # Model is significant
RsquareAdj(mod2) # 7.04% variance explained
plot(mod2)

## Inertia is another name for variation or variance in this case. Total refers to 
## total variance, Constrained refers to the amount of variance explained by the 
## explanatory variables, Unconstrained? refers to the residual variance. Constrained + 
## Unconstrained = Total. An R2 statistic can be derived simply as Constrained/Total. 
## The function RsquareAdj computes R2 and R2-adjusted. The variable Rank indicates 
## the number of variables included. The eigenvalues are displayed for both the 
## constrained and unconstrained axes. In this context, these eigenvalues indicate 
## how much variance each of the axes contribute to.

## Testing significance of the individual groups:
rda.biotic <- rda(formula = pcoa.jtu$vectors ~ biotic$Crucian + Condition(abiotic$Area))
summary(rda.biotic)
anova(rda.biotic)
## Conditioned abiotic variable explains less variation (8.6%) than crucian
## carp presence-absence (9.3%).

rda.abiotic <- rda(formula = pcoa.jtu$vectors ~ abiotic$Area + Condition(biotic$Crucian))
summary(rda.abiotic)
anova(rda.abiotic)
## Abiotic variable explains less variation (6.8%) than crucian carp 
## presence-absence (11.1%).

## Check whether any variables should be dropped from model
drop1(mod2, test="perm")
## Pond area should potentially be dropped as not significant

## Plot model result to get a sense of which variables are correlating with which 
## ponds along which axes
plot(mod2, type='n', scaling=1)
orditorp(mod2, display='sp', cex=0.5, scaling=1, col='blue')
text(mod2, display='cn', col='red')

## Plot of variables against sites
plot(mod2, type = "text")

## Perform hypothesis testing
anova(mod2, perm=1000)               # overall model fit
anova(mod2, by="margin", perm=1000)  # significance of marginal effects (Type III)
anova(mod2, by="term", perm=1000)    # significance of each term in model
anova(mod2, by="axis", perm=1000)    # significance of each axis


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
comb.bd.nest <- betadisper(comb.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
comb.nest <- with(metadata, comb.bd.nest)
comb.nest

## Compute mean distance to centroid per group
tapply(comb.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(comb.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(comb.bd.nest)

## Boxplot of turnover partition
boxplot(comb.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is no difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(comb.bd.nest)     # No significant difference between ponds
permutest(comb.bd.nest) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(comb.bd.nest)  # No significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
comb.comm.nest <- metaMDS(comb.dist$beta.jne, 
                          dist="jaccard", 
                          k=2,
                          maxit=999,
                          trymax=1000,
                          wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
comb.comm.nest$stress
stressplot(comb.comm.nest)

## plot site scores as text
ordiplot(comb.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
comb.NMDS1 <- comb.comm.nest$points[,1]
comb.NMDS2 <- comb.comm.nest$points[,2]
comb.nest.NMDS <- data.frame(NMDS1=comb.NMDS1, 
                             NMDS2=comb.NMDS2,
                             Crucian = metadata$Crucian)

## Check data
head(comb.nest.NMDS)

## Plot data frame
p15b <- ggplot(comb.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="NMDS1", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p15b

## Statistically check difference in nestedness of communities
comb.nest.anosim <- anosim(comb.dist$beta.jne, metadata$Crucian)

## Inspect results
comb.nest.anosim
summary(comb.nest.anosim)
plot(comb.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
comb.nest.adonis <- adonis(comb.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
comb.nest.adonis

## Result is not significant. The combined data indicates there is no
## substantial difference in species loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## Use RDA and variance partitioning to examine the influence of biotic and 
## abiotic variables on nestedness partition of beta diversity.
## First, compute principal coordinate decomposition (classical scaling) for 
## nestedness distance matrix using the pcoa() function from the ape package. Use 
## lingoes correction to account for negative eigenvalues.
pcoa.jne <- pcoa(comb.dist$beta.jne, correction="lingoes", rn=NULL)

## Inspect outputs
pcoa.jne$correction
pcoa.jne$note
pcoa.jne$trace

## The principal component coordinates for each diversity matrix can be written
## to a csv file for inspection
#write.csv(pcoa.jne$vectors, "../Data/Combined_data_species_nestedness_eigenvectors.csv")

## Groups of variables were previously created for RDA of turnover
## Perform forward model selection for each group of explanatory variables
## Adjusted R-squared will be used as a goodness-of-fit measure. Model selection
## needs to be performed for each group of variables, excluding those containing
## only one variable.
abiotic.rda.null <- rda(pcoa.jne$vectors ~ 1, abiotic)  # null model
abiotic.rda.all <- rda(pcoa.jne$vectors ~ ., abiotic)   # model with all variables
step.env <- ordiR2step(abiotic.rda.null, abiotic.rda.all, 
                       direction="forward", perm.max = 1000)

## Summary table and plot:
step.env$anova
plot(step.env)
## Indicates that no abiotic variables should be included in model

## Run RDA with biotic variable
## Test for significance of the whole model:
mod2 <- rda(pcoa.jne$vectors ~ biotic$Crucian)
summary(mod2)
coef(mod2)       # canonical coefficients
vif.cca(mod2)    # no collinearity present in model
mod2             # 0.002% variance is explained by variables in model
anova(mod2)      # Model is not significant
RsquareAdj(mod2) # -0.04% variance explained?
plot(mod2)

## Check whether any variables should be dropped from model
drop1(mod2, test="perm")
## Crucian carp presence-absence not significant

## Plot model result to get a sense of which variables are correlating with which 
## ponds along which axes
plot(mod2, type='n', scaling=1)
orditorp(mod2, display='sp', cex=0.5, scaling=1, col='blue')
text(mod2, display='cn', col='red')

## Plot of variables against sites
plot(mod2, type = "text")

## Perform hypothesis testing
anova(mod2, perm=1000)               # overall model fit
anova(mod2, by="margin", perm=1000)  # significance of marginal effects (Type III)
anova(mod2, by="term", perm=1000)    # significance of each term in model
anova(mod2, by="axis", perm=1000)    # significance of each axis


## 3. TOTAL BETA DIVERSITY
comb.bd.total <- betadisper(comb.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(metadata, comb.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(comb.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(comb.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(comb.bd.total)

## Boxplot of total beta diversity
boxplot(comb.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is substantial difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(comb.bd.total)     # Significant difference between ponds
permutest(comb.bd.total) # Significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(comb.bd.total)  # Significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
comb.comm.total <- metaMDS(comb.dist$beta.jac, 
                           dist="jaccard", 
                           k=2,
                           maxit=999,
                           trymax=1000,
                           wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
comb.comm.total$stress
stressplot(comb.comm.total)

## plot site scores as text
ordiplot(comb.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
comb.NMDS1 <- comb.comm.total$points[,1]
comb.NMDS2 <- comb.comm.total$points[,2]
comb.total.NMDS <- data.frame(NMDS1=comb.NMDS1,
                              NMDS2=comb.NMDS2,
                              Crucian = metadata$Crucian)

## Check data
head(comb.total.NMDS)

## Plot data frame
p15c <- ggplot(comb.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="NMDS1", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p15c

## Statistically check difference in total beta diversity of ponds
comb.total.anosim <- anosim(comb.dist$beta.jac, metadata$Crucian)

## Inspect results
comb.total.anosim
summary(comb.total.anosim)
plot(comb.total.anosim)

## There appears to be a significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
comb.total.adonis <- adonis(comb.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
comb.total.adonis

## Again result is significant. The combined data indicates there is 
## substantial variation in overall species composition of
## ponds with crucian carp and ponds without crucian carp


## Use RDA and variance partitioning to examine the influence of biotic and 
## abiotic variables on total beta diversity

## First, compute principal coordinate decomposition (classical scaling) for 
## total beta distance matrix using the pcoa() function from the ape package. Use 
## lingoes correction to account for negative eigenvalues.
pcoa.jac <- pcoa(comb.dist$beta.jac, correction="lingoes", rn=NULL)

## Inspect outputs
pcoa.jac$correction
pcoa.jac$note
pcoa.jac$trace

## The principal component coordinates for each diversity matrix can be written
## to a csv file for inspection
#write.csv(pcoa.jac$vectors, "../Data/Combined_data_species_totalbeta_eigenvectors.csv")

## Groups of variables were previously created for RDA of turnover
## Perform forward model selection for each group of explanatory variables
## Adjusted R-squared will be used as a goodness-of-fit measure. Model selection
## needs to be performed for each group of variables, excluding those containing
## only one variable.
abiotic.rda.null <- rda(pcoa.jac$vectors ~ 1, abiotic)  # null model
abiotic.rda.all <- rda(pcoa.jac$vectors ~ ., abiotic)   # model with all variables
step.env <- ordiR2step(abiotic.rda.null, abiotic.rda.all, 
                       direction="forward", perm.max = 1000)

## Summary table and plot:
step.env$anova
plot(step.env)
## Indicates that only pond area should be included in model.

## Run RDA with selected biotic and abiotic variables
## Varpart between the groups
mod <- varpart(pcoa.jac$vectors, ~ biotic$Crucian, ~ abiotic$Area)
mod
plot(mod, digits=2, cutoff=0, bg=c("skyblue2","grey60"), 
     Xnames=c("Biotic", "Abiotic"), id.size=1)

## The output indicates that crucian carp presence-absence has a larger
## individual fraction of explained variance. We can also see that there is 
## a fraction of shared variance which indicates that biotic and abiotic 
## variable effects are somewhat dependent on one another.

## Test for significance of the whole model:
mod2 <- rda(pcoa.jac$vectors ~ biotic$Crucian + abiotic$Area)
summary(mod2)
coef(mod2)       # canonical coefficients
vif.cca(mod2)    # no collinearity present in model
mod2             # 15.95% variance is explained by variables in model
anova(mod2)      # Model is significant
RsquareAdj(mod2) # 4.75% variance explained
plot(mod2)

## Testing significance of the individual groups:
rda.biotic <- rda(formula = pcoa.jac$vectors ~ biotic$Crucian + Condition(abiotic$Area))
summary(rda.biotic)
anova(rda.biotic)
## Conditioned abiotic variable explains less variation (7.7%) than crucian
## carp presence-absence (8.3%).

rda.abiotic <- rda(formula = pcoa.jac$vectors ~ abiotic$Area + Condition(biotic$Crucian))
summary(rda.abiotic)
anova(rda.abiotic)
## Abiotic variable explains less variation (6.3%) than crucian carp 
## presence-absence (9.7%).

## Check whether any variables should be dropped from model
drop1(mod2, test="perm")
## Pond area should potentially be dropped as not significant

## Plot model result to get a sense of which variables are correlating with which 
## ponds along which axes
plot(mod2, type='n', scaling=1)
orditorp(mod2, display='sp', cex=0.5, scaling=1, col='blue')
text(mod2, display='cn', col='red')

## Plot of variables against sites
plot(mod2, type = "text")

## Perform hypothesis testing
anova(mod2, perm=1000)               # overall model fit
anova(mod2, by="margin", perm=1000)  # significance of marginal effects (Type III)
anova(mod2, by="term", perm=1000)    # significance of each term in model
anova(mod2, by="axis", perm=1000)    # significance of each axis


#-----------------------#
# FAMILY-LEVEL ANALYSIS #
#-----------------------#

## BASIC SUMMARIES
## Number of ponds where each family occurs:
fam.pres <- apply(comb.fam.df > 0, 2, sum) 
sort(fam.pres, decreasing = TRUE)[1:18]

## Tranform datafrom counts (specimens or reads) to presence-absence by
## rows (ponds)
comb.fam.pa <- decostand(comb.fam.df, method = "pa", MARGIN = 1)
head(comb.fam.pa[,1:3], n = 3)


## ALPHA DIVERSITY
## Basic richness
comb.fam.richness <- specnumber(comb.fam.pa)
comb.fam.richness

## Statistically compare and plot alpha diversity
## Create data frame
comb.fam.alpha <- data.frame(comb.fam.richness)

## add metadata from external file
comb.fam.alpha <- cbind(metadata[,1:2], comb.fam.alpha)

## Reset row names of data frame for further indexing
rownames(comb.fam.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
comb.fam.regression <- glm.nb(comb.fam.richness ~ Crucian, data=comb.fam.alpha)
summary(comb.fam.regression)
anova(comb.fam.regression, test = "Chi")
drop1(comb.fam.regression, test = "Chi")
summary(glht(glm.nb(comb.fam.richness ~ Crucian, data=comb.fam.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
18.011/16
1-pchisq(18.011, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(comb.fam.alpha$comb.fam.richness ~ fitted(comb.fam.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(comb.fam.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(comb.fam.regression)
shapiro.test(sresid) # P = 0.3465

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ comb.fam.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ comb.fam.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(comb.fam.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(comb.fam.regression)
summary(influence)
CD <- cooks.distance(comb.fam.regression)
plot(CD ~ sresid)

## Plot family richness
p16 <- ggplot(comb.fam.alpha, aes(x=Crucian, y=comb.fam.richness)) + 
  geom_jitter(aes(colour=Crucian), cex=2, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,80), breaks=seq(0,80,10)) + 
  labs(subtitle="",x="Crucian carp", y="") + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p16

## Family richness appears to be the same in ponds with and without crucian
## carp, indicating this species may have negligible impact on the 
## invertebrate community.
## Was sample size enough to accurately represent invertebrate diversity?
## Plot family accumulation curve across ponds
plot(specaccum(comb.fam.df), 
     xlab = "Number of ponds", 
     ylab = "Number of families",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## Compare alpha diversity of major invertebrate groups individually in 
## ponds with and without crucian carp.
## Import data containing order information for invertebrate families
fam.order <- read.csv("../Data/Invert_family_order_data_20200818.csv", 
                      header=TRUE)

## Make columns in each dataframe characters rather than factor variables
fam.order$Family <- as.character(fam.order$Family)
fam.order$Group <- as.character(fam.order$Group)

## Transpose dataframe
order.fam.df <- data.frame(t(comb.order.fam))

## Tranform data to presence-absence
order.fam.df[order.fam.df > 0] <- 1

## Make row names a new column in dataframe
order.fam.df <- rownames_to_column(order.fam.df, "Family")

## Add order data to family data
order.fam.df <- merge(order.fam.df, fam.order, by.x="Family", by.y="Family", all.x=TRUE)
order.fam.df <- order.fam.df[,c(1,20,2:19)]

## Remove species column, and merge data by order
order.fam.df <- order.fam.df[,-1]
order.fam.df <- ddply(order.fam.df, .(Group), numcolwise(sum))

## Make Group column row names and transpose dataframe
order.fam.df <- column_to_rownames(order.fam.df, "Group")
order.fam.dat <- data.frame(t(order.fam.df))

## Make row names column in dataframe
order.fam.dat <- rownames_to_column(order.fam.dat, "Pond")

## Create new column specifying whether pond contains crucian carp
order.fam.dat$Crucian <- factor(ifelse(order.fam.dat$Pond %in% metadata$Site, 
                                       as.character(metadata$Crucian), "NA"))
order.fam.dat <- order.fam.dat[,c(1,24,2:23)]

## Melt dataframe and rename new columns
order.fam.melt <- melt(order.fam.dat, id=c("Pond","Crucian"))
colnames(order.fam.melt)[3:4] <- c("Group","Richness")

## Statistically test for differences in invertebrate species richness
order.alpha.regression <- glm.nb(Richness ~ Group/Crucian, data=order.fam.melt)
summary(order.alpha.regression)
anova(order.alpha.regression, test = "Chi")
drop1(order.alpha.regression, test = "Chi")

## Check model meets GLM assumptions
## Test for overdispersion
251.25/352
1-pchisq(251.25, df=352)  # not overdispersed

## Plot the fitted data against the observed data
plot(order.fam.melt$Richness ~ fitted(order.alpha.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(order.alpha.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(order.alpha.regression)
shapiro.test(sresid) # P = 1.82e-13

## Some deviation from normality as residuals are normally distributed
## therefore model may not be reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ order.alpha.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ order.fam.melt$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity
plot(sresid ~ order.fam.melt$Group, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
vif(order.alpha.regression)

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(order.alpha.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(order.alpha.regression)
summary(influence)
CD <- cooks.distance(order.alpha.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot Group richness in crucian carp and non-fish ponds
p17a <- ggplot(order.fam.melt[order.fam.melt$Group %in% c("Annelida","Arachnida",
                                                          "Coleoptera","Collembola",
                                                          "Crustacea","Diptera",
                                                          "Ephemeroptera","Gastrotricha",
                                                          "Hemiptera","Hirudinea"),],
               aes(x=Crucian, y=Richness)) + 
  geom_jitter(aes(colour=Crucian), width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  coord_cartesian(ylim=c(0,10)) + 
  scale_y_continuous(breaks=seq(0,15,1)) + 
  labs(title="(b) Family-level",
       x="", y=expression(alpha~diversity)) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin=unit(c(2,0,2,0), "mm")),
        axis.line.x = element_line(colour="black", size=0.5, linetype="solid"),
        axis.line.y = element_line(colour="black", size=0.5, linetype="solid"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin=unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=20),
        strip.text.x = element_text(size=12),
        legend.position = "none",
        text = element_text(size=20)) + 
  facet_grid(. ~ Group)
p17a

p17b <- ggplot(order.fam.melt[order.fam.melt$Group %in% c("Hymenoptera","Lepidoptera",
                                                          "Megaloptera","Mollusca",
                                                          "Nematoda","Odonata",
                                                          "Platyhelminthes","Psocoptera",
                                                          "Rotifera","Tardigrada",
                                                          "Thysanoptera","Trichoptera"),],
                aes(x=Crucian, y=Richness)) + 
  geom_jitter(aes(colour=Crucian), width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  coord_cartesian(ylim=c(0,10)) + 
  scale_y_continuous(breaks=seq(0,15,1)) + 
  labs(x="Crucian carp", y=expression(alpha~diversity)) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        axis.line.x = element_line(colour="black", size=0.5, linetype="solid"),
        axis.line.y = element_line(colour="black", size=0.5, linetype="solid"),
        axis.title.x = element_text(margin=unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin=unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=20),
        strip.text.x = element_text(size=12),
        legend.position = "none",
        text = element_text(size=20)) + 
  facet_grid(. ~ Group)
p17b


## BETA DIVERSITY
## Examine beta diversity across all ponds
comb.fam.multi <- beta.multi(comb.fam.pa, index.family="jaccard")
print(comb.fam.multi)

## The majority of total beta diversity arises from species turnover
## rather than nestedness.
## Now get pairwise between-site values of each component of beta diversity
comb.fam.dist <- beta.pair(comb.fam.pa, index.family="jaccard")

## Convert each partition of beta diversity to a matrix if you want to 
## write as a csv file
#comb.fam.turn <- as.matrix(dist(comb.fam.dist$beta.jtu))
#comb.fam.nest <- as.matrix(dist(comb.fam.dist$beta.jne))
#comb.fam.total <- as.matrix(dist(comb.fam.dist$beta.jac))


## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)

## First, test whether ponds with and without crucian carp have different
## communities using PERMANOVA and visualise differences with NMDS
comb.fam.bd.turn <- betadisper(comb.fam.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(metadata, comb.fam.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(comb.fam.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(comb.fam.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(comb.fam.bd.turn)

## Boxplot of turnover partition
boxplot(comb.fam.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is some difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(comb.fam.bd.turn)     # No significant difference between ponds
permutest(comb.fam.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(comb.fam.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
comb.fam.bd.turn <- metaMDS(comb.fam.dist$beta.jtu, 
                            dist="jaccard", 
                            k=2,
                            maxit=999,
                            trymax=1000,
                            wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
comb.fam.bd.turn$stress
stressplot(comb.fam.bd.turn)

## plot site scores as text
ordiplot(comb.fam.bd.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
comb.NMDS1 <- comb.fam.bd.turn$points[,1]
comb.NMDS2 <- comb.fam.bd.turn$points[,2]
comb.fam.turn.NMDS <- data.frame(NMDS1=comb.NMDS1, 
                                 NMDS2=comb.NMDS2,
                                 Crucian = metadata$Crucian)

## Check data
head(comb.fam.turn.NMDS)

## Plot data frame
p18a <- ggplot(comb.fam.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="(iv) Methods combined", 
       x="NMDS1", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p18a

## Statistically check difference in spatial turnover of communities
comb.fam.turn.anosim <- anosim(comb.fam.dist$beta.jtu, metadata$Crucian)

## Inspect results
comb.fam.turn.anosim
summary(comb.fam.turn.anosim)
plot(comb.fam.turn.anosim)

## There appears to be no significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
comb.fam.turn.adonis <- adonis(comb.fam.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
comb.fam.turn.adonis

## Again result is not significant. The combined data indicates there is no
## substantial difference in species replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## species in ponds with no fish are not substituted by species in ponds
## with crucian carp.


## Use RDA and variance partitioning to examine the influence of biotic and 
## abiotic variables on turnover partition of beta diversity

## First, compute principal coordinate decomposition (classical scaling) for 
## turnover distance matrix using the pcoa() function from the ape package. Use 
## lingoes correction to account for negative eigenvalues.
pcoa.jtu <- pcoa(comb.fam.dist$beta.jtu, correction="lingoes", rn=NULL)

## Inspect outputs
pcoa.jtu$correction
pcoa.jtu$note
pcoa.jtu$trace

## The principal component coordinates for each diversity matrix can be written
## to a csv file for inspection
#write.csv(pcoa.jtu$vectors, "../Data/Combined_data_family_turnover_eigenvectors.csv")

## Groups of variables for RDA were already created for species-level analysis
## Perform forward model selection for each group of explanatory variables
## Adjusted R-squared will be used as a goodness-of-fit measure. Model selection
## needs to be performed for each group of variables, excluding those containing
## only one variable.
abiotic.rda.null <- rda(pcoa.jtu$vectors ~ 1, abiotic)  # null model
abiotic.rda.all <- rda(pcoa.jtu$vectors ~ ., abiotic)   # model with all variables
step.env <- ordiR2step(abiotic.rda.null, abiotic.rda.all, 
                       direction="forward", perm.max = 1000)

## Summary table and plot:
step.env$anova
plot(step.env)
## Indicates that no abiotic variables should be included in model.

## Run RDA with biotic variable
## Test for significance of the whole model:
mod2 <- rda(pcoa.jtu$vectors ~ biotic$Crucian)
summary(mod2)
coef(mod2)       # canonical coefficients
vif.cca(mod2)    # no collinearity present in model
mod2             # 8.04% variance is explained by variables in model
anova(mod2)      # Model is not significant
RsquareAdj(mod2) # 2.29% variance explained
plot(mod2)

## Check whether any variables should be dropped from model
drop1(mod2, test="perm")
## Crucian carp presence-absence is borderline significant should perhaps
## be dropped

## Plot model result to get a sense of which variables are correlating with which 
## ponds along which axes
plot(mod2, type='n', scaling=1)
orditorp(mod2, display='sp', cex=0.5, scaling=1, col='blue')
text(mod2, display='cn', col='red')

## Plot of variables against sites
plot(mod2, type = "text")

## Perform hypothesis testing
anova(mod2, perm=1000)               # overall model fit
anova(mod2, by="margin", perm=1000)  # significance of marginal effects (Type III)
anova(mod2, by="term", perm=1000)    # significance of each term in model
anova(mod2, by="axis", perm=1000)    # significance of each axis


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
comb.fam.bd.nest <- betadisper(comb.fam.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
comb.fam.nest <- with(metadata, comb.fam.bd.nest)
comb.fam.nest

## Compute mean distance to centroid per group
tapply(comb.fam.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(comb.fam.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(comb.fam.bd.nest)

## Boxplot of turnover partition
boxplot(comb.fam.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is no difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(comb.fam.bd.nest)     # No significant difference between ponds
permutest(comb.fam.bd.nest) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(comb.fam.bd.nest)  # No significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
comb.fam.comm.nest <- metaMDS(comb.fam.dist$beta.jne, 
                              dist="jaccard", 
                              k=2,
                              maxit=999,
                              trymax=1000,
                              wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
comb.fam.comm.nest$stress
stressplot(comb.fam.comm.nest)

## plot site scores as text
ordiplot(comb.fam.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
comb.NMDS1 <- comb.fam.comm.nest$points[,1]
comb.NMDS2 <- comb.fam.comm.nest$points[,2]
comb.fam.nest.NMDS <- data.frame(NMDS1=comb.NMDS1, 
                                 NMDS2=comb.NMDS2,
                                 Crucian = metadata$Crucian)

## Check data
head(comb.fam.nest.NMDS)

## Plot data frame
p18b <- ggplot(comb.fam.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="NMDS1", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p18b

## Statistically check difference in nestedness of communities
comb.fam.nest.anosim <- anosim(comb.fam.dist$beta.jne, metadata$Crucian)

## Inspect results
comb.fam.nest.anosim
summary(comb.fam.nest.anosim)
plot(comb.fam.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
comb.fam.nest.adonis <- adonis(comb.fam.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
comb.fam.nest.adonis

## Result is not significant. The combined data indicates there is no
## substantial difference in species loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## Use RDA and variance partitioning to examine the influence of biotic and 
## abiotic variables on turnover partition of beta diversity

## First, compute principal coordinate decomposition (classical scaling) for 
## turnover distance matrix using the pcoa() function from the ape package. Use 
## lingoes correction to account for negative eigenvalues.
pcoa.jne <- pcoa(comb.fam.dist$beta.jne, correction="lingoes", rn=NULL)

## Inspect outputs
pcoa.jne$correction
pcoa.jne$note
pcoa.jne$trace

## The principal component coordinates for each diversity matrix can be written
## to a csv file for inspection
#write.csv(pcoa.jne$vectors, "../Data/Combined_data_family_nestedness_eigenvectors.csv")

## Groups of variables for RDA were already created for species-level analysis
## Perform forward model selection for each group of explanatory variables
## Adjusted R-squared will be used as a goodness-of-fit measure. Model selection
## needs to be performed for each group of variables, excluding those containing
## only one variable.
abiotic.rda.null <- rda(pcoa.jne$vectors ~ 1, abiotic)  # null model
abiotic.rda.all <- rda(pcoa.jne$vectors ~ ., abiotic)   # model with all variables
step.env <- ordiR2step(abiotic.rda.null, abiotic.rda.all, 
                       direction="forward", perm.max = 1000)

## Summary table and plot:
step.env$anova
plot(step.env)
## Indicates that no abiotic variables should be included in model.

## Run RDA with biotic variable
## Test for significance of the whole model:
mod2 <- rda(pcoa.jne$vectors ~ biotic$Crucian)
summary(mod2)
coef(mod2)       # canonical coefficients
vif.cca(mod2)    # no collinearity present in model
mod2             # 10.9% variance is explained by variables in model
anova(mod2)      # Model is not significant
RsquareAdj(mod2) # 5.4% variance explained
plot(mod2)

## Check whether any variables should be dropped from model
drop1(mod2, test="perm")
## Crucian carp presence-absence should potentially be dropped

## Plot model result to get a sense of which variables are correlating with which 
## ponds along which axes
plot(mod2, type='n', scaling=1)
orditorp(mod2, display='sp', cex=0.5, scaling=1, col='blue')
text(mod2, display='cn', col='red')

## Plot of variables against sites
plot(mod2, type = "text")

## Perform hypothesis testing
anova(mod2, perm=1000)               # overall model fit
anova(mod2, by="margin", perm=1000)  # significance of marginal effects (Type III)
anova(mod2, by="term", perm=1000)    # significance of each term in model
anova(mod2, by="axis", perm=1000)    # significance of each axis


## 3. TOTAL BETA DIVERSITY
comb.fam.bd.total <- betadisper(comb.fam.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(metadata, comb.fam.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(comb.fam.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(comb.fam.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(comb.fam.bd.total)

## Boxplot of total beta diversity
boxplot(comb.fam.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is substantial difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(comb.fam.bd.total)     # Significant difference between ponds
permutest(comb.fam.bd.total) # Significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(comb.fam.bd.total)  # Significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
comb.fam.comm.total <- metaMDS(comb.fam.dist$beta.jac, 
                               dist="jaccard", 
                               k=2,
                               maxit=999,
                               trymax=1000,
                               wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
comb.fam.comm.total$stress
stressplot(comb.fam.comm.total)

## plot site scores as text
ordiplot(comb.fam.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
comb.NMDS1 <- comb.fam.comm.total$points[,1]
comb.NMDS2 <- comb.fam.comm.total$points[,2]
comb.fam.total.NMDS <- data.frame(NMDS1=comb.NMDS1,
                                  NMDS2=comb.NMDS2,
                                  Crucian = metadata$Crucian)

## Check data
head(comb.fam.total.NMDS)

## Plot data frame
p18c <- ggplot(comb.fam.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="NMDS1", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p18c

## Statistically check difference in total beta diversity of ponds
comb.fam.total.anosim <- anosim(comb.fam.dist$beta.jac, metadata$Crucian)

## Inspect results
comb.fam.total.anosim
summary(comb.fam.total.anosim)
plot(comb.fam.total.anosim)

## There appears to be a significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
comb.fam.total.adonis <- adonis(comb.fam.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
comb.fam.total.adonis

## Again result is significant. The combined data indicates there is 
## substantial variation in overall species composition of
## ponds with crucian carp and ponds without crucian carp


## Use RDA and variance partitioning to examine the influence of biotic and 
## abiotic variables on turnover partition of beta diversity

## First, compute principal coordinate decomposition (classical scaling) for 
## turnover distance matrix using the pcoa() function from the ape package. Use 
## lingoes correction to account for negative eigenvalues.
pcoa.jac <- pcoa(comb.fam.dist$beta.jac, correction="lingoes", rn=NULL)

## Inspect outputs
pcoa.jac$correction
pcoa.jac$note
pcoa.jac$trace

## The principal component coordinates for each diversity matrix can be written
## to a csv file for inspection
#write.csv(pcoa.jac$vectors, "../Data/Combined_data_family_totalbeta_eigenvectors.csv")

## Groups of variables for RDA were already created for species-level analysis
## Perform forward model selection for each group of explanatory variables
## Adjusted R-squared will be used as a goodness-of-fit measure. Model selection
## needs to be performed for each group of variables, excluding those containing
## only one variable.
abiotic.rda.null <- rda(pcoa.jac$vectors ~ 1, abiotic)  # null model
abiotic.rda.all <- rda(pcoa.jac$vectors ~ ., abiotic)   # model with all variables
step.env <- ordiR2step(abiotic.rda.null, abiotic.rda.all, 
                       direction="forward", perm.max = 1000)

## Summary table and plot:
step.env$anova
plot(step.env)
## Indicates that no abiotic variables should be included in model.

## Run RDA with biotic variable
## Test for significance of the whole model:
mod2 <- rda(pcoa.jac$vectors ~ biotic$Crucian)
summary(mod2)
coef(mod2)       # canonical coefficients
vif.cca(mod2)    # no collinearity present in model
mod2             # 8.43% variance is explained by variables in model
anova(mod2)      # Model is significant
RsquareAdj(mod2) # 2.71% variance explained
plot(mod2)

## Check whether any variables should be dropped from model
drop1(mod2, test="perm")
## Crucian carp presence-absence should not be dropped as it is significant

## Plot model result to get a sense of which variables are correlating with which 
## ponds along which axes
plot(mod2, type='n', scaling=1)
orditorp(mod2, display='sp', cex=0.5, scaling=1, col='blue')
text(mod2, display='cn', col='red')

## Plot of variables against sites
plot(mod2, type = "text")

## Perform hypothesis testing
anova(mod2, perm=1000)               # overall model fit
anova(mod2, by="margin", perm=1000)  # significance of marginal effects (Type III)
anova(mod2, by="term", perm=1000)    # significance of each term in model
anova(mod2, by="axis", perm=1000)    # significance of each axis



######################################################
# INFLUENCE OF CRUCIAN CARP ON POND INVERTEBRATES    #
# (MEIOFAUNA AND TAXA WITH NO REF SEQUENCES REMOVED) #
######################################################

################################
# SWEEP-NETTING AND MICROSCOPY #
################################

#------------------------#
# SPECIES-LEVEL ANALYSIS #
#------------------------#

## BASIC SUMMARIES
## Total number of individuals per pond:
sum.of.rows <- apply(net.noref, 1, sum)
sort(sum.of.rows, decreasing = TRUE) 
sum(sum.of.rows)

## Total number of individuals per species across all sites:
sum.of.columns <- apply(net.noref, 2, sum)

## Actual number of individuals:
sort(sum.of.columns, decreasing = TRUE)
sum(sum.of.columns)

## Proportional number of individuals:
sort(round(sum.of.columns/sum(sum.of.columns)*100, 2), decreasing = TRUE)

## Number of ponds where each species occurs:
net.spec.pres <- apply(net.noref > 0, 2, sum) 
sort(net.spec.pres, decreasing = TRUE)


## ALPHA DIVERSITY
## Species richness
net.noref.sample.richness <- specnumber(net.noref.spp.pa)
net.noref.sample.richness

## Statistically compare and plot alpha diversity
## Create data frame
net.noref.alpha <- data.frame(net.noref.sample.richness)

## add metadata from external file
net.noref.alpha <- cbind(metadata[,1:2], net.noref.alpha)

## Reset row names of data frame for further indexing
rownames(net.noref.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
net.noref.alpha.regression <- glm.nb(net.noref.sample.richness ~ Crucian, data=net.noref.alpha)
summary(net.noref.alpha.regression)
anova(net.noref.alpha.regression, test = "Chi")
drop1(net.noref.alpha.regression, test = "Chi")
summary(glht(glm.nb(net.noref.sample.richness ~ Crucian, data=net.noref.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
18.052/16
1-pchisq(18.052, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(net.noref.alpha$net.noref.sample.richness ~ fitted(net.noref.alpha.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(net.noref.alpha.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(net.alpha.regression)
shapiro.test(sresid) # P = 0.3798

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ net.noref.alpha.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ net.noref.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(net.noref.alpha.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(net.noref.alpha.regression)
summary(influence)
CD <- cooks.distance(net.noref.alpha.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot species richness
p19 <- ggplot(net.noref.alpha, aes(x=Crucian, y=net.noref.sample.richness)) + 
  geom_jitter(aes(colour=Crucian), width=0.2, cex=2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,30)) + 
  labs(title="(a) Species-level",
       subtitle="(i) Netting and microscopy",
       x="", y=expression(alpha~diversity)) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p19

## Species richness does not appear to be higher in ponds without crucian
## carp, indicating this species may have a negligible impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot species accumulation curve across ponds
plot(specaccum(net.noref.spp.pa), 
     xlab = "Number of ponds", 
     ylab = "Number of species",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## BETA DIVERSITY
net.noref.multi <- beta.multi(net.noref.spp.pa, index.family="jaccard")
print(net.noref.multi)

## The majority of total beta diversity arises from species turnover
## rather than nestedness.
## Pairwise between-site values of each component of beta diversity
net.noref.dist <- beta.pair(net.noref.spp.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
net.noref.bd.turn <- betadisper(net.noref.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(metadata, net.noref.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(net.noref.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(net.noref.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(net.noref.bd.turn)

## Boxplot of turnover partition
boxplot(net.noref.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is some difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(net.noref.bd.turn)     # No significant difference between ponds
permutest(net.noref.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(net.noref.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
net.noref.comm.turn <- metaMDS(net.noref.dist$beta.jtu,
                               dist="jaccard", 
                               k=2,
                               maxit=999,
                               trymax=1000,
                               wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
net.noref.comm.turn$stress
stressplot(net.noref.comm.turn)

## plot site scores as text
ordiplot(net.noref.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
net.noref.NMDS1 <- net.noref.comm.turn$points[,1]
net.noref.NMDS2 <- net.noref.comm.turn$points[,2]
net.noref.turn.NMDS <- data.frame(NMDS1=net.noref.NMDS1, 
                                  NMDS2=net.noref.NMDS2,
                                  Crucian = metadata$Crucian)

## Check data
head(net.noref.turn.NMDS)

## Plot data frame
p20a <- ggplot(net.noref.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(x="", y="NMDS2",
       title="(a) Turnover",
       subtitle=" (i) Netting and microscopy") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position = "bottom",
        legend.key=element_blank(),
        legend.key.size = unit(2, 'lines'))
p20a

## Statistically check difference in spatial turnover of communities
net.noref.turn.anosim <- anosim(net.noref.dist$beta.jtu, metadata$Crucian)

## Inspect results
net.noref.turn.anosim
summary(net.noref.turn.anosim)
plot(net.noref.turn.anosim)

## There appears to be no significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
net.noref.turn.adonis <- adonis(net.noref.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
net.noref.turn.adonis

## Result is not significant. Netting indicates there is no substantial 
## difference in species replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## species in ponds with no fish are not substituted by species in ponds
## with crucian carp.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
net.noref.bd.nest <- betadisper(net.noref.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.nest <- with(metadata, net.noref.bd.nest)
mod.nest

## Compute mean distance to centroid per group
tapply(net.noref.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(net.noref.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of nestedness partition
plot(net.noref.bd.nest)

## Boxplot of nestedness partition
boxplot(net.noref.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is slight difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(net.noref.bd.nest)     # Significant difference between ponds
permutest(net.noref.bd.nest) # Significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(net.noref.bd.nest)  # Significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
net.noref.comm.nest <- metaMDS(net.noref.dist$beta.jne, 
                               dist="jaccard", 
                               k=2,
                               maxit=999,
                               trymax=1000,
                               wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
net.noref.comm.nest$stress
stressplot(net.noref.comm.nest)

## plot site scores as text
ordiplot(net.noref.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
net.noref.NMDS1 <- net.noref.comm.nest$points[,1]
net.noref.NMDS2 <- net.noref.comm.nest$points[,2]
net.noref.nest.NMDS <- data.frame(NMDS1=net.noref.NMDS1, 
                                  NMDS2=net.noref.NMDS2,
                                  Crucian = metadata$Crucian)

## Check data
head(net.noref.nest.NMDS)

## Plot data frame
p20b <- ggplot(net.noref.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(title="(b) Nestedness", 
       subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p20b

## Statistically check difference in nestedness of communities
net.noref.nest.anosim <- anosim(net.noref.dist$beta.jne, metadata$Crucian)

## Inspect results
net.noref.nest.anosim
summary(net.noref.nest.anosim)
plot(net.noref.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
net.noref.nest.adonis <- adonis(net.noref.dist$beta.jne ~ Crucian, metadata)
net.noref.nest.adonis

## Again result is not significant. Netting indicates there is no 
## substantial difference in species loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## 3. TOTAL BETA DIVERSITY
net.noref.bd.total <- betadisper(net.noref.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(metadata, net.noref.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(net.noref.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(net.noref.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(net.noref.bd.total)

## Boxplot of total beta diversity
boxplot(net.noref.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is little difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(net.noref.bd.total)     # No significant difference between ponds
permutest(net.noref.bd.total) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(net.noref.bd.total)  # No significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
net.noref.comm.total <- metaMDS(net.noref.dist$beta.jac, 
                                dist="jaccard", 
                                k=2,
                                maxit=999,
                                trymax=1000,
                                wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
net.noref.comm.total$stress
stressplot(net.noref.comm.total)

## plot site scores as text
ordiplot(net.noref.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
net.noref.NMDS1 <- net.noref.comm.total$points[,1]
net.noref.NMDS2 <- net.noref.comm.total$points[,2]
net.noref.total.NMDS <- data.frame(NMDS1=net.noref.NMDS1,
                                   NMDS2=net.noref.NMDS2,
                                   Crucian = metadata$Crucian)

## Check data
head(net.noref.total.NMDS)

## Plot data frame
p20c <- ggplot(net.noref.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(title=expression(bold("(c) Total"~beta~"Diversity")), 
       subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p20c

## Statistically check difference in total beta diversity of ponds
net.noref.total.anosim <- anosim(net.noref.dist$beta.jac, metadata$Crucian)

## Inspect results
net.noref.total.anosim
summary(net.noref.total.anosim)
plot(net.noref.total.anosim)

## There appears to be borderline significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
net.noref.total.adonis <- adonis(net.noref.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
net.noref.total.adonis

## Result is borderline significant. Netting indicates there is no 
## substantial variation in overall species composition of ponds with 
## crucian carp and ponds without crucian carp


#-----------------------#
# FAMILY-LEVEL ANALYSIS #
#-----------------------#

## BASIC SUMMARIES
## Total number of individuals per pond:
sum.of.rows <- apply(net.fam.noref, 1, sum)
sort(sum.of.rows, decreasing = TRUE) 
sum(sum.of.rows)

## Total number of individuals per family across all sites:
sum.of.columns <- apply(net.fam.noref, 2, sum)

## Actual number of individuals:
sort(sum.of.columns, decreasing = TRUE)
sum(sum.of.columns)

## Proportional number of individuals:
sort(round(sum.of.columns/sum(sum.of.columns)*100, 2), decreasing = TRUE)

## Number of ponds where each species occurs:
net.fam.noref.pres <- apply(net.fam.noref > 0, 2, sum) 
sort(net.fam.noref.pres, decreasing = TRUE)


## ALPHA DIVERSITY
## Family richness
net.noref.fam.richness <- specnumber(net.noref.fam.pa)
net.noref.fam.richness

## Statistically compare and plot alpha diversity
## Create data frame
net.noref.fam.alpha <- data.frame(net.noref.fam.richness)

## add metadata from external file
net.noref.fam.alpha <- cbind(metadata[,1:2], net.noref.fam.alpha)

## Reset row names of data frame for further indexing
rownames(net.noref.fam.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
net.noref.fam.regression <- glm.nb(net.noref.fam.richness ~ Crucian, data=net.noref.fam.alpha)
summary(net.noref.fam.regression)
anova(net.noref.fam.regression, test = "Chi")
drop1(net.noref.fam.regression, test = "Chi")
summary(glht(glm.nb(net.noref.fam.richness ~ Crucian, data=net.noref.fam.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
18.390/16
1-pchisq(18.390, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(net.noref.fam.alpha$net.noref.fam.richness ~ fitted(net.noref.fam.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(net.noref.fam.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(net.noref.fam.regression)
shapiro.test(sresid) # P = 0.5228

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ net.noref.fam.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ net.noref.fam.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(net.noref.fam.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(net.noref.fam.regression)
summary(influence)
CD <- cooks.distance(net.noref.fam.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot family richness
p21 <- ggplot(net.noref.fam.alpha, aes(x=Crucian, y=net.noref.fam.richness)) + 
  geom_jitter(aes(colour=Crucian), width=0.2, cex=2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,30)) + 
  labs(title="(b) Family-level", 
       subtitle="", 
       x="", y="") + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p21

## Family richness does not appear much higher in ponds without crucian
## carp, indicating this species may have a negligible impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot family accumulation curve across ponds
plot(specaccum(net.noref.fam.pa), 
     xlab = "Number of ponds", 
     ylab = "Number of families",
     ci.type="polygon", ci.col="grey")

## Number of ponds samples was adequate to fully represent invertebrate
## diversity at family-level.


## BETA DIVERSITY
net.noref.fam.multi <- beta.multi(net.noref.fam.pa, index.family="jaccard")
print(net.noref.fam.multi)

## The majority of total beta diversity arises from family turnover
## rather than nestedness.
## Pairwise between-site values of each component of beta diversity
net.noref.fam.dist <- beta.pair(net.noref.fam.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
net.noref.fam.bd.turn <- betadisper(net.noref.fam.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(metadata, net.noref.fam.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(net.noref.fam.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(net.noref.fam.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(net.noref.fam.bd.turn)

## Boxplot of turnover partition
boxplot(net.noref.fam.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is some difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(net.noref.fam.bd.turn)     # No significant difference between ponds
permutest(net.noref.fam.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(net.noref.fam.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
net.noref.fam.comm.turn <- metaMDS(net.noref.fam.dist$beta.jtu,
                                   dist="jaccard", 
                                   k=2,
                                   maxit=999,
                                   trymax=1000,
                                   wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
net.noref.fam.comm.turn$stress
stressplot(net.noref.fam.comm.turn)

## plot site scores as text
ordiplot(net.noref.fam.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
net.noref.NMDS1 <- net.noref.fam.comm.turn$points[,1]
net.noref.NMDS2 <- net.noref.fam.comm.turn$points[,2]
net.noref.fam.turn.NMDS <- data.frame(NMDS1=net.noref.NMDS1, 
                                      NMDS2=net.noref.NMDS2,
                                      Crucian = metadata$Crucian)

## Check data
head(net.noref.fam.turn.NMDS)

## Plot data frame
p22a <- ggplot(net.noref.fam.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(title="(a) Turnover",
       subtitle="(i) Netting and microscopy",
       x="", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position = "bottom",
        legend.key=element_blank(),
        legend.key.size = unit(2, 'lines'))
p22a

## Statistically check difference in spatial turnover of communities
net.noref.fam.turn.anosim <- anosim(net.noref.fam.dist$beta.jtu, metadata$Crucian)

## Inspect results
net.noref.fam.turn.anosim
summary(net.noref.fam.turn.anosim)
plot(net.noref.fam.turn.anosim)

## There appears to be no significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
net.noref.fam.turn.adonis <- adonis(net.noref.fam.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
net.noref.fam.turn.adonis

## Results are not significant. Netting indicates there is no real 
## difference in family replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## families in ponds with no fish are not substituted by families in ponds
## with crucian carp.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
net.noref.fam.bd.nest <- betadisper(net.noref.fam.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.nest <- with(metadata, net.noref.fam.bd.nest)
mod.nest

## Compute mean distance to centroid per group
tapply(net.noref.fam.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(net.noref.fam.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(net.noref.fam.bd.nest)

## Boxplot of turnover partition
boxplot(net.noref.fam.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is slight difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(net.noref.fam.bd.nest)     # No significant difference between ponds
permutest(net.noref.fam.bd.nest) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(net.noref.fam.bd.nest)  # No significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
net.noref.fam.comm.nest <- metaMDS(net.noref.fam.dist$beta.jne, 
                                   dist="jaccard", 
                                   k=2,
                                   maxit=999,
                                   trymax=1000,
                                   wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
net.noref.fam.comm.nest$stress
stressplot(net.noref.fam.comm.nest)

## plot site scores as text
ordiplot(net.noref.fam.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
net.noref.NMDS1 <- net.noref.fam.comm.nest$points[,1]
net.noref.NMDS2 <- net.noref.fam.comm.nest$points[,2]
net.noref.fam.nest.NMDS <- data.frame(NMDS1=net.noref.NMDS1, 
                                      NMDS2=net.noref.NMDS2,
                                      Crucian = metadata$Crucian)

## Check data
head(net.noref.fam.nest.NMDS)

## Plot data frame
p22b <- ggplot(net.noref.fam.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(title="(b) Nestedness", 
       subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p22b

## Statistically check difference in nestedness of communities
net.noref.fam.nest.anosim <- anosim(net.noref.fam.dist$beta.jne, metadata$Crucian)

## Inspect results
net.noref.fam.nest.anosim
summary(net.noref.fam.nest.anosim)
plot(net.noref.fam.nest.anosim)

## There appears to be a significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
net.noref.fam.nest.adonis <- adonis(net.noref.fam.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
net.noref.fam.nest.adonis

## Result is not significant. Netting indicates there is no 
## substantial difference in family loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## 3. TOTAL BETA DIVERSITY
net.noref.fam.bd.total <- betadisper(net.noref.fam.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(metadata, net.noref.fam.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(net.noref.fam.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(net.noref.fam.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(net.noref.fam.bd.total)

## Boxplot of total beta diversity
boxplot(net.noref.fam.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is little difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(net.noref.fam.bd.total)     # No significant difference between ponds
permutest(net.noref.fam.bd.total) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(net.noref.fam.bd.total)  # No significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
net.noref.fam.comm.total <- metaMDS(net.noref.fam.dist$beta.jac, 
                                    dist="jaccard", 
                                    k=2,
                                    maxit=999,
                                    trymax=1000,
                                    wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
net.noref.fam.comm.total$stress
stressplot(net.noref.fam.comm.total)

## plot site scores as text
ordiplot(net.noref.fam.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
net.noref.NMDS1 <- net.noref.fam.comm.total$points[,1]
net.noref.NMDS2 <- net.noref.fam.comm.total$points[,2]
net.noref.fam.total.NMDS <- data.frame(NMDS1=net.noref.NMDS1,
                                       NMDS2=net.noref.NMDS2,
                                       Crucian = metadata$Crucian)

## Check data
head(net.noref.fam.total.NMDS)

## Plot data frame
p22c <- ggplot(net.noref.fam.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(title=expression(bold("(c) Total"~beta~"Diversity")), 
       subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p22c

## Statistically check difference in total beta diversity of ponds
net.noref.fam.total.anosim <- anosim(net.noref.fam.dist$beta.jac, metadata$Crucian)

## Inspect results
net.noref.fam.total.anosim
summary(net.noref.fam.total.anosim)
plot(net.noref.fam.total.anosim)

## There appears to be no significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
net.noref.fam.total.adonis <- adonis(net.noref.fam.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
net.noref.fam.total.adonis

## Result is not significant. Netting indicates there is no variation in 
## overall family composition of ponds with crucian carp and ponds without 
## crucian carp.



#####################
# DNA METABARCODING #
#####################

#------------------------#
# SPECIES-LEVEL ANALYSIS #
#------------------------#

## BASIC SUMMARIES
## Total number of sequences per pond:
sum.of.rows <- apply(DNA.nomicro, 1, sum)
sort(sum.of.rows, decreasing = TRUE) 
sum(sum.of.rows)

## Total number of sequences per species across all sites
sum.of.columns <- apply(DNA.nomicro, 2, sum)

## Actual read counts:
sort(sum.of.columns, decreasing = TRUE)
sum(sum.of.columns)

## Proportional read counts:
sort(round(sum.of.columns/sum(sum.of.columns)*100, 2), decreasing = TRUE)

## Number of ponds where each species occurs
DNA.nomicro.spec.pres <- apply(DNA.nomicro > 0, 2, sum) 
sort(DNA.nomicro.spec.pres, decreasing = TRUE)


## ALPHA DIVERSITY
## Convert read counts to presence-absence
DNA.nomicro.spp.pa <- decostand(DNA.nomicro, method = "pa", MARGIN = 1)
head(DNA.nomicro.spp.pa[,1:3], n = 3)

## Basic richness
DNA.nomicro.spp.richness <- specnumber(DNA.nomicro.spp.pa)
DNA.nomicro.spp.richness

## Statistically compare and plot alpha diversity
## Create data frame
DNA.nomicro.spp.alpha <- data.frame(DNA.nomicro.spp.richness)

## add metadata from external file
DNA.nomicro.spp.alpha <- cbind(metadata[,1:2], DNA.nomicro.spp.alpha)

## Reset row names of data frame for further indexing
rownames(DNA.nomicro.spp.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
DNA.nomicro.spp.regression <- glm.nb(DNA.nomicro.spp.richness ~ Crucian, data=DNA.nomicro.spp.alpha)
summary(DNA.nomicro.spp.regression)
anova(DNA.nomicro.spp.regression, test = "Chi")
drop1(DNA.nomicro.spp.regression, test = "Chi")
summary(glht(glm(DNA.nomicro.spp.richness ~ Crucian, data=DNA.nomicro.spp.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
18.503/16
1-pchisq(18.503, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(DNA.nomicro.spp.alpha$DNA.nomicro.spp.richness ~ fitted(DNA.nomicro.spp.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(DNA.nomicro.spp.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(DNA.spp.regression)
shapiro.test(sresid) # P = 0.3337

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ DNA.nomicro.spp.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ DNA.nomicro.spp.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(DNA.nomicro.spp.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(DNA.nomicro.spp.regression)
summary(influence)
CD <- cooks.distance(DNA.nomicro.spp.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot species richness
p23 <- ggplot(DNA.nomicro.spp.alpha, aes(x=Crucian, y=DNA.nomicro.spp.richness)) + 
  geom_jitter(aes(colour=Crucian), cex=2, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  scale_y_continuous(limits=c(0,40)) + 
  labs(subtitle="(ii) DNA metabarcoding",
       x="", y=expression(alpha~diversity)) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"), 
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p23

## Species richness does not appear higher in ponds without crucian
## carp, indicating this species may have a negligible impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot species accumulation curve across ponds
plot(specaccum(DNA.nomicro.spp.pa), 
     xlab = "Number of ponds", 
     ylab = "Number of species",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## BETA DIVERSITY
DNA.nomicro.multi <- beta.multi(DNA.nomicro.spp.pa, index.family="jaccard")
print(DNA.nomicro.multi)

## The majority of total beta diversity arises from species turnover
## rather than nestedness.
## Pairwise between-site values of each component of beta diversity
DNA.nomicro.dist <- beta.pair(DNA.nomicro.spp.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
DNA.nomicro.bd.turn <- betadisper(DNA.nomicro.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(metadata, DNA.nomicro.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(DNA.nomicro.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(DNA.nomicro.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(DNA.nomicro.bd.turn)

## Boxplot of turnover partition
boxplot(DNA.nomicro.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is no real difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(DNA.nomicro.bd.turn)     # No significant difference between ponds
permutest(DNA.nomicro.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(DNA.nomicro.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
DNA.nomicro.comm.turn <- metaMDS(DNA.nomicro.dist$beta.jtu, 
                                 dist="jaccard", 
                                 k=2,
                                 maxit=999,
                                 trymax=1000,
                                 wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
DNA.nomicro.comm.turn$stress
stressplot(DNA.nomicro.comm.turn)

## plot site scores as text
ordiplot(DNA.nomicro.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
DNA.nomicro.NMDS1 <- DNA.nomicro.comm.turn$points[,1]
DNA.nomicro.NMDS2 <- DNA.nomicro.comm.turn$points[,2]
DNA.nomicro.turn.NMDS <- data.frame(NMDS1=DNA.nomicro.NMDS1, 
                                    NMDS2=DNA.nomicro.NMDS2,
                                    Crucian = metadata$Crucian)

## Check data
head(DNA.nomicro.turn.NMDS)

## Plot data frame
p24a <- ggplot(DNA.nomicro.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="(ii) DNA metabarcoding", 
       x="", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, color="black"),
        plot.subtitle = element_text(face="bold", hjust=0, color="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p24a

## Statistically check difference in spatial turnover of communities
DNA.nomicro.turn.anosim <- anosim(DNA.nomicro.dist$beta.jtu, metadata$Crucian)

## Inspect results
DNA.nomicro.turn.anosim
summary(DNA.nomicro.turn.anosim)
plot(DNA.nomicro.turn.anosim)

## There appears to be no significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
DNA.nomicro.turn.adonis <- adonis(DNA.nomicro.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
DNA.nomicro.turn.adonis

## Again result is not significant. DNA metabarcoding indicates there is no 
## substantial difference in species replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## species in ponds with no fish are not substituted by species in ponds
## with crucian carp.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
DNA.nomicro.bd.nest <- betadisper(DNA.nomicro.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
DNA.nomicro.nest <- with(metadata, DNA.nomicro.bd.nest)
DNA.nomicro.nest

## Compute mean distance to centroid per group
tapply(DNA.nomicro.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(DNA.nomicro.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(DNA.nomicro.bd.nest)

## Boxplot of turnover partition
boxplot(DNA.nomicro.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is slight difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(DNA.nomicro.bd.nest)     # No significant difference between ponds
permutest(DNA.nomicro.bd.nest) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(DNA.nomicro.bd.nest)  # No significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
DNA.nomicro.comm.nest <- metaMDS(DNA.nomicro.dist$beta.jne, 
                                 dist="jaccard", 
                                 k=2,
                                 maxit=999,
                                 trymax=1000,
                                 wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
DNA.nomicro.comm.nest$stress
stressplot(DNA.nomicro.comm.nest)

## plot site scores as text
ordiplot(DNA.nomicro.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
DNA.nomicro.NMDS1 <- DNA.nomicro.comm.nest$points[,1]
DNA.nomicro.NMDS2 <- DNA.nomicro.comm.nest$points[,2]
DNA.nomicro.nest.NMDS <- data.frame(NMDS1=DNA.nomicro.NMDS1, 
                                    NMDS2=DNA.nomicro.NMDS2,
                                    Crucian = metadata$Crucian)

## Check data
head(DNA.nomicro.nest.NMDS)

## Plot data frame
p24b <- ggplot(DNA.nomicro.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p24b

## Statistically check difference in nestedness of communities
DNA.nomicro.nest.anosim <- anosim(DNA.nomicro.dist$beta.jne, metadata$Crucian)

## Inspect results
DNA.nomicro.nest.anosim
summary(DNA.nomicro.nest.anosim)
plot(DNA.nomicro.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
DNA.nomicro.nest.adonis <- adonis(DNA.nomicro.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
DNA.nomicro.nest.adonis

## Again result is not significant. DNA metabarcoding indicates there is no 
## substantial difference in species loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## 3. TOTAL BETA DIVERSITY
DNA.nomicro.bd.total <- betadisper(DNA.nomicro.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(metadata, DNA.nomicro.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(DNA.nomicro.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(DNA.nomicro.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(DNA.nomicro.bd.total)

## Boxplot of total beta diversity
boxplot(DNA.nomicro.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is some difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(DNA.nomicro.bd.total)     # No significant difference between ponds
permutest(DNA.nomicro.bd.total) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(DNA.nomicro.bd.total)  # No significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
DNA.nomicro.comm.total <- metaMDS(DNA.nomicro.dist$beta.jac, 
                                  dist="jaccard", 
                                  k=2,
                                  maxit=999,
                                  trymax=1000,
                                  wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
DNA.nomicro.comm.total$stress
stressplot(DNA.nomicro.comm.total)

## plot site scores as text
ordiplot(DNA.nomicro.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
DNA.nomicro.NMDS1 <- DNA.nomicro.comm.total$points[,1]
DNA.nomicro.NMDS2 <- DNA.nomicro.comm.total$points[,2]
DNA.nomicro.total.NMDS <- data.frame(NMDS1=DNA.nomicro.NMDS1,
                                     NMDS2=DNA.nomicro.NMDS2,
                                     Crucian = metadata$Crucian)

## Check data
head(DNA.nomicro.total.NMDS)

## Plot data frame
p24c <- ggplot(DNA.nomicro.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p24c

## Statistically check difference in total beta diversity of ponds
DNA.nomicro.total.anosim <- anosim(DNA.nomicro.dist$beta.jac, metadata$Crucian)

## Inspect results
DNA.nomicro.total.anosim
summary(DNA.nomicro.total.anosim)
plot(DNA.nomicro.total.anosim)

## There appears to be no significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
DNA.nomicro.total.adonis <- adonis(DNA.nomicro.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
DNA.nomicro.total.adonis


#-----------------------#
# FAMILY-LEVEL ANALYSIS #
#-----------------------#

## BASIC SUMMARIES 
## Total number of sequences per pond:
sum.of.rows <- apply(DNA.fam.nomicro, 1, sum)
sort(sum.of.rows, decreasing = TRUE) 
sum(sum.of.rows)

## Total number of sequences per species across all sites:
sum.of.columns <- apply(DNA.fam.nomicro, 2, sum)

## Actual read counts:
sort(sum.of.columns, decreasing = TRUE)
sum(sum.of.columns)

## Proportional read counts:
sort(round(sum.of.columns/sum(sum.of.columns)*100, 2), decreasing = TRUE)

## Number of ponds where each species occurs
DNA.fam.nomicro.pres <- apply(DNA.fam.nomicro > 0, 2, sum) 
sort(DNA.fam.nomicro.pres, decreasing = TRUE)


## ALPHA DIVERSITY
## Convert read counts to presence-absence
DNA.nomicro.fam.pa <- decostand(DNA.fam.nomicro, method = "pa", MARGIN = 1)
head(DNA.nomicro.fam.pa[,1:3], n = 3)

## Basic richness
DNA.nomicro.fam.richness <- specnumber(DNA.nomicro.fam.pa)
DNA.nomicro.fam.richness

## Statistically compare and plot alpha diversity
## Create data frame
DNA.nomicro.fam.alpha <- data.frame(DNA.nomicro.fam.richness)

## Add sample metadata
DNA.nomicro.fam.alpha <- cbind(metadata[,1:2], DNA.nomicro.fam.alpha)

## Reset row names of data frame for further indexing
rownames(DNA.nomicro.fam.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
DNA.nomicro.fam.regression <- glm.nb(DNA.nomicro.fam.richness ~ Crucian, data=DNA.nomicro.fam.alpha)
summary(DNA.nomicro.fam.regression)
anova(DNA.nomicro.fam.regression, test = "Chi")
drop1(DNA.nomicro.fam.regression, test = "Chi")
summary(glht(glm(DNA.nomicro.fam.richness ~ Crucian, data=DNA.nomicro.fam.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
15.673/16
1-pchisq(15.673, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(DNA.nomicro.fam.alpha$DNA.nomicro.fam.richness ~ fitted(DNA.nomicro.fam.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(DNA.nomicro.fam.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(DNA.fam.regression)
shapiro.test(sresid) # P = 0.9263

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ DNA.nomicro.fam.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ DNA.nomicro.fam.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(DNA.nomicro.fam.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(DNA.nomicro.fam.regression)
summary(influence)
CD <- cooks.distance(DNA.nomicro.fam.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot family richness
p25 <- ggplot(DNA.nomicro.fam.alpha, aes(x=Crucian, y=DNA.nomicro.fam.richness)) + 
  geom_jitter(aes(colour=Crucian), cex=2, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,40)) + 
  labs(subtitle="", 
       x="", y="") + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p25

## Family richness does not appear higher in ponds without crucian
## carp, indicating this species may have a negligible impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot family accumulation curve across ponds
plot(specaccum(DNA.nomicro.fam.pa), 
     xlab = "Number of ponds", 
     ylab = "Number of species",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## BETA DIVERSITY
DNA.nomicro.fam.multi <- beta.multi(DNA.nomicro.fam.pa, index.family="jaccard")
print(DNA.nomicro.fam.multi)

## The majority of total beta diversity arises from family turnover
## rather than nestedness.
## Pairwise between-site values of each component of beta diversity
DNA.nomicro.fam.dist <- beta.pair(DNA.nomicro.fam.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
DNA.nomicro.fam.bd.turn <- betadisper(DNA.nomicro.fam.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(metadata, DNA.nomicro.fam.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(DNA.nomicro.fam.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(DNA.nomicro.fam.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(DNA.nomicro.fam.bd.turn)

## Boxplot of turnover partition
boxplot(DNA.nomicro.fam.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is no real difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(DNA.nomicro.fam.bd.turn)     # No significant difference between ponds
permutest(DNA.nomicro.fam.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(DNA.nomicro.fam.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
DNA.nomicro.fam.comm.turn <- metaMDS(DNA.nomicro.fam.dist$beta.jtu, 
                                     dist="jaccard", 
                                     k=2,
                                     maxit=999,
                                     trymax=1000,
                                     wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
DNA.nomicro.fam.comm.turn$stress
stressplot(DNA.nomicro.fam.comm.turn)

## plot site scores as text
ordiplot(DNA.nomicro.fam.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
DNA.nomicro.NMDS1 <- DNA.nomicro.fam.comm.turn$points[,1]
DNA.nomicro.NMDS2 <- DNA.nomicro.fam.comm.turn$points[,2]
DNA.nomicro.fam.turn.NMDS <- data.frame(NMDS1=DNA.nomicro.NMDS1, 
                                        NMDS2=DNA.nomicro.NMDS2,
                                        Crucian = metadata$Crucian)

## Check data
head(DNA.nomicro.fam.turn.NMDS)

## Plot data frame
p26a <- ggplot(DNA.nomicro.fam.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="(ii) DNA metabarcoding", 
       x="", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p26a

## Statistically check difference in spatial turnover of communities
DNA.nomicro.fam.turn.anosim <- anosim(DNA.nomicro.fam.dist$beta.jtu, metadata$Crucian)

## Inspect results
DNA.nomicro.fam.turn.anosim
summary(DNA.nomicro.fam.turn.anosim)
plot(DNA.nomicro.fam.turn.anosim)

## There appears to be no significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
DNA.nomicro.fam.turn.adonis <- adonis(DNA.nomicro.fam.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
DNA.nomicro.fam.turn.adonis

## Again result is not significant. DNA metabarcoding indicates there is no
## substantial difference in family replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## families in ponds with no fish are not substituted by families in ponds
## with crucian carp.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
DNA.nomicro.fam.bd.nest <- betadisper(DNA.nomicro.fam.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
DNA.nomicro.fam.nest <- with(metadata, DNA.nomicro.fam.bd.nest)
DNA.nomicro.fam.nest

## Compute mean distance to centroid per group
tapply(DNA.nomicro.fam.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(DNA.nomicro.fam.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(DNA.nomicro.fam.bd.nest)

## Boxplot of turnover partition
boxplot(DNA.nomicro.fam.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is no difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(DNA.nomicro.fam.bd.nest)     # No significant difference between ponds
permutest(DNA.nomicro.fam.bd.nest) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(DNA.nomicro.fam.bd.nest)  # No significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
DNA.nomicro.fam.comm.nest <- metaMDS(DNA.nomicro.fam.dist$beta.jne, 
                                     dist="jaccard", 
                                     k=2,
                                     maxit=999,
                                     trymax=1000,
                                     wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
DNA.nomicro.fam.comm.nest$stress
stressplot(DNA.nomicro.fam.comm.nest)

## plot site scores as text
ordiplot(DNA.nomicro.fam.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
DNA.nomicro.NMDS1 <- DNA.nomicro.fam.comm.nest$points[,1]
DNA.nomicro.NMDS2 <- DNA.nomicro.fam.comm.nest$points[,2]
DNA.nomicro.fam.nest.NMDS <- data.frame(NMDS1=DNA.nomicro.NMDS1, 
                                        NMDS2=DNA.nomicro.NMDS2,
                                        Crucian = metadata$Crucian)

## Check data
head(DNA.nomicro.fam.nest.NMDS)

## Plot data frame
p26b <- ggplot(DNA.nomicro.fam.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p26b

## Statistically check difference in nestedness of communities
DNA.nomicro.fam.nest.anosim <- anosim(DNA.nomicro.fam.dist$beta.jne, metadata$Crucian)

## Inspect results
DNA.nomicro.fam.nest.anosim
summary(DNA.nomicro.fam.nest.anosim)
plot(DNA.nomicro.fam.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
DNA.nomicro.fam.nest.adonis <- adonis(DNA.nomicro.fam.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
DNA.nomicro.fam.nest.adonis

## Again result is not significant. DNA metabarcoding indicates there is no 
## substantial difference in family loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## 3. TOTAL BETA DIVERSITY
DNA.nomicro.fam.bd.total <- betadisper(DNA.nomicro.fam.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(metadata, DNA.nomicro.fam.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(DNA.nomicro.fam.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(DNA.nomicro.fam.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(DNA.nomicro.fam.bd.total)

## Boxplot of total beta diversity
boxplot(DNA.nomicro.fam.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is no difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(DNA.nomicro.fam.bd.total)     # No significant difference between ponds
permutest(DNA.nomicro.fam.bd.total) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(DNA.nomicro.fam.bd.total)  # No significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
DNA.nomicro.fam.comm.total <- metaMDS(DNA.nomicro.fam.dist$beta.jac, 
                                      dist="jaccard", 
                                      k=2,
                                      maxit=999,
                                      trymax=1000,
                                      wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
DNA.nomicro.fam.comm.total$stress
stressplot(DNA.nomicro.fam.comm.total)

## plot site scores as text
ordiplot(DNA.nomicro.fam.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
DNA.nomicro.NMDS1 <- DNA.nomicro.fam.comm.total$points[,1]
DNA.nomicro.NMDS2 <- DNA.nomicro.fam.comm.total$points[,2]
DNA.nomicro.fam.total.NMDS <- data.frame(NMDS1=DNA.nomicro.NMDS1,
                                         NMDS2=DNA.nomicro.NMDS2,
                                         Crucian = metadata$Crucian)

## Check data
head(DNA.nomicro.fam.total.NMDS)

## Plot data frame
p26c <- ggplot(DNA.nomicro.fam.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p26c

## Statistically check difference in total beta diversity of ponds
DNA.nomicro.fam.total.anosim <- anosim(DNA.nomicro.fam.dist$beta.jac, metadata$Crucian)

## Inspect results
DNA.nomicro.fam.total.anosim
summary(DNA.nomicro.fam.total.anosim)
plot(DNA.nomicro.fam.total.anosim)

## There appears to be no significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
DNA.nomicro.fam.total.adonis <- adonis(DNA.nomicro.fam.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
DNA.nomicro.fam.total.adonis

## Again result is not significant. DNA metabarcoding indicates there is no 
## variation in overall species composition of ponds with crucian carp and 
## ponds without crucian carp.



######################
# eDNA METABARCODING #
######################

#------------------------#
# SPECIES-LEVEL ANALYSIS #
#------------------------#

## BASIC SUMMARIES
## Total number of sequences per pond:
sum.of.rows <- apply(eDNA.nomicro, 1, sum)
sort(sum.of.rows, decreasing = TRUE) 
sum(sum.of.rows)

## Total number of sequences per species across all sites:
sum.of.columns <- apply(eDNA.nomicro, 2, sum)

## Actual read counts:
sort(sum.of.columns, decreasing = TRUE)
sum(sum.of.columns)

## Proportional read counts:
sort(round(sum.of.columns/sum(sum.of.columns)*100, 2), decreasing = TRUE)

## Number of ponds where each species occurs:
eDNA.nomicro.spec.pres <- apply(eDNA.nomicro > 0, 2, sum) 
sort(eDNA.nomicro.spec.pres, decreasing = TRUE)


## ALPHA DIVERSITY
## Convert read counts to presence-absence
eDNA.nomicro.spp.pa <- decostand(eDNA.nomicro, method = "pa", MARGIN = 1)
head(eDNA.nomicro.spp.pa[,1:3], n = 3)

## Basic richness
eDNA.nomicro.spp.richness <- specnumber(eDNA.nomicro.spp.pa)
eDNA.nomicro.spp.richness

## Statistically compare and plot alpha diversity
## Create data frame
eDNA.nomicro.spp.alpha <- data.frame(eDNA.nomicro.spp.richness)

## Add sample metadata
eDNA.nomicro.spp.alpha <- cbind(metadata[,1:2], eDNA.nomicro.spp.alpha)

## Reset row names of data frame for further indexing
rownames(eDNA.nomicro.spp.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
eDNA.nomicro.spp.regression <- glm.nb(eDNA.nomicro.spp.richness ~ Crucian, data=eDNA.nomicro.spp.alpha)
summary(eDNA.nomicro.spp.regression)
anova(eDNA.nomicro.spp.regression, test = "Chi")
drop1(eDNA.nomicro.spp.regression, test = "Chi")
summary(glht(glm.nb(eDNA.nomicro.spp.richness ~ Crucian, data=eDNA.nomicro.spp.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
20.014/16
1-pchisq(20.014, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(eDNA.nomicro.spp.alpha$eDNA.nomicro.spp.richness ~ fitted(eDNA.nomicro.spp.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(eDNA.nomicro.spp.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(eDNA.nomicro.spp.regression)
shapiro.test(sresid) # P = 0.9046

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ eDNA.nomicro.spp.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ eDNA.nomicro.spp.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(eDNA.nomicro.spp.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(eDNA.nomicro.spp.regression)
summary(influence)
CD <- cooks.distance(eDNA.nomicro.spp.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot species richness
p27 <- ggplot(eDNA.nomicro.spp.alpha, aes(x=Crucian, y=eDNA.nomicro.spp.richness)) + 
  geom_jitter(aes(colour=Crucian), cex=2, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,30), breaks=seq(0,30,10)) + 
  labs(subtitle="(iii) eDNA metabarcoding", 
       x="", y=expression(alpha~diversity)) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p27

## Species richness appears to be the same in ponds without crucian
## carp, indicating this species may have little impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot species accumulation curve across ponds
plot(specaccum(eDNA.nomicro.spp.pa), 
     xlab = "Number of ponds", 
     ylab = "Number of species",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## BETA DIVERSITY
eDNA.nomicro.multi <- beta.multi(eDNA.nomicro.spp.pa, index.family="jaccard")
print(eDNA.nomicro.multi)

## The majority of total beta diversity arises from species turnover
## rather than nestedness.
## Pairwise between-site values of each component of beta diversity
eDNA.nomicro.dist <- beta.pair(eDNA.nomicro.spp.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
eDNA.nomicro.bd.turn <- betadisper(eDNA.nomicro.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(metadata, eDNA.nomicro.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(eDNA.nomicro.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(eDNA.nomicro.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(eDNA.nomicro.bd.turn)

## Boxplot of turnover partition
boxplot(eDNA.nomicro.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is substantial difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(eDNA.nomicro.bd.turn)     # No significant difference between ponds
permutest(eDNA.nomicro.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(eDNA.nomicro.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
eDNA.nomicro.comm.turn <- metaMDS(eDNA.nomicro.dist$beta.jtu, 
                                  dist="jaccard", 
                                  k=2,
                                  maxit=999,
                                  trymax=1000,
                                  wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
eDNA.nomicro.comm.turn$stress
stressplot(eDNA.nomicro.comm.turn)

## plot site scores as text
ordiplot(eDNA.nomicro.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
eDNA.nomicro.NMDS1 <- eDNA.nomicro.comm.turn$points[,1]
eDNA.nomicro.NMDS2 <- eDNA.nomicro.comm.turn$points[,2]
eDNA.nomicro.turn.NMDS <- data.frame(NMDS1=eDNA.nomicro.NMDS1, 
                                     NMDS2=eDNA.nomicro.NMDS2,
                                     Crucian = metadata$Crucian)

## Check data
head(eDNA.nomicro.turn.NMDS)

## Plot data frame
p28a <- ggplot(eDNA.nomicro.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="(iii) eDNA metabarcoding", 
       x="", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p28a

## Statistically check difference in spatial turnover of communities
eDNA.nomicro.turn.anosim <- anosim(eDNA.nomicro.dist$beta.jtu, metadata$Crucian)

## Inspect results
eDNA.nomicro.turn.anosim
summary(eDNA.nomicro.turn.anosim)
plot(eDNA.nomicro.turn.anosim)

## There appears to be no significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
eDNA.nomicro.turn.adonis <- adonis(eDNA.nomicro.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
eDNA.nomicro.turn.adonis

## Again result is not significant. eDNA metabarcoding indicates there is no 
## substantial difference in species replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## species in ponds with no fish are not substituted by species in ponds
## with crucian carp.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
eDNA.nomicro.bd.nest <- betadisper(eDNA.nomicro.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
eDNA.nomicro.nest <- with(metadata, eDNA.nomicro.bd.nest)
eDNA.nomicro.nest

## Compute mean distance to centroid per group
tapply(eDNA.nomicro.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(eDNA.nomicro.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(eDNA.nomicro.bd.nest)

## Boxplot of turnover partition
boxplot(eDNA.nomicro.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is no difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(eDNA.nomicro.bd.nest)     # No significant difference between ponds
permutest(eDNA.nomicro.bd.nest) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(eDNA.nomicro.bd.nest)  # No significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
eDNA.nomicro.comm.nest <- metaMDS(eDNA.nomicro.dist$beta.jne, 
                                  dist="jaccard", 
                                  k=2,
                                  maxit=999,
                                  trymax=1000,
                                  wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
eDNA.nomicro.comm.nest$stress
stressplot(eDNA.nomicro.comm.nest)

## plot site scores as text
ordiplot(eDNA.nomicro.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
eDNA.nomicro.NMDS1 <- eDNA.nomicro.comm.nest$points[,1]
eDNA.nomicro.NMDS2 <- eDNA.nomicro.comm.nest$points[,2]
eDNA.nomicro.nest.NMDS <- data.frame(NMDS1=eDNA.nomicro.NMDS1, 
                                     NMDS2=eDNA.nomicro.NMDS2,
                                     Crucian = metadata$Crucian)

## Check data
head(eDNA.nomicro.nest.NMDS)

## Plot data frame
p28b <- ggplot(eDNA.nomicro.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p28b

## Statistically check difference in nestedness of communities
eDNA.nomicro.nest.anosim <- anosim(eDNA.nomicro.dist$beta.jne, metadata$Crucian)

## Inspect results
eDNA.nomicro.nest.anosim
summary(eDNA.nomicro.nest.anosim)
plot(eDNA.nomicro.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
eDNA.nomicro.nest.adonis <- adonis(eDNA.nomicro.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
eDNA.nomicro.nest.adonis

## Result is not significant. eDNA metabarcoding indicates there is no 
## substantial difference in species loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## 3. TOTAL BETA DIVERSITY
eDNA.nomicro.bd.total <- betadisper(eDNA.nomicro.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
eDNA.nomicro.total <- with(metadata, eDNA.nomicro.bd.total)
eDNA.nomicro.total

## Compute mean distance to centroid per group
tapply(eDNA.nomicro.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(eDNA.nomicro.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(eDNA.nomicro.bd.total)

## Boxplot of total beta diversity
boxplot(eDNA.nomicro.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is substantial difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(eDNA.nomicro.bd.total)     # Significant difference between ponds
permutest(eDNA.nomicro.bd.total) # Significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(eDNA.nomicro.bd.total)  # Significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
eDNA.nomicro.comm.total <- metaMDS(eDNA.nomicro.dist$beta.jac, 
                                   dist="jaccard", 
                                   k=2,
                                   maxit=999,
                                   trymax=1000,
                                   wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
eDNA.nomicro.comm.total$stress
stressplot(eDNA.nomicro.comm.total)

## plot site scores as text
ordiplot(eDNA.nomicro.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
eDNA.nomicro.NMDS1 <- eDNA.nomicro.comm.total$points[,1]
eDNA.nomicro.NMDS2 <- eDNA.nomicro.comm.total$points[,2]
eDNA.nomicro.total.NMDS <- data.frame(NMDS1=eDNA.nomicro.NMDS1,
                                      NMDS2=eDNA.nomicro.NMDS2,
                                      Crucian = metadata$Crucian)

## Check data
head(eDNA.nomicro.total.NMDS)

## Plot data frame
p28c <- ggplot(eDNA.nomicro.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p28c

## Statistically check difference in total beta diversity of ponds
eDNA.nomicro.total.anosim <- anosim(eDNA.nomicro.dist$beta.jac, metadata$Crucian)

## Inspect results
eDNA.nomicro.total.anosim
summary(eDNA.nomicro.total.anosim)
plot(eDNA.nomicro.total.anosim)

## There appears to be a significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
eDNA.nomicro.total.adonis <- adonis(eDNA.nomicro.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
eDNA.nomicro.total.adonis

## Again result is significant. eDNA metabarcoding indicates there is 
## substantial variation in overall species composition of
## ponds with crucian carp and ponds without crucian carp


#-----------------------#
# FAMILY-LEVEL ANALYSIS #
#-----------------------#

## BASIC SUMMARIES
## Total number of sequences per pond:
sum.of.rows <- apply(eDNA.fam.nomicro, 1, sum)
sort(sum.of.rows, decreasing = TRUE) 
sum(sum.of.rows)

## Total number of sequences per species across all sites:
sum.of.columns <- apply(eDNA.fam.nomicro, 2, sum)

## Actual read counts:
sort(sum.of.columns, decreasing = TRUE)
sum(sum.of.columns)

## Proportional read counts:
sort(round(sum.of.columns/sum(sum.of.columns)*100, 2), decreasing = TRUE)

## Number of ponds where each species occurs:
eDNA.fam.nomicro.pres <- apply(eDNA.fam.nomicro > 0, 2, sum) 
sort(eDNA.fam.nomicro.pres, decreasing = TRUE)


## ALPHA DIVERSITY
## Convert read counts to presence-absence
eDNA.fam.nomicro.pa <- decostand(eDNA.fam.nomicro, method = "pa", MARGIN = 1)
head(eDNA.fam.nomicro.pa)

## Basic richness
eDNA.nomicro.fam.richness <- specnumber(eDNA.fam.nomicro.pa)
eDNA.nomicro.fam.richness

## Statistically compare and plot alpha diversity
## Create data frame
eDNA.nomicro.fam.alpha <- data.frame(eDNA.nomicro.fam.richness)

## Add sample metadata
eDNA.nomicro.fam.alpha <- cbind(metadata[,1:2], eDNA.nomicro.fam.alpha)

## Reset row names of data frame for further indexing
rownames(eDNA.nomicro.fam.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
eDNA.nomicro.fam.regression <- glm.nb(eDNA.nomicro.fam.richness ~ Crucian, data=eDNA.nomicro.fam.alpha)
summary(eDNA.nomicro.fam.regression)
anova(eDNA.nomicro.fam.regression, test = "Chi")
drop1(eDNA.nomicro.fam.regression, test = "Chi")
summary(glht(glm.nb(eDNA.nomicro.fam.richness ~ Crucian, data=eDNA.nomicro.fam.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
20.138/16
1-pchisq(20.138, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(eDNA.nomicro.fam.alpha$eDNA.nomicro.fam.richness ~ fitted(eDNA.nomicro.fam.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(eDNA.nomicro.fam.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(eDNA.nomicro.fam.regression)
shapiro.test(sresid) # P = 0.6574

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ eDNA.nomicro.fam.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ eDNA.nomicro.fam.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(eDNA.nomicro.fam.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(eDNA.nomicro.fam.regression)
summary(influence)
CD <- cooks.distance(eDNA.nomicro.fam.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot family richness
p29 <- ggplot(eDNA.nomicro.fam.alpha, aes(x=Crucian, y=eDNA.nomicro.fam.richness)) + 
  geom_jitter(aes(colour=Crucian), cex=2, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,30), breaks=seq(0,30,10)) + 
  labs(subtitle="",
       x="", y="") + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p29

## Family richness does not appear higher in ponds without crucian
## carp, indicating this species may have a negligible impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot family accumulation curve across ponds
plot(specaccum(eDNA.fam.nomicro.pa), 
     xlab = "Number of ponds", 
     ylab = "Number of species",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## BETA DIVERSITY
eDNA.nomicro.fam.multi <- beta.multi(eDNA.fam.nomicro.pa, index.family="jaccard")
print(eDNA.nomicro.fam.multi)

## The majority of total beta diversity arises from family turnover
## rather than nestedness.
## Pairwise between-site values of each component of beta diversity
eDNA.nomicro.fam.dist <- beta.pair(eDNA.fam.nomicro.pa, index.family="jaccard")

## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
eDNA.nomicro.fam.bd.turn <- betadisper(eDNA.nomicro.fam.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
eDNA.nomicro.turn <- with(metadata, eDNA.nomicro.fam.bd.turn)
eDNA.nomicro.turn

## Compute mean distance to centroid per group
tapply(eDNA.nomicro.fam.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(eDNA.nomicro.fam.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(eDNA.nomicro.fam.bd.turn)

## Boxplot of turnover partition
boxplot(eDNA.nomicro.fam.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is no real difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(eDNA.nomicro.fam.bd.turn)     # No significant difference between ponds
permutest(eDNA.nomicro.fam.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(eDNA.nomicro.fam.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
eDNA.nomicro.fam.comm.turn <- metaMDS(eDNA.nomicro.fam.dist$beta.jtu, 
                                      dist="jaccard", 
                                      k=2,
                                      maxit=999,
                                      trymax=1000,
                                      wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
eDNA.nomicro.fam.comm.turn$stress
stressplot(eDNA.nomicro.fam.comm.turn)

## plot site scores as text
ordiplot(eDNA.nomicro.fam.comm.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
eDNA.nomicro.NMDS1 <- eDNA.nomicro.fam.comm.turn$points[,1]
eDNA.nomicro.NMDS2 <- eDNA.nomicro.fam.comm.turn$points[,2]
eDNA.nomicro.fam.turn.NMDS <- data.frame(NMDS1=eDNA.nomicro.NMDS1, 
                                         NMDS2=eDNA.nomicro.NMDS2,
                                         Crucian = metadata$Crucian)

## Check data
head(eDNA.nomicro.fam.turn.NMDS)

## Plot data frame
p30a <- ggplot(eDNA.nomicro.fam.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="(iii) eDNA metabarcoding", 
       x="", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p30a

## Statistically check difference in spatial turnover of communities
eDNA.nomicro.fam.turn.anosim <- anosim(eDNA.nomicro.fam.dist$beta.jtu, metadata$Crucian)

## Inspect results
eDNA.nomicro.fam.turn.anosim
summary(eDNA.nomicro.fam.turn.anosim)
plot(eDNA.nomicro.fam.turn.anosim)

## There appears to be a significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
eDNA.nomicro.fam.turn.adonis <- adonis(eDNA.nomicro.fam.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
eDNA.nomicro.fam.turn.adonis

## Result is significant. eDNA metabarcoding indicates there is
## some difference in family replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## families in ponds with no fish are substituted by families in ponds
## with crucian carp.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
eDNA.nomicro.fam.bd.nest <- betadisper(eDNA.nomicro.fam.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
eDNA.nomicro.fam.nest <- with(metadata, eDNA.nomicro.fam.bd.nest)
eDNA.nomicro.fam.nest

## Compute mean distance to centroid per group
tapply(eDNA.nomicro.fam.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(eDNA.nomicro.fam.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(eDNA.nomicro.fam.bd.nest)

## Boxplot of turnover partition
boxplot(eDNA.nomicro.fam.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is slight difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(eDNA.nomicro.fam.bd.nest)     # No significant difference between ponds
permutest(eDNA.nomicro.fam.bd.nest) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(eDNA.nomicro.fam.bd.nest)  # No significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
eDNA.nomicro.fam.comm.nest <- metaMDS(eDNA.nomicro.fam.dist$beta.jne, 
                                      dist="jaccard", 
                                      k=2,
                                      maxit=999,
                                      trymax=1000,
                                      wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
eDNA.nomicro.fam.comm.nest$stress
stressplot(eDNA.nomicro.fam.comm.nest)

## plot site scores as text
ordiplot(eDNA.nomicro.fam.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
eDNA.nomicro.NMDS1 <- eDNA.nomicro.fam.comm.nest$points[,1]
eDNA.nomicro.NMDS2 <- eDNA.nomicro.fam.comm.nest$points[,2]
eDNA.nomicro.fam.nest.NMDS <- data.frame(NMDS1=eDNA.nomicro.NMDS1, 
                                         NMDS2=eDNA.nomicro.NMDS2,
                                         Crucian = metadata$Crucian)

## Check data
head(eDNA.nomicro.fam.nest.NMDS)

## Plot data frame
p30b <- ggplot(eDNA.nomicro.fam.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p30b

## Statistically check difference in nestedness of communities
eDNA.nomicro.fam.nest.anosim <- anosim(eDNA.nomicro.fam.dist$beta.jne, metadata$Crucian)

## Inspect results
eDNA.nomicro.fam.nest.anosim
summary(eDNA.nomicro.fam.nest.anosim)
plot(eDNA.nomicro.fam.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
eDNA.nomicro.fam.nest.adonis <- adonis(eDNA.nomicro.fam.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
eDNA.nomicro.fam.nest.adonis

## Again result is not significant. eDNA metabarcoding indicates there is no 
## substantial difference in family loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## 3. TOTAL BETA DIVERSITY
eDNA.nomicro.fam.bd.total <- betadisper(eDNA.nomicro.fam.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
eDNA.nomicro.total <- with(metadata, eDNA.nomicro.fam.bd.total)
eDNA.nomicro.total

## Compute mean distance to centroid per group
tapply(eDNA.nomicro.fam.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(eDNA.nomicro.fam.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(eDNA.nomicro.fam.bd.total)

## Boxplot of total beta diversity
boxplot(eDNA.nomicro.fam.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is large difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(eDNA.nomicro.fam.bd.total)     # No significant difference between ponds
permutest(eDNA.nomicro.fam.bd.total) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(eDNA.nomicro.fam.bd.total)  # No significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
eDNA.nomicro.fam.comm.total <- metaMDS(eDNA.nomicro.fam.dist$beta.jac, 
                                       dist="jaccard", 
                                       k=2,
                                       maxit=999,
                                       trymax=1000,
                                       wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
eDNA.nomicro.fam.comm.total$stress
stressplot(eDNA.nomicro.fam.comm.total)

## plot site scores as text
ordiplot(eDNA.nomicro.fam.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
eDNA.nomicro.NMDS1 <- eDNA.nomicro.fam.comm.total$points[,1]
eDNA.nomicro.NMDS2 <- eDNA.nomicro.fam.comm.total$points[,2]
eDNA.nomicro.fam.total.NMDS <- data.frame(NMDS1=eDNA.nomicro.NMDS1,
                                          NMDS2=eDNA.nomicro.NMDS2,
                                          Crucian = metadata$Crucian)

## Check data
head(eDNA.nomicro.fam.total.NMDS)

## Plot data frame
p30c <- ggplot(eDNA.nomicro.fam.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p30c

## Statistically check difference in total beta diversity of ponds
eDNA.nomicro.fam.total.anosim <- anosim(eDNA.nomicro.fam.dist$beta.jac, metadata$Crucian)

## Inspect results
eDNA.nomicro.fam.total.anosim
summary(eDNA.nomicro.fam.total.anosim)
plot(eDNA.nomicro.fam.total.anosim)

## There appears to be a significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
eDNA.nomicro.fam.total.adonis <- adonis(eDNA.nomicro.fam.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
eDNA.nomicro.fam.total.adonis

## Again result is significant. eDNA metabarcoding indicates there 
## is some variation in overall family composition of ponds with 
## crucian carp and ponds without crucian carp.



####################
# METHODS COMBINED #
####################

#------------------------#
# SPECIES-LEVEL ANALYSIS #
#------------------------#

## BASIC SUMMARIES
## Number of ponds where each species occurs
comb.nomicro.spec.pres <- apply(comb.nomicro > 0, 2, sum)
sort(comb.nomicro.spec.pres, decreasing = TRUE)

## Tranform data from counts (specimens or reads) to presence-absence by
## rows (ponds)
comb.nomicro.pa <- decostand(comb.nomicro, method = "pa", MARGIN = 1)
head(comb.nomicro.pa[,1:3], n = 3)


## ALPHA DIVERSITY
## Basic richness
comb.nomicro.richness <- specnumber(comb.nomicro.pa)
comb.nomicro.richness

## Statistically compare and plot alpha diversity
## Create data frame
comb.nomicro.alpha <- data.frame(comb.nomicro.richness)

## Add sample metadata
comb.nomicro.alpha <- cbind(metadata[,1:2], comb.nomicro.alpha)

## Reset row names of data frame for further indexing
rownames(comb.nomicro.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
comb.nomicro.regression <- glm.nb(comb.nomicro.richness ~ Crucian, data=comb.nomicro.alpha)
summary(comb.nomicro.regression)
anova(comb.nomicro.regression, test = "Chi")
drop1(comb.nomicro.regression, test = "Chi")
summary(glht(glm.nb(comb.nomicro.richness ~ Crucian, data=comb.nomicro.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
17.892/16
1-pchisq(17.892, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(comb.nomicro.alpha$comb.nomicro.richness ~ fitted(comb.nomicro.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(comb.nomicro.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(comb.nomicro.regression)
shapiro.test(sresid) # P = 0.4579

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ comb.nomicro.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ comb.nomicro.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(comb.nomicro.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(comb.nomicro.regression)
summary(influence)
CD <- cooks.distance(comb.nomicro.regression)
plot(CD ~ sresid)

## Plot species richness
p31 <- ggplot(comb.nomicro.alpha, aes(x=Crucian, y=comb.nomicro.richness)) + 
  geom_jitter(aes(colour=Crucian), cex=2, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,60), breaks=seq(0,60,10)) + 
  labs(subtitle="(iv) Methods combined", 
       x="Crucian carp", y=expression(alpha~diversity)) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p31

## Species richness appears to be the same in ponds with and without crucian
## carp, indicating this species may have negligible impact on the 
## invertebrate community.

## Was sample size enough to accurately represent invertebrate diversity?
## Plot species accumulation curve across ponds
plot(specaccum(comb.nomicro.pa), 
     xlab = "Number of ponds", 
     ylab = "Number of species",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## Examine alpha diversity of the major invertebrate groups individually 
## Transpose dataframe
nomicro.order.df <- data.frame(t(comb.nomicro.order.spp))

## Tranform data to presence-absence
nomicro.order.df[nomicro.order.df > 0] <- 1

## Make row names a new column in dataframe
nomicro.order.df <- rownames_to_column(nomicro.order.df, "Species")

## Add order data to netting data
nomicro.order.df <- merge(nomicro.order.df, spp.order, by.x="Species", by.y="Species", all.x=TRUE)
nomicro.order.df <- nomicro.order.df[,c(1,20,2:19)]

## Remove species column, and merge data by order
nomicro.order.df <- nomicro.order.df[,-1]
nomicro.order.df <- ddply(nomicro.order.df, .(Group), numcolwise(sum))

## Make Group column row names and transpose dataframe
nomicro.order.df <- column_to_rownames(nomicro.order.df, "Group")
nomicro.order.dat <- data.frame(t(nomicro.order.df))

## Make row names column in dataframe
nomicro.order.dat <- rownames_to_column(nomicro.order.dat, "Pond")

## Create new column specifying whether pond contains crucian carp
nomicro.order.dat$Crucian <- factor(ifelse(nomicro.order.dat$Pond %in% metadata$Site, as.character(metadata$Crucian), "NA"))
nomicro.order.dat <- nomicro.order.dat[,c(1,20,2:19)]

## Melt dataframe and rename new columns
nomicro.order.melt <- melt(nomicro.order.dat, id=c("Pond","Crucian"))
colnames(nomicro.order.melt)[3:4] <- c("Group","Richness")

## Statistically test for differences in invertebrate species richness
nomicro.order.alpha.regression <- glm.nb(Richness ~ Group/Crucian, data=nomicro.order.melt)
summary(nomicro.order.alpha.regression)
anova(nomicro.order.alpha.regression, test = "Chi")
drop1(nomicro.order.alpha.regression, test = "Chi")

## Check model meets GLM assumptions
## Test for overdispersion
251.29/288
1-pchisq(251.29, df=288)  # not overdispersed

## Plot the fitted data against the observed data
plot(nomicro.order.melt$Richness ~ fitted(nomicro.order.alpha.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(nomicro.order.alpha.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(nomicro.order.alpha.regression)
shapiro.test(sresid) # P = 6.183e-14

## Some deviation from normality as residuals are normally distributed
## therefore model may not be reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ nomicro.order.alpha.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ nomicro.order.melt$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity
plot(sresid ~ nomicro.order.melt$Group, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
vif(nomicro.order.alpha.regression)

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(nomicro.order.alpha.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(nomicro.order.alpha.regression)
summary(influence)
CD <- cooks.distance(nomicro.order.alpha.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Create text to be added to plot indicating significance
text1 <- data.frame(label=c("","",".","","","","","",""),
                    Group=c("Annelida","Arachnida",
                            "Coleoptera","Collembola",
                            "Crustacea","Diptera",
                            "Ephemeroptera","Hemiptera",
                            "Hirudinea"))

text2 <- data.frame(label=c("","","","**","","","","",""),
                    Group=c("Hymenoptera","Lepidoptera",
                            "Megaloptera","Mollusca",
                            "Odonata","Platyhelminthes",
                            "Psocoptera","Rotifera",
                            "Trichoptera"))

## Plot Group richness in crucian carp and non-fish ponds
p32a <- ggplot(nomicro.order.melt[nomicro.order.melt$Group %in% c("Annelida","Arachnida",
                                                                  "Coleoptera","Collembola",
                                                                  "Crustacea","Diptera",
                                                                  "Ephemeroptera","Hemiptera",
                                                                  "Hirudinea"),],
               aes(x=Crucian, y=Richness)) + 
  geom_jitter(aes(colour=Crucian), width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  coord_cartesian(ylim=c(0,15)) + 
  scale_y_continuous(breaks=seq(0,15,1)) + 
  labs(title="(a) Species-level",
       x="", y=expression(alpha~diversity)) + 
  geom_text(data=text1, mapping=aes(x=1.5, y=15, label=label), size=10) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin=unit(c(2,0,2,0), "mm")),
        axis.line.x = element_line(colour="black", size=0.5, linetype="solid"),
        axis.line.y = element_line(colour="black", size=0.5, linetype="solid"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin=unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=20),
        strip.text.x = element_text(size=12),
        legend.position = "none",
        text = element_text(size=20)) + 
  facet_grid(. ~ Group)
p32a

p32b <- ggplot(nomicro.order.melt[nomicro.order.melt$Group %in% c("Hymenoptera","Lepidoptera",
                                                                  "Megaloptera","Mollusca",
                                                                  "Odonata","Platyhelminthes",
                                                                  "Psocoptera","Rotifera",
                                                                  "Trichoptera"),],
               aes(x=Crucian, y=Richness)) + 
  geom_jitter(aes(colour=Crucian), width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  coord_cartesian(ylim=c(0,15)) + 
  scale_y_continuous(breaks=seq(0,15,1)) + 
  labs(x="Crucian carp", y=expression(alpha~diversity)) + 
  geom_text(data=text2, mapping=aes(x=1.5, y=15, label=label), size=10) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        axis.line.x = element_line(colour="black", size=0.5, linetype="solid"),
        axis.line.y = element_line(colour="black", size=0.5, linetype="solid"),
        axis.title.x = element_text(margin=unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin=unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=20),
        strip.text.x = element_text(size=12),
        legend.position = "none",
        text = element_text(size=20)) + 
  facet_grid(. ~ Group)
p32b

## Plot together
p32 <- ggarrange(p32a,p32b, ncol=1, nrow=2)


## BETA DIVERSITY
## Examine beta diversity across all ponds
comb.nomicro.multi <- beta.multi(comb.nomicro.pa, index.family="jaccard")
print(comb.nomicro.multi)

## The majority of total beta diversity arises from species turnover
## rather than nestedness.
## Now get pairwise between-site values of each component of beta diversity
comb.nomicro.dist <- beta.pair(comb.nomicro.pa, index.family="jaccard")

## Convert each partition of beta diversity to a matrix if you want to 
## write as a csv file
#comb.turn <- as.matrix(dist(net.beta$beta.jtu))
#comb.nest <- as.matrix(dist(net.beta$beta.jne))
#comb.total <- as.matrix(dist(net.beta$beta.jac))


## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)

## First, test whether ponds with and without crucian carp have different
## communities using PERMANOVA and visualise differences with NMDS
comb.nomicro.bd.turn <- betadisper(comb.nomicro.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
comb.nomicro.turn <- with(metadata, comb.nomicro.bd.turn)
comb.nomicro.turn

## Compute mean distance to centroid per group
tapply(comb.nomicro.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(comb.nomicro.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(comb.nomicro.bd.turn)

## Boxplot of turnover partition
boxplot(comb.nomicro.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is substantial difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(comb.nomicro.bd.turn)     # Significant difference between ponds
permutest(comb.nomicro.bd.turn) # Significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(comb.nomicro.bd.turn)  # Significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
comb.nomicro.bd.turn <- metaMDS(comb.nomicro.dist$beta.jtu, 
                                dist="jaccard",
                                k=2,
                                maxit=999,
                                trymax=1000,
                                wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
comb.nomicro.bd.turn$stress
stressplot(comb.nomicro.bd.turn)

## plot site scores as text
ordiplot(comb.nomicro.bd.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
comb.nomicro.NMDS1 <- comb.nomicro.bd.turn$points[,1]
comb.nomicro.NMDS2 <- comb.nomicro.bd.turn$points[,2]
comb.nomicro.turn.NMDS <- data.frame(NMDS1=comb.nomicro.NMDS1, 
                                     NMDS2=comb.nomicro.NMDS2,
                                     Crucian = metadata$Crucian)

## Check data
head(comb.nomicro.turn.NMDS)

## Plot data frame
p33a <- ggplot(comb.nomicro.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="(iv) Methods combined", 
       x="NMDS1", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p33a

## Statistically check difference in spatial turnover of communities
comb.nomicro.turn.anosim <- anosim(comb.nomicro.dist$beta.jtu, metadata$Crucian)

## Inspect results
comb.nomicro.turn.anosim
summary(comb.nomicro.turn.anosim)
plot(comb.nomicro.turn.anosim)

## There appears to be a significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
comb.nomicro.turn.adonis <- adonis(comb.nomicro.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
comb.nomicro.turn.adonis

## Again result is significant. The combined data indicates there is a
## substantial difference in species replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## species in ponds with no fish are substituted by species in ponds
## with crucian carp.

## Use RDA and variance partitioning to examine the influence of biotic and 
## abiotic variables on turnover partition of beta diversity.
## First, compute principal coordinate decomposition (classical scaling) for 
## turnover distance matrix using the pcoa() function from the ape package. Use 
## lingoes correction to account for negative eigenvalues.
pcoa.jtu <- pcoa(comb.nomicro.dist$beta.jtu, correction="lingoes", rn=NULL)

## Inspect outputs
pcoa.jtu$correction
pcoa.jtu$note
pcoa.jtu$trace

## The principal component coordinates for each diversity matrix can be written
## to a csv file for inspection
#write.csv(pcoa.jtu$vectors, "../Data/Combined_data_species_turnover_eigenvectors.csv")

## Create groups of variables for RDA
## Biotic:
biotic <- data.frame(env.data$Crucian)
colnames(biotic) <- "Crucian"

## Abiotic:
abiotic <- env.data[,c(2:6)]

## log10 transform abiotic data to remove units
abiotic <- log10(abiotic+1)

## Perform forward model selection for each group of explanatory variables
## Adjusted R-squared will be used as a goodness-of-fit measure. Model selection
## needs to be performed for each group of variables, excluding those containing
## only one variable.
abiotic.rda.null <- rda(pcoa.jtu$vectors ~ 1, abiotic)  # null model
abiotic.rda.all <- rda(pcoa.jtu$vectors ~ ., abiotic)   # model with all variables
step.env <- ordiR2step(abiotic.rda.null, abiotic.rda.all, 
                       direction="forward", perm.max = 1000)

## Summary table and plot:
step.env$anova
plot(step.env)
## Indicates that no abiotic variables should be included in model

## Run RDA with biotic variable
## Test for significance of the whole model:
mod2 <- rda(pcoa.jtu$vectors ~ biotic$Crucian)
summary(mod2)
coef(mod2)       # canonical coefficients
vif.cca(mod2)    # no collinearity present in model
mod2             # 9.6% variance is explained by variables in model
anova(mod2)      # Model is significant
RsquareAdj(mod2) # 3.9% variance explained
plot(mod2)

## Check whether any variables should be dropped from model
drop1(mod2, test="perm")
## Crucian carp presence-absence significant

## Plot model result to get a sense of which variables are correlating with which 
## ponds along which axes
plot(mod2, type='n', scaling=1)
orditorp(mod2, display='sp', cex=0.5, scaling=1, col='blue')
text(mod2, display='cn', col='red')

## Plot of variables against sites
plot(mod2, type = "text")

## Perform hypothesis testing
anova(mod2, perm=1000)               # overall model fit
anova(mod2, by="margin", perm=1000)  # significance of marginal effects (Type III)
anova(mod2, by="term", perm=1000)    # significance of each term in model
anova(mod2, by="axis", perm=1000)    # significance of each axis


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
comb.nomicro.bd.nest <- betadisper(comb.nomicro.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
comb.nomicro.nest <- with(metadata, comb.nomicro.bd.nest)
comb.nomicro.nest

## Compute mean distance to centroid per group
tapply(comb.nomicro.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(comb.nomicro.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(comb.nomicro.bd.nest)

## Boxplot of turnover partition
boxplot(comb.nomicro.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is no difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(comb.nomicro.bd.nest)     # No significant difference between ponds
permutest(comb.nomicro.bd.nest) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(comb.nomicro.bd.nest)  # No significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
comb.nomicro.comm.nest <- metaMDS(comb.nomicro.dist$beta.jne, 
                                  dist="jaccard", 
                                  k=2,
                                  maxit=999,
                                  trymax=1000,
                                  wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
comb.nomicro.comm.nest$stress
stressplot(comb.nomicro.comm.nest)

## plot site scores as text
ordiplot(comb.nomicro.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
comb.nomicro.NMDS1 <- comb.nomicro.comm.nest$points[,1]
comb.nomicro.NMDS2 <- comb.nomicro.comm.nest$points[,2]
comb.nomicro.nest.NMDS <- data.frame(NMDS1=comb.nomicro.NMDS1, 
                                     NMDS2=comb.nomicro.NMDS2,
                                     Crucian = metadata$Crucian)

## Check data
head(comb.nomicro.nest.NMDS)

## Plot data frame
p33b <- ggplot(comb.nomicro.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", x="NMDS1", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p33b

## Statistically check difference in nestedness of communities
comb.nomicro.nest.anosim <- anosim(comb.nomicro.dist$beta.jne, metadata$Crucian)

## Inspect results
comb.nomicro.nest.anosim
summary(comb.nomicro.nest.anosim)
plot(comb.nomicro.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
comb.nomicro.nest.adonis <- adonis(comb.nomicro.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
comb.nomicro.nest.adonis

## Result is not significant. The combined data indicates there is no
## substantial difference in species loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## Use RDA and variance partitioning to examine the influence of biotic and 
## abiotic variables on nestedness partition of beta diversity.
## First, compute principal coordinate decomposition (classical scaling) for 
## nestedness distance matrix using the pcoa() function from the ape package. Use 
## lingoes correction to account for negative eigenvalues.
pcoa.jne <- pcoa(comb.nomicro.dist$beta.jne, correction="lingoes", rn=NULL)

## Inspect outputs
pcoa.jne$correction
pcoa.jne$note
pcoa.jne$trace

## The principal component coordinates for each diversity matrix can be written
## to a csv file for inspection
#write.csv(pcoa.jne$vectors, "../Data/Combined_data_species_nestedness_eigenvectors.csv")

## Groups of variables were previously created for RDA of turnover
## Perform forward model selection for each group of explanatory variables
## Adjusted R-squared will be used as a goodness-of-fit measure. Model selection
## needs to be performed for each group of variables, excluding those containing
## only one variable.
abiotic.rda.null <- rda(pcoa.jne$vectors ~ 1, abiotic)  # null model
abiotic.rda.all <- rda(pcoa.jne$vectors ~ ., abiotic)   # model with all variables
step.env <- ordiR2step(abiotic.rda.null, abiotic.rda.all, 
                       direction="forward", perm.max = 1000)

## Summary table and plot:
step.env$anova
plot(step.env)
## Indicates that no abiotic variables should be included in model

## Run RDA with biotic variable
## Test for significance of the whole model:
mod2 <- rda(pcoa.jne$vectors ~ biotic$Crucian)
summary(mod2)
coef(mod2)       # canonical coefficients
vif.cca(mod2)    # no collinearity present in model
mod2             # 3% variance is explained by variables in model
anova(mod2)      # Model is not significant
RsquareAdj(mod2) # -3.1% variance explained?
plot(mod2)

## Check whether any variables should be dropped from model
drop1(mod2, test="perm")
## Crucian carp presence-absence not significant

## Plot model result to get a sense of which variables are correlating with which 
## ponds along which axes
plot(mod2, type='n', scaling=1)
orditorp(mod2, display='sp', cex=0.5, scaling=1, col='blue')
text(mod2, display='cn', col='red')

## Plot of variables against sites
plot(mod2, type = "text")

## Perform hypothesis testing
anova(mod2, perm=1000)               # overall model fit
anova(mod2, by="margin", perm=1000)  # significance of marginal effects (Type III)
anova(mod2, by="term", perm=1000)    # significance of each term in model
anova(mod2, by="axis", perm=1000)    # significance of each axis


## 3. TOTAL BETA DIVERSITY
comb.nomicro.bd.total <- betadisper(comb.nomicro.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
comb.nomicro.total <- with(metadata, comb.nomicro.bd.total)
comb.nomicro.total

## Compute mean distance to centroid per group
tapply(comb.nomicro.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(comb.nomicro.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(comb.nomicro.bd.total)

## Boxplot of total beta diversity
boxplot(comb.nomicro.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is substantial difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(comb.nomicro.bd.total)     # Significant difference between ponds
permutest(comb.nomicro.bd.total) # Significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(comb.nomicro.bd.total)  # Significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
comb.nomicro.comm.total <- metaMDS(comb.nomicro.dist$beta.jac, 
                                   dist="jaccard", 
                                   k=2,
                                   maxit=999,
                                   trymax=1000,
                                   wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
comb.nomicro.comm.total$stress
stressplot(comb.nomicro.comm.total)

## plot site scores as text
ordiplot(comb.nomicro.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
comb.nomicro.NMDS1 <- comb.nomicro.comm.total$points[,1]
comb.nomicro.NMDS2 <- comb.nomicro.comm.total$points[,2]
comb.nomicro.total.NMDS <- data.frame(NMDS1=comb.nomicro.NMDS1,
                                      NMDS2=comb.nomicro.NMDS2,
                                      Crucian = metadata$Crucian)

## Check data
head(comb.nomicro.total.NMDS)

## Plot data frame
p33c <- ggplot(comb.nomicro.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="NMDS1", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p33c

## Statistically check difference in total beta diversity of ponds
comb.nomicro.total.anosim <- anosim(comb.nomicro.dist$beta.jac, metadata$Crucian)

## Inspect results
comb.nomicro.total.anosim
summary(comb.nomicro.total.anosim)
plot(comb.nomicro.total.anosim)

## There appears to be a significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
comb.nomicro.total.adonis <- adonis(comb.nomicro.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
comb.nomicro.total.adonis

## Again result is significant. The combined data indicates there is 
## substantial variation in overall species composition of
## ponds with crucian carp and ponds without crucian carp


## Use RDA and variance partitioning to examine the influence of biotic and 
## abiotic variables on total beta diversity

## First, compute principal coordinate decomposition (classical scaling) for 
## total beta distance matrix using the pcoa() function from the ape package. Use 
## lingoes correction to account for negative eigenvalues.
pcoa.jac <- pcoa(comb.nomicro.dist$beta.jac, correction="lingoes", rn=NULL)

## Inspect outputs
pcoa.jac$correction
pcoa.jac$note
pcoa.jac$trace

## The principal component coordinates for each diversity matrix can be written
## to a csv file for inspection
#write.csv(pcoa.jac$vectors, "../Data/Combined_data_species_totalbeta_eigenvectors.csv")

## Groups of variables were previously created for RDA of turnover
## Perform forward model selection for each group of explanatory variables
## Adjusted R-squared will be used as a goodness-of-fit measure. Model selection
## needs to be performed for each group of variables, excluding those containing
## only one variable.
abiotic.rda.null <- rda(pcoa.jac$vectors ~ 1, abiotic)  # null model
abiotic.rda.all <- rda(pcoa.jac$vectors ~ ., abiotic)   # model with all variables
step.env <- ordiR2step(abiotic.rda.null, abiotic.rda.all, 
                       direction="forward", perm.max = 1000)

## Summary table and plot:
step.env$anova
plot(step.env)
## Indicates that no abiotic variables should be included in model

## Run RDA with biotic variable
## Test for significance of the whole model:
mod2 <- rda(pcoa.jac$vectors ~ biotic$Crucian)
summary(mod2)
coef(mod2)       # canonical coefficients
vif.cca(mod2)    # no collinearity present in model
mod2             # 8.8% variance is explained by variables in model
anova(mod2)      # Model is not significant
RsquareAdj(mod2) # 3.1% variance explained
plot(mod2)

## Check whether any variables should be dropped from model
drop1(mod2, test="perm")
## Crucian carp presence-absence not significant

## Plot model result to get a sense of which variables are correlating with which 
## ponds along which axes
plot(mod2, type='n', scaling=1)
orditorp(mod2, display='sp', cex=0.5, scaling=1, col='blue')
text(mod2, display='cn', col='red')

## Plot of variables against sites
plot(mod2, type = "text")

## Perform hypothesis testing
anova(mod2, perm=1000)               # overall model fit
anova(mod2, by="margin", perm=1000)  # significance of marginal effects (Type III)
anova(mod2, by="term", perm=1000)    # significance of each term in model
anova(mod2, by="axis", perm=1000)    # significance of each axis


#-----------------------#
# FAMILY-LEVEL ANALYSIS #
#-----------------------#

## BASIC SUMMARIES
## Number of ponds where each family occurs:
fam.nomicro.pres <- apply(comb.fam.nomicro > 0, 2, sum) 
sort(fam.nomicro.pres, decreasing = TRUE)[1:18]

## Tranform datafrom counts (specimens or reads) to presence-absence by
## rows (ponds)
comb.nomicro.fam.pa <- decostand(comb.fam.nomicro, method = "pa", MARGIN = 1)
head(comb.nomicro.fam.pa[,1:3], n = 3)


## ALPHA DIVERSITY
## Basic richness
comb.nomicro.fam.richness <- specnumber(comb.nomicro.fam.pa)
comb.nomicro.fam.richness

## Statistically compare and plot alpha diversity
## Create data frame
comb.nomicro.fam.alpha <- data.frame(comb.nomicro.fam.richness)

## add metadata from external file
comb.nomicro.fam.alpha <- cbind(metadata[,1:2], comb.nomicro.fam.alpha)

## Reset row names of data frame for further indexing
rownames(comb.nomicro.fam.alpha) <- NULL

## Statistically compare whether crucian carp influences alpha diversity
comb.nomicro.fam.regression <- glm.nb(comb.nomicro.fam.richness ~ Crucian, data=comb.nomicro.fam.alpha)
summary(comb.nomicro.fam.regression)
anova(comb.nomicro.fam.regression, test = "Chi")
drop1(comb.nomicro.fam.regression, test = "Chi")
summary(glht(glm.nb(comb.nomicro.fam.richness ~ Crucian, data=comb.nomicro.fam.alpha), 
             linfct=mcp(Crucian="Tukey")))

## Check model meets GLM assumptions
## Test for overdispersion
16.087/16
1-pchisq(16.087, df=16)  # not overdispersed

## Plot the fitted data against the observed data
plot(comb.nomicro.fam.alpha$comb.nomicro.fam.richness ~ fitted(comb.nomicro.fam.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(comb.nomicro.fam.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(comb.nomicro.fam.regression)
shapiro.test(sresid) # P = 0.1542

## No deviation from normality as residuals are normally distributed
## therefore model is reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ comb.nomicro.fam.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ comb.nomicro.fam.alpha$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(comb.nomicro.fam.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(comb.nomicro.fam.regression)
summary(influence)
CD <- cooks.distance(comb.nomicro.fam.regression)
plot(CD ~ sresid)

## Plot family richness
p34 <- ggplot(comb.nomicro.fam.alpha, aes(x=Crucian, y=comb.nomicro.fam.richness)) + 
  geom_jitter(aes(colour=Crucian), cex=2, width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape = NA) + 
  scale_y_continuous(limits=c(0,60), breaks=seq(0,60,10)) + 
  labs(subtitle="",
       x="Crucian carp", y="") + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        legend.position = "none",
        text = element_text(size=20))
p34

## Family richness appears to be the same in ponds with and without crucian
## carp, indicating this species may have negligible impact on the 
## invertebrate community.
## Was sample size enough to accurately represent invertebrate diversity?
## Plot family accumulation curve across ponds
plot(specaccum(comb.nomicro.fam.pa), 
     xlab = "Number of ponds", 
     ylab = "Number of families",
     ci.type="polygon", ci.col="grey")

## More ponds would have been required to fully represent invertebrate
## diversity.


## Compare alpha diversity of major invertebrate groups individually in 
## ponds with and without crucian carp.
## Transpose dataframe
nomicro.order.fam.df <- data.frame(t(comb.nomicro.order.fam))

## Tranform data to presence-absence
nomicro.order.fam.df[nomicro.order.fam.df > 0] <- 1

## Make row names a new column in dataframe
nomicro.order.fam.df <- rownames_to_column(nomicro.order.fam.df, "Family")

## Add order data to family data
nomicro.order.fam.df <- merge(nomicro.order.fam.df, fam.order, by.x="Family", by.y="Family", all.x=TRUE)
nomicro.order.fam.df <- nomicro.order.fam.df[,c(1,20,2:19)]

## Remove species column, and merge data by order
nomicro.order.fam.df <- nomicro.order.fam.df[,-1]
nomicro.order.fam.df <- ddply(nomicro.order.fam.df, .(Group), numcolwise(sum))

## Make Group column row names and transpose dataframe
nomicro.order.fam.df <- column_to_rownames(nomicro.order.fam.df, "Group")
nomicro.order.fam.dat <- data.frame(t(nomicro.order.fam.df))

## Make row names column in dataframe
nomicro.order.fam.dat <- rownames_to_column(nomicro.order.fam.dat, "Pond")

## Create new column specifying whether pond contains crucian carp
nomicro.order.fam.dat$Crucian <- factor(ifelse(nomicro.order.fam.dat$Pond %in% metadata$Site, as.character(metadata$Crucian), "NA"))
nomicro.order.fam.dat <- nomicro.order.fam.dat[,c(1,20,2:19)]

## Melt dataframe and rename new columns
nomicro.order.fam.melt <- melt(nomicro.order.fam.dat, id=c("Pond","Crucian"))
colnames(nomicro.order.fam.melt)[3:4] <- c("Group","Richness")

## Statistically test for differences in invertebrate species richness
nomicro.order.alpha.regression <- glm.nb(Richness ~ Group/Crucian, data=nomicro.order.fam.melt)
summary(nomicro.order.alpha.regression)
anova(nomicro.order.alpha.regression, test = "Chi")
drop1(nomicro.order.alpha.regression, test = "Chi")

## Check model meets GLM assumptions
## Test for overdispersion
200.55/288
1-pchisq(200.55, df=288)  # not overdispersed

## Plot the fitted data against the observed data
plot(nomicro.order.fam.melt$Richness ~ fitted(nomicro.order.alpha.regression))

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(nomicro.order.alpha.regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
chkres(nomicro.order.alpha.regression)
shapiro.test(sresid) # P = 1.33e-13

## Some deviation from normality as residuals are normally distributed
## therefore model may not be reliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ nomicro.order.alpha.regression$fitted.values, pch=20, cex=2, cex.lab=1.5)
plot(sresid ~ nomicro.order.fam.melt$Crucian, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity
plot(sresid ~ nomicro.order.fam.melt$Group, pch=20, cex=2, cex.lab=1.5) # No heteroscedascity

## Assumption 3: no collinearity
vif(nomicro.order.alpha.regression)

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(nomicro.order.alpha.regression)

## No significant serial auto-correlation as p-value not significant
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Some autocorrelation at first time lag

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(nomicro.order.alpha.regression)
summary(influence)
CD <- cooks.distance(nomicro.order.alpha.regression)
plot(CD ~ sresid)

## Cook's distance < 1 so observations not exerting strong influence on model
## parameters

## Plot Group richness in crucian carp and non-fish ponds
p35a <- ggplot(nomicro.order.fam.melt[nomicro.order.fam.melt$Group %in% c("Annelida","Arachnida",
                                                                          "Coleoptera","Collembola",
                                                                          "Crustacea","Diptera",
                                                                          "Ephemeroptera","Hemiptera",
                                                                          "Hirudinea"),],
               aes(x=Crucian, y=Richness)) + 
  geom_jitter(aes(colour=Crucian), width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  coord_cartesian(ylim=c(0,10)) + 
  scale_y_continuous(breaks=seq(0,16,1)) + 
  labs(title="(b) Family-level",
       x="", y=expression(alpha~diversity)) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) + 
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin=unit(c(2,0,2,0), "mm")),
        axis.line.x = element_line(colour="black", size=0.5, linetype="solid"),
        axis.line.y = element_line(colour="black", size=0.5, linetype="solid"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin=unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=20),
        strip.text.x = element_text(size=12),
        legend.position = "none",
        text = element_text(size=20)) + 
  facet_grid(. ~ Group)
p35a

p35b <- ggplot(nomicro.order.fam.melt[nomicro.order.fam.melt$Group %in% c("Hymenoptera","Lepidoptera",
                                                                          "Megaloptera","Mollusca",
                                                                          "Odonata","Platyhelminthes",
                                                                          "Psocoptera","Thysanoptera",
                                                                          "Trichoptera"),],
               aes(x=Crucian, y=Richness)) + 
  geom_jitter(aes(colour=Crucian), width=0.2, show.legend=FALSE) + 
  geom_boxplot(alpha=0.7, outlier.shape=NA) + 
  coord_cartesian(ylim=c(0,10)) + 
  scale_y_continuous(breaks=seq(0,16,1)) + 
  labs(x="Crucian carp", y=expression(alpha~diversity)) + 
  scale_colour_manual(values=c("grey40","deepskyblue2")) + 
  scale_x_discrete(breaks=c("N","Y"),
                   labels=c("Absent","Present")) +
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        axis.line.x = element_line(colour="black", size=0.5, linetype="solid"),
        axis.line.y = element_line(colour="black", size=0.5, linetype="solid"),
        axis.title.x = element_text(margin=unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin=unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y = element_text(colour="black", size=20),
        strip.text.x = element_text(size=12),
        legend.position = "none",
        text = element_text(size=20)) + 
  facet_grid(. ~ Group)
p35b


## BETA DIVERSITY
## Examine beta diversity across all ponds
comb.nomicro.fam.multi <- beta.multi(comb.nomicro.fam.pa, index.family="jaccard")
print(comb.nomicro.fam.multi)

## The majority of total beta diversity arises from species turnover
## rather than nestedness.
## Now get pairwise between-site values of each component of beta diversity
comb.nomicro.fam.dist <- beta.pair(comb.nomicro.fam.pa, index.family="jaccard")

## Convert each partition of beta diversity to a matrix if you want to 
## write as a csv file
#comb.nomicro.fam.turn <- as.matrix(dist(comb.nomicro.fam.dist$beta.jtu))
#comb.nomicro.fam.nest <- as.matrix(dist(comb.nomicro.fam.dist$beta.jne))
#comb.nomicro.fam.total <- as.matrix(dist(comb.nomicro.fam.dist$beta.jac))


## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)

## First, test whether ponds with and without crucian carp have different
## communities using PERMANOVA and visualise differences with NMDS
comb.nomicro.fam.bd.turn <- betadisper(comb.nomicro.fam.dist$beta.jtu, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.turn <- with(metadata, comb.nomicro.fam.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(comb.nomicro.fam.bd.turn$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(comb.nomicro.fam.bd.turn$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(comb.nomicro.fam.bd.turn)

## Boxplot of turnover partition
boxplot(comb.nomicro.fam.bd.turn, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicate that there is some difference in turnover partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether turnover is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(comb.nomicro.fam.bd.turn)     # No significant difference between ponds
permutest(comb.nomicro.fam.bd.turn) # No significant difference between ponds

## Analyse pairwise differences between groups (crucian/non-crucian) using 
## parametric Tukey's HSD test.
TukeyHSD(comb.nomicro.fam.bd.turn)  # No significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
comb.nomicro.fam.bd.turn <- metaMDS(comb.nomicro.fam.dist$beta.jtu, 
                                    dist="jaccard", 
                                    k=2,
                                    maxit=999,
                                    trymax=1000,
                                    wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
comb.nomicro.fam.bd.turn$stress
stressplot(comb.nomicro.fam.bd.turn)

## plot site scores as text
ordiplot(comb.nomicro.fam.bd.turn, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
comb.nomicro.NMDS1 <- comb.nomicro.fam.bd.turn$points[,1]
comb.nomicro.NMDS2 <- comb.nomicro.fam.bd.turn$points[,2]
comb.nomicro.fam.turn.NMDS <- data.frame(NMDS1=comb.nomicro.NMDS1, 
                                         NMDS2=comb.nomicro.NMDS2,
                                         Crucian = metadata$Crucian)

## Check data
head(comb.nomicro.fam.turn.NMDS)

## Plot data frame
p36a <- ggplot(comb.nomicro.fam.turn.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="(iv) Methods combined", 
       x="NMDS1", y="NMDS2") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p36a

## Statistically check difference in spatial turnover of communities
comb.nomicro.fam.turn.anosim <- anosim(comb.nomicro.fam.dist$beta.jtu, metadata$Crucian)

## Inspect results
comb.nomicro.fam.turn.anosim
summary(comb.nomicro.fam.turn.anosim)
plot(comb.nomicro.fam.turn.anosim)

## There appears to be no significant difference in spatial turnover
## Look at PermANOVA using adonis(), considered more robust than anosim()
comb.nomicro.fam.turn.adonis <- adonis(comb.nomicro.fam.dist$beta.jtu ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
comb.nomicro.fam.turn.adonis

## Again result is not significant. The combined data indicates there is no
## substantial difference in species replacement (i.e. turnover) between
## ponds with crucian carp and ponds without crucian carp. Therefore,
## species in ponds with no fish are not substituted by species in ponds
## with crucian carp.


## Use RDA and variance partitioning to examine the influence of biotic and 
## abiotic variables on turnover partition of beta diversity

## First, compute principal coordinate decomposition (classical scaling) for 
## turnover distance matrix using the pcoa() function from the ape package. Use 
## lingoes correction to account for negative eigenvalues.
pcoa.jtu <- pcoa(comb.nomicro.fam.dist$beta.jtu, correction="lingoes", rn=NULL)

## Inspect outputs
pcoa.jtu$correction
pcoa.jtu$note
pcoa.jtu$trace

## The principal component coordinates for each diversity matrix can be written
## to a csv file for inspection
#write.csv(pcoa.jtu$vectors, "../Data/Combined_data_family_turnover_eigenvectors.csv")

## Groups of variables for RDA were already created for species-level analysis
## Perform forward model selection for each group of explanatory variables
## Adjusted R-squared will be used as a goodness-of-fit measure. Model selection
## needs to be performed for each group of variables, excluding those containing
## only one variable.
abiotic.rda.null <- rda(pcoa.jtu$vectors ~ 1, abiotic)  # null model
abiotic.rda.all <- rda(pcoa.jtu$vectors ~ ., abiotic)   # model with all variables
step.env <- ordiR2step(abiotic.rda.null, abiotic.rda.all, 
                       direction="forward", perm.max = 1000)

## Summary table and plot:
step.env$anova
plot(step.env)
## Indicates that no abiotic variables should be included in model.

## Run RDA with biotic variable
## Test for significance of the whole model:
mod2 <- rda(pcoa.jtu$vectors ~ biotic$Crucian)
summary(mod2)
coef(mod2)       # canonical coefficients
vif.cca(mod2)    # no collinearity present in model
mod2             # 7% variance is explained by variables in model
anova(mod2)      # Model is not significant
RsquareAdj(mod2) # 1.2% variance explained
plot(mod2)

## Check whether any variables should be dropped from model
drop1(mod2, test="perm")

## Plot model result to get a sense of which variables are correlating with which 
## ponds along which axes
plot(mod2, type='n', scaling=1)
orditorp(mod2, display='sp', cex=0.5, scaling=1, col='blue')
text(mod2, display='cn', col='red')

## Plot of variables against sites
plot(mod2, type = "text")

## Perform hypothesis testing
anova(mod2, perm=1000)               # overall model fit
anova(mod2, by="margin", perm=1000)  # significance of marginal effects (Type III)
anova(mod2, by="term", perm=1000)    # significance of each term in model
anova(mod2, by="axis", perm=1000)    # significance of each axis


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
comb.nomicro.fam.bd.nest <- betadisper(comb.nomicro.fam.dist$beta.jne, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
comb.nomicro.fam.nest <- with(metadata, comb.nomicro.fam.bd.nest)
comb.nomicro.fam.nest

## Compute mean distance to centroid per group
tapply(comb.nomicro.fam.bd.nest$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(comb.nomicro.fam.bd.nest$distances, metadata$Crucian, var)

## Ordination plot of turnover partition
plot(comb.nomicro.fam.bd.nest)

## Boxplot of turnover partition
boxplot(comb.nomicro.fam.bd.nest, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is no difference in nestedness partition
## of beta diversity between crucian and non-crucian ponds. 
## Statistically check whether nestedness is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(comb.nomicro.fam.bd.nest)     # No significant difference between ponds
permutest(comb.nomicro.fam.bd.nest) # No significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(comb.nomicro.fam.bd.nest)  # No significant difference between ponds

## Ordination of beta diversity partitioned by nestedness:
comb.nomicro.fam.comm.nest <- metaMDS(comb.nomicro.fam.dist$beta.jne, 
                                      dist="jaccard", 
                                      k=2,
                                      maxit=999,
                                      trymax=1000,
                                      wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
comb.nomicro.fam.comm.nest$stress
stressplot(comb.nomicro.fam.comm.nest)

## plot site scores as text
ordiplot(comb.nomicro.fam.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
comb.nomicro.NMDS1 <- comb.nomicro.fam.comm.nest$points[,1]
comb.nomicro.NMDS2 <- comb.nomicro.fam.comm.nest$points[,2]
comb.nomicro.fam.nest.NMDS <- data.frame(NMDS1=comb.nomicro.NMDS1, 
                                         NMDS2=comb.nomicro.NMDS2,
                                         Crucian = metadata$Crucian)

## Check data
head(comb.nomicro.fam.nest.NMDS)

## Plot data frame
p36b <- ggplot(comb.nomicro.fam.nest.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", x="NMDS1", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p36b

## Statistically check difference in nestedness of communities
comb.nomicro.fam.nest.anosim <- anosim(comb.nomicro.fam.dist$beta.jne, metadata$Crucian)

## Inspect results
comb.nomicro.fam.nest.anosim
summary(comb.nomicro.fam.nest.anosim)
plot(comb.nomicro.fam.nest.anosim)

## There appears to be no significant difference in nestedness
## Look at PermANOVA using adonis(), considered more robust than anosim()
comb.nomicro.fam.nest.adonis <- adonis(comb.nomicro.fam.dist$beta.jne ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
comb.nomicro.fam.nest.adonis

## Result is not significant. The combined data indicates there is no
## substantial difference in species loss or gain (i.e. nestedness) 
## between ponds with crucian carp and ponds without crucian carp.


## Use RDA and variance partitioning to examine the influence of biotic and 
## abiotic variables on turnover partition of beta diversity

## First, compute principal coordinate decomposition (classical scaling) for 
## turnover distance matrix using the pcoa() function from the ape package. Use 
## lingoes correction to account for negative eigenvalues.
pcoa.jne <- pcoa(comb.nomicro.fam.dist$beta.jne, correction="lingoes", rn=NULL)

## Inspect outputs
pcoa.jne$correction
pcoa.jne$note
pcoa.jne$trace

## The principal component coordinates for each diversity matrix can be written
## to a csv file for inspection
#write.csv(pcoa.jne$vectors, "../Data/Combined_data_family_nestedness_eigenvectors.csv")

## Groups of variables for RDA were already created for species-level analysis
## Perform forward model selection for each group of explanatory variables
## Adjusted R-squared will be used as a goodness-of-fit measure. Model selection
## needs to be performed for each group of variables, excluding those containing
## only one variable.
abiotic.rda.null <- rda(pcoa.jne$vectors ~ 1, abiotic)  # null model
abiotic.rda.all <- rda(pcoa.jne$vectors ~ ., abiotic)   # model with all variables
step.env <- ordiR2step(abiotic.rda.null, abiotic.rda.all, 
                       direction="forward", perm.max = 1000)

## Summary table and plot:
step.env$anova
plot(step.env)
## Indicates that no abiotic variables should be included in model.

## Run RDA with biotic variable
## Test for significance of the whole model:
mod2 <- rda(pcoa.jne$vectors ~ biotic$Crucian)
summary(mod2)
coef(mod2)       # canonical coefficients
vif.cca(mod2)    # no collinearity present in model
mod2             # 14.6% variance is explained by variables in model
anova(mod2)      # Model is not significant
RsquareAdj(mod2) # 9.4% variance explained
plot(mod2)

## Check whether any variables should be dropped from model
drop1(mod2, test="perm")

## Plot model result to get a sense of which variables are correlating with which 
## ponds along which axes
plot(mod2, type='n', scaling=1)
orditorp(mod2, display='sp', cex=0.5, scaling=1, col='blue')
text(mod2, display='cn', col='red')

## Plot of variables against sites
plot(mod2, type = "text")

## Perform hypothesis testing
anova(mod2, perm=1000)               # overall model fit
anova(mod2, by="margin", perm=1000)  # significance of marginal effects (Type III)
anova(mod2, by="term", perm=1000)    # significance of each term in model
anova(mod2, by="axis", perm=1000)    # significance of each axis


## 3. TOTAL BETA DIVERSITY
comb.nomicro.fam.bd.total <- betadisper(comb.nomicro.fam.dist$beta.jac, metadata$Crucian)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(metadata, comb.nomicro.fam.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(comb.nomicro.fam.bd.total$distances, metadata$Crucian, mean)

## Compute variance per group
tapply(comb.nomicro.fam.bd.total$distances, metadata$Crucian, var)

## Ordination plot of total beta diversity
plot(comb.nomicro.fam.bd.total)

## Boxplot of total beta diversity
boxplot(comb.nomicro.fam.bd.total, xlab="Crucian carp", xaxt="n", bty="n")
axis(side=1, at=c(1,2), labels=c("Absent","Present"))

## Plots indicates that there is substantial difference in total beta diversity 
## between crucian and non-crucian ponds. 
## Statistically check whether total beta is different in crucian and
## non-crucian ponds using standard parametric anova or permutation
## tests.
anova(comb.nomicro.fam.bd.total)     # Significant difference between ponds
permutest(comb.nomicro.fam.bd.total) # Significant difference between ponds

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(comb.nomicro.fam.bd.total)  # Significant difference between ponds

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
comb.nomicro.fam.comm.total <- metaMDS(comb.nomicro.fam.dist$beta.jac, 
                                       dist="jaccard", 
                                       k=2,
                                       maxit=999,
                                       trymax=1000,
                                       wascores=TRUE)

## Assess goodness of ordination fit (stress plot)
comb.nomicro.fam.comm.total$stress
stressplot(comb.nomicro.fam.comm.total)

## plot site scores as text
ordiplot(comb.nomicro.fam.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
comb.nomicro.NMDS1 <- comb.nomicro.fam.comm.total$points[,1]
comb.nomicro.NMDS2 <- comb.nomicro.fam.comm.total$points[,2]
comb.nomicro.fam.total.NMDS <- data.frame(NMDS1=comb.nomicro.NMDS1,
                                          NMDS2=comb.nomicro.NMDS2,
                                          Crucian = metadata$Crucian)

## Check data
head(comb.nomicro.fam.total.NMDS)

## Plot data frame
p36c <- ggplot(comb.nomicro.fam.total.NMDS, aes(x=NMDS1, y=NMDS2, colour=Crucian)) + 
  geom_point() + 
  stat_ellipse() + 
  labs(subtitle="", 
       x="NMDS1", y="") + 
  coord_cartesian(xlim=c(-1,1), ylim=c(-1,1)) + 
  scale_colour_manual(name="Crucian carp",
                      breaks=c("N","Y"),
                      labels=c("Absent","Present"),
                      values=c("grey40","deepskyblue2")) + 
  theme(panel.background = element_rect(fill = 'white'),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        plot.title = element_text(face="bold", hjust=0, colour="black"),
        plot.subtitle = element_text(face="bold", hjust=0, colour="black", margin = unit(c(2, 0, 0, 0), "mm")),
        text = element_text(size=20),
        legend.position="none")
p36c

## Statistically check difference in total beta diversity of ponds
comb.nomicro.fam.total.anosim <- anosim(comb.nomicro.fam.dist$beta.jac, metadata$Crucian)

## Inspect results
comb.nomicro.fam.total.anosim
summary(comb.nomicro.fam.total.anosim)
plot(comb.nomicro.fam.total.anosim)

## There appears to be no significant difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
comb.nomicro.fam.total.adonis <- adonis(comb.nomicro.fam.dist$beta.jac ~ Crucian, metadata)

## Inspect results
## no summary() or plot() methods included
comb.nomicro.fam.total.adonis

## Again result is not significant. The combined data indicates there is 
## no substantial variation in overall species composition of
## ponds with crucian carp and ponds without crucian carp


## Use RDA and variance partitioning to examine the influence of biotic and 
## abiotic variables on turnover partition of beta diversity

## First, compute principal coordinate decomposition (classical scaling) for 
## turnover distance matrix using the pcoa() function from the ape package. Use 
## lingoes correction to account for negative eigenvalues.
pcoa.jac <- pcoa(comb.nomicro.fam.dist$beta.jac, correction="lingoes", rn=NULL)

## Inspect outputs
pcoa.jac$correction
pcoa.jac$note
pcoa.jac$trace

## The principal component coordinates for each diversity matrix can be written
## to a csv file for inspection
#write.csv(pcoa.jac$vectors, "../Data/Combined_data_family_totalbeta_eigenvectors.csv")

## Groups of variables for RDA were already created for species-level analysis
## Perform forward model selection for each group of explanatory variables
## Adjusted R-squared will be used as a goodness-of-fit measure. Model selection
## needs to be performed for each group of variables, excluding those containing
## only one variable.
abiotic.rda.null <- rda(pcoa.jac$vectors ~ 1, abiotic)  # null model
abiotic.rda.all <- rda(pcoa.jac$vectors ~ ., abiotic)   # model with all variables
step.env <- ordiR2step(abiotic.rda.null, abiotic.rda.all, 
                       direction="forward", perm.max = 1000)

## Summary table and plot:
step.env$anova
plot(step.env)
## Indicates that no abiotic variables should be included in model.

## Run RDA with biotic variable
## Test for significance of the whole model:
mod2 <- rda(pcoa.jac$vectors ~ biotic$Crucian)
summary(mod2)
coef(mod2)       # canonical coefficients
vif.cca(mod2)    # no collinearity present in model
mod2             # 7.9% variance is explained by variables in model
anova(mod2)      # Model is significant
RsquareAdj(mod2) # 2.2% variance explained
plot(mod2)

## Check whether any variables should be dropped from model
drop1(mod2, test="perm")
## Crucian carp presence-absence should not be dropped as it is significant

## Plot model result to get a sense of which variables are correlating with which 
## ponds along which axes
plot(mod2, type='n', scaling=1)
orditorp(mod2, display='sp', cex=0.5, scaling=1, col='blue')
text(mod2, display='cn', col='red')

## Plot of variables against sites
plot(mod2, type = "text")

## Perform hypothesis testing
anova(mod2, perm=1000)               # overall model fit
anova(mod2, by="margin", perm=1000)  # significance of marginal effects (Type III)
anova(mod2, by="term", perm=1000)    # significance of each term in model
anova(mod2, by="axis", perm=1000)    # significance of each axis



#################
# SUMMARY PLOTS #
#################

## Difference in alpha diversity of ponds between standard methods
## and unbiased methods
g1 <- ggarrange(m1,m5,m3,m7, ncol=2, nrow=2)
ggsave(filename="Figures/Fig2.png", 
       plot=g1, width=12, height=10, dpi=300, units="in")

## Difference in beta diversity of ponds between standard methods
g2 <- ggarrange(m2a,m2b,m2c,
                m4a,m4b,m4c,
                m6a,m6b,m6c,
                m8a,m8b,m8c,
                ncol=3, nrow=4,
                common.legend=TRUE,
                legend="bottom")
ggsave(filename="Figures/Fig3.png", 
       plot=g2, width=18, height=18, dpi=300, units="in")

## Difference in alpha diversity of ponds in relation to crucian carp
## by standard methods individually and combined
g3 <- ggarrange(p1,p3,p5,p7,p9,p11,p13,p16,
                ncol=2, nrow=4,
                common.legend=TRUE,
                legend="bottom")
ggsave(filename="Figures/Fig4.png", 
       plot=g3, width=12, height=15, dpi=300, units="in")

## Difference in beta diversity of ponds in relation to crucian carp
## by standard methods individually and combined at species-level
g4 <- ggarrange(p2a,p2b,p2c,
                p6a,p6b,p6c,
                p10a,p10b,p10c,
                p15a,p15b,p15c,
                nrow=4, ncol=3,
                common.legend=TRUE,
                legend="bottom")
ggsave(filename="Figures/Fig5.png", 
       plot=g4, width=15, height=15, dpi=300, units="in")

## Difference in alpha diversity of major invertebrate groups in ponds
## in relation to crucian carp by standard methods
g5 <- ggarrange(p14a,p14b,p17a,p17b, 
                ncol=1, nrow=4)
ggsave(filename="Figures/FigS11.png", 
       plot=g5, width=17, height=21, dpi=300, units="in")

## Difference in beta diversity of ponds in relation to crucian carp
## by standard methods individually and combined at family-level
g6 <- ggarrange(p4a,p4b,p4c,
                p8a,p8b,p8c,
                p12a,p12b,p12c,
                p18a,p18b,p18c,
                nrow=4, ncol=3,
                common.legend=TRUE,
                legend="bottom")
ggsave(filename="Figures/FigS12.png", 
       plot=g6, width=15, height=15, dpi=300, units="in")

## Difference in alpha diversity of ponds in relation to crucian carp
## by unbiased methods individually and combined
g7 <- ggarrange(p19,p21,p23,p25,p27,p29,p31,p34,
                ncol=2, nrow=4,
                common.legend=TRUE,
                legend="bottom")
ggsave(filename="Figures/FigS13.png", 
       plot=g7, width=12, height=15, dpi=300, units="in")

## Difference in alpha diversity of major invertebrate groups in ponds
## in relation to crucian carp by unbiased methods
g8 <- ggarrange(p32a,p32b,p35a,p35b, 
                ncol=1, nrow=4)
ggsave(filename="Figures/FigS14.png", 
       plot=g8, width=15, height=21, dpi=300, units="in")

## Difference in beta diversity of ponds in relation to crucian carp
## by unbiased methods individually and combined at species-level
g9 <- ggarrange(p20a,p20b,p20c,
                p24a,p24b,p24c,
                p28a,p28b,p28c,
                p33a,p33b,p33c,
                nrow=4, ncol=3,
                common.legend=TRUE,
                legend="bottom")
ggsave(filename="Figures/FigS15.png", 
       plot=g9, width=15, height=15, dpi=300, units="in")

## Difference in beta diversity of ponds in relation to crucian carp
## by unbiased methods individually and combined at family-level
g10 <- ggarrange(p22a,p22b,p22c,
                 p26a,p26b,p26c,
                 p30a,p30b,p30c,
                 p36a,p36b,p36c,
                 nrow=4, ncol=3,
                 common.legend=TRUE,
                 legend="bottom")
ggsave(filename="Figures/FigS16.png", 
       plot=g10, width=15, height=15, dpi=300, units="in")




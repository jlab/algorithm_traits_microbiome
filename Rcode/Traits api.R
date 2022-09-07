#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
## OPTION 1 - USING THE R-PACKAGE BACDIVE
#+++++++++++++++ BACDIVE API WEB SERVICE +++++++++++++++++++++++++++++++++#
#+
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
## AIM: COMPARE MY ACTUAL TAXA DATAFRAME WITH ALL TAXA IN BACDIVE;
## IF ANY NEWLY ADDED SPECIE ,NA,SORT IT OUT;
## AT THE END CREATE A NEW DATABASE WITH ALL SPECIES.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

## Call up all strain-type Species in BACDIVE
library(readxl)
advsearch_bacdive_2022_06_02 <- read_excel("C:/Users/Bew/Dropbox/PC/Downloads/advsearch_bacdive_2022-06-02.xlsx")
#View(advsearch_bacdive_2022_06_02)

## Rename database
total_strain_bacdive <- advsearch_bacdive_2022_06_02

# compare species column of two databases
new_taxa_bacdive <- total_strain_bacdive$species[!total_strain_bacdive$species %in% search_taxon$search_string]
new_taxa_bacdive <- as.data.frame(new_taxa_bacdive)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++ RECONSTRUCTING TRAIT TABLE +++++++++++++++++++#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#
#curated_trait_data <- readRDS("curated_trait_data.rds")
curated_trait_data <- read.csv("C:/Users/Anwender/Downloads/curated_trait_data.csv")

# build search string 
search_taxon <- curated_trait_data %>%
  transmute(GENUS = Genus, SPECIES = Species) %>% 
  select(GENUS,SPECIES)%>%
  mutate(search_string = paste0(GENUS ," ",SPECIES))

## Registration for BacDive is required but free and easy to accomplish.
## In real applications username and password could of course also be stored
## within the R code itself,or read from a file.

##+++++++++++++++++++++++++ BACDIVE +++++++++++++++++++++++++++++++++++++#
################################################################################
# Run custom R file to activate bacdive open access object
source("C:/Users/Anwender/Documents/bacdive_custom_func.R") # run custom func script to activate,
                                                            # refresh the bacdive access object.
# create an empty list for bacID                                                           
# ID_list <- list()

# Sending inquiries to bacdive - request func
# The extraction of the bacdive ids was done in 2 runs.
t1 <- Sys.time()
#ID_list <- list()
for(taxon in search_taxon[1277:10906,3]){
  result <- request(bacdive,taxon)
  ID_list <- append(ID_list,result$results)
}
t2<- Sys.time()
t2-t1

# convert large list of bacdive ids to a dataframe
bacid <- t(as.data.frame(ID_list))
rownames(bacid) <- NULL

##++++++++++++++++++++ Extracting the traits from Bacdive ++++++++++++++++++++#
t3 <- Sys.time()
traits_tbl <-lapply(bacid[1:38611,1],bacdat)
Traits <- do.call(rbind,traits_tbl)
t4 <- Sys.time()
t4-t3 # Laufzeit des Prozessen

# Add the bacid dataframe to Trait table
Traits$ID <- bacid[,1]

# reorder colnames in dataframe
Traits <- select(Traits, 12,1,2,3,4,5,6,7,8,9,10,11)

# clean aggregration score column
Traits$aggregation_score <- ifelse(Traits$aggregation_score %in% c("aggregates in clumps","aggregates in chains"),
                                   Traits$aggregation_score,NA)
refs <- list(
  Oxygen_tolerance = c('obligate aerobe'      = 5,
                       'aerobe'               = 5,
                       'facultative aerobe'   = 4,
                       'aerotolerant'         = 4,
                       'microaerophile'       = 3,
                       'microaerotolerant'    = 3,
                       'facultative anaerobe' = 2,
                       'anaerobe'             = 1,
                       'obligate anaerobe'    = 1),
  Gram_positive = c('positive' = 1,
                    'variable' = 0.5,
                    'negative' = 0),
  Spore = c('yes'  = 1,
            'no' = 0),
  Motility = c('yes'  = 1,
               'no' = 0),
  Aggregation_score = c('aggregates in chains'  = 1,
                        'aggregates in clumps' = 1)
)

Traits <- Traits %>%
  mutate_all(funs(ifelse(. == '', NA, .))) %>%
  mutate(
    oxygen_tolerance = refs$Oxygen_tolerance[match(oxygen_tolerance, names(refs$Oxygen_tolerance))],
    gram = refs$Gram_positive[match(gram, names(refs$Gram_positive))],
    spore = refs$Spore[match(spore, names(refs$Spore))],
    motility = refs$Motility[match(motility, names(refs$Motility))],
    aggregation_score = refs$Aggregation_score[match(aggregation_score, names(refs$Aggregation_score))],
    length = ifelse(grepl('<|>', length), NA, length),
    length = gsub("µm","",length),
    length = ifelse(grepl('-', length), 
                    rowMeans(cbind(as.numeric(sub('-.*', '', length)),
                                   as.numeric(sub('.*-', '', length)))), length),
    width = ifelse(grepl('<|>', width), NA, width),
    width = gsub("µm","",width),
    width = ifelse(grepl('-', width), rowMeans(cbind(as.numeric(sub('-.*', '', width)),
                                                     as.numeric(sub('.*-', '', width)))), width),
    GC_content = ifelse(grepl('-', GC_content), rowMeans(cbind(as.numeric(sub('-.*', '', GC_content)),
                                                          as.numeric(sub('.*-', '', GC_content))),na.rm = T), GC_content),
    GC_content = ifelse(grepl('±',GC_content),sapply(sapply(GC_content, strsplit, '±'), function(x) x[[1]]),GC_content)
  )
# rename & restructure trait db
bac <- Traits %>% 
  transmute(Genus = Genus,
            Species = Species,
            Oxygen_tolerance = oxygen_tolerance,
            Gram_positive = gram,
            Spore = spore,
            Motility = motility,
            Aggregation_score = aggregation_score,
            Length = length,
            Width = width,
            GC_content)

#convert to narrow table
bac <- bac %>%
  gather(trait, val, -Genus, -Species) %>%
  filter(!is.na(val)) %>%
  mutate(val = as.numeric(val))

# split the specie column into Genus and Species 
bac <- separate(bac, Species, c("Genus","Species"),fill = "left") %>% 
  select(2,3,1,4)

################################################################################
## OPTIONAL: Retreive bacdive ID
#+ Getting the number of strains per species - using request_bac func
#ID_test <- list()
#t3 <- Sys.time()
#for (taxon in search_taxon[1:10906,3]){
  #id <- request_bac(taxon)
 #ID_test <- append(ID_test,id)}
#t4<- Sys.time()
#t4-t3

###############################################################################
##++++++++++++++++++++++++ NCBI +++++++++++++++++++++++++++++++++++++++++++++##
# get Genome size, GC content, and Gene Number from NCBI ftp site.
# https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt
# download procaryote traits from link above.

library(tidyverse)
library(data.table)

genos <- fread("C:/Users/Anwender/Documents/prokaryotes.txt")
genos <- genos[order(genos$`#Organism/Name`),]
genos <- genos %>% 
  filter(Status == 'Complete Genome') %>%
  transmute(
    sp = gsub("'|\\[|\\]", "", `#Organism/Name`),
    Genus = ifelse(grepl("Candidatus",sp),sapply(sapply(sp, strsplit, ' '), function(x) paste0(x[1]," ",x[2])),
                   sapply(sapply(sp, strsplit, ' '), function(x) x[1])),
    Species = ifelse(grepl("Candidatus",sp),sapply(sapply(sp, strsplit, ' '), function(x) x[3]),sapply(sapply(sp, strsplit, ' '), function(x) x[2])),
    Genome_Mb = as.numeric(ifelse(`Size (Mb)` == '-', NA, `Size (Mb)`)),
    GC_content = `GC%`, 
    Gene_number = as.numeric(ifelse(Genes == '-', NA, Genes))) %>%
  mutate(Species = gsub("endosymbiont","", Species)) %>% 
  gather(trait, val, Genome_Mb, GC_content, Gene_number) %>%
  select(-sp) %>%
  filter(!(is.na(val)|is.na(Species)|is.na(Genus)))

# convert val column in genos to numeric - before merging to genos1 db
genos$val <- as.numeric(genos$val)

# Downloaded another NCBI database from the link below
#https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/

genos1 <- fread("C:/Users/Anwender/Downloads/prokaryotes.csv")
genos1 <- genos1[order(genos1$`#Organism Name`),]
genos1 <- genos1[-c(1:2,30083:30121),]%>% 
  transmute(
    sp = gsub("'|\\[|\\]", "", `#Organism Name`),
    Genus = ifelse(grepl("Candidatus",sp),sapply(sapply(sp, strsplit, ' '), function(x) paste0(x[1]," ",x[2])),
                                                 sapply(sapply(sp, strsplit, ' '), function(x) x[1])),
    Species = ifelse(grepl("Candidatus",sp),sapply(sapply(sp, strsplit, ' '), function(x) x[3]),sapply(sapply(sp, strsplit, ' '), function(x) x[2])),
    Genome_Mb = as.numeric(ifelse(`Size(Mb)` == '-', NA, `Size(Mb)`)),
    GC_content = `GC%`) %>% 
  mutate(Species = gsub("endosymbiont","", Species)) %>% 
  gather(trait, val, Genome_Mb, GC_content) %>%
  select(-sp) %>%
  filter(!(is.na(val))|is.na(Species)|is.na(Genus))

## Still downloaded from ncbi
# also found this ncbi database from the website
#https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/
genos2 <- fread('C:/Users/Anwender/Downloads/genomes.csv') 

# Select just Bacteria using the Organism groups column
genos2 <- subset(genos2, grepl("Bacteria.*",genos2$`Organism Groups`))

# Descending order using genus names
genos2 <- genos2[order(genos2$`#Organism Name`),]

genos2 <- genos2 %>% 
  transmute(sp = `#Organism Name`, val = `Size(Mb)`) %>%
  mutate(sp = gsub(" bacterium.*|\\[|\\]", '', sp),
         sp = gsub("'.*", "", sub("'?", "", sp)),
         sp = gsub("Type-.*|unicellular.*|unidentified.*|uncultured.*|SAR.*|SCandidatus.*|Vaccinium.*|endosymbiont.*|subdivision.*",
                   "",sp)) %>%
  filter(!grepl('\\d', sp)) %>%
  mutate(Genus = ifelse(grepl("Candidatus",sp),sapply(sapply(sp, strsplit, ' '), function(x) paste0(x[1]," ",x[2])),
                 sapply(sapply(sp, strsplit, ' '), function(x) x[1])),
        Species = ifelse(grepl("Candidatus",sp),sapply(sapply(sp, strsplit, ' '), 
                   function(x) x[3]),sapply(sapply(sp, strsplit, ' '), function(x) x[2])))%>%
  transmute(Genus, Species, trait = 'Genome_Mb', val)
genos2 <- filter(genos2, !(is.na(Genus) | is.na(Species)| Species == "x"))

# merge all 03 databases
genos <- bind_rows(genos,genos1,genos2)
# remove from the Genus column - Candidatus
genos <- genos %>% mutate(Genus = gsub("Candidatus","",Genus))

#  mutate(sp = ifelse(str_count(sp, "\\S+") > 1, word(sp, 1, 2), sp)) %>%
#separate(sp, c("Genus", "Species"), " ", fill = 'right')
##########################################################################
#+++++++++++++++++ GOLD DATABASE ++++++++++++++++++++++++++++++++++++++++#
##########################################################################
library(readxl)
library(tidyverse)
gold_DB <- read_excel('C:/Users/Anwender/Downloads/gold_DB.xlsx')

# Keep just bacteria and rename column names of DB
gold_DB <- gold_DB %>% 
  filter(`ORGANISM NCBI SUPERKINGDOM`=="Bacteria") %>% 
  transmute(Species = `ORGANISM NCBI SPECIES`,
            Width = `ORGANISM CELL DIAMETER`,
            Length = `ORGANISM CELL LENGTH`,
            Gram_positive = `ORGANISM GRAM STAIN`,
            Motility = `ORGANISM MOTILITY`,
            Oxygen_tolerance = `ORGANISM OXYGEN REQUIREMENT`,
            pH_optimum = `ORGANISM PH`,
            Spore = `ORGANISM SPORULATION`)
# Fix taxonomy and filter out unclassified species
gold_DB <- gold_DB %>% 
  mutate(Species = gsub("\\[|\\]|'", "",Species),
        Species = gsub(' bacterium.*|uncultured.*|candidate.*', '',Species)) %>%
  mutate(Genus = ifelse(grepl("Candidatus",Species),sapply(sapply(Species, strsplit, ' '), function(x) paste0(x[1]," ",x[2])),
                        sapply(sapply(Species, strsplit, ' '), function(x) x[1])),
         Species = ifelse(grepl("Candidatus",Species),sapply(sapply(Species, strsplit, ' '), 
                                                        function(x) x[3]),sapply(sapply(Species, strsplit, ' '), function(x) x[2])))%>%
  mutate(Species = gsub("sp.?|endosymbiont", "", Species)) %>% 
  filter(!(is.na(Genus)|is.na(Species)))

# remove empty spaces in species column
gold_DB <- gold_DB %>%  filter(Species!= "")

# Fix the width - take averages of ranges when necessary.
gold_DB <- gold_DB %>%
  mutate(
    Width = ifelse(Width == '803nm', 0.803, Width),
    Width = ifelse(Width == '500 nm', 0.5, Width),
    Width = ifelse(Width == 'less than 1 ?min diameter', 0.9, Width),
    Width = ifelse(grepl(".*nm",Width),gsub("nm","mm",Width),Width),
    Width = sapply(Width,mm_to_micrometre_width),
    Width = gsub("\\s","", Width),
    Width = gsub("\\+/|\\??|¿¿", "",Width),
    Width = ifelse(Width == "06-08", 0.6-0.8, Width),
    Width = sapply(Width,function(Width) ifelse(grepl("\\d+\\.?\\d?.?[to|?|¿].?\\d+\\.?\\d?.?",Width),
                   gsub("to|\\?|\\¿","-",Width),Width)),
    Width = gsub("¿|[[:alpha:]]","",Width),
    Width = gsub("4.9±2.6","4.9", Width)) %>%
  separate(Width, c("Width0","Width1"), fill = 'right', sep = '[-|–| ]')%>%
  mutate(Width = ifelse(is.na(Width0), Width0, 
                        ifelse(is.na(Width1), as.numeric(Width0),
                               (as.numeric(Width1) + as.numeric(Width0)) / 2))) %>%
  mutate(Width = as.numeric(Width)) %>%
  select(-Width0, -Width1)

# Fix the length 
gold_DB <- gold_DB %>%
  mutate(
    Length = gsub('variable', NA, Length),
    Length = gsub("I .8-2.5 μm", "0.8-2.5", Length),
    Length = gsub(".65um", "0.65", Length),
    Length = gsub("[()]","",Length),
    Length = ifelse(grepl(".*nm",Length),gsub("nm","mm",Length),Length),
    Length = sapply(Length,mm_to_micrometre_length),
    Length = gsub("\\s","",Length),
    Length = gsub("\\+/|\\??|¿¿", "",Length),
    Length = sapply(Length,function(Length) ifelse(grepl("\\d+\\.?\\d?.?[to|?|¿].?\\d+\\.?\\d?.?",Length),
                                                gsub("to|\\?|\\¿","-",Length),Length)),
    Length = gsub("¿|[[:alpha:]]","",Length),
    Length = gsub("1.5-.5.0","1.5-5.0",Length),
    Length = gsub("1.42.0","1.4-2.0",Length),
    Length = gsub("1.22.4","1.2-2.4",Length),
    Length = gsub("1.52.5","1.5-2.5",Length),
    Length = gsub("1.55.0","1.5-5.0",Length)) %>%
  separate(Length, c("Length0","Length1"), fill = 'right', sep = '[-|–| ]')%>%
  mutate(Length = ifelse(is.na(Length0), Length0, 
                        ifelse(is.na(Length1), as.numeric(Length0),
                               (as.numeric(Length1) + as.numeric(Length0)) / 2))) %>%
  mutate(Length = as.numeric(Length)) %>%
  select(-Length0, -Length1)

# gram positive, motility, oxygen tolerance, sporulation
gold_DB <- gold_DB %>%
  mutate(
    Gram_positive = c(1,0)[match(Gram_positive, c('Gram+','Gram-'))],
    Motility = c(1,1,0,0)[match(Motility, c('Motile','Chemotactic','Non-motile','Nonmotile'))],
    Oxygen_tolerance = c(5,5,4,3,2,1,1)[match(Oxygen_tolerance, 
                                              c('Obligate aerobe','Aerobe','Microaerophilic','Facultative','Facultative anaerobe','Anaerobe','Obligate anaerobe'))],
    Spore = c(1,1,0)[match(Spore, c('Non-sporulating','Nonsporulating','Sporulating'))])

#ph optimum
gold_DB <- gold_DB %>%
  mutate(
    pH_optimum = gsub(" |~", "", pH_optimum),
    pH_optimum = ifelse(pH_optimum %in% c('acido-sensible','Notknown',"5557296000"), NA, pH_optimum)) %>%
  separate(pH_optimum, c('pH0','pH1'), fill = 'right', sep = '[-|-]') %>%
  mutate(pH_optimum = ifelse(is.na(pH0), pH0, 
                             ifelse(is.na(pH1), as.numeric(pH0),
                                    (as.numeric(pH1) + as.numeric(pH0)) / 2))) %>%
  mutate(pH_optimum = as.numeric(pH_optimum)) %>%
  select(-pH0, -pH1)

jgi <- gold_DB %>%
  gather(trait, val, -Genus, -Species) %>%
  filter(!is.na(val)) %>% 
  select(Genus,Species,trait,val) 
 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
##################################################################################
#+++++++++++++++++++++++++++++++++++++++ IJSEM ++++++++++++++++++++++++++++++++++#
#read ijsem table
ijsem<-read.delim('C:/Users/Anwender/Downloads/IJSEM_pheno_db_v1.0.txt', sep="\t", header=T, check.names=F, fill=T,
                  na.strings=c("NA", "", "Not indicated", " Not indicated","not indicated", "Not Indicated", "n/a", "N/A", "Na", "Not given", "not given","Not given for yeasts", "not indicated, available in the online version", "Not indicated for yeasts", "Not Stated", "Not described for yeasts", "Not determined", "Not determined for yeasts"))

#simplify column names
colnames(ijsem)<-c("Habitat", "Year", "DOI", "rRNA16S", "GC", "Oxygen",
                   "Length", "Width", "Motility", "Spore", "MetabAssays", "Genus", "Species", "Strain", "pH_optimum", "pH_range", "Temp_optimum", "Temp_range", "Salt_optimum", "Salt_range", "Pigment", "Shape", "Aggregation", "FirstPage", "CultureCollection", "CarbonSubstrate", "Genome", "Gram", "Subhabitat", "Biolog")

#clean Habitat column
levels(ijsem$Habitat)[levels(ijsem$Habitat)=="freshwater (river, lake, pond)"]<-"freshwater"
levels(ijsem$Habitat)[levels(ijsem$Habitat)=="freshwater sediment (river, lake, pond"]<-"freshwater sediment"

#clean Oxygen column
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="aerobic"]<-"obligate aerobe"
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="anerobic"]<-"obligate anerobe"
levels(ijsem$Oxygen)[levels(ijsem$Oxygen)=="microerophile"]<-"microaerophile"

#clean pH_optimum column
#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
ijsem$pH_optimum<-as.character(ijsem$pH_optimum)
ijsem$pH_optimum<-sapply(ijsem$pH_optimum, simplify=T, function(x){mean(swan(unlist(strsplit(x, split="-", fixed=T))), na.rm = T)})

#remove pH values <0 and >10
ijsem$pH_optimum[ijsem$pH_optimum<0 | ijsem$pH_optimum>10]<-NA

#clean Temp_optimum column
ijsem$Temp_optimum<-as.character(ijsem$Temp_optimum)

#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
ijsem$Temp_optimum<-sapply(ijsem$Temp_optimum, simplify=T, function(x){mean(swan(unlist(strsplit(x, split="-", fixed=T))))})

#clean Salt_optimum column
ijsem$Salt_optimum<-as.character(ijsem$Salt_optimum)

#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
ijsem$Salt_optimum<-sapply(ijsem$Salt_optimum, simplify=T, function(x){mean(swan(unlist(strsplit(x, split="-", fixed=T))))})

#there are some formatting issues that should be solved
############## Now my additional edits to the IJSEM data
#Assign oxygen score
ijsem$Oxygen_score <- c(5,4,3,2,1)[match(ijsem$Oxygen, 
                                         c('obligate aerobe','microaerophile','facultative aerobe','facultative anaerobe','anaerobic'))]

#turn Aggregation into a binary
ijsem$Aggregation <- c(0,1,1)[match(ijsem$Aggregation, c('none','chain','clump'))]

#clean Length and Width
#this step splits the range values and takes the mean value
#values that are not numeric are transformed to NAs
j <- as.character(ijsem$Length)
j <- gsub(' in diameter| \\(diameter\\)|mm|Indicated|diameter: | ', '', j)
j <- gsub(',', '.', j, fixed = TRUE)
j[grepl(">", j)] <- NA
j[j == ''] <- NA
j[j == c('0.61.6')] <- '0.6-1.6'
j[j == c('0.4 0.9 ')] <- '0.4-0.9'
j <- gsub('äóñ *', '-', j)

j <- ifelse(grepl("-[a-zA-Z]+", j) & !is.na(j), 
            as.character(format(as.Date(j, "%d-%b"), "%d.%m")), 
            j)
j <- ifelse(grepl("[a-zA-Z]+-", j) & !is.na(j), 
            as.character(format(as.Date(paste(1, j), "%d %b-%y"), "%m.%y")), 
            j)

#if it has a dash, assume it was a range and take the mean
filt <- grepl("-", j) & !is.na(j)
j[filt] <- rowMeans(cbind(
  as.numeric(sapply(strsplit(j[filt], '-'), function(x) x[[1]])),
  as.numeric(sapply(strsplit(j[filt], '-'), function(x) x[[2]]))))

#if it has a backslash, assume it was a decimal point
filt <- grepl("\\/", j) & !is.na(j)
j[filt] <- paste(sapply(strsplit(j[filt], '\\/'), function(x) x[[1]]),
                 sapply(strsplit(j[filt], '\\/'), function(x) x[[2]]), sep = '.')

ijsem$Length <- swan(j)


#now width
j <- as.character(ijsem$Width)
j <- gsub(' |\\(|\\)|not.+|indiameter|Filament.+', '', j)
j <- gsub(',', '.', j, fixed = TRUE)
j[grepl("<", j)] <- NA
j[j == ''] <- NA
j <- gsub('äóñ *', '-', j)
j <- ifelse(grepl("-[a-zA-Z]+", j) & !is.na(j), 
            as.character(format(as.Date(j, "%d-%b"), "%d.%m")), 
            j)
j <- ifelse(grepl("[a-zA-Z]+-", j) & !is.na(j), 
            as.character(format(as.Date(paste(1, j), "%d %b-%y"), "%m.%y")), 
            j)

#if it has a dash, assume it was a range and take the mean
filt <- grepl("-", j) & !is.na(j)
j[filt] <- rowMeans(cbind(
  as.numeric(sapply(strsplit(j[filt], '-'), function(x) x[[1]])),
  as.numeric(sapply(strsplit(j[filt], '-'), function(x) x[[2]]))))

#if it has a backsclash, assume it was a decimal point
filt <- grepl("\\/", j) & !is.na(j)
j[filt] <- paste(sapply(strsplit(j[filt], '\\/'), function(x) x[[1]]),
                 sapply(strsplit(j[filt], '\\/'), function(x) x[[2]]), sep = '.')

ijsem$Width <- as.numeric(j)

#Turn motility into a binary
ijsem$Motility <- c(0,1,1,1,1)[match(ijsem$Motility, c('non-motile', 'flagella',
                                                       'motile, but unspecified structure', 'gliding', 'axial filament'))]

#turn spore-forming into binary
ijsem$Spore <- c(0,1)[match(ijsem$Spore, c('no','yes'))]

#turn gram status into a binary
ijsem$Gram <- c(0, 0.5, 1)[match(ijsem$Gram, c('negative','variable','positive'))]

# fix pH
j <- as.character(ijsem$pH_range)
j <- gsub(' |\\(|\\)|at 6\\.0 but not 5\\.5|>|<|not inidcated', '', j)
j <- gsub('MinimumpHof5|t|4\\.0\\+|5\\.0\\+;notat4\\.5|14\\.08', '', j)
j[j %in% c('13.05','4.3','7.1','7.2','5.5','7.5','30.04','5.1')] <- NA
j[j == '5.7and6.8'] <- '5.7-6.8'
j[j == '4.0and8.5'] <- '4-8.5'
j[j == '5,6,10,12,'] <- '5-12'
j[j == '5and9.5'] <- '5-9.5'
j[j == '5.59.5'] <- '5.5-9.5'

# remove anything with an inequality, or other special characters
j[grepl("<|>", j)] <- NA
j[j == ''] <- NA
j <- gsub('äóñ *', '-', j)
j <- ifelse(grepl("-[a-zA-Z]+", j) & !is.na(j), 
            as.character(format(as.Date(j, "%d-%b"), "%d.%m")), 
            j)

filt <- grepl("\\/", j) & !is.na(j)
j[filt] <- paste(sapply(strsplit(j[filt], '\\/'), function(x) x[[1]]),
                 sapply(strsplit(j[filt], '\\/'), function(x) x[[2]]), sep = '.')

j[!grepl("-", j) & !is.na(j) & !grepl("\\.", j)] <- NA
filt <- !grepl("-", j) & !is.na(j) & grepl("\\.", j)
j[filt] <- gsub("\\.", "-", j[filt])
lows <- swan(sapply(strsplit(j, '-'), function(x) x[[1]]))
highs <- swan(sapply(strsplit(j, '-'), function(x) if (length(x) > 1) x[[2]] else NA))

ijsem$pH_lows <- lows
ijsem$pH_highs <- highs

#if there isn't a pH optimum, but there are ph high and lows, we take the midpoint (same as done in Barberan 2016)
ijsem$pH_optimum <- ifelse(is.na(ijsem$pH_optimum), (lows + highs) / 2, ijsem$pH_optimum)

####GC
j <- as.character(ijsem$GC)
j <- gsub('äóñ *|\xe4\xf3\xf1', '-', j)
j <- gsub(',', '\\.', j)
j[j %in% c('DQ987877','marine organism','BAOL01000001','40.50%','GU323338')] <- NA
j <- gsub('O±.+|\\+/-.+|%|o|\x8c\xb10.5|\x8c\xb10.4', '', j)
filt <- grepl("-", j) & !is.na(j)
j[filt] <- rowMeans(cbind(
  as.numeric(sapply(strsplit(j[filt], '-'), function(x) x[[1]])),
  as.numeric(sapply(strsplit(j[filt], '-'), function(x) x[[2]]))))

ijsem$GC <- swan(j)

# select traits for analysis...
ijsem <- select(ijsem, Genus, Species, GC_content = GC, 
                Oxygen_tolerance = Oxygen_score, Length, Width, 
                Motility = Motility, Spore, 
                pH_optimum, Temp_optimum, Salt_optimum, Aggregation_score = Aggregation, 
                Gram_positive = Gram) %>%
  gather(trait, val, -Genus, -Species) %>%
  filter(!is.na(val)) %>%
  filter(!Genus %in% c("10.1099/ijs.0.006320-0","10.1099/ijs.0.02424-0")) %>%
  mutate_if(is.factor, as.character)

# Fix the Species column - checking the species names
ijsem <- ijsem %>% 
  mutate(Species = trimws(Species),
         Species = ifelse(str_count(Species,"\\w+") > 1,sapply(sapply(Species, strsplit, '\\s+'), function(x) x[2]),Species))

###########################################################
# now, sporulation data from Browne et al 2016
# https://www.nature.com/articles/nature17645#supplementary-information 
# [data edited in excel for easier processing)

spo <- read.csv('C:/Users/Anwender/Documents/microbiome_trait_succession/data/Browne2016_sporulationTable.csv', header = T, skip = 1)

# fix long colnames
colnames(spo) <- substr(colnames(spo), 1, 2)

spo <- spo %>%
  transmute(
    Genus = unlist(lapply(strsplit(as.character(Sp), " |_"), function(x) x[[1]])),
    Species = unlist(lapply(strsplit(as.character(Sp), " |_"), function(x) x[[2]])),
    trait = 'Spore_score',
    val = si)

################ extract IgA data from Palm et al 2014
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4174347/
#i modified it a bit in excel before processing

iga <- read.csv('C:/Users/Anwender/Documents/microbiome_trait_succession/data/Palm2014_IGA.csv', stringsAsFactors = FALSE)
iga <- t(iga)
colnames(iga) <- iga[1, ]
iga <- as.data.frame(iga[c(2:nrow(iga)), ])
iga[, c(4:ncol(iga))] <- sapply(iga[, c(4:ncol(iga))], function(x) as.numeric(as.character(x)))

iga <- iga %>%
  gather(tax, val, -var, -status, -subj) %>%
  spread(var, val) %>%
  mutate(tax = gsub("\\s|.__|\\[|\\]|Other", "", tax)) %>%
  separate(tax, sep = ';', c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  filter(!(Genus == 'unclassified' | Genus == '')) %>%
  mutate(Species = ifelse(Species == '', 'unclassified', Species)) %>%
  group_by(Genus, Species, trait = 'IgA') %>%
  summarise(val = log(mean(ici) + 1)) %>%
  as.data.frame()

#################
#extract number of b-vitamins in genomes from Magnusdottir 2015
#https://www.frontiersin.org/articles/10.3389/fgene.2015.00148/full
bvit <- read.csv('C:/Users/Anwender/Documents/microbiome_trait_succession/data/Bvitamins_Magnusdottir2015.csv', stringsAsFactors = FALSE) %>%
  mutate(Genus = sapply(sapply(tax, strsplit, ' '), function(x) x[[1]]),
         Species = sapply(sapply(tax, strsplit, ' '), function(x) x[[2]])) %>%
  mutate(Species = ifelse(Species == 'sp.', 'unclassified', Species)) %>%
  select(-NCBI.Taxonomy.ID, -tax, -Body.Site) %>%
  gather(Bvit, val, -Genus, -Species) %>%
  group_by(Genus, Species, trait = 'B_vitamins') %>%
  summarise(val = length(unique(Bvit[val > 0]))) %>%
  as.data.frame()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++ RRNDB (Ribosomal RNA operon copy number DATABASE) +++++++++++++

rrnDB.5.7 <- read.delim("C:/Users/Anwender/Documents/rrnDB-5.7.tsv") %>% 
  mutate(
    sp = gsub("'|\\[|\\]", "", NCBI.scientific.name),
    Genus = sapply(sapply(sp, strsplit, ' '), function(x) x[[1]]),
    Species = sapply(sapply(sp, strsplit, ' '), function(x) ifelse(length(x) > 1, x[[2]], NA)),
    trait = 'Copies_16S') %>%
  transmute(Genus, Species, trait, val = X16S.gene.count) %>%
  filter(!is.na(val)) %>%
  mutate(Species = ifelse(is.na(Species), 'unclassified', Species))

# put it all together
x <- bind_rows(
  mutate(ijsem, source = 'IJSEM'),
  mutate(spo, source = 'Browne2016'),
  mutate(bac, source = 'BacDive'),
  mutate(genos, source = 'NCBI'),
  mutate(rrnDB.5.7, source = 'rrnDB'),
  mutate(iga, source = 'Palm2014'),
  mutate(bvit, source = 'Mag2015'),
  mutate(jgi, source = 'JGI')
)

#manual entry
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3346390/
x$val[x$Genus == 'Bifidobacterium' & x$Species != 'unclassified' & x$trait == 'Aggregation_score'] <- 1 


x_unclean <- nrow(x)
### I looked at all the points plotted using 
c <- ggplot(x, aes(y = val)) + geom_boxplot() + facet_wrap(~trait, scales = 'free')
# identifying and removing outliers
# remove inordinately large values...

x <- x %>% 
  filter(Genus!="")
x <- x %>%
  filter(!(trait == 'Copies_16S' & val > 20.0)) %>%
  filter(!(trait == 'GC_content' & (val > 80.0 | val < 20.0))) %>%
  filter(!(trait == 'Gene_number' & val > 11000.0)) %>%
  filter(!(trait == 'Genome_Mb' & val > 14.0)) %>%
  filter(!(trait == 'Length' & val > 30.0)) %>%
  filter(!(trait == 'pH_optimum' & val < 2.5)) %>%
  filter(!(trait == 'Salt_optimum' & val > 25.0)) %>%
  filter(!(trait == 'Temp_optimum' & val > 80.0)) %>%
  filter(!(trait == 'Width' & val > 8.0))

print(paste(x_unclean - nrow(x), "outliers of", nrow(x), "trait value data points were removed"))

#Genus/Species fixes
#any species with a special character, including 'sp.', becomes 'unclassified'
x <- x %>% 
  mutate(
    Genus = trimws(Genus),
    Genus = gsub('\\[|\\]', '', Genus),
    Genus = gsub('Sul<a5><e5><81>tobacter', 'Sulfitobacter', Genus),
    Genus = ifelse(grepl('Type-', Genus), NA, Genus),
    Genus = gsub('-.*', '', Genus),
    Species = ifelse(grepl('_| ', Genus), gsub('.*_', '', Genus), Species),
    Genus = gsub('_.*| .*', '', Genus),
    Species = trimws(Species),
    Species = ifelse(grepl('(?!")[[:punct:]]|sp.?', Species, perl = TRUE), 'unclassified', Species))

#remove unclassified and NAs
x <- x[x$Genus != "",] %>%
  filter(!(is.na(Genus) | Genus == 'unclassified')) %>%
  filter(!(is.na(Species) | Species == 'unclassified'))
x <- x %>% filter(Species!="")
x <- x[x$Species!="AB1",]

#for plotting comparisons by source
x_by_source <- x %>%
  group_by(Genus, Species, source, trait) %>%
  summarise(val = ifelse(trait[1] %in% c('Length','Width'), log(mean(val)), mean(val))) %>%
  spread(trait, val)

study_names <- c(Mag2015 = "Magnúsdóttir et al. 2015", 
                 Palm2014 = "Palm et al. 2014", 
                 Browne2016 = "Browne et al. 2016")

trait_sources <- x_by_source %>% 
  ungroup() %>%
  mutate(source = ifelse(source %in% names(study_names), study_names[match(source, names(study_names))], source)) %>%
  gather(trait, val, -Genus, -Species, -source) %>% 
  filter(!is.na(val)) %>% 
  mutate(trait = ifelse(trait %in% c('Spore','Spore_score'), 'Sporulation', trait)) %>%
  group_by(Trait = trait) %>% 
  summarise(Sources = paste(sort(unique(source)), collapse = '; '))

#calculate means among species. Log length/width data. trim whitespace
x <- x %>%
  group_by(Genus, Species, trait) %>%
  summarise(val = ifelse(trait[1] %in% c('Length','Width'), log(mean(val)), mean(val))) %>%
  spread(trait, val) 
  View(gather(x,trait, val, -Genus, -Species))

# Merge spore scores (if no score and ijsem says 0, then 0; 
# otherwise the median of spore score when we know for ijsem spore == 1
x <- x %>% 
  mutate(
    Sporulation = ifelse(is.na(Spore_score), 
                         ifelse(Spore > 0, median(Spore[Spore > 0], na.rm = T), 0), Spore_score)) %>%
  ungroup() %>%
  select(-Spore, -Spore_score)

#remove unclassified and NAs
x <- x %>%
  filter(!(is.na(Genus) | Genus == 'unclassified')) %>%
  filter(!(is.na(Species) | Species == 'unclassified'))

#remove taxa not present in the Living tree Project or the Silva-derived taxonomy file from our usearch pipeline
tmp <- bind_rows(select(tax_succ, Genus, Species), select(tax_LTP, Genus, Species)) %>%
  filter(!(Genus == 'unclassified' | Species == 'unclassified'))
x <- filter(x, paste(Genus, Species) %in% paste(tmp$Genus, tmp$Species))

###save RDS
saveRDS(x, file = 'data\\traits_sparse.RDS')
saveRDS(trait_sources, file = 'data\\trait_sources.RDS')
#++++++++++++++++++++++
#+ Import and read tree from LTP
library(phytools)
library(ape)
library(maps)
ltp_tree <- read.newick('C:/Users/Anwender/Downloads/LTPs132_SSU_tree.newick')

#ltp meta table
ltp_meta <- read.csv('C:/Users/Anwender/Downloads/LTPs132_SSU.csv',sep = "\t",header = F)

# trying to find out which species exist both,
# in my trait table as in the LTP
x <- x %>% 
  mutate(Organism = paste0(Genus," ",Species))

LTP_headers <- read.delim("C:/Users/Anwender/Downloads/LTP_headers", header=FALSE)
# attribute colnames to LTP
colnames(LTP_headers)<- c("Species","Family")
# reorder LTP
LTP_headers <- LTP_headers %>% arrange(Species)
# match ltp to trait table
z <- merge(x,LTP_headers, by.x = 'Organism',by.y = 'V1')
#########################################################
# ++ First Trait version
# arrange format in val column
# int col
c <- c('Oxygen_tolerance','Spore','Gram_positive','Gene_number','Motility','Aggregation_score','Copies_16S','Spore_score','B_vitamins')
traitdata_tresor$val[traitdata_tresor$trait %in% c] <- as.integer(traitdata_tresor$val[traitdata_tresor$trait %in% c])
# double col
d <- c('GC_content','Genome_Mb','IgA','Length','pH_optimum','Salt_optimum','Temp_optimum','Width')
traitdata_tresor$val[traitdata_tresor$trait %in% d] <- as.double(traitdata_tresor$val[traitdata_tresor$trait %in% d])

# checking if d vector are all doubles i.e contains atleast one decimal number
traitdata_tresor$val[traitdata_tresor$trait %in% d]<-ifelse(!grepl('\\.',traitdata_tresor$val[traitdata_tresor$trait %in% d]),
                                                            paste0(traitdata_tresor$val[traitdata_tresor$trait %in% d],".",0),
                                                            traitdata_tresor$val[traitdata_tresor$trait %in% d])

# round decimal places in pH_opt col to 2 decimal places
traitdata_tresor$val[traitdata_tresor$trait=="pH_optimum"] <- round(as.numeric(traitdata_tresor$val[traitdata_tresor$trait=="pH_optimum"]),digits = 2)
# check presence of decimal point, adding . and 2 decimal places
traitdata_tresor$val[traitdata_tresor$trait=="pH_optimum"] <- ifelse(!grepl('\\.',traitdata_tresor$val[traitdata_tresor$trait=="pH_optimum"]),
                                                                     paste0(traitdata_tresor$val[traitdata_tresor$trait=="pH_optimum"],".",00),
                                                                     traitdata_tresor$val[traitdata_tresor$trait=="pH_optimum"])

traitdata_tresor$val[traitdata_tresor$trait=="pH_optimum"] <- ifelse(grepl('\\.\\d{1}$',traitdata_tresor$val[traitdata_tresor$trait=="pH_optimum"]),
                                                                     paste0(traitdata_tresor$val[traitdata_tresor$trait=="pH_optimum"],0),
                                                                     traitdata_tresor$val[traitdata_tresor$trait=="pH_optimum"])
# sort
traitdata_tresor <- traitdata_tresor[order(traitdata_tresor$Genus,traitdata_tresor$Species,traitdata_tresor$trait,traitdata_tresor$source),]

# +++ Latest trait version
# arrange format in val column
# int col
c <- c('Oxygen_tolerance','Spore','Gram_positive','Gene_number','Motility','Aggregation_score','Copies_16S','Spore_score','B_vitamins')
y$val[y$trait %in% c] <- as.integer(y$val[y$trait %in% c])
# double col
d <- c('GC_content','Genome_Mb','IgA','Length','pH_optimum','Salt_optimum','Temp_optimum','Width')
y$val[y$trait %in% d] <- as.double(y$val[y$trait %in% d])

# checking if d vector are all doubles i.e contains atleast one decimal number
y$val[y$trait %in% d]<-ifelse(!grepl('\\.',y$val[y$trait %in% d]),
                                                            paste0(y$val[y$trait %in% d],".",0),
                                                            y$val[y$trait %in% d])

# sort table using genus,species,trait and finally source
y <- y[order(y$Genus,y$Species,y$trait,y$source),]

# arrange traitdata_tresor specie column
traitdata_tresor <- separate(traitdata_tresor, Species, c('Genus1','Species'),fill = "left") %>% 
  select(1,3,4,5,6)

# write file in a tsv file
library(readr)
write_csv(traitdata_tresor,'C:\\Users\\Anwender\\Downloads\\traitdata_tresor',append = F)
write_tsv(y,'C:\\Users\\Anwender\\Downloads\\traitdata_tresor_2.tsv',append = F)

# write in xlsx file
library(writexl) #export data from R format to Excel format
write_xlsx(traitdata_tresor, 'C:\\Users\\Anwender\\Downloads\\traitdata_tresor')
x <- select(x,1,2,3,5,4)

write_tsv(x,'C:\\Users\\Anwender\\Downloads\\traitdata_tresor.tsv',append = F)




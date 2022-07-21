# Custom functions for microbiome analyses.
# Created by John Guittar - modifications by Tresor Guiraud Tonou

scaleFUN <- function(x) sprintf("%.1f", x)

loadpax <- function(pkg){
  # (1) checks package installation, (2) installs them if not, then (3) loads them
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

swan <- function (x) {
  # suppress warnings when converting to numeric
  suppressWarnings(as.numeric(x))
}
## INSTALL BAC_DIVE PACKAGE
#install.packages("BacDive", repos="http://R-Forge.R-project.org")

## IMPORT PACKAGE
library(BacDive) 
library(tidyverse)
## create a bac_dive_open access object
username <- "guiraud.tresor.tonou-fotso@bioinfsys.uni-giessen.de"
password <- "VDw6iDX_pQV.!N6"
bacdive <- open_bacdive(username = username, password = password)
print(bacdive)
# it would be frustrating if the object was already expired on creation
stopifnot(
  inherits(bacdive, "bacdive_access"),
  !summary(bacdive)[c("expired", "refresh_expired")]
)

## run search + fetch in one step
# (a) get BacDive IDs using request func

# build a request function
request_bac <- function(taxon){
  taxon_id <- request(object = bacdive, # request returns just the bacdive IDs 
                      query = taxon,     
                      search = "taxon", page= 0L, handler = NULL, sleep = 0.1)
  output <- ifelse(taxon_id$count == 1,taxon_id$results, 
                   paste0(taxon," has ", taxon_id$count, " Strains!"))
  return(output)
}

# Extract Traits of Taxon/Bacteria species from BacDive,
# and returns a data frame.

bacdat <- function(ID){
  x <- fetch(bacdive,ID)
  y <- x$results
  z <- unlist(y)
  df <- data.frame(
    Genus = ifelse(paste0(ID,".Name and taxonomic classification.genus")%in% names(z),z[[paste0(ID,".Name and taxonomic classification.genus")]],"unclassified"),
    
    Species = ifelse(paste0(ID,".Name and taxonomic classification.species")%in% names(z),z[[paste0(ID,".Name and taxonomic classification.species")]],"unclassified"),
    
    Synonyms = ifelse(paste0(ID,".Name and taxonomic classification.LPSN.synonyms.synonym")%in% names(z),z[[paste0(ID,".Name and taxonomic classification.LPSN.synonyms.synonym")]],NA),
    
    gram = ifelse(paste0(ID,".Morphology.cell morphology.gram stain")%in% names(z),z[[paste0(ID,".Morphology.cell morphology.gram stain")]],NA),
    
    length = ifelse(paste0(ID,".Morphology.cell morphology.cell length")%in% names(z),z[[paste0(ID,".Morphology.cell morphology.cell length")]],NA),
    
    width = ifelse(paste0(ID,".Morphology.cell morphology.cell width")%in% names(z),z[[paste0(ID,".Morphology.cell morphology.cell width")]],NA),
    
    motility = ifelse(paste0(ID,".Morphology.cell morphology.motility")%in% names(z),z[[paste0(ID,".Morphology.cell morphology.motility")]],NA),
    
    oxygen_tolerance = ifelse(paste0(ID,".Physiology and metabolism.oxygen tolerance.oxygen tolerance")%in% names(z),z[[paste0(ID,".Physiology and metabolism.oxygen tolerance.oxygen tolerance")]],NA),
    
    spore = ifelse(paste0(ID,".Physiology and metabolism.spore formation.spore formation")%in% names(z),z[[paste0(ID,".Physiology and metabolism.spore formation.spore formation")]],NA),
    
    aggregation_score = ifelse(paste0(ID,".Physiology and metabolism.observation.observation")%in% names(z),z[[paste0(ID,".Physiology and metabolism.observation.observation")]],NA),
    
    GC_content = ifelse(paste0(ID,".Sequence information.GC content.GC-content")%in% names(z),z[[paste0(ID,".Sequence information.GC content.GC-content")]],NA),
    
    row.names = NULL,stringsAsFactors = FALSE,check.names = T)
  
  return(df)
}

mm_to_micrometre_width <- function(string){
  if(grepl("\\d+\\.?\\d?.*\\-.*\\d+\\.?\\d?.?mm",string)){ #grabs interval ranges
    string <- gsub("mm","",string)
    string <- str_split(string,"-")
    lower <- string[[1]][1]
    upper <- string[[1]][2]
    range <- paste0(as.numeric(lower)*1000,"-",as.numeric(upper)*1000)
    return(range)
  }else if(grepl("\\d+\\.?\\d?.?mm",string)){ #grabs non-interval ranges
    string <- gsub("mm","",string)
    return(as.numeric(string)*1000)
  }else if(grepl("\\d+\\.?\\d?.?\\-.?\\d+\\.?\\d?.?[×|X|x]\\d+\\.?\\d?.?\\-.?\\d+\\.?\\d?.?",string)){# grabs "1.2-0.4×0.4-0.1"
    width <- str_split(string,"×|x|X")
    return(width[[1]][1])
  }else if(grepl("\\d+\\.?\\d?.?[×|±].?\\d+\\.?\\d?.?",string)){ #split &returns just the width
    width <- str_split(string,"×|±")
    return(width[[1]][1])
  }
  else{ 
      return(string)
    }
}

mm_to_micrometre_length <- function(string){
  if(grepl("\\d+\\.?\\d?.*\\-.*\\d+\\.?\\d?.?mm",string)){ #grabs interval ranges
    string <- gsub("mm","",string)
    string <- str_split(string,"-")
    lower <- string[[1]][1]
    upper <- string[[1]][2]
    range <- paste0(as.numeric(lower)*1000,"-",as.numeric(upper)*1000)
    return(range)
  }else if(grepl("\\d+\\.?\\d?.?mm",string)){ #grabs non-interval ranges
    string <- gsub("mm","",string)
    return(as.numeric(string)*1000)
  }else if(grepl("\\d+\\.?\\d?.?\\-.?\\d+\\.?\\d?.?[×|X|x]\\d+\\.?\\d?.?\\-.?\\d+\\.?\\d?.?",string)){# grabs "1.2-0.4×0.4-0.1"
    width <- str_split(string,"×|x|X")
    return(width[[1]][2])
  }else if(grepl("\\d+\\.?\\d?.?[×|±].?\\d+\\.?\\d?.?",string)){ #split &returns just the width
    width <- str_split(string,"×|±")
    return(width[[1]][2])
  }
  else{ 
    return(string)
  }
}


# Code to concatenate and remove duplicates
# from a systematic search of multiple databases

# search details----------------------------------------------------------------

# Search performed 2/8/2018

# search strings used: ("emerg* infect*" OR spillover OR zoonotic OR zoonos*) AND disease

# Web of Science:
# Refined by: DOCUMENT TYPES: ( REVIEW )
# Timespan: 1992-2018. Indexes: SCI-EXPANDED, SSCI, A&HCI.

# Exporting instructions: Add to Marked List; select: Author(s)/Editor(s), Abstract, Accesssion Number, Title, Source; Save to other file formats; Tab-delimited (Win)

# PubMed:
# Refined by: Review, English, 1992-2018.

# Anthropology Plus, EconLit, Sociological Collection
# Boolean/Phrase search
# Language: English
# Date: 1992-2018
# Share -> Email a link to download exported files -> XML format
# https://conversiontools.io/convert_xml_to_excel/

# Annual Reviews papers transferred manually (no export available)

# prep workspace---------------------------------------------------------------

# packages
library(metagear)
library(stringr)

# load in data files
ann <- read.csv("./Data/ARresult.csv", header = T, as.is = T, encoding = "UTF-8")
pub <- read.csv("./Data/pubmedresult.csv", header = T, as.is = T, 
                encoding = "UTF-8")
soc <- read.csv("./Data/socresult.csv", header = T, as.is = T, encoding = "UTF-8")
wos <- read.csv("./Data/WOSresult.csv", header = T, as.is = T, encoding = "UTF-8")


cleannames <- c("authors", "title", "year", "journal", "search", "abstract")

# Annual Reviews----------------------------------------------------------------
# paper info extracted manually
# 278 studies

# PubMed----------------------------------------------------------------
# choose relevant columns to keep
pub$AB <- NA
pub$year <- str_sub(pub$ShortDetails, start= -4)
pub$journal <- str_sub(pub$ShortDetails, 1, str_length(pub$ShortDetails)-6)
keeppub <- c("Description", "Title", "year", "journal", "Resource", "AB")
pub <- pub[keeppub]
# rename to sensible titles
names(pub) <- cleannames
# 3063 records

# social science abstracts------------------------------------------------------

# choose relevant columns to keep
keepsoc <- c("au", "atl", "year", "jtl", "longDbName", "ab")
soc <- soc[keepsoc]
# rename to sensible titles
names(soc) <- cleannames
# remove all blank or NA rows
soc <- soc[!apply(is.na(soc) | soc == "", 1, all), ]
# for now, remove all authors beyond 1st
soc <- na.omit(soc)
# 137 records

# Web of Science----------------------------------------------------------------
# choose relevant columns to keep
wos$search <- "Web of Science"
keepwos <- c("AU", "TI", "PY", "SO", "search", "AB")
wos <- wos[keepwos]
# rename to sensible titles
names(wos) <- cleannames
# 2264 records

# concatenating all searches ---------------------------------------------------
# rbind searches together
allResults <- rbind.data.frame(ann, pub, soc, wos)
# sort by title
allResults <- allResults[order(allResults$title, decreasing = F), ]
# convert to character
allResults$title <- as.character(allResults$title)
# remove periods in the titles
allResults$title <- gsub("[.]", "", allResults$title)
# remove commas in the titles
allResults$title <- gsub("[,]", "", allResults$title)
# remove all non-graphical characters
allResults$title <- str_replace_all(allResults$title, "[^[:graph:]]", " ") 
# trim whitespace
allResults$title <- str_trim(allResults$title)
# convert titles to lower case
allResults$title <- tolower(allResults$title)
# remove duplicate titles
trimmed <- allResults[!duplicated(allResults$title), ]

# number of duplicates
nrow(allResults) - nrow(trimmed)

# export to folder
#write.csv(trimmed, "./Data/trimmed_noduplicates.csv", row.names = FALSE)

# use metagear to divide papers randomly for 1st round of screening
#refs <- effort_initialize(trimmed)
#team <- c("Cecilia", "Joy")
#refs_unscreened <- effort_distribute(refs, reviewers = team, save_split = TRUE)

# note that during manual screening, more duplicate papers were found and removed
# therefore, Fig. 1 shows 4517 papers as having their title/abstract screened
# rather than the 4524 of the "trimmed_noduplicates"

# the final results of screening can be found in the "Data/effort_all.csv"
## Combining data for John Casey

require(readr)
library(lattice)
library(dplyr)
library(tidyr)

setwd("c:/Users/Angela/Documents/UW lab/KM1513/Diel/")

## Load data ----------------
HILIC.metabdata <- read.csv('HILIC.allIS.BlankFilled.ISremoved.csv', header=TRUE, row.names = 1)
HILIC.metabdata <- arrange(HILIC.metabdata,SampID,replicate)
cyano.aq.metabdata <- read.csv('Cyano.aq.allIS.BlankFilled.csv', header=TRUE, row.names = 1)
cyano.aq.metabdata <- arrange(cyano.aq.metabdata,SampID,replicate)
cyano.DCM.metabdata <- read.csv('Cyano.DCM.allIS.BlankFilled.csv', header=TRUE, row.names = 1)
cyano.DCM.metabdata <- arrange(cyano.DCM.metabdata,SampID,replicate)

c(names(cyano.aq.metabdata),names(cyano.DCM.metabdata),
  names(HILIC.metabdata))[duplicated(c(names(cyano.aq.metabdata),
                                       names(cyano.DCM.metabdata),
                                       names(HILIC.metabdata)))]

HILIC.metabdata <- HILIC.metabdata %>%
     select(-Adenosyl.Homocysteine, -Adenosyl.Methionine,
            -Chitobiose, -runDate, -type)
cyano.aq.metabdata <- cyano.aq.metabdata %>%
     select(-runDate, -type, -Betaine, -Choline,
            -Vitamin.B3)
cyano.DCM.metabdata <- cyano.DCM.metabdata %>%
     select(-runDate, -type)

allMetabData <- full_join(HILIC.metabdata, cyano.aq.metabdata)
allMetabData <- full_join(allMetabData, cyano.DCM.metabdata)

sampLog <- read.csv("Clean Metabolomics Sample Log for SCOPE Cruise_DIEL.csv",
                    header = TRUE)
sampLog <- sampLog %>% filter(Sample.ID<21) %>%
     select(Sample.ID, Station.Number, Replicate, Time) %>%
     rename(SampID = Sample.ID, replicate = Replicate)
combData <- full_join(sampLog, allMetabData)

## Output with comment -------------------------
comment.text <- "# Targeted Metabolite Data from KM1513 Diel 1. Values are normalized to the best-matched internal standard and scaled to represent their relative peak area."
con <- file("Diel1TargetedMetabData.csv", open="wt")
writeLines(paste(comment.text), con)
write.csv(combData, con)
close(con)

require(readr)
library(lattice)
library(dplyr)
library(tidyr)

setwd("c:/Users/Angela/Documents/UW lab/KM1513/Diel/")

## Load data ----------------
HILIC.metabdata <- read.csv('Normalized_data_QC_output_wBlkSubDiel_both_HILIC_allData_011617new.csv', header=TRUE, row.names = 1)
HILIC.metabdata <- arrange(HILIC.metabdata,SampID,replicate)
HILIC.metabdata <- HILIC.metabdata %>% 
     mutate(Compound.Name = as.character(Compound.Name))
cyano.aq.metabdata <- read.csv('Normalized_data_QC_output_wBlkSubDiel_all_Cyano_Aq_011017.csv', header=TRUE, row.names = 1)
cyano.aq.metabdata <- arrange(cyano.aq.metabdata,SampID,replicate)
cyano.aq.metabdata <- cyano.aq.metabdata %>% 
     mutate(Compound.Name = as.character(Compound.Name))
cyano.DCM.metabdata <- read.csv('Normalized_data_QC_output_wBlkSubDiel_both_Cyano_DCM_all.csv', header=TRUE, row.names = 1)
cyano.DCM.metabdata <- arrange(cyano.DCM.metabdata,SampID,replicate)
cyano.DCM.metabdata <- cyano.DCM.metabdata %>% 
     mutate(Compound.Name = as.character(Compound.Name))
# c(names(cyano.aq.metabdata),names(cyano.DCM.metabdata),
#   names(HILIC.metabdata))[duplicated(c(names(cyano.aq.metabdata),
#                                        names(cyano.DCM.metabdata),
#                                        names(HILIC.metabdata)))]

HILIC.sub <- HILIC.metabdata %>%
     filter(Compound.Name!="Adenosyl Homocysteine") %>%
     filter(Compound.Name!="Adenosyl Methionine") %>%
     filter(Compound.Name!="Chitobiose") %>%
     filter(type != "Blk")

cyano.aq.sub <- cyano.aq.metabdata %>%
     filter(Compound.Name!="Vitamin B1") %>%
     filter(Compound.Name!="4-hydroxybenzaldehyde") %>%
     filter(Compound.Name!="Caffeine") %>%
     filter(Compound.Name!="Betaine") %>%
     filter(Compound.Name!="Choline") %>%
     filter(Compound.Name!="Vitamin B3") %>%
     filter(type != "Blk")

c(unique(cyano.aq.sub$Compound.Name),unique(cyano.DCM.metabdata$Compound.Name))[
     duplicated(c(unique(cyano.aq.sub$Compound.Name),
                  unique(cyano.DCM.metabdata$Compound.Name)))]

c(unique(HILIC.sub$Compound.Name),unique(cyano.DCM.metabdata$Compound.Name))[
     duplicated(c(unique(HILIC.sub$Compound.Name),
                  unique(cyano.DCM.metabdata$Compound.Name)))]

cyano.DCM.sub <- cyano.DCM.metabdata %>%
     filter(Compound.Name!="Chitobiose") %>%
     filter(Compound.Name!="Glutathione") %>%
     filter(Compound.Name!="Caffeine") %>%
     filter(Compound.Name!="Pyridoxal") %>%
     filter(Compound.Name!="Adenosyl Methionine") %>%
     filter(Compound.Name!="Phenylalanine") %>%
     filter(Compound.Name!="Pyridoxal Phosphate") %>%
     filter(Compound.Name!="Tryptophan") %>%
     filter(Compound.Name!="Choline") %>%
     filter(Compound.Name!="Vitamin B3") %>%
     filter(type != "Blk")

allMetabData <- full_join(cyano.DCM.sub, cyano.aq.sub)
allMetabData <- full_join(allMetabData, HILIC.sub)

sampLog <- read.csv("Clean Metabolomics Sample Log for SCOPE Cruise_DIEL.csv",
                    header = TRUE)
sampLog <- sampLog %>% 
     select(Sample.ID, Station.Number, Replicate, Time) %>%
     rename(SampID = Sample.ID, replicate = Replicate) %>%
     mutate(replicate = as.character(replicate),
            SampID = as.character(SampID))
sampLog$Time[sampLog$Time==700] <- 600
combData <- full_join(sampLog, allMetabData)

## Output with comment -------------------------
comment.text <- "# Targeted Metabolite Data from KM1513, both Diel periods, Values are normalized to the best-matched internal standard, and scaled to represent their relative peak area."
con <- file("AllDielTargetedMetabData.csv", open="wt")
writeLines(paste(comment.text), con)
write.csv(combData, con)
close(con)

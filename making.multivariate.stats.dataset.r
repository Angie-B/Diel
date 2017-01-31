##Making multivariate Stats data set

## observations are in the rows, so each sample is an observation
## variables are in the columns, so each compound is a variable
library(dplyr)
library(readr)
library(tidyr)

dat.join.aq.allIS <- read_csv("Cyano_Aq_diel_1_Samples_Normalized_BlkSub_backfill_allIS_noblk_new0927.csv")
dat.join.aq.allIS <- dat.join.aq.allIS[,-1] %>%
     mutate(column = "aq") %>%
     mutate(SampID = ifelse(SampID=="AllCyanoAqExtracts3",300,
                            ifelse(SampID=="AllCyanoAqExtracts2",200,
                                   ifelse(SampID=="AllCyanoAqExtracts1",100,
                                          as.numeric(SampID)))))

organized.data <- dat.join.aq.allIS %>%
     select(runDate, type, SampID, replicate, 
            Compound.Name, PooPlusModel) %>%
     filter(type=="Smp") %>%
     spread(Compound.Name, PooPlusModel)
     

other.data.data <- dat.join.aq.allIS %>%
     select(runDate, type, SampID,replicate, Compound.Name,
            PooPlusModel.IS, Notes) %>%
     filter(type=="Smp")

other.sample.data <- dat.join.aq.allIS %>%
     filter(type=="Smp") %>%
     select(runDate, type, SampID,replicate,
            Vol, Date) %>%
     distinct() %>%
     mutate(HoursSinceStart = ifelse(SampID == 1, 1,
                                     (SampID-1)*4) ,
            ToD = ifelse(SampID == 1, 7, 
                         (SampID*4)+2),
            ToD = ToD %% 24)

Notes <- dat.join.aq.allIS %>%
     select(runDate, type, SampID, replicate, 
            Compound.Name, Notes) %>%
     filter(type=="Smp") %>%
     filter(grepl("overloaded",Notes)) %>%
     spread(Compound.Name, Notes) 
names(Notes)[-1:-4]

tf.data <- dat.join.aq.allIS %>%
     select(runDate, type, SampID, replicate, 
            Compound.Name, Area) %>%
     filter(type=="Smp") %>%
     filter(!grepl(paste(names(Notes)[-1:-4], collapse="|"),Compound.Name)) %>%
     mutate(tf.area = !is.na(Area)) %>% select(-Area) %>%
     spread(Compound.Name, tf.area)
     
write.csv(organized.data, "Cyano.aq.allIS.BlankFilled.csv")

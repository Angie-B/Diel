## Let's make some NMDS plots of Diel Data!

library(tidyr)
library(stringr)
library(readr)
library(dplyr)
library(vegan)
source("C:/Users/Angela/Documents/UW School work/Fall 2016/multivariate stats/files/biostats.R")
source("C:/Users/Angela/Documents/UW School work/Fall 2016/multivariate stats/files/coldiss.R")

### Get Good Data --------

## Read in Data -----
dat <- read.csv("021617_AllDielTargetedMetabData.csv", comment.char = "#")
dat <- dat[,-1]

##get only useful data --------
all.stations.long <- dat %>%
     filter(type=="Smp") %>%
     filter(Compound.Name!="d-IAA")%>%
     filter(Compound.Name!="D-Phenylalanine")%>%
     filter(Compound.Name!="D-Tryptophan")%>%
     filter(Compound.Name!="Vitamin B2_IS")%>%
     filter(Compound.Name!="D8-Arachidonic Acid")%>%
     filter(Compound.Name!="d4-Tryptamine")%>%
     filter(Compound.Name!="Vitamin B1_IS")%>%
     filter(Compound.Name!="13C Acetyl CoA")%>%
     filter(Compound.Name!="13C-Isethionic Acid")%>%
     filter(Compound.Name!="d4 Taurine")%>%
     filter(Compound.Name!="13C-Sulfoacetic Acid")%>%
     filter(Compound.Name!="13C-Sulfolactic Acid")%>%
     filter(Compound.Name!="d3-Cysteic Acid")%>%
     filter(Compound.Name!="d4 Succinic Acid")%>%
     filter(Compound.Name!="D-Methionine")%>%
     filter(Compound.Name!="Heavy Alanine")%>%
     filter(Compound.Name!="Heavy Histidine")%>%
     filter(Compound.Name!="Heavy Isoleucine") %>%
     select(SampID, replicate, Time, Date, Compound.Name, PooPlusNormed) %>%
     mutate(Day = as.character(Date) %>%
                 str_replace("26-Jul","0")%>%
                 str_replace("27-Jul","1")%>%
                 str_replace("28-Jul","2")%>%
                 str_replace("29-Jul","3")%>%
                 str_replace("30-Jul","4")%>%
                 str_replace("31-Jul","5")%>%
                 str_replace("1-Aug","6")%>%
                 str_replace("2-Aug","7")%>%
                 str_replace("3-Aug","8"),
            Day = as.numeric(Day),
            HourOfExperiment = ((Day*2400)+Time)/100 )

overloaded.cmpds <- dat %>%
     mutate(Flag = ifelse(test = grepl("overloaded",
                                       Notes),
                          yes = "over",
                          no = "ok")) %>%
     filter(Flag == "over") %>%
     mutate(Compound.Name = as.character(Compound.Name))
overloaded.cmpds <- unique(overloaded.cmpds$Compound.Name)
overloaded.cmpds

all.stations.long <- all.stations.long %>%
     filter(Compound.Name!="DMSP")%>%
     filter(Compound.Name!="Homarine")%>%
     filter(Compound.Name!="X5P")

wide.all.stations <- all.stations.long %>%
     select(SampID, replicate, HourOfExperiment, Compound.Name, PooPlusNormed) %>%
     spread(Compound.Name, PooPlusNormed) %>%
     mutate(SampID = as.numeric(as.character(SampID)))
## fill in NAs with mean of the other replicates
NAs <- all.stations.long %>%
     group_by(HourOfExperiment, Compound.Name) %>%
     summarize(mean = mean(PooPlusNormed, na.rm=T))

for(i in 4:ncol(wide.all.stations)){
     if(any(is.na(wide.all.stations[,i]))){
          IDs <- wide.all.stations[is.na(wide.all.stations[,i]),"HourOfExperiment"]
          name <- colnames(wide.all.stations)[i]
          for(j in 1:length(IDs)){
               val = as.numeric(NAs[(NAs$HourOfExperiment==IDs[j] & 
                                          NAs$Compound.Name==name),
                                    "mean"])
               wide.all.stations[(is.na(wide.all.stations[,i]) & 
                                       wide.all.stations$HourOfExperiment==IDs[j]),i] <-
                    val
          }
     }
}

## split diel1 and 2 --------
wide.all.stations <- wide.all.stations %>% arrange(SampID)
diel1 <- wide.all.stations %>% 
     filter(SampID<25)
diel2 <- wide.all.stations %>%
     filter(SampID>25)

## NMDS ------
nmds.data <- decostand(diel2[,-1:-3], method = 'standardize', 2)
rownames(nmds.data) <- paste(diel2$SampID, diel2$replicate,sep="_")

spe.nmds<-metaMDS(nmds.data, distance='euclidean', k=2, try = 100,
                  autotransform=FALSE, noshare = FALSE, wascores = FALSE,
                  trymax=1000)
spe.nmds
## stress 0.1739219
## extra culled 0.16844 
spe.nmds
stressplot(spe.nmds)

## scree plot vs nmber of dimentions
nmds.scree(nmds.data, distance='euclidean', k=15, 
           autotransform=FALSE, trymax=20) 

## monte carlo randomization
nmds.monte(nmds.data, distance = 'euclidean', k=2,
           autotransform= FALSE,
           trymax = 20)
## p val of 0.0099
colorkey1 <- rep(rainbow(n=6), each = 3,3)
colorkey <- c(colorkey1, colorkey1[1:6])
colorkey2 <- c(colorkey1, colorkey1[1:3])

colorkey.legend <- rainbow(n=6)
# ringkey <-c(rep("grey",10),rep("black",10))
plot(spe.nmds,type='n')
text(spe.nmds,labels=row.names(nmds.data))
points(spe.nmds,  bg = colorkey,  pch=21, cex=1.2)
# legend(x=20, y=10, legend = c(c("6","10","14","18","22","2")),
       # fill=c(colorkey.legend))

legend(x=20, y=10, legend = c(c("18","22","2","6","10","14")),
fill=c(colorkey.legend))


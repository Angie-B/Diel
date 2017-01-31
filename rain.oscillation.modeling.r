## RAIN example
require(readr)
library(lattice)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stats)
library(cowplot)
library(graphics)
require(rain)
setwd("c:/Users/Angela/Documents/UW lab/KM1513/Diel/")

## Load data ----------------
HILIC.metabdata <- read.csv('HILIC.allIS.BlankFilled.ISremoved.csv', header=TRUE, row.names = 1)
HILIC.metabdata <- arrange(HILIC.metabdata,SampID,replicate)
IDs <- paste(HILIC.metabdata$SampID,HILIC.metabdata$replicate,sep="_")
row.names(HILIC.metabdata) <- IDs
HILIC.metabdata <- HILIC.metabdata[,5:ncol(HILIC.metabdata)]
HILIC.tdata <- t(HILIC.metabdata)

cyano.aq.metabdata <- read.csv('Cyano.aq.allIS.BlankFilled.csv', header=TRUE, row.names = 1)
cyano.aq.metabdata <- arrange(cyano.aq.metabdata,SampID,replicate)
IDs <- paste(cyano.aq.metabdata$SampID,cyano.aq.metabdata$replicate,sep="_")
row.names(cyano.aq.metabdata) <- IDs
cyano.aq.metabdata <- cyano.aq.metabdata[,5:ncol(cyano.aq.metabdata)]
cyano.aq.tdata <- t(cyano.aq.metabdata)

cyano.DCM.metabdata <- read.csv('Cyano.DCM.allIS.BlankFilled.csv', header=TRUE, row.names = 1)
cyano.DCM.metabdata <- arrange(cyano.DCM.metabdata,SampID,replicate)
IDs <- paste(cyano.DCM.metabdata$SampID,cyano.DCM.metabdata$replicate,sep="_")
row.names(cyano.DCM.metabdata) <- IDs
cyano.DCM.metabdata <- cyano.DCM.metabdata[,5:ncol(cyano.DCM.metabdata)]
cyano.DCM.tdata <- t(cyano.DCM.metabdata)

## Number of Compounds measured -------
cmpd.names <- unique(c(names(cyano.aq.metabdata), 
                       names(cyano.DCM.metabdata), names(HILIC.metabdata)))
cmpd.names <- cmpd.names[-1:-4]
cmpd.names[order(cmpd.names)]
cmpd.names
length(cmpd.names)

## rain --------
timestamp()
HILIC.results <- rain(HILIC.metabdata, deltat = 4, period = 24, nr.series = 3,
                peak.border = c(0.3, 0.7), method = "independent",
                verbose = TRUE)
timestamp()
best <- order(HILIC.results$pVal)[1:10]

HILIC.significant <- HILIC.results[HILIC.results$pVal<0.05,]

xyplot(as.matrix(HILIC.tdata[best, 0:19 * 3 +  rep(c(1:3), each = 20)])
        ~rep(0:19 * 4 + 6, each = 10) |rownames(HILIC.tdata)[best],
        scales = list(y = list(relation = 'free')),
        layout = c(2, 5), type = 'b', pch = 16, xlab = 'time',
        ylab = 'expression value', cex.lab = 1, title = "HILIC cmpds")

write.csv(HILIC.significant, "significantOscillations.HILIC.diel1.csv")
write.csv(HILIC.results, "Oscillation.Modeling.HILIC.diel1.csv")

cyano.aq.results$column <- "Cyano.Aq"
cyano.aq.results$Compound <- row.names(cyano.aq.results)
cyano.DCM.results$column <- "Cyano.DCM"
cyano.DCM.results$Compound <- row.names(cyano.DCM.results)
HILIC.results$Compound <- row.names(HILIC.results)

HILIC.results$column <- "HILIC"

all.results <- full_join(cyano.aq.results, cyano.DCM.results)
all.results <- full_join(HILIC.results, all.results)

# write.csv(all.results, "oscillation.modeling.all.diel1.csv")

## plots ----------
all.results <- read.csv("oscillation.modeling.all.diel1.csv")
cyano.master <- read.csv("../CYANO_MasterList_Summer2016.csv")
HILIC.master <- read.csv("../HILIC_MasterList_Summer2016.csv")
require(stringr)
cyano.master <- cyano.master %>% select(Compound.Name, Group) %>% 
     group_by(Compound.Name) %>% summarize(Group = Group[1]) %>%
     filter(Group!="Internal Std") %>%
     mutate(Compound.Name = Compound.Name %>%
                 str_replace(" ",
                             ".") %>%
                 str_replace("-",
                             ".") %>%
                 str_replace(",",
                             ".")%>%
                 str_replace("Tocopherol.(Vit E)",
                             "Tocopherol.Vit.E")%>%
                 str_replace("Indole.3.methyl acetate",
                             "Indole.3.methyl.acetate")%>%
                 str_replace("7.dehydrocholesterol",
                             "X7.dehydrocholesterol")%>%
                 str_replace("4.hydroxybenzaldehyde",
                             "X4.hydroxybenzaldehyde")%>%
                 str_replace("Methyl.indole.3 carboxylate",
                             "Methyl.indole.3.carboxylate")%>%
                 str_replace("Indole.3.carboxylic acid",
                             "Indole.3.carboxylic.acid")%>%
                 str_replace("Indole.3 acetamide","Indole.3.acetamide"))
HILIC.master <- HILIC.master %>% select(Compound.Name, Group) %>% 
     group_by(Compound.Name) %>% summarize(Group = Group[1]) %>%
     filter(Group!="Internal Std") %>%
     mutate(Compound.Name = Compound.Name %>%
                 str_replace(" ",
                             ".") %>%
                 str_replace("-",
                             ".") %>%
                 str_replace(",",
                             ".") %>%
                 str_replace("Fructose.6 phosphate",
                             "Fructose.6.phosphate")%>%
                 str_replace("Glucose.6 phosphate",
                             "Glucose.6.phosphate")%>%
                 str_replace("glycerol.3 phosphate",
                             "glycerol.3.phosphate")%>%
                 str_replace("Ribose.5 phosphate",
                             "Ribose.5.phosphate") %>%
                 str_replace("trans.Hydroxyl proline",
                             "trans.Hydroxyl.proline")) %>%
     mutate(column = "HILIC")
test <- left_join(all.results, HILIC.master, by = c("Compound"="Compound.Name","column"))
test2 <- left_join(all.results, cyano.master, by=c("Compound"="Compound.Name"))
complete.results <- rbind(test[which(test$column=="HILIC"),],
                          test2[which(test2$column!="HILIC"),])
complete.results <- complete.results %>% mutate(Group = as.character(Group)) %>%
     mutate(Group = ifelse(Compound=="Tocopherol.Vit.E","Vitamin",Group)) %>%
     mutate(Group = ifelse(grepl("Vitamin",Group),"Vitamin",
                         ifelse(grepl("Amino Acid -",Group),"Amino Acid related",
                         ifelse(grepl("Amino Acid synthesis",Group),"Amino Acid related",
                         ifelse(grepl("Amino Acid derivative",Group),"Amino Acid related",
                         ifelse(grepl("Amino sugar",Group),"Amino Acid related",  
                         ifelse(grepl("Antioxidant",Group),"Antioxidant",
                                   Group)))))))
complete.results$sig <- complete.results$pVal<0.05
complete.results$peak.time <- complete.results$phase+2
complete.results <- complete.results %>% 
     mutate(peak.time = ifelse(peak.time == 26, 2, peak.time))
complete.results$cmpd.class <- complete.results$Group

complete.results <- read.csv("Oscillation.Peaktime.Plot.Data.csv")
alpha <- ifelse(complete.results$sig, 0.8, 0.25)
levels(complete.results$cmpd.class)
complete.results$cmpd.class = factor(complete.results$cmpd.class,
          levels(complete.results$cmpd.class)[c(5,4,17,12,11,18,8,19,3,9,10,
                                                20,16,6,15,14,2,1,7,13)])

ggplot(complete.results, aes(x = peak.time, y = cmpd.class,
                        text = Compound)) + 
     geom_jitter(width = 0.5, height = 0.05, size = 2, alpha = alpha)  


## test -----
HILIC.test <- HILIC.metabdata %>% 
     mutate(ToD = c(rep(c(6,10,14,18,22,2),each=3,times=3),6,6,6,10,10,10))
ggplot(HILIC.test, aes(x=ToD, y=DMSP)) + geom_point()

## FDR ------
## Currently broken -- 
## "Error Error in plot.new() : figure margins too large
## In addition: Warning message:
##      In fdrtool(complete.results$pVal, statistic = "pvalue") :
##      There may be too few input test statistics for reliable FDR calculations!"

# require(fdrtool)
# fdrtool(complete.results$pVal, statistic = "pvalue")

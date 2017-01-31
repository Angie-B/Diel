## Cluster Metabolites
require(readr)
library(lattice)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stats)
library(cowplot)
library(graphics)
require(rain)
library(vegan)
setwd("c:/Users/Angela/Documents/UW lab/KM1513/Diel/")

## Load data ----------------
HILIC.metabdata <- read.csv('C:/Users/Angela/Documents/UW School work/Fall 2016/multivariate stats/project.files/HILIC.allIS.BlankFilled.badcmpdsremoved.csv', header=TRUE, row.names = 1)
HILIC.metabdata <- arrange(HILIC.metabdata,SampID,replicate)
IDs <- paste(HILIC.metabdata$SampID,HILIC.metabdata$replicate,sep="_")
row.names(HILIC.metabdata) <- IDs
HILIC.metabdata <- HILIC.metabdata[,5:(ncol(HILIC.metabdata))]
HILIC.metabdata.trans <- decostand(HILIC.metabdata, 
                                   method = "standardize", MARGIN = 2,
                                   na.rm=T)
Metab.data <- t(HILIC.metabdata.trans)
## Euclidean distance between samples
## based on metabolite characteristics
samp.eucd <- vegdist(HILIC.metabdata.trans, method = 'euclidean')
metab.eucd <- vegdist(Metab.data, method = "euclidean")

## dendrograms of samples -----
#space conserving
sampcl.ave <- hclust(samp.eucd, method='average')
plot(sampcl.ave)
#space conserving with Ward
sampcl.ward <- hclust(samp.eucd, method = 'ward.D')
plot(sampcl.ward)

## dendrograms of metabs ----

#space conserving
metabcl.ave <- hclust(metab.eucd, method='average')
plot(metabcl.ave)
#space conserving with Ward
metabcl.ward <- hclust(metab.eucd, method = 'ward.D')
plot(metabcl.ward)

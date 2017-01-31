## plot FCM data


library(ggplot2)
library(tidyr)
library(stringr)
library(readr)
library(dplyr)
library(vegan)
library(cowplot)
# source("../files/biostats.R")
# source("../files/coldiss.R")
setwd("~/UW School work/Fall 2016/multivariate stats/project.files/")

samp.data <- read.csv("All.Samp.Data.Clean.csv", row.names = 1)
FCM.data <- samp.data %>%
     filter(SampID<21) %>%
     group_by(SampID) %>%
     summarise(FCM.PRO = mean(FCM.PRO),
               FCM.SYN2 = mean(FCM.SYN2),
               FCM.HET = mean(FCM.HET),
               FCM.PEUK = mean(FCM.PEUK),
               TimeOfDay = mean(Timepoint..HST)) %>%
     gather(Organism, FCM.Count, -SampID, -TimeOfDay)

ggplot(FCM.data, aes(x=SampID, y=FCM.Count, color = Organism)) + geom_line() +
     geom_vline(xintercept = c(5.5,11.5,17.5,23.5),linetype=2) 
ggplot(FCM.data, aes(x=TimeOfDay, y=FCM.Count, color = Organism)) + geom_point()

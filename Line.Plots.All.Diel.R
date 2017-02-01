library(ggplot2)
library(tidyr)
library(tibble)
require(graphics); require(grDevices)
library(Hmisc)
library(gtools)
library(cowplot)
require(RColorBrewer)
library(readr)
library(plotly)
library(stringr)
library(ggbiplot)
library(stats)
library(dplyr)

## Read in Data -----
dat <- read.csv("AllDielTargetedMetabData.csv", comment.char = "#")
dat <- dat[,-1]

##get useful data --------
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
     select(SampID, replicate, Time, Date, Compound.Name, PooPlusModel) %>%
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
     select(SampID, replicate, HourOfExperiment, Compound.Name, PooPlusModel) %>%
     spread(Compound.Name, PooPlusModel)
## fill in NAs with mean of the other replicates
NAs <- all.stations.long %>%
     group_by(HourOfExperiment, Compound.Name) %>%
     summarize(mean = mean(PooPlusModel, na.rm=T))

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

## Standardize ------
## zero mean an dunit variance for each cmpd

std.all.stations <- wide.all.stations
std.all.stations[,-1:-3] <- decostand(std.all.stations[,-1:-3], 
                                      method = 'standardize', 2)
std.all.stations.long <- std.all.stations %>%
     gather(Compound.Name, Value, -SampID, -replicate, -HourOfExperiment )

## Mean and SD -------
summary.stations.long <- std.all.stations.long %>% ungroup() %>%
     group_by(SampID, Compound.Name, HourOfExperiment) %>%
     summarise(mean = mean(Value, na.rm=T),
               sd = sd(Value),
               min = mean-sd,
               max = mean+sd,
               median = median(Value, na.rm=T)) %>%
     ungroup() %>%
     filter(SampID!=25)%>%
     mutate(SampID = as.numeric(as.character(SampID)),
            b = ifelse(SampID < 25,"a","b"))

## plot all --------
p1 <- ggplot(summary.stations.long, aes(x=HourOfExperiment, y=mean, 
                                        group = Compound.Name)) +
     scale_x_continuous(breaks=seq(0,216,24))+
     geom_line() + facet_wrap(~b, scales = "free_x")

## read in significance ------
signif <- read.csv("significantOscillations.allCol.diel1.csv")
signif <- signif %>% select(Compound.Name, Significant, phase)
merged <- full_join(signif,summary.stations.long)
merged$Significant[is.na(merged$Significant)] <- FALSE
     

## plot only sig ------
plot.data <- merged %>% filter(!is.na(b))
p2 <- ggplot(plot.data, aes(x=HourOfExperiment, y=mean, 
                            group = Compound.Name)) +
     scale_x_continuous(breaks=seq(0,216,24)) +
     geom_line(aes(alpha = Significant)) + facet_wrap(~b, scales = "free_x")

plot.data <- plot.data %>% filter(Significant==TRUE)

p3 <- ggplot(plot.data, aes(x=HourOfExperiment, y=mean, 
                         group = Compound.Name)) +
     geom_line(aes(color = as.factor(phase))) + 
     geom_ribbon(aes(ymin = min,
                     ymax= max),alpha=0.1) +
     scale_x_continuous(breaks=seq(0,216,24)) +
     # geom_vline(xintercept = c(24,48,72,96,120,144,168,192),linetype=2) +
     facet_wrap(phase~b, scales = "free_x")

## get cmpd data -----
cmpd.dat <- read.csv("C:/Users/Angela/Documents/UW lab/MRM_Methods_Table_noIS.csv")
cmpd.dat.sub <- cmpd.dat %>%
     select(Compound.Name, C, H, O, N, P, S)
cmpd.dat.sub <- cmpd.dat.sub[!duplicated(cmpd.dat.sub), ]
ratio.merge <- full_join(merged, cmpd.dat.sub)

plot.data <- ratio.merge %>% 
     filter(Significant==TRUE) %>%
     filter(!is.na(b)) %>%
     filter(S>0)

p4 <- ggplot(plot.data, aes(x=HourOfExperiment, y=mean, 
                            group = Compound.Name)) +
     geom_ribbon(aes(ymin = min,
                     ymax= max),alpha=0.1) +
     scale_x_continuous(breaks=seq(0,216,12))+
     geom_line(aes(color = as.factor(phase))) + facet_grid(S~b, scales="free_x")

## multiplot ---
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
     require(grid)
     
     plots <- c(list(...), plotlist)
     
     numPlots = length(plots)
     
     if (is.null(layout)) {
          layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                           ncol = cols, nrow = ceiling(numPlots/cols))
     }
     
     if (numPlots == 1) {
          print(plots[[1]])
          
     } else {
          grid.newpage()
          pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
          
          for (i in 1:numPlots) {
               matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
               
               print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                               layout.pos.col = matchidx$col))
          }
     }
}
multiplot(p1,p2, cols=1)


par(mfrow=c(1,1))

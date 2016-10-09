##Angie messing with diel Ottesen data ---------------------------------------
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stats)
library(cowplot)
library(graphics)

diel.metabolites <- read_csv("~/UW lab/KM1513/Diel/Cyano_Aq_diel_1_Samples_Normalized_BlkSub_backfill_2IS_noblk_new1008.csv")
diel.sub1 <- diel.metabolites[,c(-1:-2)] 
diel.sub2 <- diel.sub1 %>%
     filter(type!="Poo") %>%
     filter(SampID!="FilterBlk") %>%
     select(Compound.Name, PooPlusModel, SampID) %>%
     mutate(SampID = as.numeric(SampID)) %>%
     mutate(Time = ifelse(SampID==1, 7, SampID*4+2)) %>%
     select(-SampID) %>%
     rename(Relative.Concentration.Data = PooPlusModel)

times <- c()
all.model.data <- c()
compounds <- unique(diel.sub2$Compound.Name)
compounds.for.modeling <- compounds[c(1,2,3,4,5,8,11,12,13,14,15,16,17,
                                      18,19,20,21,22,23,24,25,26,28,
                                      29,32)]
pdf(paste(Sys.Date(), "Aq_cyano_cmpds_and_models.pdf"), 6,4)

for (k in 1:length(compounds.for.modeling)) {
     
     start.time <- Sys.time()
     cmpd.name <- compounds.for.modeling[k]
     
     diel.cmpdspecific <- diel.sub2 %>%
          filter(Compound.Name ==cmpd.name)
     
     ggplot(diel.cmpdspecific, aes(x=Time, y=Relative.Concentration.Data)) + geom_point() +
          ggtitle(cmpd.name)
     
     glm.1 <- glm(diel.cmpdspecific$Relative.Concentration.Data ~ cos(2*pi/24*diel.cmpdspecific$Time) + 
                       sin(2*pi/24*diel.cmpdspecific$Time))
     
     
     glm.1.wlinear <- glm(diel.cmpdspecific$Relative.Concentration.Data ~ cos(2*pi/24*diel.cmpdspecific$Time) + 
                       sin(2*pi/24*diel.cmpdspecific$Time) + diel.cmpdspecific$Time)
     
     
     amplitude <- glm.1$coefficients[2]/cos(atan(-glm.1$coefficients[3]/glm.1$coefficients[2]))
     w.phase <- atan(-glm.1$coefficients[3]/glm.1$coefficients[2])
     pseudoR2 <- 1-(glm.1$deviance/glm.1$null.deviance)
     
     predicted.peak <- w.phase * 24/(2*pi)+12
     
     dev.diff <- glm.1$null.deviance-glm.1$deviance
     
     diel.cmpdspecific$fitted <- fitted(glm.1)
     diel.cmpdspecific.fitted <- diel.cmpdspecific %>% 
          gather(dat.or.mod, value, -Compound.Name, -Time)
     
     diel.cmpdspecific.fitted.aves <- diel.cmpdspecific.fitted %>%
          group_by(Time,dat.or.mod) %>%
          summarise(mean = mean(value),
                    sd = sd(value))
     
     just.data <- diel.cmpdspecific.fitted.aves %>%
          filter(dat.or.mod!="fitted")
     p <- ggplot(just.data, aes(x=Time, y = mean, 
                                color = dat.or.mod, group = dat.or.mod)) +
          geom_point() + geom_line(size = 1) + xlab("Time (hours)") +
          ylab("Concentration (arbitrary units)") +
          scale_color_manual(values=c("black")) +
          geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd),alpha=0.1) + 
          ggtitle(cmpd.name) +
          scale_x_continuous(breaks=c(12,24,36,48,60,72,84,96))
     print(p 
           + annotate("rect",xmin=7+12,xmax=6+24,ymin=-Inf,ymax=Inf,
                      alpha=0.1,fill="black")
           + annotate("rect",xmin=7+12+24,xmax=6+48,ymin=-Inf,ymax=Inf,
                      alpha=0.1,fill="black")
           + annotate("rect",xmin=7+12+48,xmax=6+24+48,ymin=-Inf,ymax=Inf,
                      alpha=0.1,fill="black")
           + annotate("rect",xmin=7+12+72,xmax=Inf,ymin=-Inf,ymax=Inf,
                      alpha=0.1,fill="black"))
     
     g <- ggplot(diel.cmpdspecific.fitted.aves, aes(x=Time, y = mean, 
                                          color = dat.or.mod, group = dat.or.mod)) +
          geom_point() + geom_line(size = 1) + xlab("Time (hours)") +
          ylab("Concentration (arbitrary units)") +
          scale_color_manual(values=c("red","black")) +
               geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd),alpha=0.1) + 
               ggtitle(cmpd.name) +
          scale_x_continuous(breaks=c(12,24,36,48,60,72,84,96))
     
     print(g 
           + annotate("rect",xmin=7+12,xmax=6+24,ymin=-Inf,ymax=Inf,
                      alpha=0.1,fill="black")
           + annotate("rect",xmin=7+12+24,xmax=6+48,ymin=-Inf,ymax=Inf,
                      alpha=0.1,fill="black")
           + annotate("rect",xmin=7+12+48,xmax=6+24+48,ymin=-Inf,ymax=Inf,
                      alpha=0.1,fill="black")
           + annotate("rect",xmin=7+12+72,xmax=Inf,ymin=-Inf,ymax=Inf,
                      alpha=0.1,fill="black")
           # +
                # geom_vline(xintercept = c(24,48,72,96),linetype=2) +
                # geom_vline(xintercept = c(12,36,60,84),linetype=3)
           # + annotate("text", x = -Inf, y = Inf, 
                        # hjust=-0.2, vjust =1,
                        # label = paste("Predicted peak time: ~", round(predicted.peak)*100)
                        # parse=TRUE,
                        # label = "alpha * cos(frac(2 * pi, 24) * t) -
                        # beta * sin(frac(2 * pi, 24) * t)"
                        # )
           )
     
     chi.sq <- anova(glm.1, test = "Chisq")
     chi.sq$`Pr(>Chi)`
     
     plot(glm.1$fitted.values~diel.cmpdspecific$Relative.Concentration.Data)
     abline(0,1,col="darkgray")
     ## Making models with randomized data --------------
     simple.Relative.Concentration.Data <- diel.cmpdspecific %>% select(Relative.Concentration.Data)
     simple.time <- diel.cmpdspecific %>% select(Time)
     
     # set.seed(NULL)

     i=0
     n=0

     while (n<10 | i<500) {
          randorder <- runif(nrow(simple.time), 0, 1)
          simple.time1 <- cbind(simple.time,randorder)
          simple.time2 <- simple.time1[order(randorder),]
          model.data <- cbind(simple.time2,simple.Relative.Concentration.Data)

          # ggplot(model.data, aes(x=Time, y=Relative.Concentration.Data)) + geom_point()

          glm.2 <- glm(model.data$Relative.Concentration.Data ~ cos(2*pi/24*model.data$Time) +
                            sin(2*pi/24*model.data$Time))

          temp.pseudoR2 <- 1-(glm.2$deviance/glm.2$null.deviance)

          temp.dev.diff <- glm.2$null.deviance-glm.2$deviance

          if (temp.dev.diff>dev.diff){
               n=n+1
          }

          if (i > 10000) break

          i=i+1

     }

     perm.pval <- n/i

     this.cmpd.data <- c(amplitude, dev.diff,w.phase, predicted.peak, pseudoR2,chi.sq$`Pr(>Chi)`,
                         i,n,perm.pval)
     names(this.cmpd.data) <- c("amplitude","dev.diff", "w.phase",
                                "predicted.peak",
                                "pseudoR2","chi.sq.Pr.Null",
                                "chi.sq.Pr.cos","chi.sq.Pr.Sin",
                                "iterations","n.better","perm.pval")

     all.model.data <- cbind(all.model.data, this.cmpd.data)
     colnames(all.model.data)[ncol(all.model.data)] <- cmpd.name

     end.time <- Sys.time()
     time.diff <- end.time - start.time
     times <- c(times, time.diff)
}
dev.off()

write.csv(all.model.data, "Cyano_aq_models_noextraNorm.csv")


##From Jacob------------- not currently working
library(mgcv)
diel.cmpdspecific$hod <- diel.cmpdspecific$Time %% 24
######
#knots are how long our time interval
mod <- gamm(Relative.Concentration.Data ~ s(hod, bs = "cc") + s(Time, bs = "ts"), 
            data = diel.cmpdspecific,
            correlation = corCAR1(form = ~Time), 
            knots = list(hod = c(0, 24)), 
            control =lmeControl( opt = "optim", msMaxIter = 1000, msVerbose = T))
########
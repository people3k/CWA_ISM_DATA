##SI Argentina Population Expansion Global Taphonomic Adjustment

library(rcarbon)
library(ggplot2)
library(nlme)
library(cowplot)
library(plyr)
library(tidyverse)
#library(tidyr)

#####Running KDEs by regions--North, Center and South
###Load data
FullData<-read.csv(file="S1File.csv", header=T) ##(Supplemental Table 1)
SPD<-subset(FullData, Area=="North")

##calibrate dates from the North region of CW Argentina
CalMz <- calibrate(x = SPD$Age,  errors = SPD$Error, calCurves = "shcal20",  normalised = FALSE)
boxbins <- binPrep(sites = SPD$SiteName, ages = SPD$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(4000,100))
plot(spd.mz, runm=200, xlim=c(4000,100), type="simple")

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##KDE
####Attempt KDE
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(4000,100),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

Ntaph<-transformSPD(D.ckde)
#Ntaph2<-transformSPD(D.ckde, correction = expression(PrDens/(21149.57*(calBP+1788.03)^-1.26)))

plot(Ntaph,type='multiline')
#plot(Ntaph2,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)
Check2<-as.data.frame(Ntaph$res.matrix)

Check3 <- replace(Check2, is.na(Check2), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Check3)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Check3, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Check3, #select columns
            1, # row-wise calcs
            quantile, probs=0.95) # give `quantile`

##Cbind spd and KDEs and write
dd<-cbind(calBP,PrDens, Check3, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd4<-dd %>%  filter(MKDE >0)
write.table(dd4, file = "SI/NorthTaph.csv", sep = ",", col.names=NA)

##Import taphonomically adjusted data for the Northern area
dd2n<- read.csv("SI/NorthTaph.csv")

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKDn <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dd2n,control=nls.control(maxiter=700),start=list(A=.00005, B=0.00015, xmid=2000, scale=-100))
summary(logFitKDn)
predictn<-predict(logFitKDn)

ntc<-ggplot(dd2n,aes(x=(calBP), y=(MKDE*10000))) +
  geom_ribbon(aes(ymin = lo*10000, ymax = hi*10000), fill = "grey80") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF" ) +
  geom_line(aes(y=predictn*10000), color="Blue", size=2, linetype = "dashed" ) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Adj. KDE of radiocarbon ages", title = "B. Northern Taph. Adj. KDE and Logistic Model")+
  geom_vline(xintercept = 2400)+
  geom_vline(xintercept = 1300)+
  geom_vline(xintercept = 500)
  #annotate("text", x =3500, y = .075, label = "Phase 1", size = 6 )+
  #annotate("text", x =2000, y = .075, label = "Phase 2", size = 6)+
  #annotate("text", x =900, y = .075, label = "Phase 3", size = 6)+
  #annotate("text", x =310, y = .075, label = "Phase 4", size = 6)
ntc

##################CENTRAL AREA--------------------------------------------------------------
SPD<-subset(FullData, Area=="Center")

##calibrate dates from the North region of CW Argentina
CalMz <- calibrate(x = SPD$Age,  errors = SPD$Error, calCurves = "shcal20",  normalised = FALSE)
boxbins <- binPrep(sites = SPD$SiteName, ages = SPD$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(4000,100))
plot(spd.mz, runm=200, xlim=c(4000,100), type="simple")

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##KDE
####Attempt KDE
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(4000,100),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

Ctaph<-transformSPD(D.ckde)
#Ntaph2<-transformSPD(D.ckde, correction = expression(PrDens/(21149.57*(calBP+1788.03)^-1.26)))

plot(Ctaph,type='multiline')
#plot(Ntaph2,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)
Check2<-as.data.frame(Ctaph$res.matrix)

Check3 <- replace(Check2, is.na(Check2), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Check3)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Check3, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Check3, #select columns
            1, # row-wise calcs
            quantile, probs=0.95) # give `quantile`

##Cbind spd and KDEs and write
dd<-cbind(calBP,PrDens, Check3, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd3<-dd %>%  filter(MKDE >0)
write.table(dd3, file = "SI/CenterTaph.csv", sep = ",", col.names=NA)

dd2c2<- read.csv("SI/CenterTaph.csv")

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKDc2 <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dd2c2,control=nls.control(maxiter=700),start=list(A=.000005, B=0.000015, xmid=2000, scale=-25))
summary(logFitKDc2)
predictc2<-predict(logFitKDc2)

ctc<-ggplot(dd2c2,aes(x=(calBP), y=(MKDE*10000))) +
  geom_ribbon(aes(ymin = lo*10000, ymax = hi*10000), fill = "grey80") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF" ) +
  geom_line(aes(y=predictc2*10000), color="Blue", size=2, linetype = "dashed" ) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Adj. KDE of radiocarbon ages", title = "B. Central Taph. Adj. KDE and Logistic Model")+
  geom_vline(xintercept = 2400)+
  geom_vline(xintercept = 1300)+
  geom_vline(xintercept = 500)
#annotate("text", x =3500, y = .06, label = "Phase 1", size = 6 )+
#annotate("text", x =2000, y = .06, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .06, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .06, label = "Phase 4", size = 6)
ctc

###

##################South AREA--------------------------------------------------------------
SPD<-subset(FullData, Area=="South")

##calibrate dates from the North region of CW Argentina
CalMz <- calibrate(x = SPD$Age,  errors = SPD$Error, calCurves = "shcal20",  normalised = FALSE)
boxbins <- binPrep(sites = SPD$SiteName, ages = SPD$Age, h = 100)

####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(4000,100))
plot(spd.mz, runm=200, xlim=c(4000,100), type="simple")

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##KDE
####Attempt KDE
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(4000,100),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')

Ntaph<-transformSPD(D.ckde)
#Ntaph2<-transformSPD(D.ckde, correction = expression(PrDens/(21149.57*(calBP+1788.03)^-1.26)))

plot(Ntaph,type='multiline')
#plot(Ntaph2,type='multiline')

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)
Check2<-as.data.frame(Ntaph$res.matrix)

Check3 <- replace(Check2, is.na(Check2), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Check3)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Check3, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Check3, #select columns
            1, # row-wise calcs
            quantile, probs=0.95) # give `quantile`

##Cbind spd and KDEs and write
dd<-cbind(calBP,PrDens, Check3, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd5<-dd %>%  filter(MKDE >0)
write.table(dd, file = "SI/SouthTaph.csv", sep = ",", col.names=NA)

###Load south taphonomically adjusted data
dds2<- read.csv("SI/SouthTaph.csv")

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKDs2 <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dds2,control=nls.control(maxiter=700),start=list(A=.000005, B=0.000015, xmid=2000, scale=-50))
summary(logFitKDs2)
predicts2<-predict(logFitKDs2)

stc<-ggplot(dds2,aes(x=(calBP), y=(MKDE*10000))) +
  geom_ribbon(aes(ymin = lo*10000, ymax = hi*10000), fill = "grey80") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF" ) +
  geom_line(aes(y=predicts2*10000), color="Blue", size=2, linetype = "dashed" ) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Adj. KDE of radiocarbon ages", title = "B. South Taph. Adj. KDE and Logistic Model")+
  geom_vline(xintercept = 2400)+
  geom_vline(xintercept = 1300)+
  geom_vline(xintercept = 500)
#annotate("text", x =3500, y = .06, label = "Phase 1", size = 6 )+
#annotate("text", x =2000, y = .06, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .06, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .06, label = "Phase 4", size = 6)
stc


###At this point, one can combine plots from the un-adjusted and taphonomically adjusted KDE 
#Graph the unadjusted KDE in S3File_MendozaV1.2. Then combine the plots. One can do the same for the
#Northern and Central areas as well.
library(cowplot)
Fig5Rev<-plot_grid(p1south,stc, ncol=1, align="hv", axis = "rl")
Fig5Rev

##Export the plot.
pdf("SI/TaphSouth.pdf", width=12, height=14.55)
Fig5Rev
dev.off()





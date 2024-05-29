library(rcarbon)
library(ggplot2)
library(nlme)
library(cowplot)
library(tidyr)
library(conover.test)
library(viridis)
library(dplyr)
library(tidyverse)
library(purrr)
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

##Check the effect of h function on SPD if you desire
#binsense(x=CalMz,y=SPD$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####make KDEs
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(4000,100),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
D.ckde$timeRange

##Write matrix of KDEs as a data frame
Check<-as.data.frame(D.ckde$res.matrix)
#I then convert NAs at the end due to KDE bandwidth to zeros to easily cbind the data frames below
Check2 <- replace(Check, is.na(Check), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Check2)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Check2, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Check2, #select columns
         1, # row-wise calcs
         quantile, probs=0.95) # give `quantile`

calBP<-spd.mz$grid$calBP
PrDens<-spd.mz$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
dd<-cbind(calBP,PrDens, Check, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2<-dd %>%  filter(MKDE >0)
##Write the table
write.table(dd2, file = "Updated/NorthKDE50bin.csv", sep = ",", col.names=NA)

#load North KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("Updated/NorthKDE50bin.csv") %>%
dplyr::select(-X,-calBP,-PrDens, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps..........
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
##Add in the 30 year bin dates
calBP<-c(3970, 3940,3910,3880, 3850, 3820, 3790, 3760, 3730, 3700, 3670, 3640, 3610, 3580, 3550, 3520, 3490, 3460, 3430,
         3400, 3370, 3340, 3310, 3280, 3250, 3220, 3190, 3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950, 2920, 2890,
         2860, 2830, 2800, 2770, 2740, 2710, 2680, 2650, 2620, 2590, 2560, 2530, 2500, 2470, 2440, 2410, 2380, 2350, 2320,
         2290, 2260, 2230, 2200, 2170, 2140, 2110, 2080, 2050, 2020, 1990, 1960, 1930, 1900, 1870, 1840, 1810, 1780, 1750, 1720,
         1690, 1660, 1630, 1600, 1570, 1540, 1510, 1480, 1450,
         1420, 1390, 1360, 1330, 1300, 1270, 1240, 1210, 1180, 1150, 1120, 1090, 1060, 1030, 1000, 970, 940, 910,
         880, 850, 820, 790, 760, 730, 700, 670, 640, 610, 580, 550, 520, 490, 460, 430, 400, 370, 340, 310, 280, 250, 220)

sums<-cbind(calBP, MKDE, out50)
write.table(sums, file = "Updated/North30Sumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("Updated/North30Sumbin.csv") %>%
 dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
write.table(pcgrowth, file = "Updated/North30PerCap.csv", sep = ",", col.names=NA)


#Load KDE data frame with calculated mean KDE and per capita growth
dd2<- read.csv("Updated/NorthKDE50bin.csv")
#dd2<-na.omit(dd2)

###transform dd2 into long table format to plot all of the simulated KDEs
#dd3<-dd2 %>% pivot_longer(cols=c(V1:V200),names_to='Runs', values_to='KDE')
#write.table(dd3, file = "Mendoza/KDELong.csv", sep = ",")
#dd3<- read.csv("Mendoza/KDELong.csv")

##Now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKD <- nls(MKDE~SSfpl(calBP, A, B, xmid, scale), data=dd2,control=nls.control(maxiter=700),start=list(A=.00005, B=0.00015, xmid=2000, scale=-100))
summary(logFitKD)
predict<-predict(logFitKD)

##Combine predicted logistic values and KDE table for plotting
dd4<-cbind(dd2, predict)

###Simple plot
plot(dd2$MKDE~dd2$calBP)
lines(dd2$calBP, dd2$MKDE,col="green",lty=2,lwd=3)
lines(dd2$calBP,predict(logFitKD),col="red",lty=2,lwd=3)

###Plot mean KDE and the 95th and 5th percentiles of KDE models at each time step using ggplot

#p1north<-ggplot(data = dd3) +
 # geom_line(aes(x = calBP, y = (KDE*100), group = Runs), color = "purple", alpha = 0.5) +
#  geom_line(data = dd4, aes(x = calBP, y = (predict*100)), linetype = "dashed", color="green", size=1.5) +
 # geom_point(data = dd2, aes(x = calBP, y = (MKDE*100)), color="black", size=2) +
  #scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(4000,200))+
  #theme_bw() +
  #theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
   #     axis.title.y=element_text(size=24), axis.text.y = element_text(
    #      size=28), plot.title = element_text(size=18, face = "bold"))+
  #labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "A. North Mendoza KDE and Logistic Model")+
 # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  #geom_vline(xintercept = 2400)+
  #geom_vline(xintercept = 1300)+
  #geom_vline(xintercept = 500)
#p1north


p1north <- ggplot(dd2,aes(x=(calBP), y=(MKDE*100))) +
  geom_ribbon(aes(ymin = lo*100, ymax = hi*100), fill = "grey80") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF" ) +
  geom_line(aes(y=predict*100), color="Blue", size=2, linetype = "dashed" ) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "A. Northern Area KDE and Logistic Model")+
  geom_vline(xintercept = 2400)+
  geom_vline(xintercept = 1300)+
  geom_vline(xintercept = 500)+
  annotate("text", x =3500, y = .06, label = "Phase 1", size = 6 )+
  annotate("text", x =2000, y = .06, label = "Phase 2", size = 6)+
  annotate("text", x =900, y = .06, label = "Phase 3", size = 6)+
  annotate("text", x =310, y = .06, label = "Phase 4", size = 6)
p1north


###Plot per capita growth of mean KDE against time
pc30ct<- read.csv("Updated/North30PerCap.csv")

Npc <- ggplot(pc30ct,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
 scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "B. Northern Area KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 2400)+
  geom_vline(xintercept = 1300)+
  geom_vline(xintercept = 500)+
  annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
  annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
  annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
  annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Npc

###Plot mean KDE against the per capita growth rate in the North

Npc2 <- ggplot(pc30ct,aes(x=(MKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=factor(PeriodID)), size=2.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
 # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
 # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Mean KDE (density)", y="KDE per capita growth", title = "A. Northern Area KDE Per Capita Growth vs. Density")
# geom_hline(yintercept = 0)+
  #geom_vline(xintercept = 2400)+
  #geom_vline(xintercept = 1300)+
  #geom_vline(xintercept = 500)+
  #annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
  #annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
  #annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
  #annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Npc2


####Now that we have made plots for the North region, let's do the same for the Center region
#First, subset the data to the Center Region. I am going to add a ``c" to all named objects to denote
#the center region.

SPDc<-subset(FullData, Area=="Center")

##calibrate dates from lower 48 states
CalMzc <- calibrate(x = SPDc$Age,  errors = SPDc$Error, calCurves = "shcal20",  normalised = FALSE)
boxbins <- binPrep(sites = SPDc$SiteName, ages = SPDc$Age, h = 100)
####Run SPD
spd.mzc <- spd(CalMzc, bins=boxbins, runm=200, timeRange=c(4000,100))
plot(spd.mzc, runm=200, xlim=c(4000,100), type="simple")

#binsense(x=CalMzc,y=SPDc$SiteName,h=seq(0,500,100),timeRange=c(4000,100))

##KDE
####Make KDEs
US.randatesc = sampleDates(CalMzc, bins=boxbins, nsim=200,verbose=FALSE)
D.ckdec = ckde(US.randatesc,timeRange=c(4000,100),bw=50, normalised = FALSE)
plot(D.ckdec,type='multiline')
D.ckdec$timeRange

##Write matrix of KDEs as a data frame
Checkc<-as.data.frame(D.ckdec$res.matrix)
#I then convert NAs at the end due to KDE bandwidth to zeros to easily cbind the data frames below
Checkc2 <- replace(Checkc, is.na(Checkc), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Checkc2)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Checkc2, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Checkc2, #select columns
            1, # row-wise calcs
            quantile, probs=0.95) # give `quantile`
#save cal BP and prdens of the SPD
calBPc<-spd.mzc$grid$calBP
PrDensc<-spd.mzc$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
ddc<-cbind(calBPc,PrDensc, Checkc, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2c<-ddc %>%  filter(MKDE >0)
##Write the table
write.table(dd2c, file = "Updated/CenterKDE50bin.csv", sep = ",", col.names=NA)

#load Center KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("Updated/CenterKDE50bin.csv") %>%
  dplyr::select(-X,-calBPc,-PrDensc, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps.
###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
##Add in the 30 year bin dates
calBP<-c(3970, 3940,3910,3880, 3850, 3820, 3790, 3760, 3730, 3700, 3670, 3640, 3610, 3580, 3550, 3520, 3490, 3460, 3430,
         3400, 3370, 3340, 3310, 3280, 3250, 3220, 3190, 3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950, 2920, 2890,
         2860, 2830, 2800, 2770, 2740, 2710, 2680, 2650, 2620, 2590, 2560, 2530, 2500, 2470, 2440, 2410, 2380, 2350, 2320,
         2290, 2260, 2230, 2200, 2170, 2140, 2110, 2080, 2050, 2020, 1990, 1960, 1930, 1900, 1870, 1840, 1810, 1780, 1750, 1720,
         1690, 1660, 1630, 1600, 1570, 1540, 1510, 1480, 1450,
         1420, 1390, 1360, 1330, 1300, 1270, 1240, 1210, 1180, 1150, 1120, 1090, 1060, 1030, 1000, 970, 940, 910,
         880, 850, 820, 790, 760, 730, 700, 670, 640, 610, 580, 550, 520, 490, 460, 430, 400, 370, 340, 310, 280, 250, 220)

sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "Updated/Center30Sumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("Updated/Center30Sumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
write.table(pcgrowth, file = "Updated/Center30PerCap.csv", sep = ",", col.names=NA)

#Load KDE data frame with calculated mean KDE at annual resolution 
dd2c<- read.csv("Updated/CenterKDE50bin.csv")

###transform dd2 into long table format if you want to plot all of the simulated KDEs
#dd3c<-dd2c %>% pivot_longer(cols=c(V1:V200),names_to='Runs', values_to='KDE')
#write.table(dd3, file = "Mendoza/KDELong.csv", sep = ",")
#dd3<- read.csv("Mendoza/KDELong.csv")

##Now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKDc <- nls(MKDE~SSfpl(calBPc, A, B, xmid, scale), data=dd2c,control=nls.control(maxiter=700),start=list(A=.000001, B=0.000035, xmid=2000, scale=-50))
summary(logFitKDc)
predictc<-predict(logFitKDc)

##Combine predicted logistic values and KDE table for plotting
dd4c<-cbind(dd2c, predictc)

###Simple plot
plot(dd2c$MKDE~dd2c$calBP)
lines(dd2c$calBPc, dd2c$MKDE,col="green",lty=2,lwd=3)
lines(dd2c$calBPc,predict(logFitKDc),col="red",lty=2,lwd=3)

###Plot mean KDE and 95 and 5th percentile confidence envelope

p1center <- ggplot(dd2c,aes(x=(calBPc), y=(MKDE*100))) +
  geom_ribbon(aes(ymin = lo*100, ymax = hi*100), fill = "grey80") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF" ) +
  geom_line(aes(y=predictc*100), color="Blue", size=2, linetype = "dashed" ) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "A. Central Area KDE and Logistic Model")+
  geom_vline(xintercept = 2400)+
  geom_vline(xintercept = 1300)+
  geom_vline(xintercept = 500)+
  annotate("text", x =3500, y = .06, label = "Phase 1", size = 6 )+
  annotate("text", x =2000, y = .06, label = "Phase 2", size = 6)+
  annotate("text", x =900, y = .06, label = "Phase 3", size = 6)+
  annotate("text", x =310, y = .06, label = "Phase 4", size = 6)
p1center


###Plot per capita growth of mean KDE against time
pc30ct2<- read.csv("Updated/Center30PerCap.csv")

Cpc <- ggplot(pc30ct2,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  #scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "B. Central Area KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 2400)+
  geom_vline(xintercept = 1300)+
  geom_vline(xintercept = 500)+
  annotate("text", x =3500, y = .5, label = "Phase 1", size = 6)+
  annotate("text", x =2000, y = .5, label = "Phase 2", size = 6)+
  annotate("text", x =900, y = .5, label = "Phase 3", size = 6)+
  annotate("text", x =310, y = .5, label = "Phase 4", size = 6)
Cpc

###Plot the mean KDE against the per capita growth to create a phase space.
Cpc2 <- ggplot(pc30ct2,aes(x=(MKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=factor(PeriodID)), size=2.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Mean KDE (density)", y="KDE per capita growth", title = "B. Central Area KDE Per Capita Growth vs. Density")
# geom_hline(yintercept = 0)+
#geom_vline(xintercept = 2400)+
#geom_vline(xintercept = 1300)+
#geom_vline(xintercept = 500)+
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Cpc2

###South Region===================================================================
SPD<-subset(FullData, Area=="South")

##calibrate dates Southern Mendoza
CalMz <- calibrate(x = SPD$Age,  errors = SPD$Error, calCurves = "shcal20",  normalised = FALSE)
boxbins <- binPrep(sites = SPD$SiteName, ages = SPD$Age, h = 100)
####Run SPD
spd.mz <- spd(CalMz, bins=boxbins, runm=200, timeRange=c(4000,100))
plot(spd.mz, runm=200, xlim=c(4000,100), type="simple")

##KDE
####Attempt KDE
US.randates = sampleDates(CalMz, bins=boxbins, nsim=200,verbose=FALSE)
D.ckde = ckde(US.randates,timeRange=c(4000,100),bw=50, normalised = FALSE)
plot(D.ckde,type='multiline')
D.ckde$timeRange

##Write matrix of KDEs as a data frame
Checks<-as.data.frame(D.ckde$res.matrix)
#I then convert NAs at the end due to KDE bandwidth to zeros to easily cbind the data frames below
Checks2 <- replace(Checks, is.na(Checks), 0)
#Calculate the mean KDE of each time step of the 200 simulations
MKDE<-rowMeans(Checks2)
#calculate the 5th and 95th percentile of each time step of the 200 KDE simulations
lo<-apply( Checks2, #select columns
           1, # row-wise calcs
           quantile, probs=0.05) # give `quantile`
hi <-apply( Checks2, #select columns
            1, # row-wise calcs
            quantile, probs=0.95) # give `quantile`
calBPs<-spd.mzc$grid$calBP
PrDenss<-spd.mzc$grid$PrDens

##Cbind spd and KDEs, mean KDE, percentiles and write
ddc<-cbind(calBPs,PrDenss, Checks, MKDE, hi, lo)
##Remove end rows with zeros due to undefined KDE values at a 50 bandwidth at the end of the sequence
dd2c<-ddc %>%  filter(MKDE >0)
##Write the table
write.table(dd2c, file = "Updated/SouthKDE50bin.csv", sep = ",", col.names=NA)

#load Center KDE data set and select columns for removal that we do not want to sum 
dd2c<- read.csv("Updated/SouthKDE50bin.csv") %>%
  dplyr::select(-X,-calBPs,-PrDenss, -MKDE,-hi,-lo)

### Sum into 30 year generation time steps.
###Sum SPD by 30 year intervals
library(zoo)
# sum and save new csvs.
out50 <- rollapply(dd2c,30,(sum),by=30,by.column=TRUE,align='right')

###Calculate the mean KDE of the 30 year sums of the 200 KDEs
MKDE<-rowMeans(out50)
##Add in the 30 year bin dates
calBP<-c(3970, 3940,3910,3880, 3850, 3820, 3790, 3760, 3730, 3700, 3670, 3640, 3610, 3580, 3550, 3520, 3490, 3460, 3430,
         3400, 3370, 3340, 3310, 3280, 3250, 3220, 3190, 3160, 3130, 3100, 3070, 3040, 3010, 2980, 2950, 2920, 2890,
         2860, 2830, 2800, 2770, 2740, 2710, 2680, 2650, 2620, 2590, 2560, 2530, 2500, 2470, 2440, 2410, 2380, 2350, 2320,
         2290, 2260, 2230, 2200, 2170, 2140, 2110, 2080, 2050, 2020, 1990, 1960, 1930, 1900, 1870, 1840, 1810, 1780, 1750, 1720,
         1690, 1660, 1630, 1600, 1570, 1540, 1510, 1480, 1450,
         1420, 1390, 1360, 1330, 1300, 1270, 1240, 1210, 1180, 1150, 1120, 1090, 1060, 1030, 1000, 970, 940, 910,
         880, 850, 820, 790, 760, 730, 700, 670, 640, 610, 580, 550, 520, 490, 460, 430, 400, 370, 340, 310, 280, 250, 220)

sums<-cbind(calBP, MKDE, out50)

write.table(sums, file = "Updated/South30Sumbin.csv", sep = ",", col.names=NA)

###Calculate per capita growth rate of 30 year time steps for each of the 200 simulations.

d <-read_csv("Updated/South30Sumbin.csv") %>%
  dplyr::select(-...1)
d2<-arrange(d, calBP)
### We calculate capita growth as LN(MKDE at t+1/ MKDE at time t).
pcgrowth<-d2 %>% mutate(PerCap = (log(lag(MKDE)/MKDE)))
write.table(pcgrowth, file = "Updated/South30PerCap.csv", sep = ",", col.names=NA)


#Load KDE data frame with calculated mean KDE and per capita growth
dd2s<- read.csv("Updated/SouthKDE50bin.csv")

###transform dd2 into long table format to plot all of the simulated KDEs
#dd3s<-dd2s %>% pivot_longer(cols=c(V1:V200),names_to='Runs', values_to='KDE')
#write.table(dd3, file = "Mendoza/KDELong.csv", sep = ",")
#dd3<- read.csv("Mendoza/KDELong.csv")

##OK, now we fit a logistic model to the mean KDE. This can be a bit tricky due to the start parameters.
#These need to be adjusted based on the MKDE's values. In this case, the start values work well.
#fit logistic to KDE
logFitKD <- nls(MKDE~SSfpl(calBPs, A, B, xmid, scale), data=dd2s,control=nls.control(maxiter=700),start=list(A=.00005, B=0.00015, xmid=2000, scale=-100))
summary(logFitKD)
predict<-predict(logFitKD)

##Combine predicted logistic values and KDE table for plotting
dd4s<-cbind(dd2s, predict)

###Simple plot
plot(dd2s$MKDE~dd2s$calBP)
lines(dd2s$calBP, dd2$MKDE,col="green",lty=2,lwd=3)
lines(dd2s$calBP,predict(logFitKD),col="red",lty=2,lwd=3)

###Plot mean KDE and 5th and 95th percentiles of confidence envelope

p1south <- ggplot(dd2s,aes(x=(calBPs), y=(MKDE*100))) +
  geom_ribbon(aes(ymin = lo*100, ymax = hi*100), fill = "grey80") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF" ) +
  geom_line(aes(y=predict*100), color="Blue", size=2, linetype = "dashed" ) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE of radiocarbon ages", title = "A. Southern Area KDE and Logistic Model")+
  geom_vline(xintercept = 2400)+
  geom_vline(xintercept = 1300)+
  geom_vline(xintercept = 500)+
  annotate("text", x =3500, y = .06, label = "Phase 1", size = 6 )+
  annotate("text", x =2000, y = .06, label = "Phase 2", size = 6)+
  annotate("text", x =900, y = .06, label = "Phase 3", size = 6)+
  annotate("text", x =310, y = .06, label = "Phase 4", size = 6)
p1south

###Plot per capita growth of mean KDE against time
pc30ct3<- read.csv("Updated/South30PerCap.csv")

Spc <- ggplot(pc30ct3,aes(x=(calBP), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(), size=1) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  #scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="KDE per capita growth", title = "B. Southern Area KDE Per Capita Growth")+
  geom_line(size=1.25)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 2400)+
  geom_vline(xintercept = 1300)+
  geom_vline(xintercept = 500)+
  annotate("text", x =3500, y = .5, label = "Phase 1", size = 6)+
  annotate("text", x =2000, y = .5, label = "Phase 2", size = 6)+
  annotate("text", x =900, y = .5, label = "Phase 3", size = 6)+
  annotate("text", x =310, y = .5, label = "Phase 4", size = 6)
Spc

##Plot the mean KDE vs. per capita growth in the South area
Spc2 <- ggplot(pc30ct3,aes(x=(MKDE), y=(PerCap))) +
  #geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70") +
  geom_point(aes(y=PerCap, color=factor(PeriodID)), size=2.5) +
  geom_path(aes(),size=1)+
  #scale_color_gradient(low ="#F8766D", high = "#619CFF") +
  #scale_color_manual(values=c("#619CFF", "#00BA38", "#F8766D"))+
  #geom_line(aes(y=logFit3), color="blue", size=1) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3700,300))+
  # scale_y_continuous(limits=c(-.75,0.5))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Mean KDE (density)", y="KDE per capita growth", title = "C. Southern Area KDE Per Capita Growth vs. Density")
# geom_hline(yintercept = 0)+
#geom_vline(xintercept = 2400)+
#geom_vline(xintercept = 1300)+
#geom_vline(xintercept = 500)+
#annotate("text", x =3500, y = .25, label = "Phase 1", size = 6)+
#annotate("text", x =2000, y = .25, label = "Phase 2", size = 6)+
#annotate("text", x =900, y = .25, label = "Phase 3", size = 6)+
#annotate("text", x =310, y = .25, label = "Phase 4", size = 6)
Spc2

###Export plots. This example exports the mean KDE vs. per capita growth plots
Fig5Rev<-plot_grid(Npc2,Cpc2,Spc2, ncol=1, align="v", axis = "rl")
Fig5Rev

pdf("Updated/Phase.pdf", width=13, height=11.55)
Fig5Rev
dev.off()

###Isotopes=======================================================================================
d3 <- read.csv("S2File.csv")
d3c<-subset(d3, PeriodID>0  & d13Cca<0)

###Protein lines for graphs
C3Pro2<-0.555*d3c$d13Cca +-12.7
C4Pro2<-0.50*d3c$d13Cca + -7.61

##Graph carbonate against collagen for all regions for the S1 Appendix
p1 <- ggplot(d3c,aes(x=(d13Cca), y=(d13Ccol))) +
  geom_point(aes(color=(d15Ncorr)),size=4.5) +
  # scale_color_gradient(low ="#F8766D", high = "#619CFF" ) +
  scale_color_viridis()+
  theme_bw() +
  geom_line(aes(d13Cca, C3Pro2), color="Red", size=2, linetype = "dashed" )+
  geom_line(aes(d13Cca, C4Pro2), color="Green", size=2)+
  scale_x_continuous(breaks=c(-14,-11, -9, -7, -5))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Delta 13C carbonate", y="Delta 13C collagen", title = "")+
  geom_smooth(method="lm", se=FALSE)+
  facet_wrap(factor(Region)~.)
p1

p2 <- ggplot(d3c,aes(x=(d13Cca), y=(d13Ccol))) +
  geom_point(aes(color=factor(PeriodID)), size=4.5) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF" ) +
  # geom_line(aes(y=residuals), color="blue", size=.5) +
  theme_bw() +
  geom_line(aes(d13Cca, C3Pro2), color="Red", size=2, linetype = "dashed" )+
  #geom_line(aes(d13Ccol, Marine), color="Blue", size=2, linetype = "dashed")+
  geom_line(aes(d13Cca, C4Pro2), color="Green", size=2)+
  scale_x_continuous(breaks=c(-14,-11, -9, -7, -5))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Delta 13C carbonate", y="Delta 13C collagen", title = "")+
  geom_smooth(method="lm", se=FALSE)+
  facet_wrap(factor(Region)~.)
p2

pdf("TimeoWholeDiet3.pdf", width=13, height=11.55)
p1
dev.off()


###Violin plots by North region===================================
d3c<-subset(d3, PeriodID>0 & Region=="North" & d13Cca<0)
library(dplyr)
Means <- d3c %>% group_by(PeriodID) %>%
  summarize(Avg = median(d13Cca))

kruskal.test(d13Cca~ PeriodID, data = d3c)
conover.test(d3c$d13Cca, d3c$PeriodID, kw=TRUE, method="bonferroni")
#conover.test(d3c$d13Cca, d3c$PeriodID, kw=TRUE, method="by")

pn1 <- ggplot(d3c, aes(factor(PeriodID), (d13Cca), fill=factor(PeriodID)))+
  geom_violin()+
  #geom_boxplot(notch = FALSE, width=0.2)+
  stat_boxplot(geom ='errorbar')+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  geom_line(data = Means,mapping = aes(x = factor(PeriodID), y = Avg,group=0),
            color="black", size=1.2)+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Demographic phase", y="Delta 13C carbonate", title = "C. Delta 13C Carbonate by Demographic Phase" )+
  theme(legend.position="none")+
  annotate("text", x =1.5, y = -6, label = "-5.47**", size = 6)+
  annotate("text", x =2.5, y = -6, label = "1.05", size = 6)
#annotate("text", x =310, y = .06, label = "Phase 4", size = 6)
pn1


##Extra graph==================================
p2 <- ggplot(d3c,aes(x=(d13Cca), y=(d13Ccol))) +
  geom_point(aes(color=factor(PeriodID)), size=4.5) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF" ) +
  # geom_line(aes(y=residuals), color="blue", size=.5) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  #labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "A. Lower 48 Population Growth")+
  geom_smooth(method="lm", se=FALSE)+
  facet_wrap(factor(Region)~.)
p2

###Violin plots by Central region======================================================
d3cc<-subset(d3, PeriodID>0 & Region=="Center" & d13Cca<0)

kruskal.test(d13Cca~ PeriodID, data = d3cc)
conover.test(d3cc$d13Cca, d3cc$PeriodID, kw=TRUE, method="bonferroni")
#conover.test(d3c$d13Cca, d3c$PeriodID, kw=TRUE, method="by")

library(dplyr)
Meansc <- d3cc %>% group_by(PeriodID) %>%
  summarize(Avg = median(d13Cca))

pn1c <- ggplot(d3cc, aes(factor(PeriodID), (d13Cca), fill=factor(PeriodID)))+
  geom_violin()+
  #geom_boxplot(notch = FALSE, width=0.2)+
  stat_boxplot(geom ='errorbar')+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  geom_line(data = Meansc,mapping = aes(x = factor(PeriodID), y = Avg,group=0),
            color="black", size=1.2)+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Demographic phase", y="Delta 13C carbonate", title = "C. Delta 13C Carbonate by Demographic Phase" )+
  theme(legend.position="none")+
  annotate("text", x =1.5, y = -5, label = "-1.08", size = 6)+
  annotate("text", x =2.45, y = -5, label = "-3.92**", size = 6)+
  annotate("text", x =3.5, y = -5, label = "3.77**", size = 6)
pn1c

##Extra graph
p2c <- ggplot(d3cc,aes(x=(d13Cca), y=(d13Ccol))) +
  geom_point(aes(color=factor(PeriodID)), size=4.5) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF" ) +
  # geom_line(aes(y=residuals), color="blue", size=.5) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  #labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "A. Lower 48 Population Growth")+
  geom_smooth(method="lm", se=FALSE)+
  facet_wrap(factor(Region)~.)
p2c

###Violin plots by South region======================================================
d3cs<-subset(d3, PeriodID>0 & Region=="South" & d13Cca<0)

kruskal.test(d13Cca~ PeriodID, data = d3cs)
conover.test(d3cs$d13Cca, d3cs$PeriodID, kw=TRUE, method="bonferroni")
#conover.test(d3c$d13Cca, d3c$PeriodID, kw=TRUE, method="by")

library(dplyr)
Meanss <- d3cs %>% group_by(PeriodID) %>%
  summarize(Avg = median(d13Cca))

pn1s <- ggplot(d3cs, aes(factor(PeriodID), (d13Cca), fill=factor(PeriodID)))+
  geom_violin()+
  #geom_boxplot(notch = FALSE, width=0.2)+
  stat_boxplot(geom ='errorbar')+
  stat_summary(fun=median, geom="point", size=4, color="black")+
  #scale_fill_manual(values=c("#FC4E07", "#00A4CCFF"))+
  geom_line(data = Meanss,mapping = aes(x = factor(PeriodID), y = Avg,group=0),
            color="black", size=1.2)+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28),  plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Demographic phase", y="Delta 13C carbonate", title = "C. Delta 13C Carbonate by Demographic Phase" )+
  theme(legend.position="none")+
  annotate("text", x =1.5, y = -7, label = "1.35", size = 6)+
  annotate("text", x =2.5, y = -7, label = "0.65", size = 6)+
  annotate("text", x =3.5, y = -7, label = "1.12", size = 6)
pn1s


##Extra graph
p2s <- ggplot(d3cs,aes(x=(d13Cca), y=(d13Ccol))) +
  geom_point(aes(color=factor(PeriodID)), size=4.5) +
  #scale_color_gradient(low ="#F8766D", high = "#619CFF" ) +
  # geom_line(aes(y=residuals), color="blue", size=.5) +
  theme_bw() +
  # scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(3500,400))+
  # scale_y_continuous(limits=c(0,0.42))+
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  #labs(x = "Years Cal BP", y="SPD of radiocarbon ages", title = "A. Lower 48 Population Growth")+
  geom_smooth(method="lm", se=FALSE)+
  facet_wrap(factor(Region)~.)
p2s


###code to combine plots of area population, per capita growth, and isotopes
library(cowplot)

Fig5Rev<-plot_grid(p1north, Npc, pn1, ncol=1, align="hv", axis = "rl")
Fig5Rev

pdf("Updated/Northbin2.pdf", width=14, height=12.55)
Fig5Rev
dev.off()





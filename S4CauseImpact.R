###Causal Impact analysis=========================================================================
#1. Load libraries
library(ggplot2)
library(nlme)
library(cowplot)
library(CausalImpact)

#2. Load data
dd2sum<- read.csv("CausalImpact/North30Sumbin.csv")
dd2sum<-subset(dd2sum, calBP< 3730 & calBP>280)
dd2sum2<- read.csv("CausalImpact/Center30Sumbin.csv")
dd2sum2<-subset(dd2sum2, calBP< 3730 & calBP>280)
dd2sum3<- read.csv("CausalImpact/South30Sumbin.csv")
dd2sum3<-subset(dd2sum3,calBP< 3730 & calBP>280)



#3. Conduct CI analysis for each region. There are two periods of demographic transitions in each 
#region. We begin with the North area.

##North Area First period of population expansion evidenced by multiple potential waves of DTs 
#beginning at 2950 cal BP.
#Bind the mean KDE and modeled rainfall columns
CTXImpact<-cbind((dd2sum$MKDE), dd2sum$CRR)

#Set the training period, intervention point, and post period. In the North area
#the first training set is from 3700 to 2950 cal BP when we observe consecutive 
# potential demographic transitions and population expansion (see Fig. 5 in the main text).
pre.period1 <- c(1, 35)
post.period1 <- c(36,77)

##Run causal impact
impact1 <- CausalImpact(CTXImpact, pre.period1, post.period1)

##Basic results of the causal impact function
plot(impact1)
summary(impact1)
summary(impact1, "report")

####Plot the results of the counterfactual in ggplot

#dd2sumplot<-cbind(dd2sum,impact2$series$point.effect, hi, low )

pNa<-ggplot(data=dd2sum) +
  geom_ribbon(aes(x=calBP, ymin =impact1$series$point.pred.lower , ymax = impact1$series$point.pred.upper), fill = "grey70" )+
  geom_line(aes(x = calBP, y = (impact1$series$point.pred)), color="green", size=1.5) +
  geom_line(aes(x = calBP, y = (MKDE)), size=2) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(4000,200))+
  #scale_y_continuous(limits=c(-.4,0.4))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Actual and predicted MKDE", title = "A. Northern Area Population Counterfactual")+
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 2950)+
  geom_vline(xintercept = 1410)
#geom_vline(xintercept = 500)
# geom_hline(yintercept=0)
pNa

pNb<-ggplot(data=dd2sum) +
  geom_ribbon(aes(x=calBP, ymin =impact1$series$point.effect.lower , ymax = impact1$series$point.effect.upper), fill = "grey70" )+
  geom_line(aes(x = calBP, y = (impact1$series$point.effect)), color="black", size=1.5) +
  #geom_line(aes(x = calBP, y = (MKDE)), size=2) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(4000,200))+
  #scale_y_continuous(limits=c(-.4,0.4))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Pointwise effect (actual-predicted)", title = "B. Northern Area Pointwise Effects")+
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 2950)+
  geom_vline(xintercept = 1410)+
  geom_hline(yintercept=0)
pNb

###Run the same analysis for the potential demographic transition in the North
pre.period2 <- c(60, 77)
post.period2 <- c(78,114)

impact2 <- CausalImpact(CTXImpact, pre.period2, post.period2)

plot(impact2)
summary(impact2)
summary(impact2, "report")

pNa2<-ggplot(data=dd2sum) +
  geom_ribbon(aes(x=calBP, ymin =impact2$series$point.pred.lower , ymax = impact2$series$point.pred.upper), fill = "grey70" )+
  geom_line(aes(x = calBP, y = (impact2$series$point.pred)), color="green", size=1.5) +
  geom_line(aes(x = calBP, y = (MKDE)), size=2) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(4000,200))+
  #scale_y_continuous(limits=c(-.4,0.4))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Actual and predicted MKDE", title = "A. Northern Area Population Counterfactual")+
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 1930)+
  geom_vline(xintercept = 1320)
#geom_vline(xintercept = 500)
# geom_hline(yintercept=0)
pNa2

pNb2<-ggplot(data=dd2sum) +
  geom_ribbon(aes(x=calBP, ymin =impact2$series$point.effect.lower , ymax = impact2$series$point.effect.upper), fill = "grey70" )+
  geom_line(aes(x = calBP, y = (impact2$series$point.effect)), color="black", size=1.5) +
  #geom_line(aes(x = calBP, y = (MKDE)), size=2) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(4000,200))+
  #scale_y_continuous(limits=c(-.4,0.4))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Pointwise effect (actual-predicted)", title = "B. Northern Area Pointwise Effects")+
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 1930)+
  geom_vline(xintercept = 1320)+
  geom_hline(yintercept=0)
pNb2

###Central Region---------------------------------------------------
CImpact<-cbind((dd2sum2$MKDE), dd2sum2$CRR)


pre.period3 <- c(1, 45)
post.period3<- c(46,84)

impact3 <- CausalImpact(CImpact, pre.period3, post.period3)

plot(impact3)
summary(impact3)
summary(impact3, "report")

pCa<-ggplot(data=dd2sum2) +
  geom_ribbon(aes(x=calBP, ymin =impact3$series$point.pred.lower , ymax = impact3$series$point.pred.upper), fill = "grey70" )+
  geom_line(aes(x = calBP, y = (impact3$series$point.pred)), color="green", size=1.5) +
  geom_line(aes(x = calBP, y = (MKDE)), size=2) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(4000,200))+
  #scale_y_continuous(limits=c(-.4,0.4))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Actual and predicted MKDE", title = "C. Central Area Population Counterfactual")+
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 2400)+
  geom_vline(xintercept = 1210)
# geom_hline(yintercept=0)
pCa

pCb<-ggplot(data=dd2sum2) +
  geom_ribbon(aes(x=calBP, ymin =impact3$series$point.effect.lower , ymax = impact3$series$point.effect.upper), fill = "grey70" )+
  geom_line(aes(x = calBP, y = (impact3$series$point.effect)), color="black", size=1.5) +
  #geom_line(aes(x = calBP, y = (MKDE)), size=2) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(4000,200))+
  #scale_y_continuous(limits=c(-.4,0.4))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Pointwise effect (actual-predicted)", title = "D. Central Area Pointwise Effects")+
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 2400)+
  geom_vline(xintercept = 1210)+
  geom_hline(yintercept=0)
pCb

###Run CI analysis on the second suspected demographic transition between cycles
pre.period4 <- c(60, 84)
post.period4 <- c(85,114)

impact4 <- CausalImpact(CImpact, pre.period4, post.period4)

plot(impact4)
summary(impact4)
summary(impact4, "report")

pCa2<-ggplot(data=dd2sum2) +
  geom_ribbon(aes(x=calBP, ymin =impact4$series$point.pred.lower , ymax = impact4$series$point.pred.upper), fill = "grey70" )+
  geom_line(aes(x = calBP, y = (impact4$series$point.pred)), color="green", size=1.5) +
  geom_line(aes(x = calBP, y = (MKDE)), size=2) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(4000,200))+
  #scale_y_continuous(limits=c(-.4,0.4))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Actual and predicted MKDE", title = "C. Central Area Population Counterfactual")+
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 1930)+
  geom_vline(xintercept = 1210)
# geom_hline(yintercept=0)
pCa2

pCb2<-ggplot(data=dd2sum2) +
  geom_ribbon(aes(x=calBP, ymin =impact4$series$point.effect.lower , ymax = impact4$series$point.effect.upper), fill = "grey70" )+
  geom_line(aes(x = calBP, y = (impact4$series$point.effect)), color="black", size=1.5) +
  #geom_line(aes(x = calBP, y = (MKDE)), size=2) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(4000,200))+
  #scale_y_continuous(limits=c(-.4,0.4))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Pointwise effect (actual-predicted)", title = "D. Central Area Pointwise Effects")+
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 1930)+
  geom_vline(xintercept = 1210)+
  geom_hline(yintercept=0)
pCb2

##Southern Mendoza Causal Impact-------------------------------------------------
SImpact<-cbind((dd2sum3$MKDE), dd2sum3$CRR)

pre.period5 <- c(1, 55)
post.period5 <- c(56,91)

impact5 <- CausalImpact(SImpact, pre.period5, post.period5)

plot(impact5)
summary(impact5)
summary(impact5, "report")

pSa<-ggplot(data=dd2sum3) +
  geom_ribbon(aes(x=calBP, ymin =impact5$series$point.pred.lower , ymax = impact5$series$point.pred.upper), fill = "grey70" )+
  geom_line(aes(x = calBP, y = (impact5$series$point.pred)), color="green", size=1.5) +
  geom_line(aes(x = calBP, y = (MKDE)), size=2) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(4000,200))+
  #scale_y_continuous(limits=c(-.4,0.4))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Actual and predicted MKDE", title = "E. Southern Area Population Counterfactual")+
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 2030)+
  geom_vline(xintercept = 1000)
# geom_hline(yintercept=0)
pSa

pSb<-ggplot(data=dd2sum3) +
  geom_ribbon(aes(x=calBP, ymin =impact5$series$point.effect.lower , ymax = impact5$series$point.effect.upper), fill = "grey70" )+
  geom_line(aes(x = calBP, y = (impact5$series$point.effect)), color="black", size=1.5) +
  #geom_line(aes(x = calBP, y = (MKDE)), size=2) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(4000,200))+
  #scale_y_continuous(limits=c(-.4,0.4))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Pointwise effect (actual-predicted)", title = "F. Southern Area Pointwise Effects")+
  #annotate("text", x =3500, y = .01, label = "Pre-Intervention", size = 7)+
  #annotate("text", x =1500, y = .01, label = "Post-Intervention", size = 7)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 2030)+
  geom_vline(xintercept = 1000)+
  geom_hline(yintercept=0)
pSb

###run CI analysis for the second potential demographic transition
pre.period6 <- c(57, 91)
post.period6 <- c(92,114)

impact6 <- CausalImpact(SImpact, pre.period6, post.period6)

plot(impact6)
summary(impact6)
summary(impact6, "report")

pSa2<-ggplot(data=dd2sum3) +
  geom_ribbon(aes(x=calBP, ymin =impact6$series$point.pred.lower , ymax = impact6$series$point.pred.upper), fill = "grey70" )+
  geom_line(aes(x = calBP, y = (impact6$series$point.pred)), color="green", size=1.5) +
  geom_line(aes(x = calBP, y = (MKDE)), size=2) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(4000,200))+
  #scale_y_continuous(limits=c(-.4,0.4))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Actual and predicted MKDE", title = "E. Southern Area Population Counterfactual")+
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 2030)+
  geom_vline(xintercept = 1000)
# geom_hline(yintercept=0)
pSa2

pSb2<-ggplot(data=dd2sum3) +
  geom_ribbon(aes(x=calBP, ymin =impact6$series$point.effect.lower , ymax = impact6$series$point.effect.upper), fill = "grey70" )+
  geom_line(aes(x = calBP, y = (impact6$series$point.effect)), color="black", size=1.5) +
  #geom_line(aes(x = calBP, y = (MKDE)), size=2) +
  scale_x_reverse(breaks=c(3500, 2500, 1500, 500), limits=c(4000,200))+
  #scale_y_continuous(limits=c(-.4,0.4))+
  theme_bw() +
  theme(axis.text.x = element_text(size=28, colour = "black"), axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24), axis.text.y = element_text(
          size=28), plot.title = element_text(size=18, face = "bold"))+
  labs(x = "Years Cal BP", y="Pointwise effect (actual-predicted)", title = "F. Southern Area Pointwise Effects")+
  # annotate("text", x =3500, y = .047, label = "Phase 1", size = 7, angle=90)+
  #annotate("text", x =2000, y = .047, label = "Phase 2", size = 7, angle=90)+
  #annotate("text", x =1020, y = .047, label = "Phase 3", size = 7, angle=90)+
  #annotate("text", x =210, y = .047, label = "Phase 4", size = 7, angle=90)+
  #geom_vline(xintercept = 1746)+
  #geom_vline(xintercept = 1904)+
  # facet_wrap(factor(Region)~.)
  geom_vline(xintercept = 2030)+
  geom_vline(xintercept = 1000)+
  geom_hline(yintercept=0)
pSb2

####Export lots of causal impact as an array
Fig5Rev<-plot_grid(pNa,pNb,pCa, pCb, pSa, pSb, ncol=2, align="v", axis = "rl")
Fig5Rev

pdf("Updated/Cause1Updated.pdf", width=18, height=16.55)
Fig5Rev
dev.off()

FigCI2<-plot_grid(pNa2,pNb2,pCa2, pCb2, pSa2, pSb2, ncol=2, align="v", axis = "rl")
FigCI2

pdf("Updated/Cause2Updated.pdf", width=18, height=16.55)
FigCI2
dev.off()


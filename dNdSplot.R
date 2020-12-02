library(ggplot2)
library(cowplot)
library(plyr)
library(gridExtra)
library(gtable)
library(scales)
library(grid)
library(RColorBrewer)
library(tidyr)
library(egg)
library(reshape2)
library(gridGraphics)

### To produce Dos/DnDs plot ### Detail to create the input file are indicated in Dos_DnDS_Analysis.README  
xx=read.table("~/Paper/SupergeneDegeneration/Hmel2_NuSsBR.VariantSsIndels.annotate.Sorted.windowed.CountDnDs", header=TRUE, stringsAsFactors = FALSE)
P1=read.table("~/Paper/SupergeneDegeneration/P1_HomoP1.annotate.CountDnDs", header=TRUE, stringsAsFactors = FALSE)
P2=read.table("~/Paper/SupergeneDegeneration/P2_HomoP23_SNP.annotate.CountDnDs", header=TRUE, stringsAsFactors = FALSE)
P3=read.table("~/Paper/SupergeneDegeneration/P3_HomoP23_SNP.annotate.CountDnDs", header=TRUE, stringsAsFactors = FALSE)

xx=subset(xx, ((xx$Fin-xx$Deb)> 400000))
xx$type="Whole"
MeanP1=colSums(P1[,4:7])
MeanP1=c("NA", "NA", "NA", MeanP1, "P1")
MeanP2=colSums(P2[,4:7])
MeanP2=c("NA", "NA", "NA", MeanP2, "P2")
MeanP3=colSums(P3[,4:7])
MeanP3=c("NA", "NA", "NA", MeanP3, "P3")


xx=rbind(xx, MeanP1, MeanP2, MeanP3)
xx$FixMissense=as.numeric(xx$FixMissense)
xx$FixSynonymous=as.numeric(xx$FixSynonymous)
xx$HeteroMissense=as.numeric(xx$HeteroMissense)
xx$HeteroSynonymous=as.numeric(xx$HeteroSynonymous)
xx$dN.dS=xx$FixMissense/xx$FixSynonymous
xx$pN.pS=xx$HeteroMissense/xx$HeteroSynonymous
xx$MKTRatio=xx$pN.pS/xx$dN.dS
xx$DoS=xx$FixMissense/(xx$FixMissense+xx$FixSynonymous) - xx$HeteroMissense/(xx$HeteroMissense + xx$HeteroSynonymous)
xx$RatioFixed=(xx$HeteroMissense+xx$HeteroSynonymous)/(xx$HeteroMissense+xx$HeteroSynonymous+xx$FixMissense+xx$FixSynonymous)

my_thm = list(theme_bw(),
              theme(legend.position = "none",
                    axis.title.y = element_blank(),
                    axis.title.x = element_blank(),
                    axis.text.y=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks=element_blank(),
                    panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_blank(),
                    panel.border =  element_blank()
              )
)
my_thm2 = list(theme_bw(),
               theme(legend.position = "none",
                     axis.text = element_text(size = 10),
                     axis.title = element_text(size = 15),
                     panel.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     
                     axis.line = element_blank(),
                     panel.border = element_rect(colour = "black", fill=NA, size=0.8))
)


marg = theme(plot.margin=unit(rep(-0.5,4),"cm"))
marg2 = theme(plot.margin=unit(c(0,-1,0,0),"lines"))

#################Dos~pnps####################"
### RightPlot
probs <- c(0.025,0.05)
quantiles=quantile(xx$DoS, probs = probs)
densDoS <- density(xx$DoS)
xxDens=data.frame(x=densDoS$x, y=densDoS$y)
xxDens$quant<- factor(findInterval(xxDens$x,quantiles))
RightPlot=ggplot(xxDens, aes(x,y))+geom_line() +
  geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) +
  scale_x_continuous(breaks=quantiles) +
  scale_fill_brewer(guide="none")+my_thm+
  geom_vline(xintercept = xx[xx$type=="P1",12], colour="red")+
  geom_vline(xintercept = xx[xx$type=="P2",12], colour="orange")+
  geom_vline(xintercept = xx[xx$type=="P3",12], colour="blue")+
  coord_flip()+
  theme(plot.margin=unit(c(0,0,0,-0.5),"cm"))
####
### Top Plot ###
probs <- c(0.95,0.975)
quantiles=quantile(xx$pN.pS, probs = probs)
densPnPs <- density(xx$pN.pS)
xxDensPnPs=data.frame(x=densPnPs$x, y=densPnPs$y)
xxDensPnPs$quant<- factor(findInterval(xxDensPnPs$x,quantiles))
TopPlot=ggplot(xxDensPnPs, aes(x,y))+geom_line() +
  geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) +
  scale_x_continuous(breaks=quantiles) +
  scale_fill_brewer(guide="none")+my_thm+
  geom_vline(xintercept = xx[xx$type=="P1",10], colour="red")+
  geom_vline(xintercept = xx[xx$type=="P2",10], colour="orange")+
  geom_vline(xintercept = xx[xx$type=="P3",10], colour="blue")+
  marg


### Main Plot
DotPlot=ggplot(xx)+geom_point(aes(x=xx$pN.pS, y=xx$DoS, color=xx$type), shape=17,size=3, alpha=0.6)+
  ylab("Direction of Selection")+
  xlab("pN/pS")+
  scale_colour_manual(values = c("red", "orange", "blue", "black"))+
  my_thm2+ylim(min(densDoS$x),max(densDoS$x))+xlim(min(densPnPs$x),max(densPnPs$x))+
  theme(plot.margin=unit(c(0,0,0,0),"cm"))
DotPlot
####

empty = ggplot() + geom_blank() + marg
plotChoices = plot_grid(TopPlot, empty, DotPlot, RightPlot, align="hv",
                         rel_widths=c(3,1.07), rel_heights=c(1,3))

plotChoices=plotChoices + draw_grob(legend, 2.8/3.3, 0.38, 0.3/3.3, 1)

plotChoices


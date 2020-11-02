DataChr=data.frame(FP0=double(),FP1=double(),FP2=double(), Parameter=character(), Chrom=character(), Value=double(), stringsAsFactors = F)

########Parameter #####
###Mate Choice### Probability of mating from Chouteau et al., 2017
AcceptMbicFsil=0.66 #Male Bicoloratus (Hn1), Female Silvana (Hn0)
AcceptMbicFbic=0.12 #Male Bicoloratus (Hn1), Female Bicoloratus (Hn1)
AcceptMbicFtar=0.83 #Male Bicoloratus (Hn1), Female tarapotensis (Hn123)
AcceptMsilFsil=0.42 #
AcceptMsilFbic=0.80
AcceptMsilFtar=0.91
AcceptMtarFsil=0.70
AcceptMtarFbic=0.80
AcceptMtarFtar=0.385

### Survival adult ### Estimated from data of Chouteau et al., 2016, but extended for a 1 month period.
SurvABic=0.84 #Bicoloratus
SurvASil=0.27 #Silvana
SurvATar=0.65 #Tarapotensis 

### Survival Larvae ### Estimate from this study
WLP0P0=0.79069767  #Hn0/Hn0
WLP0P1=0.75770925  #Hn0/Hn1
WLP0P2=0.70422535  #Hn0/Hn123
WLP1P1=0.07142857  #Hn1/Hn1
WLP1P2=0.8055555  #Hn1/Hn123
WLP2P2=0.31  #Hn123/Hn123

### Variation of chromososome frequencies
Pas=8 #Change the frequency by 0.8 each loop
   for (P1 in seq(Pas, 1000-2*Pas, by=Pas)) #Variation of Hn1 frequency
   {
     for (P2 in seq(Pas, 1000-P1-Pas, by=Pas)) #Variation of Hn123 frequency
     {
  
 P0=1000-P1-P2
 FP0=P0/1000
 FP1=P1/1000
 FP2=P2/1000

### Freq Geno # Frequency of genotype considering hardy-weinberg
FP0P0=FP0*FP0
FP1P1=FP1*FP1
FP2P2=FP2*FP2
FP0P1=2*FP1*FP0
FP0P2=2*FP0*FP2
FP1P2=2*FP1*FP2

# Freq Pheno #Frequency of phenotype
FSil=FP0P0 #Silvana
FTar=FP2P2+FP0P2 #Tarapotensis
FBic=FP1P1+FP1P2+FP0P1 #Bicoloratus

#Relative freq of chrom in geno #Relative frequency of each haplotype (Hn1/Hn0/Hn123) in genotypes 
FP0inP0P0=FP0P0/(FP0P0+FP0P1+FP0P2) 
FP0inP0P1=FP0P1/(FP0P0+FP0P1+FP0P2)
FP0inP0P2=FP0P2/(FP0P0+FP0P1+FP0P2)
  
FP1inP1P1=FP1P1/(FP1P1+FP0P1+FP1P2)
FP1inP0P1=FP0P1/(FP1P1+FP0P1+FP1P2)
FP1inP1P2=FP1P2/(FP1P1+FP0P1+FP1P2)
  
FP2inP2P2=FP2P2/(FP2P2+FP0P2+FP1P2)
FP2inP0P2=FP0P2/(FP2P2+FP0P2+FP1P2)
FP2inP1P2=FP1P2/(FP2P2+FP0P2+FP1P2)

### Mate choice Fitness #Considering 50 % of male, 50% of female
WmcSil=(0.5*FSil*AcceptMsilFsil)+(0.5*FBic*AcceptMsilFbic)+(0.5*FTar*AcceptMsilFtar)+
  (0.5*FSil*AcceptMsilFsil)+(0.5*FBic*AcceptMbicFsil)+(0.5*FTar*AcceptMtarFsil) #Mean mate choice fitness for silvana individuals
WmcBic=(0.5*FSil*AcceptMbicFsil)+(0.5*FBic*AcceptMbicFbic)+(0.5*FTar*AcceptMbicFtar)+
  (0.5*FSil*AcceptMsilFbic)+(0.5*FBic*AcceptMbicFbic)+(0.5*FTar*AcceptMtarFbic) #Mean mate choice fitness for bicoloratus individuals
WmcTar=(0.5*FSil*AcceptMtarFsil)+(0.5*FBic*AcceptMtarFbic)+(0.5*FTar*AcceptMtarFtar)+
  (0.5*FSil*AcceptMsilFtar)+(0.5*FBic*AcceptMbicFtar)+(0.5*FTar*AcceptMtarFtar) #Mean mate choice fitness for tarapotensis individuals

### Genotype adult Fitness  ## Genotype adult fitness only depend of phenotype fitness
WAdultP0P0=SurvASil
WAdultP0P1=SurvABic
WAdultP0P2=SurvATar
WAdultP1P1=SurvABic
WAdultP1P2=SurvABic
WAdultP2P2=SurvATar

### Genotype mate choice Fitness ### Only depend on phenotype
WmcP0P0=WmcSil
WmcP0P1=WmcBic
WmcP0P2=WmcTar
WmcP1P1=WmcBic
WmcP1P2=WmcBic
WmcP2P2=WmcTar

#Fitness mesure must be relative. Mean fitness calculated here
MeanWAdult=WAdultP0P0*FP0P0+WAdultP0P1*FP0P1+WAdultP0P2*FP0P2+WAdultP1P1*FP1P1+WAdultP1P2*FP1P2+WAdultP2P2*FP2P2 #Adult fitness
MeanWmc=WmcP0P0*FP0P0+WmcP0P1*FP0P1+WmcP0P2*FP0P2+WmcP1P1*FP1P1+WmcP1P2*FP1P2+WmcP2P2*FP2P2 #Mate choice fitness
MeanWL=WLP0P0*FP0P0+WLP0P1*FP0P1+WLP0P2*FP0P2+WLP1P1*FP1P1+WLP1P2*FP1P2+WLP2P2*FP2P2 #Larvae fitness

#Relative fitness #Adult
RWAdultP0P0=WAdultP0P0/MeanWAdult
RWAdultP0P1=WAdultP0P1/MeanWAdult
RWAdultP0P2=WAdultP0P2/MeanWAdult
RWAdultP1P1=WAdultP1P1/MeanWAdult
RWAdultP1P2=WAdultP1P2/MeanWAdult
RWAdultP2P2=WAdultP2P2/MeanWAdult

#Relative fitness #Mate choice
RWmcP0P0=WmcP0P0/MeanWmc
RWmcP0P1=WmcP0P1/MeanWmc
RWmcP0P2=WmcP0P2/MeanWmc
RWmcP1P1=WmcP1P1/MeanWmc
RWmcP1P2=WmcP1P2/MeanWmc
RWmcP2P2=WmcP2P2/MeanWmc

#Relative fitness #Larvae
RWLP0P0=WLP0P0/MeanWL
RWLP0P1=WLP0P1/MeanWL
RWLP0P2=WLP0P2/MeanWL
RWLP1P1=WLP1P1/MeanWL
RWLP1P2=WLP1P2/MeanWL
RWLP2P2=WLP2P2/MeanWL

### Selection apply on genotype, but we want the fitness of haplotype, so we want to consider the frequency of haplotype in genotype
## Adult fitness
WAdultP0=FP0inP0P0*WAdultP0P0+FP0inP0P1*WAdultP0P1+FP0inP0P2*WAdultP0P2
WAdultP1=FP1inP1P1*WAdultP1P1+FP1inP0P1*WAdultP0P1+FP1inP1P2*WAdultP1P2
WAdultP2=FP2inP2P2*WAdultP2P2+FP2inP0P2*WAdultP0P2+FP2inP1P2*WAdultP1P2

## mate choice fitness
WmcP0=FP0inP0P0*WmcP0P0+FP0inP0P1*WmcP0P1+FP0inP0P2*WmcP0P2
WmcP1=FP1inP1P1*WmcP1P1+FP1inP0P1*WmcP0P1+FP1inP1P2*WmcP1P2
WmcP2=FP2inP2P2*WmcP2P2+FP2inP0P2*WmcP0P2+FP2inP1P2*WmcP1P2

## Larvae fitness
WLP0=FP0inP0P0*WLP0P0+FP0inP0P1*WLP0P1+FP0inP0P2*WLP0P2
WLP1=FP1inP1P1*WLP1P1+FP1inP0P1*WLP0P1+FP1inP1P2*WLP1P2
WLP2=FP2inP2P2*WLP2P2+FP2inP0P2*WLP0P2+FP2inP1P2*WLP1P2

## Calculating mean fitness 
MeanWAdultChr=WAdultP0*FP0+WAdultP1*FP1+WAdultP2*FP2
MeanWmcChr=WmcP0*FP0+WmcP1*FP1+WmcP2*FP2
MeanWLChr=WLP0*FP0+WLP1*FP1+WLP2*FP2

## Relative fitness
RWAdultP0=WAdultP0/MeanWAdultChr
RWAdultP1=WAdultP1/MeanWAdultChr
RWAdultP2=WAdultP2/MeanWAdultChr

RWmcP0=WmcP0/MeanWmcChr
RWmcP1=WmcP1/MeanWmcChr
RWmcP2=WmcP2/MeanWmcChr

RWLP0=WLP0/MeanWLChr
RWLP1=WLP1/MeanWLChr
RWLP2=WLP2/MeanWLChr

DataChr[nrow(DataChr)+1,]= c(FP0, FP1, FP2, "Larval survival", "Hn0", RWLP0)
DataChr[nrow(DataChr)+1,]= c(FP0, FP1, FP2, "Larval survival", "Hn1", RWLP1)
DataChr[nrow(DataChr)+1,]= c(FP0, FP1, FP2, "Larval survival", "Hn2", RWLP2)
DataChr[nrow(DataChr)+1,]= c(FP0, FP1, FP2, "Mate Choice", "Hn0", RWmcP0)
DataChr[nrow(DataChr)+1,]= c(FP0, FP1, FP2, "Mate Choice", "Hn1", RWmcP1)
DataChr[nrow(DataChr)+1,]= c(FP0, FP1, FP2, "Mate Choice", "Hn2", RWmcP2)
DataChr[nrow(DataChr)+1,]= c(FP0, FP1, FP2, "Adult Survival", "Hn0", RWAdultP0)
DataChr[nrow(DataChr)+1,]= c(FP0, FP1, FP2, "Adult Survival", "Hn1", RWAdultP1)
DataChr[nrow(DataChr)+1,]= c(FP0, FP1, FP2, "Adult Survival", "Hn2", RWAdultP2)

 }
}


DataChr$Value=as.numeric(DataChr$Value)
DataChr$FP1=as.numeric(DataChr$FP1)
DataChr$FP2=as.numeric(DataChr$FP2)
DataChr$FP0=as.numeric(DataChr$FP0)
DataChr$Log10=log10(DataChr$Value)


###### Plot with standard distribution ###
library(cowplot)
library(plyr)

col=c("#E5446D","#31BD97","#F19143")
col=c("grey5","#F19143","#E5446D")
BackTop=geom_rect(xmin=-Inf, xmax=Inf, ymin=0, ymax=Inf, fill="grey98")
BackDown= geom_rect(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0.0, fill="grey59")
TextUp=annotate("text",label="Frequency +",
                x=0.5, y=0.4,
                family="DejaVu Sans",
                colour="grey58", fontface="bold")
TextDown=annotate("text",label="Frequency -",
                  x=0.5, y=-0.13,
                  family="DejaVu Sans",
                  colour="grey98", fontface="bold")

DataChrMC=subset(DataChr, DataChr$Parameter=="Mate Choice")

SumDataMateChoiceFP1 = ddply(DataChrMC[DataChrMC$Chrom=="Hn1",], .(FP1), summarize, 
                             Mean=mean(Log10), Max=max(Log10), Min=min(Log10),
                             quant05=quantile(Log10, 0.1, na.rm=T),
                             quant95=quantile(Log10, 0.9, na.rm=T), sd=sd(Log10), Chrom="Hn1")
colnames(SumDataMateChoiceFP1)[1]="frequency"

SumDataMateChoiceFP0 = ddply(DataChrMC[DataChrMC$Chrom=="Hn0",], .(FP0), summarize, 
                             Mean=mean(Log10), Max=max(Log10), Min=min(Log10),
                             quant05=quantile(Log10, 0.1, na.rm=T),
                             quant95=quantile(Log10, 0.9, na.rm=T), sd=sd(Log10), Chrom="Hn0")
colnames(SumDataMateChoiceFP0)[1]="frequency"

SumDataMateChoiceFP2 = ddply(DataChrMC[DataChrMC$Chrom=="Hn2",], .(FP2), summarize,
                             Mean=mean(Log10), Max=max(Log10), Min=min(Log10),
                             quant05=quantile(Log10, 0.1, na.rm=T), 
                             quant95=quantile(Log10, 0.9, na.rm=T), sd=sd(Log10), Chrom="Hn2")
colnames(SumDataMateChoiceFP2)[1]="frequency"

DataBind=rbind(SumDataMateChoiceFP1,SumDataMateChoiceFP0,SumDataMateChoiceFP2)

Format=theme(axis.title=element_text(face="bold"))

base=ggplot(DataBind,aes(x=frequency, y=Mean))
PlotMC=base+
  BackTop+BackDown+TextUp+TextDown+
  geom_line(aes(color=Chrom), size=2)+ 
  geom_ribbon(aes(ymin=quant05, ymax=quant95, fill=Chrom), alpha=0.2)+
  scale_color_manual(values=col)+
  scale_fill_manual(values=col)+
  ggtitle("Mating success")+
  ylim(-0.15,0.5)+ 
  Format+
  theme(text = element_text(family = "Helvetica"))+
  geom_hline(yintercept = 0, linetype="dashed")+labs(x="Frequency", y="")

DataChrA=subset(DataChr, DataChr$Parameter=="Adult Survival")
SumDataAFP1 = ddply(DataChrA[DataChrA$Chrom=="Hn1",], .(FP1), summarize,
                    quant05=quantile(Log10, 0.1, na.rm=T),
                    quant95=quantile(Log10, 0.9, na.rm=T),
                    Mean=mean(Log10), Max=max(Log10), Min=min(Log10), Chrom="Hn1")
colnames(SumDataAFP1)[1]="frequency"

SumDataAFP0 = ddply(DataChrA[DataChrA$Chrom=="Hn0",], .(FP0), summarize,
                    quant05=quantile(Log10, 0.1, na.rm=T),
                    quant95=quantile(Log10, 0.9, na.rm=T),
                    Mean=mean(Log10), Max=max(Log10), Min=min(Log10), Chrom="Hn0")
colnames(SumDataAFP0)[1]="frequency"

SumDataAFP2 = ddply(DataChrA[DataChrA$Chrom=="Hn2",], .(FP2), summarize, 
                    quant05=quantile(Log10, 0.1, na.rm=T),
                    quant95=quantile(Log10, 0.9, na.rm=T),
                    Mean=mean(Log10), Max=max(Log10), Min=min(Log10), Chrom="Hn2")
colnames(SumDataAFP2)[1]="frequency"

DataBindA=rbind(SumDataAFP1,SumDataAFP0,SumDataAFP2)

base=ggplot(DataBindA,aes(x=frequency, y=Mean))
PlotA=base+
  BackTop+BackDown+TextUp+TextDown+
  geom_line(aes(color=Chrom), size=2)+ 
  geom_ribbon(aes(ymin=quant05, ymax=quant95, fill=Chrom), alpha=0.2)+
  ggtitle("Adult Survival")+
  Format+
  scale_color_manual(values=col)+
  scale_fill_manual(values=col)+
  ylim(-0.15,0.5)+
  theme(text = element_text(family = "Helvetica"))+
  geom_hline(yintercept = 0, linetype="dashed")+labs(x="Frequency", y="")

DataChrL=subset(DataChr, DataChr$Parameter=="Larval survival")


SumDataLFP1 = ddply(DataChrLHn1[DataChrL$Chrom=="Hn1",], .(FP1), summarize, 
                    quant05=quantile(Log10, 0.1, na.rm=T),
                    quant95=quantile(Log10, 0.9, na.rm=T),
                    Mean=mean(Log10), Max=max(Log10), Min=min(Log10), Chrom="Hn1")
colnames(SumDataLFP1)[1]="frequency"

SumDataLFP0 = ddply(DataChrLHn0[DataChrL$Chrom=="Hn0",], .(FP0), summarize,
                    quant05=quantile(Log10, 0.1, na.rm=T),
                    quant95=quantile(Log10, 0.9, na.rm=T),
                    Mean=mean(Log10), Max=max(Log10), Min=min(Log10), Chrom="Hn0")
colnames(SumDataLFP0)[1]="frequency"

SumDataLFP2 = ddply(DataChrLHn2[DataChrL$Chrom=="Hn2",], .(FP2), summarize,
                    quant05=quantile(Log10, 0.1, na.rm=T),
                    quant95=quantile(Log10, 0.9, na.rm=T),
                    Mean=mean(Log10), Max=max(Log10), Min=min(Log10), Chrom="Hn2")
colnames(SumDataLFP2)[1]="frequency"
DataBindL=rbind(SumDataLFP1,SumDataLFP0,SumDataLFP2)

base=ggplot(DataBindL,aes(x=frequency, y=Mean))
PlotL=base+
  BackTop+BackDown+TextUp+TextDown+
  geom_line(aes(color=Chrom), size=2)+ 
  geom_ribbon(aes(ymin=quant05, ymax=quant95, fill=Chrom), alpha=0.2)+
  ggtitle("Larval Survival")+
  ylim(-0.15,0.5)+
  scale_color_manual(values=col)+
  scale_fill_manual(values=col)+
  Format+
  theme(text = element_text(family = "Helvetica"))+
  geom_hline(yintercept = 0, linetype="dashed")+labs(x="Frequency", y="Selection coefficient")


legend <- cowplot::get_legend(PlotL)


Plot=plot_grid(PlotL+theme(legend.position = "none"),
               PlotMC+theme(legend.position = "none"), 
               PlotA+theme(legend.position = "none"),
               legend, 
               rel_widths = c(1, 1,1,0.75),
               ncol = 4)
Plot

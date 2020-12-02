
library(ggplot2)
library(cowplot)
library(plyr)
library(RColorBrewer)


x=read.table("Date_P1homo_P2Homo.txt", stringsAsFactors = FALSE, header = TRUE)  #Data last common ancestor (TMRCA)
x$type="TMRCA Hn23"
y=read.table("Date_P1homoP2homo_P1homoSs.txt", stringsAsFactors = FALSE, header = TRUE) #Data first common ancestor
y$type="TMRCA Hn23-Hn0"
x=rbind(x,y)

for (i in 1:nrow(x))
{
  if (x$Scaffold[i]  == "Hmel215006_P1"){
    x$Pos[i]=x$Pos[i]+1139000
    x$Scaffold[i]="Hmel215006"
  }
  else if (x$Scaffold[i]  == "Hmel215006_P2"){
    x$Pos[i]=x$Pos[i]+1548000
    x$Scaffold[i]="Hmel215006"
  }
  else if  (x$Scaffold[i]  == "Hmel215006_BefP1"){
    x$Scaffold[i]="Hmel215006"
  }
  else if  (x$Scaffold[i]  == "Hmel215025_P2"){
    x$Scaffold[i]="Hmel215025"
  }
  else if  (x$Scaffold[i]  == "Hmel215025_AftP2"){
    x$Pos[i]=x$Pos[i]+139000
    x$Scaffold[i]="Hmel215025"
  }
}

chr=read.table("Hmel2_chromosomes.agp", fill=TRUE, stringsAsFactors = FALSE) # read the agp file, fill the empty case with blank

for (i in 1:nrow(x)) # for each line
{if (x$Scaffold[i] %in% chr$V6) # if the scaffold name exist in the agp file
{
  y <- which(chr$V6 == x$Scaffold[i]) # grep de line of the scaffold in the agp file
  x$Pos[i]=x$Pos[i]+chr[y,2] # change the position of the window in the wholefile, to move from its position in the scaffold to its position on the chromosome
}
}

x$Mean =x$Mean * 0.00384
x$X95inf=x$X95inf * 0.00384
x$X95Sup=x$X95Sup * 0.00384

for (i in 1:nrow(x))
{
  if (x$Pos[i]<1420000 | x$Pos[i] > 3130000)
  {
    x$Inv[i]="WholeGenome"
  }
  else if (x$Pos[i]>1420000 & x$Pos[i] < 1826000)
  {
    x$Inv[i]="P1"
  }
  else if (x$Pos[i]<2010000 & x$Pos[i] > 1826000)
  {
    x$Inv[i]="P2"
  }
  else if (x$Pos[i]>2010000 & x$Pos[i] < 3130000)
  {
    x$Inv[i]="P3"
  }
}
x$TypeCol=paste(x$type, x$Inv)

xSub=subset(x, x$Inv!="WholeGenome")
plot1=ggplot(xSub, aes(x=Inv, y=Mean)) +
  geom_point(aes(fill=as.factor(type), fill=as.factor(type)), size=1, alpha=0.6, shape=21,
             position = position_jitterdodge(jitter.width = 0.20, jitter.height = 0, dodge.width=0.6))+
  geom_boxplot(aes(fill = as.factor(type)),position = position_dodge(0.7), alpha=0.8, outlier.shape = NA)+
  scale_y_continuous(breaks=pretty(xSub$Mean,n=5),labels = scales::number_format(accuracy = 0.1,decimal.mark = '.'))+
  scale_fill_manual(values = c("orange", "chocolate4"))+
  ylab("TMRCA")+xlab("Region")+
  theme_cowplot(14)+theme(legend.position = "none")
plot1

save_plot("TMRCAHn1.svg", plot1, base_aspect_ratio = 0.85)

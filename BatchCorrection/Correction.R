#Correction between two batches

#packages
install.packages("ggsignif")
library("ggsignif")
library(readr)
library(ggplot2)
library(readxl)

source("~/Documents/DARPA_AIM/Code/Andy_y_axis.R")

#Read in file Control for Correction
C <- read_csv("Data/Correction.csv")

S<- subset(C, Titer == "S")
R<- subset(C, Titer == "RBD")

summary(lm((S$New) ~ (S$Old)+0))
summary(lm(R$New ~ R$Old+0))

pred.func.RBD = function(x){
  0.418355*x
}

pred.func.S = function(x){
  0.4736*x
}

new.dat <- data.frame(x = seq(from = 0, to = 10^6))
new.dat$RBD <- pred.func.RBD(new.dat$x)
new.dat$S <- pred.func.S(new.dat$x)

G<- ggplot()+ geom_line(data=new.dat, aes(x=x, y=RBD), lty=2, col=1)+
  geom_line(data=new.dat, aes(x=x, y=S), lty=2, col=2)+
  theme_bw()+ geom_point(data=C, aes(x=Old, y=New, col=Titer, shape=Titer, size=Titer))+ xlab("Old antibody titers")+ ylab("New antibody titers")+ 
  scale_y_continuous(limits=c(0, 2e5), trans="log1p",breaks = as.vector(c(0,1) %o% 10^(0:8)),
                     minor_breaks = as.vector((1:9) %o% 10^(0:8)), labels = ggplot_scientific_notation_axes_labels)+
  scale_x_continuous(limits=c(0, 2e5), trans="log1p",breaks = as.vector(c(0,1) %o% 10^(0:8)),
                     minor_breaks = as.vector((1:9) %o% 10^(0:8)), labels = ggplot_scientific_notation_axes_labels)+
  scale_size_manual(values=c(3,3))+scale_shape_manual(values=c(1,20))+scale_color_manual(values=c(1,2))

save_plot("Figures/New/Correction.png", G, nrow=1,ncol=1 ,base_aspect_ratio = 1.3)

###ANTIBODY CORRECTION
#Read in files 

SARS_COV2 <- read_csv("Data/SARS_Cov2_5.csv")

#Batch correction
SARS_COV2$S <- ifelse(SARS_COV2$batch==1, SARS_COV2$s_igg * 0.4736, SARS_COV2$s_igg)
SARS_COV2$RBD <- ifelse(SARS_COV2$batch==1, SARS_COV2$rbd_igg * 0.418355, SARS_COV2$rbd_igg)

G<-ggplot(SARS_COV2, aes(as.factor(batch), log1p(s_igg)))+geom_violin(fill=2, alpha=0.5)+theme_bw()+geom_signif(comparisons=list(c("1","2")),map_signif_level = TRUE)+
  ylab("log(S-specific antibody titer)")+ xlab("Batch")+scale_x_discrete(labels=c("1" = "Old", "2" = "New"))+theme(axis.text = element_text(size=13),axis.title=element_text(size=13))+
  scale_y_continuous(limits = c(0,12))

save_plot("Figures/New/S_NoCorrection.png", G, nrow=1,ncol=1 ,base_aspect_ratio = 1.3)

G<-ggplot(SARS_COV2, aes(as.factor(batch), log1p(S)))+geom_violin(fill=2, alpha=0.5)+theme_bw()+geom_signif(comparisons=list(c("1","2")),map_signif_level = TRUE)+
  ylab("log(S-specific antibody titer)")+ xlab("corrected Batch")+scale_x_discrete(labels=c("1" = "Old", "2" = "New"))+theme(axis.text = element_text(size=13),axis.title=element_text(size=13))+
  scale_y_continuous(limits = c(0,12))

save_plot("Figures/New/S_Correction.png", G, nrow=1,ncol=1 ,base_aspect_ratio = 1.3)


G<- ggplot(SARS_COV2, aes(as.factor(batch), log1p(rbd_igg)))+geom_violin(fill=1, alpha=0.5)+theme_bw()+geom_signif(comparisons=list(c("1","2")),map_signif_level = TRUE)+
  ylab("log(RBD-specific antibody titer)")+ xlab("Batch")+scale_x_discrete(labels=c("1" = "Old", "2" = "New"))+theme(axis.text = element_text(size=13),axis.title=element_text(size=13))+
  scale_y_continuous(limits = c(0,12))

save_plot("Figures/New/RBD_NoCorrection.png", G, nrow=1,ncol=1 ,base_aspect_ratio = 1.3)

G<- ggplot(SARS_COV2, aes(as.factor(batch), log1p(RBD)))+geom_violin(fill=1, alpha=0.5)+theme_bw()+geom_signif(comparisons=list(c("1","2")),map_signif_level = TRUE)+
  ylab("log(RBD-specific antibody titer)")+ xlab("corrected Batch")+scale_x_discrete(labels=c("1" = "Old", "2" = "New"))+theme(axis.text = element_text(size=13),axis.title=element_text(size=13))+
  scale_y_continuous(limits = c(0,12))

save_plot("Figures/New/RBD_Correction.png", G, nrow=1,ncol=1 ,base_aspect_ratio = 1.3)


### NEUTRALIZATION CORRECTION
C <- read_csv("Data/Neutr_Correction.csv")

W<- subset(C, Titer == "WA-1")
D<- subset(C, Titer == "Delta")
O<- subset(C, Titer == "Omicron")

summary(lm(W$New ~ W$Old+0))
summary(lm(D$New ~ D$Old+0))
summary(lm(O$New ~ O$Old+0))

pred.func.W = function(x){
  1.13395*x
}

pred.func.D = function(x){
  1.6886*x
}

pred.func.O = function(x){
  1.61622*x
}

new.dat <- data.frame(x = seq(from = 0, to = 10^5))
new.dat$W <- pred.func.W(new.dat$x)
new.dat$D <- pred.func.D(new.dat$x)
new.dat$O <- pred.func.O(new.dat$x)

G<- ggplot()+ geom_line(data=new.dat, aes(x=x, y=W), lty=2, col=4)+
  geom_line(data=new.dat, aes(x=x, y=D), lty=2, col=5)+
  geom_line(data=new.dat, aes(x=x, y=O), lty=2, col=3)+
  theme_bw()+ geom_point(data=C, aes(x=Old, y=New, col=Titer, shape=Titer, size=Titer))+ xlab("Old antibody titers (Batch2)")+ ylab("New antibody titers (Batch3)")+ 
  scale_y_continuous(limits=c(0, 1e4), trans="log1p", breaks = as.vector(c(0,1) %o% 10^(0:8)),
                     minor_breaks = as.vector((1:9) %o% 10^(0:8)), labels = ggplot_scientific_notation_axes_labels)+
  scale_x_continuous(limits=c(0, 1e4), trans="log1p", breaks = as.vector(c(0,1) %o% 10^(0:8)),
                     minor_breaks = as.vector((1:9) %o% 10^(0:8)), labels = ggplot_scientific_notation_axes_labels)+
  scale_size_manual(values=c(3,3,4))+scale_shape_manual(values=c(1,20,4))+scale_color_manual(values=c(5,3,4))

save_plot("Figures/New/Neutr_Correction.png", G, nrow=1,ncol=1 ,base_aspect_ratio = 1.3)

#Read in file
##Batch correction
Neutr<- read_excel("Data/Julia/AllNeutralizationData_PRNT50Values.xlsx")

names(Neutr)[3]<-"Batch"
Neutr[Neutr == "<40"] <- "0"
Neutr[Neutr == "<50"] <- "0"
Neutr[Neutr == ">6400"] <- "6400"

Neutr$`WA-1 Original`<- as.numeric(Neutr$`WA-1 Original`)
Neutr$`WA-1 Repeat`<- as.numeric(Neutr$`WA-1 Repeat`)

Neutr$WA_1 <-ifelse(Neutr$Batch==1, Neutr$`WA-1 Repeat`, Neutr$`WA-1 Original` * 1.13395)
Neutr$WA_1_old <-ifelse(Neutr$Batch==1, Neutr$`WA-1 Repeat`, Neutr$`WA-1 Original`)

WA_1 <- data.frame(Batch=Neutr$Batch, Titer=Neutr$WA_1_old)
WA_1<- rbind(WA_1, cbind(Batch=3,Titer=Neutr$`WA-1 Original`[Neutr$Batch == "1"]))


Neutr$`Delta Original`<- as.numeric(Neutr$`Delta Original`)
Neutr$`Delta Repeat`<- as.numeric(Neutr$`Delta Repeat`)

Neutr$Delta <-ifelse(Neutr$Batch==1, Neutr$`Delta Repeat`, Neutr$`Delta Original` * 1.6886)

Neutr$`Omicron Original`<- as.numeric(Neutr$`Omicron Original`)
Neutr$`Omicron Repeat`<- as.numeric(Neutr$`Omicron Repeat`)

Neutr$Omicron <-ifelse(Neutr$Batch==1, Neutr$`Omicron Repeat`, Neutr$`Omicron Original` * 1.61622)
#threshold
Neutr$WA_1[Neutr$WA_1 > 6400] <- 6400
Neutr$Delta[Neutr$Delta > 6400] <- 6400
Neutr$Omicron[Neutr$Omicron > 6400] <- 6400


G<-ggplot(Neutr, aes(as.factor(Batch), log1p(WA_1)))+geom_violin(fill=4, alpha=0.5)+theme_bw()+geom_signif(comparisons=list(c("1","2")),map_signif_level = TRUE)+
  ylab("log(WA-1 neutrlization titer)")+ xlab("corrected Batch")+scale_x_discrete(labels=c("1" = "New (3)", "2" = "Old (2)"))+theme(axis.text = element_text(size=13),axis.title=element_text(size=13))+
  scale_y_continuous(limits = c(0,12))+ geom_hline(yintercept=log1p(6400), lty=2, col="black")

save_plot("Figures/New/WA1_Correction.png", G, nrow=1,ncol=1 ,base_aspect_ratio = 1.3)

G<-ggplot(WA_1, aes(as.factor(Batch), log1p(Titer)))+geom_violin(fill=4, alpha=0.5)+theme_bw()+geom_signif(comparisons=list(c("1","2"), c("1","3"), c("2","3")), y_position = c(10, 11, 10),map_signif_level = TRUE)+
  ylab("log(WA-1 neutralization titer)")+ xlab("Batch")+scale_x_discrete(labels=c("1" = "New (3)", "2" = "Old (2)", "3"="Old (1)"))+theme(axis.text = element_text(size=13),axis.title=element_text(size=13))+
  scale_y_continuous(limits = c(0,12))+ geom_hline(yintercept=log1p(6400), lty=2, col="black")

save_plot("Figures/New/WA1_NoCorrection.png", G, nrow=1,ncol=1 ,base_aspect_ratio = 1.3)

##############################
###Extra Analysis

John<- data_frame(Batch= c(Neutr$Batch[1:45], Neutr$`Repeated Batch`[1:45]), Titer=c(Neutr$`WA-1 Original`[1:45], Neutr$`WA-1 Repeat`[1:45]))

A<- ggplot(John, aes(as.factor(Batch), log1p(Titer)))+geom_hline(yintercept=log1p(6400), lty=2, col="darkgrey")+geom_violin(fill=4, alpha=0.5)+theme_bw()+geom_signif(comparisons=list(c("1","3")), map_signif_level = TRUE)+
  ylab("log(WA-1 neutralization titer)")+ xlab("Batch")+scale_x_discrete(labels=c("1" = "PU", "3"="FFU"))+theme(axis.text = element_text(size=13),axis.title=element_text(size=13))+
  ggtitle("Plaque vs Focus forming assay (N=44)")+scale_y_continuous(limits = c(0,10))

B<- ggplot(Neutr[1:45,], aes(log1p(`WA-1 Original`), log1p(`WA-1 Repeat`)))+geom_point(alpha=0.5)+theme_bw()+theme(axis.text = element_text(size=13),axis.title=element_text(size=13))+
  xlab("log(WA-1 neutralization titer) PU")+ylab("log(WA-1 neutralization titer) FFU")+
  ggtitle("Plaque vs Focus forming assay (N=44)")


Fig<- plot_grid(A, B, align="hv", nrow=1, ncol=2)
save_plot("Figures/New/Assay_Comparison.png", Fig, nrow=1,ncol=2 ,base_aspect_ratio = 1.3)

##########################

#Merge files, update corrected titers
SARS_COV2$s_igg <- SARS_COV2$S
SARS_COV2$rbd_igg<- SARS_COV2$RBD

SARS_COV2<- SARS_COV2[,-c(144,145)]

ID<- c()
for (i in SARS_COV2$donor){
  j<- which(Neutr$Donor == i)
  ID<- cbind(ID,j)
}

SARS_COV2$wa1<- Neutr$WA_1[ID]
SARS_COV2$delta<- Neutr$Delta[ID]
SARS_COV2$omicron_ba2.12.1<- Neutr$Omicron[ID]

write.csv(SARS_COV2, "Data/Data_DARPA_1.csv",row.names = FALSE)

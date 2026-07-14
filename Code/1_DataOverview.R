###FIGURE 1

#Packages
source("Code/0_Packages.R")

#Load dataset
load("Data/Dataset_DARPA.RData")
Dataset_DARPA<- as.data.frame(Dataset_DARPA)

# Titer variables
outcome <- c("s_igg", "rbd_igg", "wa1", "delta", "omicron_ba2.12.1")

##Tissue variable
Tissue <- c("bld", "bm", "spl", "lng", "lln", "mln", "iln")

#Main cell types
Cell <- c("smbc", "scd8", "scd4", "stfh", "streg")
Predictors <- c()
for (i in Tissue){
  for(j in Cell){
    V<- paste(i, j, sep="_")
    Predictors <- c(Predictors, V)
  }
}

#Subsets
Subset_4 <- c("scd4", "scd4_naive", "scd4_tem", "scd4_tcm", "scd4_temra","scd4_cd49a", "scd4_cd103", "scd4_cxcr6", "stfh", "streg")
Predictors_4 <- c()
for (i in Tissue){
  for(j in Subset_4){
    V<- paste(i, j, sep="_")
    Predictors_4 <- c(Predictors_4, V)
  }
}

Subset_8 <- c("scd8", "scd8_naive", "scd8_tem", "scd8_tcm", "scd8_temra","scd8_cd49a", "scd8_cd103", "scd8_cxcr6")
Predictors_8 <- c()
for (i in c("spl", "lln")){
  for(j in Subset_8){
    V<- paste(i, j, sep="_")
    Predictors_8 <- c(Predictors_8, V)
  }
}

Subset_B <- c("smbc", "smbc_cd69", "smbc_igm", "smbc_igg", "smbc_iga")
Predictors_B <- c()
for (i in Tissue){
  for(j in Subset_B){
    V<- paste(i, j, sep="_")
    Predictors_B <- c(Predictors_B, V)
  }
}


### TISSUE AGE GENDER DISTRIBUTION - Figure1B
Data <- data_frame(donor=Dataset_DARPA$donor, sex=Dataset_DARPA$sex, age=Dataset_DARPA$age, Dataset_DARPA[,c(Predictors_4, Predictors_8, Predictors_B)])

for (i in Tissue){
 Data[i] <-  rowSums(dplyr::select(Data, starts_with(paste(i, "_", sep=""))), na.rm=T)
}

Data[Tissue]<- ifelse(Data[Tissue]==0, NA, 1)

set.seed(777)
Fig1B <- ggplot(Data)+ theme_bw()+ ylab("Tissue")+ geom_jitter(aes(x=age, y=bld+6, shape=sex), col="red", height=0.25, width=0, alpha=0.6)+
  geom_jitter(aes(x=age, y=bm+5, shape=sex), col="darkred",  height=0.25, width=0, alpha=0.6)+
  geom_jitter(aes(x=age, y=spl+4, shape=sex), col="purple",  height=0.25, width=0, alpha=0.6)+
  geom_jitter(aes(x=age, y=lng+3, shape=sex), col="blue",  height=0.25, width=0, alpha=0.6)+
  geom_jitter(aes(x=age, y=lln+2, shape=sex), col="deepskyblue", height=0.25, width=0, alpha=0.6)+
  geom_jitter(aes(x=age, y=mln+1, shape=sex), col="forestgreen", height=0.25, width=0, alpha=0.6)+
  geom_jitter(aes(x=age, y=iln, shape=sex), col="orange", height=0.25, width=0, alpha=0.6)+
  theme(legend.position = "right")+
  xlab("Age [years]")+ 
  scale_y_continuous(breaks = c(1,2,3,4,5,6,7), labels = c("ILN", "MLN", "LLN", "LNG", "SPL", "BM", "BLD"))+
  scale_x_continuous(breaks = c(20,25,30,35,40,45,50,55,60,65,70,75,80,85))+
  theme(axis.text = element_text(size=11),axis.title=element_text(size=11), axis.text.y = element_text(color=rev(c("red", "darkred","purple","blue","deepskyblue","forestgreen", "orange"))))+
  scale_shape_manual("sex", values=c(19,1), guide = guide_legend(override.aes =list(col = "black")))+theme(legend.position = "top")

#save_plot("Figures/Poster/Tissue.png", P1, nrow=1,ncol=2 ,base_aspect_ratio = 1.3)  


###TISSUE FREQUENCY DISTRIBUTION - Figure 1C

Freq <- Dataset_DARPA[ ,Predictors] %>% 
  pivot_longer(col= everything(),
               names_to = "Cells", 
               values_to = "values", values_drop_na = TRUE)

Freq$Tissue<-str_split_i(Freq$Cells, "_",1)
Freq$Cells <- factor(Freq$Cells, levels = rev(Predictors), ordered = TRUE)
#logit transformed
Freq$Trans_value <- logit(ifelse(Freq$values/100 <= 0 | Freq$values/100 == 0.001/100, 0.001/100, Freq$values/100))

#Colour panel
CC<- c("bld" = "red", "bm" = "darkred", "spl"= "purple","lng"="blue", "lln"="deepskyblue", "mln"="forestgreen", "iln"= "orange")

Fig1C<- ggplot(Freq, aes(x = Trans_value, y = Cells, fill=Tissue)) +theme_bw()+
  geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(width = 0.0, height = 0),
                      point_size = 0.75, point_alpha=0.5, alpha = 0.8,  bandwidth=0.3)+
  scale_fill_manual(values = CC)+
  xlab("% of/within Spike-specific cells (S+)") + ylab("Cell populations") +  scale_y_discrete(expand = expand_scale(add = c(0.5, 1.5)), label= rev(c(rep(c("mBc", "CD8", "CD4", "TFH", "Treg"),7))))+
  geom_vline(xintercept = logit(1e-3/100), lty=2)+theme(axis.text = element_text(size=11),axis.title=element_text(size=11), axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank(), 
                                                        legend.position="none",  axis.text.y = element_text(color=rev(c(rep("red",5), rep("darkred",5), rep("purple",5),rep("blue",5),rep("deepskyblue",5),rep("forestgreen",5),rep("orange" ,5)))))+
  scale_x_continuous(limits=c(logit(5e-4/100), logit(0.9)), breaks = logit(c(as.vector(c(1) %o% 10^(-3:1)),50, 90)/100),
                     minor_breaks = logit(c(as.vector(c(1,2.5,5,7.5) %o% 10^(-3:1)),90)/100), labels = c("Threshold","0.01","0.1", "1", "10", "50", "90"))



#Supplementary Figure Violin plots to compare the distributions with each other -Figure S1
Freq$Type <-str_split_i(Freq$Cells, "_",2)
Freq$Tissue <- factor(Freq$Tissue, levels = Tissue, ordered = TRUE)
Freq$Type <- factor(Freq$Type, levels = Cell, ordered = TRUE)

cell_populations <- unique(Freq$Type)
tissues_to_test <- setdiff(unique(Freq$Tissue), "bld")

stat_results <- expand.grid(type = cell_populations, tissue = tissues_to_test) %>%
  group_split(type, tissue) %>%
    map_df(~{
    curr_type   <- as.character(.x$type[[1]])
    curr_tissue <- as.character(.x$tissue[[1]])
    # Subset data for the test
    blood_vals <- Freq %>% filter(Type == curr_type, Tissue == "bld") %>% pull(Trans_value)
    tissue_vals <- Freq %>% filter(Type == curr_type, Tissue == curr_tissue) %>% pull(Trans_value)
    
    # Standard K-S test
    if(length(blood_vals) > 1 && length(tissue_vals) > 1) {
    ks_res <- ks.test(blood_vals, tissue_vals)
    
    data.frame(
      type = curr_type,
      group1 = "bld",
      group2 = curr_tissue,
      p = ks_res$p.value
    )
    } else {
      NULL # Skip if data is missing for a specific combination
    }
  }) %>%
  # Set the height for the p-value label (manually adjust if needed)
  group_by(type) %>%
  mutate(y.position = max(Freq$Trans_value, na.rm = TRUE) * 1.05)

stat_results_plot <- stat_results %>%
  ungroup() %>%                          # This fixes the "size 6 or 1" error
  rstatix::add_significance("p") %>%     # Creates 'p.signif'
  mutate(
    Type = type,  
    Tissue = group2,
    y.position = logit(0.9)              # Adjust height to be within your logit(0.9) limit
  )

FigS1A<- ggplot(Freq, aes(x=Type, y=Trans_value)) + theme_bw()+
  geom_violin(aes(fill=Tissue))+ 
  stat_pvalue_manual(
    stat_results_plot, 
    label = "p.signif", 
    x = "Type",              # Pass column name as string to avoid aesthetic conflicts
    inherit.aes = FALSE, 
    tip.length = 0.1,
    position = position_dodge(0.9)         # Aligns with the dodged violins
  )+
  scale_y_continuous(limits=c(logit(5e-4/100), logit(0.9)), breaks = logit(c(as.vector(c(1) %o% 10^(-3:1)),50, 90)/100),
                                    minor_breaks = logit(c(as.vector(c(1,2.5,5,7.5) %o% 10^(-3:1)),90)/100), labels = c("Threshold","0.01","0.1", "1", "10", "50", "90"))+
  scale_fill_manual(values = CC, labels = c("Blood", "Bone marrow", "Spleen" ,"Lung", "Lung draining lymph nodes", "Mesenteric lymph nodes", "Inguinal lymph nodes")) + theme(axis.text = element_text(size=11),axis.title=element_text(size=11))+
  geom_hline(yintercept = logit(1e-3/100), lty=2)+ ggtitle("Main cell populations") +
  geom_point(aes(fill=Tissue),position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),size=0.25, alpha=0.5)+
  ylab("% of/within Spike-specific cells (S+)") + xlab("Cell populations")+
  scale_x_discrete(label= c("mBc", "CD8", "CD4", "TFH", "Treg"))
  

#Subsets
Subset_4 <- c("scd4_naive", "scd4_tem", "scd4_tcm", "scd4_temra","scd4_cd49a", "scd4_cd103", "scd4_cxcr6")
Predictors_4 <- c()
for (i in Tissue){
  for(j in Subset_4){
    V<- paste(i, j, sep="_")
    Predictors_4 <- c(Predictors_4, V)
  }
}

Subset_8 <- c("scd8_naive", "scd8_tem", "scd8_tcm", "scd8_temra","scd8_cd49a", "scd8_cd103", "scd8_cxcr6")
Predictors_8 <- c()
for (i in c("spl", "lln")){
  for(j in Subset_8){
    V<- paste(i, j, sep="_")
    Predictors_8 <- c(Predictors_8, V)
  }
}

Subset_B <- c("smbc_cd69", "smbc_igm", "smbc_igg", "smbc_iga")
Predictors_B <- c()
for (i in Tissue){
  for(j in Subset_B){
    V<- paste(i, j, sep="_")
    Predictors_B <- c(Predictors_B, V)
  }
}


#Figure S1B
FreqB <- Dataset_DARPA[ ,Predictors_B] %>% 
  pivot_longer(col= everything(),
               names_to = "Cells", 
               values_to = "values", values_drop_na = TRUE)


FreqB$Tissue<-str_split_i(FreqB$Cells, "_",1)
FreqB$Cells <- factor(FreqB$Cells, levels = rev(Predictors_B), ordered = TRUE)
#logit transformed
FreqB$Trans_value <- logit(ifelse(FreqB$values/100 <= 0 | FreqB$values/100 == 0.0025/100, 0.001/100, FreqB$values/100))

FreqB$Type <- sapply(str_split(FreqB$Cells, "_"), function(x) paste(x[2], x[3], sep="_"))
FreqB$Tissue <- factor(FreqB$Tissue, levels = Tissue, ordered = TRUE)
FreqB$Type <- factor(FreqB$Type, levels = Subset_B, ordered = TRUE)

FigS1B<- ggplot(FreqB, aes(x=Type, y=Trans_value, fill=Tissue)) + theme_bw()+
  geom_violin()+ scale_y_continuous(limits=c(logit(5e-4/100), logit(0.9)), breaks = logit(c(as.vector(c(1) %o% 10^(-3:1)),50, 90)/100),
                                    minor_breaks = logit(c(as.vector(c(1,2.5,5,7.5) %o% 10^(-3:1)),90)/100), labels = c("Threshold","0.01","0.1", "1", "10", "50", "90"))+
  scale_fill_manual(values = CC) + theme(axis.text = element_text(size=11),axis.title=element_text(size=11),  legend.position="none")+
  geom_hline(yintercept = logit(1e-3/100), lty=2)+ggtitle("B-cell subsets") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),size=0.25, alpha=0.5)+
  ylab("% within Spike-specific cells (S+)") + xlab("Cell populations")+
  scale_x_discrete(label= c("CD69", "IgM", "IgG", "IgA"))


#Figure S1C
Freq8 <- Dataset_DARPA[ ,Predictors_8] %>% 
  pivot_longer(col= everything(),
               names_to = "Cells", 
               values_to = "values", values_drop_na = TRUE)


Freq8$Tissue<-str_split_i(Freq8$Cells, "_",1)
Freq8$Cells <- factor(Freq8$Cells, levels = rev(Predictors_8), ordered = TRUE)
#logit transformed
Freq8$Trans_value <- logit(ifelse(Freq8$values/100 <= 0 | Freq8$values/100 == 0.0025/100, 0.001/100, Freq8$values/100))

Freq8$Type <- sapply(str_split(Freq8$Cells, "_"), function(x) paste(x[2], x[3], sep="_"))
Freq8$Tissue <- factor(Freq8$Tissue, levels = Tissue, ordered = TRUE)
Freq8$Type <- factor(Freq8$Type, levels = Subset_8, ordered = TRUE)

FigS1C<- ggplot(Freq8, aes(x=Type, y=Trans_value, fill=Tissue)) + theme_bw()+
  geom_violin()+ scale_y_continuous(limits=c(logit(5e-4/100), logit(0.9)), breaks = logit(c(as.vector(c(1) %o% 10^(-3:1)),50, 90)/100),
                                    minor_breaks = logit(c(as.vector(c(1,2.5,5,7.5) %o% 10^(-3:1)),90)/100), labels = c("Threshold","0.01","0.1", "1", "10", "50", "90"))+
  scale_fill_manual(values = CC) + theme(axis.text = element_text(size=11),axis.title=element_text(size=11),  legend.position="none")+
  geom_hline(yintercept = logit(1e-3/100), lty=2)+ggtitle("CD8+ T-cell subsets") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),size=0.25, alpha=0.5)+
  ylab("% within Spike-specific cells (S+)") + xlab("Cell populations")+
  scale_x_discrete(label= c("Naive", "TEM", "TCM", "TEMRA", "CD49a", "CD103", "CXCR6"))


#Figure S1D
Freq4 <- Dataset_DARPA[ ,Predictors_4] %>% 
  pivot_longer(col= everything(),
               names_to = "Cells", 
               values_to = "values", values_drop_na = TRUE)


Freq4$Tissue<-str_split_i(Freq4$Cells, "_",1)
Freq4$Cells <- factor(Freq4$Cells, levels = rev(Predictors_4), ordered = TRUE)
#logit transformed
Freq4$Trans_value <- logit(ifelse(Freq4$values/100 <= 0 | Freq4$values/100 == 0.0025/100, 0.001/100, Freq4$values/100))

Freq4$Type <- sapply(str_split(Freq4$Cells, "_"), function(x) paste(x[2], x[3], sep="_"))
Freq4$Tissue <- factor(Freq4$Tissue, levels = Tissue, ordered = TRUE)
Freq4$Type <- factor(Freq4$Type, levels = Subset_4, ordered = TRUE)

FigS1D<- ggplot(Freq4, aes(x=Type, y=Trans_value, fill=Tissue)) + theme_bw()+
  geom_violin()+ scale_y_continuous(limits=c(logit(5e-4/100), logit(0.9)), breaks = logit(c(as.vector(c(1) %o% 10^(-3:1)),50, 90)/100),
                                    minor_breaks = logit(c(as.vector(c(1,2.5,5,7.5) %o% 10^(-3:1)),90)/100), labels = c("Threshold","0.01","0.1", "1", "10", "50", "90"))+
  scale_fill_manual(values = CC) + theme(axis.text = element_text(size=11),axis.title=element_text(size=11),  legend.position="none")+
  geom_hline(yintercept = logit(1e-3/100), lty=2)+ggtitle("CD4+ T-cell subsets") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9),size=0.25, alpha=0.5)+
  ylab("% within Spike-specific cells (S+)") + xlab("Cell populations")+
  scale_x_discrete(label= c("Naive", "TEM", "TCM", "TEMRA", "CD49a", "CD103", "CXCR6"))


### DISTRIBUTION FOR bAB - Figure 1D

Titers_bAb <- Dataset_DARPA[, c("s_igg", "rbd_igg")]%>% 
  pivot_longer(col= everything(),
               names_to = "Titers", 
               values_to = "values")

Titers_bAb$Titers <- factor(Titers_bAb$Titers, levels = outcome, ordered = TRUE)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

Fig1D<- ggplot(Titers_bAb, aes(y = values, x = Titers)) +theme_bw()+
  geom_vridgeline(aes(width = after_stat(density)), stat="ydensity", trim=FALSE, alpha = 0.85, scale = 2)+
  #scale_y_continuous(trans = "log1p")+
  scale_y_continuous(limits=c(0, 5e5), trans="log1p", labels = ggplot_scientific_notation_axes_labels ,breaks = 10^seq(1, 6, by = 1),
                     # as.vector(1 %o% 10^(1:6)), 
                     minor_breaks = as.vector((1:9) %o% 10^(1:6)))+
  scale_x_discrete(expand = expansion(add = c(0.5, 1)), guide = guide_axis(angle=45), label= addline_format(c("S-IgG", "RBD-IgG")))+
  geom_jitter(size=0.5, alpha=0.5, width = 0.0, height = 0)+ylab("AUC")+xlab("")+
  geom_hline(yintercept = 0, lty=2) + theme(axis.text = element_text(size=11),axis.title=element_text(size=11))


### DISTRIBUTION FOR nAB - Figure 1E

Titers_nAb <- Dataset_DARPA[,c("wa1", "delta", "omicron_ba2.12.1")] %>% 
  pivot_longer(col= everything(),
               names_to = "Titers", 
               values_to = "values")

Titers_nAb$Titers <- factor(Titers_nAb$Titers, levels = outcome, ordered = TRUE)

Fig1E<- ggplot(Titers_nAb, aes(y = values, x = Titers)) +theme_bw()+
  geom_vridgeline(aes(width = after_stat(density)), stat="ydensity", trim=FALSE, alpha = 0.85, scale = 2)+
  #scale_y_continuous(trans = "log1p")+  
  scale_x_discrete(expand = expansion(add = c(0.5, 1)), guide = guide_axis(angle=45), label= addline_format(c("WA1", "Delta", "Omicron")))+
  geom_jitter(size=0.5, alpha=0.5, width = 0.0, height = 0)+ylab("PRNT50 [FFU/mL]")+xlab("")+
  geom_hline(yintercept = 0, lty=2)+
  scale_y_continuous(limits=c(0.0, 5e5), trans="log1p" , labels = ggplot_scientific_notation_axes_labels ,breaks = 10^seq(1, 6, by = 1),
                     #                     as.vector(1 %o% 10^(1:6)), 
                     minor_breaks =as.vector((1:9) %o% 10^(1:6))) + theme(axis.text = element_text(size=11),axis.title=element_text(size=11))



###TITER CORRELATIONS - Figure 1F

TiterCorr_NA<- data.frame()
TiterCorrP_NA<- data.frame()

i<-1
j<-1
m<-1
n<-1
for (i in outcome){
  for (j in outcome){
    if(dim(Dataset_DARPA[!is.na(Dataset_DARPA[[i]]) & !is.na(Dataset_DARPA[[j]]), ])[1] >1){
      Spear_NA<- cor.test(as.vector(Dataset_DARPA[[i]]), as.vector(Dataset_DARPA[[j]]), method="spearman", exact = FALSE)
      TiterCorr_NA[m,n] <- Spear_NA$estimate
      TiterCorrP_NA[m,n] <- Spear_NA$p.value
    }else{     
      TiterCorr_NA[m,n] <- NA
      TiterCorrP_NA[m,n] <- NA
    }
    m <- m+1
  }
  n <- n+1
  m <- 1
}

#Benjamin-Hochberg Correction
k<- length(c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron"))
Titeradjust <- as.matrix(sapply(TiterCorrP_NA, as.numeric))
TiterCorrP_NA_BH <-  matrix(p.adjust(as.vector(Titeradjust), method = "BH"), nrow = k, ncol = k)


rownames(TiterCorr_NA) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
colnames(TiterCorr_NA)<- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
rownames(TiterCorrP_NA) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
colnames(TiterCorrP_NA)<- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
rownames(TiterCorrP_NA_BH) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
colnames(TiterCorrP_NA_BH)<- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")


CC <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(200))


#png("Manuscript/Figures/A2.png", width = 12 , height =12, units = "in", res = 300)

pdf("Manuscript/Figures/Fig1F.pdf", width = 8 , height =8)
#par(mar = c(5.1, 4.1, 4.1, 2.1))
par(mar = c(10, 2, 4, 6), ps=12) 

#par(mar = c(10, 2, 2, 5), pty="s") 
corrplot(as.matrix(TiterCorr_NA), p.mat = as.matrix(TiterCorrP_NA_BH), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=1, cl.cex=1, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         number.cex = 1, pch.cex=2, method = 'circle', col=CC, type="upper",  na.label.col = "black" , diag = FALSE, cl.pos = "b")
mtext("Spearman coefficient", side = 1, line = 9, cex = 1, xpd = TRUE)

dev.off()


### NUMBERS of NEUTRALIZATIONS - Figure 1G/H
Dataset_DARPA$number_neutral <- as.numeric(Dataset_DARPA$wa1>0) + as.numeric(Dataset_DARPA$delta>0) + as.numeric(Dataset_DARPA$omicron_ba2.12.1>0)

Dataset_DARPA$ID_neutral <- ifelse(Dataset_DARPA$number_neutral==0, 0, ifelse(Dataset_DARPA$number_neutral==1 | Dataset_DARPA$number_neutral==3, 2, 
                                                                      ifelse(Dataset_DARPA$delta>0, 3, 4)))

Dataset_DARPA <- Dataset_DARPA[-20,]

set.seed(10)
Fig1G <- ggplot(Dataset_DARPA, aes(x=as.factor(number_neutral), y=rbd_igg))+geom_boxplot()+geom_jitter(aes(shape=as.factor(ID_neutral)), alpha=0.5, width=0.25, size=2)+theme_bw()+xlab("Number of SARS-CoV-2 strain neutralization")+ylab("RBD-IgG AUC")+scale_y_log10()+
  theme(legend.position = "none")+theme(axis.text = element_text(size=11),axis.title=element_text(size=11))+
  scale_shape_manual(values=c(1,19, 15,17))


set.seed(10)
Fig1H <- ggplot(Dataset_DARPA, aes(x=as.factor(number_neutral), y=s_igg))+geom_boxplot()+geom_jitter(aes(shape=as.factor(ID_neutral)), alpha=0.5, width=0.25, size=2)+theme_bw()+xlab("Number of SARS-CoV-2 strain neutralization")+ylab("S-IgG AUC")+scale_y_log10()+
  theme(legend.position = "none")+theme(axis.text = element_text(size=11),axis.title=element_text(size=11))+
  scale_shape_manual(values=c(1,19, 15,17))


#Figure1
#Fig1A <- ggdraw() + draw_image("Figures/New/Sample.png", scale = 0.95)

#Fig1AB<- plot_grid(NULL, Fig1B, labels = c("A", "B"), nrow=2)

#Fig1DE <- plot_grid(Fig1D, Fig1E, nrow=1, ncol=2, labels = c("D", "E"), align="hv", rel_widths = c(3/7, 4/7))
#Fig1CDE <- aplot::plot_list(Fig1C, Fig1DE, ncol = 2, widths = c(2, 1), tag_levels = list(c("C", ""))) &
 # theme(
 #  plot.tag = element_text(
 #     size = 14,
 #     face = "bold"
 #   ),
 #   plot.tag.position = c(0, 0.99)
 # )


#Fig1FGH<- plot_grid(NULL, Fig1G, Fig1H, labels = c("F", "G", "H"), nrow=1, ncol=3, rel_widths = c(1/4,3/8,3/8))

#Fig1ABCDE <- plot_grid(Fig1AB,Fig1CDE, nrow=1, ncol=2, rel_widths = c(2/3, 4/3))

#Figure1<- plot_grid(Fig1ABCDE, Fig1FGH, nrow=2, ncol=1, rel_heights = c(2, 1)) 


#########Alternative arrangement

Fig1AB<- plot_grid(NULL, Fig1B, labels = c("A", "B"), nrow=1, rel_widths = c(0.6, 0.4))

Fig1DE <- plot_grid(Fig1D, Fig1E, nrow=1, ncol=2, labels = c("D", "E"), align="hv", rel_widths = c(3/7, 4/7))

Fig1CDE <- plot_grid(Fig1C, Fig1DE, labels = c("C", ""), nrow=1, ncol=2, rel_widths = c(2/3, 1/3))


#Fig1CDE <- aplot::plot_list(Fig1C, Fig1DE, ncol = 2, widths = c(2/3, 1/3), tag_levels = list(c("C", ""))) &
 # theme(
  # plot.tag = element_text(
  #    size = 14,
  #    face = "bold"
 # ),
  # plot.tag.position = c(0, 0.99)
  #)


Fig1FGH<- plot_grid(NULL, Fig1G, Fig1H, labels = c("F", "G", "H"), nrow=1, ncol=3, rel_widths = c(0.3,0.35,0.35))
                                                                                                

Figure1<- plot_grid(Fig1AB, Fig1CDE, Fig1FGH, nrow=3, ncol=1, rel_heights = c(1.2, 2.2, 1))

save_plot("Manuscript/Figures/Figure1_Final.pdf", Figure1, nrow=3,ncol=3 ,base_aspect_ratio = 1)


###FigureS1

FigS1BC<- plot_grid(FigS1B, FigS1C, labels = c("B", "C"), nrow=1, rel_widths = c(0.6, 0.4))

FigureS1<- plot_grid(FigS1A, FigS1BC, FigS1D, labels = c("A", " ", "D"), nrow=3, ncol=1)

save_plot("Manuscript/Figures/FigureS1_Final.pdf", FigureS1, nrow=3,ncol=4 ,base_aspect_ratio = 1)


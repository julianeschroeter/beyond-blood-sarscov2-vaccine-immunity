### Figure 2
#Packages
source("Code/0_Packages.R")

## load data
#Raw dataset
load("Data/Dataset_DARPA.RData")
data<- as.data.frame(Dataset_DARPA)

#Imputation
load("imp/imp.Rdata")
imp_data<- mi.res
rm(mi.res)

#Variables
outcome <-  c("s_igg", "rbd_igg", "wa1", "delta", "omicron_ba2.12.1")

Tissue <- c("bld", "bm", "spl", "lng", "lln", "mln", "iln")
Cell <- c("smbc", "scd8", "scd4", "stfh", "streg")

Predictors <- c()
for (i in Tissue){
  for(j in Cell){
    V<- paste(i, j, sep="_")
    Predictors <- c(Predictors, V)
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

Subset_8 <- c("scd8_naive", "scd8_tem", "scd8_tcm", "scd8_temra","scd8_cd49a", "scd8_cd103", "scd8_cxcr6")
Predictors_8 <- c()
for (i in c("spl", "lln")){
  for(j in Subset_8){
    V<- paste(i, j, sep="_")
    Predictors_8 <- c(Predictors_8, V)
  }
}

Subset_4 <- c("scd4_naive", "scd4_tem", "scd4_tcm", "scd4_temra","scd4_cd49a", "scd4_cd103", "scd4_cxcr6")
Predictors_4 <- c()
for (i in Tissue){
  for(j in Subset_4){
    V<- paste(i, j, sep="_")
    Predictors_4 <- c(Predictors_4, V)
  }
}

#Correlations before and after imputations
#Correlation between predictors and titers of imputed and raw dataset

#Before Imputation
TVal_out_NA<- data.frame()
TPVal_out_NA<- data.frame()
#after Imputation
TVal_out<- data.frame()
TPVal_out<- data.frame()

i<-1
j<-1
m<-1
n<-1
for (i in c(outcome, Predictors, Predictors_B, Predictors_8, Predictors_4)){
  for (j in c(outcome, Predictors, Predictors_B, Predictors_8, Predictors_4)){
    if(i!=j){
      Spear<- miceadds::micombine.cor(mi.res=imp_data, variables=c(i, j), method="spearman")
      TVal_out[m,n] <- Spear$r[1]
      TPVal_out[m,n] <- Spear$p[1]
      
      if(dim(data[!is.na(data[[i]]) & !is.na(data[[j]]), ])[1] >1){
        Spear_NA<- cor.test(as.vector(data[[i]]), as.vector(data[[j]]), method="spearman", exact = FALSE)
        TVal_out_NA[m,n] <- Spear_NA$estimate
        TPVal_out_NA[m,n] <- Spear_NA$p.value
      }else{     
        TVal_out_NA[m,n] <- NA
        TPVal_out_NA[m,n] <- NA
      }
    }else{
      TVal_out[m,n] <- 1
      TPVal_out[m,n] <- 2.2e-16
      
      TVal_out_NA[m,n] <- 1
      TPVal_out_NA[m,n] <-  2.2e-16
    }
    m <- m+1
  }
  n <- n+1
  m <- 1
}

rownames(TVal_out) <- c(outcome,  Predictors, Predictors_B, Predictors_8, Predictors_4)
colnames(TVal_out)<- c(outcome, Predictors, Predictors_B, Predictors_8, Predictors_4)
rownames(TPVal_out) <- c(outcome,  Predictors, Predictors_B, Predictors_8, Predictors_4)
colnames(TPVal_out)<- c(outcome ,Predictors, Predictors_B, Predictors_8, Predictors_4)

rownames(TVal_out_NA) <-  c(outcome, Predictors, Predictors_B, Predictors_8, Predictors_4)
colnames(TVal_out_NA)<- c(outcome, Predictors, Predictors_B, Predictors_8, Predictors_4)
rownames(TPVal_out_NA) <-  c(outcome, Predictors, Predictors_B, Predictors_8, Predictors_4)
colnames(TPVal_out_NA)<- c(outcome, Predictors, Predictors_B, Predictors_8, Predictors_4)

save(TVal_out, file = "Data/Correlation.RData")
save(TPVal_out, file = "Data/Correlation_P.RData")
save(TVal_out_NA, file = "Data/Correlation_NA.RData")
save(TPVal_out_NA, file = "Data/Correlation_P_NA.RData")


load("Data/Correlation.RData")
load("Data/Correlation_P.RData")
load("Data/Correlation_NA.RData")
load("Data/Correlation_P_NA.RData")

#Benjamin-Hochberg correction
k<- length(c(outcome, Predictors, Predictors_B, Predictors_8, Predictors_4))
TPVal_out_1 <- as.matrix(sapply(TPVal_out, as.numeric))
TPVal_out_BH <-  matrix(p.adjust(as.vector(TPVal_out_1), method = "BH"), nrow = k, ncol = k)
TPVal_out_NA_1 <- as.matrix(sapply(TPVal_out_NA, as.numeric))
TPVal_out_NA_BH <-  matrix(p.adjust(as.vector(TPVal_out_NA_1), method = "BH"), nrow = k, ncol = k)

rownames(TPVal_out_BH) <- c(outcome,  Predictors, Predictors_B, Predictors_8, Predictors_4)
colnames(TPVal_out_BH)<- c(outcome,Predictors, Predictors_B, Predictors_8, Predictors_4)

rownames(TPVal_out_NA_BH) <-  c(outcome, Predictors, Predictors_B, Predictors_8, Predictors_4)
colnames(TPVal_out_NA_BH)<- c(outcome, Predictors, Predictors_B, Predictors_8, Predictors_4)

##Figure 1F
CC <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(200))


sub_data <- TVal_out_NA[outcome, outcome] 
sub_data_P <-  TPVal_out_NA_BH[outcome, outcome]
sub_data_P_notcorrected <-  TPVal_out_NA[outcome, outcome]

colnames(sub_data) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
rownames(sub_data) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
colnames(sub_data_P) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
rownames(sub_data_P) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
colnames(sub_data_P_notcorrected) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
rownames(sub_data_P_notcorrected) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")


#png("Manuscript/Figures/A2.png", width = 12 , height =12, units = "in", res = 300)
pdf("Manuscript/Figures/Fig1F.pdf", width = 8 , height =8)
par(mar = c(10, 2, 4, 6), ps=12) 

#par(mar = c(10, 2, 2, 5), pty="s") 
corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=1, cl.cex=1, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         number.cex = 1, pch.cex=2, method = 'circle', col=CC, type="upper",  na.label.col = "black" , diag = FALSE, cl.pos = "b")
mtext("Spearman coefficient", side = 1, line = 9, cex = 1, xpd = TRUE)

dev.off()


pdf("Manuscript/Figures/Fig1F_notcorrected.pdf", width = 8 , height =8)
par(mar = c(10, 2, 4, 6), ps=12) 

#par(mar = c(10, 2, 2, 5), pty="s") 
corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P_notcorrected), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=1, cl.cex=1, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         number.cex = 1, pch.cex=2, method = 'circle', col=CC, type="upper",  na.label.col = "black" , diag = FALSE, cl.pos = "b")
mtext("Spearman coefficient", side = 1, line = 9, cex = 1, xpd = TRUE)

dev.off()




#Supplementary Figure S3A

sub_data <- TVal_out[outcome, outcome] 
sub_data_P <-  TPVal_out_BH[outcome, outcome]
sub_data_P_notcorrected <-  TPVal_out[outcome, outcome]

colnames(sub_data) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
rownames(sub_data) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
colnames(sub_data_P) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
rownames(sub_data_P) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
colnames(sub_data_P_notcorrected) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
rownames(sub_data_P_notcorrected) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")

pdf("Manuscript/Figures/FigS3A.pdf", width = 8 , height =8)
par(mar = c(10, 2, 4, 6), ps=12) 

#par(mar = c(10, 2, 2, 5), pty="s") 
corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=1, cl.cex=1, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         number.cex = 1, pch.cex=2, method = 'circle', col=CC, type="upper",  na.label.col = "black" , diag = FALSE, cl.pos = "b")
mtext("Spearman coefficient", side = 1, line = 9, cex = 1, xpd = TRUE)

dev.off()


pdf("Manuscript/Figures/FigS3A_notcorrected.pdf", width = 8 , height =8)
par(mar = c(10, 2, 4, 6), ps=12) 

#par(mar = c(10, 2, 2, 5), pty="s") 
corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P_notcorrected), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=1, cl.cex=1, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         number.cex = 1, pch.cex=2, method = 'circle', col=CC, type="upper",  na.label.col = "black" , diag = FALSE, cl.pos = "b")
mtext("Spearman coefficient", side = 1, line = 9, cex = 1, xpd = TRUE)

dev.off()


##Figure 2A
my_labels <- c(rep(c("mBc", "CD8", "CD4", "TFH", "Treg"),7))
my_custom_colours <- c(rep("red",5), rep("darkred",5), rep("purple",5),rep("blue",5),rep("deepskyblue",5),rep("forestgreen",5),rep("orange" ,5))

N<- as.vector(colSums(!is.na(data[,names(TVal_out_NA[Predictors])])))

sub_data <- TVal_out_NA[outcome, Predictors] 
sub_data_P <-  TPVal_out_NA_BH[outcome, Predictors]
sub_data_P_notcorrected <-  TPVal_out_NA[outcome, Predictors]

colnames(sub_data) <- rep("", ncol(sub_data))
rownames(sub_data) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")

pdf("Manuscript/Figures/Fig2A.pdf", width = 35, height = 9)
# Increase bottom margin to make space for sample size labels
par(mar = c(8, 2, 2, 5), xpd=TRUE, ps=12) 

corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "n")
n_col <- ncol(as.matrix(sub_data))
for(i in 1:n_col){
  text(x=i+0.25, y=nrow(as.matrix(sub_data))+1.25, labels=my_labels[i], srt=45, adj=1, xpd=TRUE, cex=2, col=my_custom_colours[i])
}
text(c(0:37), -0.5, c("n=",N, " ", " "), cex=2, xpd = TRUE)
mtext("Correlations of raw dataset", at=-1, side = 3, line = -2, cex = 3, font = 2, adj=0)

dev.off()


pdf("Manuscript/Figures/Fig2A_notcorrected.pdf", width = 35, height = 9)
# Increase bottom margin to make space for sample size labels
par(mar = c(8, 2, 2, 5), xpd=TRUE, ps=12) 

corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P_notcorrected), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "n")
n_col <- ncol(as.matrix(sub_data))
for(i in 1:n_col){
  text(x=i+0.25, y=nrow(as.matrix(sub_data))+1.25, labels=my_labels[i], srt=45, adj=1, xpd=TRUE, cex=2, col=my_custom_colours[i])
}
text(c(0:37), -0.5, c("n=",N, " ", " "), cex=2, xpd = TRUE)
mtext("Correlations of raw dataset", at=-1, side = 3, line = -2, cex = 3, font = 2, adj=0)

dev.off()


##Figure 2D

sub_data <- TVal_out[outcome, Predictors] 
sub_data_P <-  TPVal_out_BH[outcome, Predictors]
sub_data_P_notcorrected <-  TPVal_out[outcome, Predictors]

colnames(sub_data) <- rep("", ncol(sub_data))
rownames(sub_data) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")

pdf("Manuscript/Figures/Fig2D.pdf", width = 35, height = 9)
# Increase bottom margin to make space for sample size labels
par(mar = c(8, 2, 2, 5), xpd=TRUE, ps=12) 

corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "b")
n_col <- ncol(as.matrix(sub_data))
for(i in 1:n_col){
  text(x=i+0.25, y=nrow(as.matrix(sub_data))+1.25, labels=my_labels[i], srt=45, adj=1, xpd=TRUE, cex=2, col=my_custom_colours[i])
}
mtext("Pooled correlations of imputed datasets", at=-1, side = 3, line = -1, cex = 3, font = 2, adj=0)
text(18, -0.5, "Spearman coefficient", cex=2, xpd = TRUE)
#mtext("Spearman coefficient", side = 1, line = 9, cex = 2, xpd = TRUE)

dev.off()

pdf("Manuscript/Figures/Fig2D_notcorrected.pdf", width = 35, height = 9)
# Increase bottom margin to make space for sample size labels
par(mar = c(8, 2, 2, 5), xpd=TRUE, ps=12) 

corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P_notcorrected), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "b")
n_col <- ncol(as.matrix(sub_data))
for(i in 1:n_col){
  text(x=i+0.25, y=nrow(as.matrix(sub_data))+1.25, labels=my_labels[i], srt=45, adj=1, xpd=TRUE, cex=2, col=my_custom_colours[i])
}
mtext("Pooled correlations of imputed datasets", at=-1, side = 3, line = -1, cex = 3, font = 2, adj=0)
text(18, -0.5, "Spearman coefficient", cex=2, xpd = TRUE)
#mtext("Spearman coefficient", side = 1, line = 9, cex = 2, xpd = TRUE)

dev.off()


##Supplementary Figure S3B-D

###B cells
sub_data <- TVal_out[outcome, Predictors_B] 
sub_data_P <-  TPVal_out_BH[outcome, Predictors_B]
sub_data_P_notcorrected <-  TPVal_out[outcome, Predictors_B]

colnames(sub_data) <- rep("", ncol(sub_data))
rownames(sub_data) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")

my_labels <- c(rep(c("CD69", "IgM", "IgG", "IgA"),7))
my_custom_colours <- c(rep("red",4), rep("darkred",4), rep("purple",4),rep("blue",4),rep("deepskyblue",4),rep("forestgreen",4),rep("orange" ,4))

pdf("Manuscript/Figures/FigS3B.pdf", width = 25, height = 9)
# Increase bottom margin to make space for sample size labels
par(mar = c(8, 2, 2, 5), xpd=TRUE, ps=12) 

corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "n", title="B-cell subsets")
n_col <- ncol(as.matrix(sub_data))
for(i in 1:n_col){
  text(x=i-0.25, y=nrow(as.matrix(sub_data))+0.75, labels=my_labels[i], srt=45, adj=0, xpd=TRUE, cex=2, col=my_custom_colours[i])
}
dev.off()

pdf("Manuscript/Figures/FigS3B_notcorrected.pdf", width = 25, height = 9)
# Increase bottom margin to make space for sample size labels
par(mar = c(8, 2, 2, 5), xpd=TRUE, ps=12) 

corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P_notcorrected), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "n")
n_col <- ncol(as.matrix(sub_data))
for(i in 1:n_col){
  text(x=i-0.25, y=nrow(as.matrix(sub_data))+0.75, labels=my_labels[i], srt=45, adj=0, xpd=TRUE, cex=2, col=my_custom_colours[i])
}
dev.off()


###CD8 T cells
sub_data <- TVal_out[outcome, Predictors_8] 
sub_data_P <-  TPVal_out_BH[outcome, Predictors_8]
sub_data_P_notcorrected <-  TPVal_out[outcome, Predictors_8]

colnames(sub_data) <- rep("", ncol(sub_data))
rownames(sub_data) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")

my_labels <- c(rep(c("Naive", "TEM", "TCM", "TEMRA", "CD49a", "CD103", "CXCR6"),2))
my_custom_colours <- c(rep("purple",7),rep("deepskyblue",7))

pdf("Manuscript/Figures/FigS3C.pdf", width = 17.5, height = 9)
# Increase bottom margin to make space for sample size labels
par(mar = c(8, 2, 2, 5), xpd=TRUE, ps=12) 

corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "n", title="CD8 T-cell subsets")
n_col <- ncol(as.matrix(sub_data))
for(i in 1:n_col){
  text(x=i-0.25, y=nrow(as.matrix(sub_data))+0.75, labels=my_labels[i], srt=45, adj=0, xpd=TRUE, cex=2, col=my_custom_colours[i])
}
dev.off()


pdf("Manuscript/Figures/FigS3C_notcorrected.pdf", width = 17.5, height = 9)
# Increase bottom margin to make space for sample size labels
par(mar = c(8, 2, 2, 5), xpd=TRUE, ps=12) 

corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P_notcorrected), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "n")
n_col <- ncol(as.matrix(sub_data))
for(i in 1:n_col){
  text(x=i-0.25, y=nrow(as.matrix(sub_data))+0.75, labels=my_labels[i], srt=45, adj=0, xpd=TRUE, cex=2, col=my_custom_colours[i])
}
dev.off()


#CD4 T cells

sub_data <- TVal_out[outcome, Predictors_4] 
sub_data_P <-  TPVal_out_BH[outcome, Predictors_4]
sub_data_P_notcorrected <-  TPVal_out[outcome, Predictors_4]

colnames(sub_data) <- rep("", ncol(sub_data))
rownames(sub_data) <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")

my_labels <- c(rep(c("Naive", "TEM", "TCM", "TEMRA", "CD49a", "CD103", "CXCR6"),7))
my_custom_colours <- c(rep("red",7), rep("darkred",7), rep("purple",7),rep("blue",7),rep("deepskyblue",7),rep("forestgreen",7),rep("orange" ,7))

pdf("Manuscript/Figures/FigS3D.pdf", width = 41, height = 9)
# Increase bottom margin to make space for sample size labels
par(mar = c(8, 2, 2, 5), xpd=TRUE, ps=12) 

corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "n", title="CD4 T-cell subsets")
n_col <- ncol(as.matrix(sub_data))
for(i in 1:n_col){
  text(x=i-0.25, y=nrow(as.matrix(sub_data))+0.75, labels=my_labels[i], srt=45, adj=0, xpd=TRUE, cex=2, col=my_custom_colours[i])
}
dev.off()

pdf("Manuscript/Figures/FigS3D_notcorrected.pdf", width = 41, height = 9)
# Increase bottom margin to make space for sample size labels
par(mar = c(8, 2, 2, 5), xpd=TRUE, ps=12) 

corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P_notcorrected), is.corr = TRUE, col.lim = c(-1, 1), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "n")
n_col <- ncol(as.matrix(sub_data))
for(i in 1:n_col){
  text(x=i-0.25, y=nrow(as.matrix(sub_data))+0.75, labels=my_labels[i], srt=45, adj=0, xpd=TRUE, cex=2, col=my_custom_colours[i])
}
dev.off()

#Correlations of significance - Figure 2C

##BH corrected p-values in imputed datasets
idx <- which(TPVal_out_BH <= 0.05,  arr.ind = TRUE)

significant_pairs <- data.frame(
  row  = rownames(TVal_out)[idx[,1]],
  col  = colnames(TVal_out)[idx[,2]],
  r    = TVal_out[idx],
  p_BH = TPVal_out_BH[idx]
)

significant_pairs[which(significant_pairs$r != 1),]

##BH corrected p-values in raw dataset
idx <- which(TPVal_out_NA_BH <= 0.05 ,  arr.ind = TRUE)

significant_pairs_NA <- data.frame(
  row  = rownames(TVal_out_NA)[idx[,1]],
  col  = colnames(TVal_out_NA)[idx[,2]],
  r    = TVal_out_NA[idx],
  p_BH = TPVal_out_NA_BH[idx]
)

significant_pairs_NA[which(significant_pairs_NA$r != 1),]


####Correlation scatterplott before after imputation

ut <- upper.tri(TVal_out_NA, diag = FALSE)

df_corr <- data.frame(
  row = rownames(TVal_out_NA)[row(TVal_out_NA)[ut]],
  col = colnames(TVal_out_NA)[col(TVal_out_NA)[ut]],
  value_before = TVal_out_NA[ut],
  p_before= TPVal_out_NA_BH[ut],
  value_imp = TVal_out[ut],
  p_imp= TPVal_out_BH[ut]
)

#lower limit
df_corr$p_before[df_corr$p_before <= 2.2e-16] <- 2.2e-16
df_corr$name <- paste(df_corr$row, df_corr$col, sep = "\n")

LL<- subset(df_corr, (p_imp <= 0.05 ))

LL$name1<- c("S-IgG<br>RBD-IgG", "S-IgG<br>WA1", "RBD-IgG<br>WA1", "S-IgG<br>Delta", "RBD-IgG<br>Delta", "WA1<br>Delta",
             "S-IgG<br>Omicron", "RBD-IgG<br>Omicron", "WA1<br>Omicron", "Delta<br>Omicron",
             "Omicron<br><span style='color:deepskyblue'>CD8</span>","mBc<br><span style='color:deepskyblue'>IgM mBc</span>",
             "CD8 Naive<br>CD8 EM", "CD4 Naive<br>CD4 EM", "CD4 EM<br>CD4 CM", "CD4 EM<br>CD4 CM")

repelled <- ggplot_build(
  ggplot(LL, aes(p_before, p_imp, label = name1)) +
    geom_text_repel(seed = 144, max.overlaps = 100, label.padding=.1)
)$data[[1]]

LL$x_rep <- repelled$x
LL$y_rep <- repelled$y

LL$Tissue <- toupper(sub("_.*", "", LL$name))     
LL$Tissue[1:11] <- "Titer"

LL$adjust<- 0
LL$adjust[c(7,8,15)] <- 1
LL$adjust[c(4,14)] <- 2
LL$adjust[13] <- 3

LL$hjust <-0.5
LL$hjust[LL$adjust==0] <-0
LL$hjust[LL$adjust==3] <-1

LL$vjust <- 0.5
LL$vjust[LL$adjust==1] <-0
LL$vjust[LL$adjust==2] <-1

LL$hjust[8] <- 1
LL$hjust[14] <- 1
LL$hjust[15] <- 1


#Figure 2C
mat_colors <- c("BLD" = "red", "BM" = "darkred", "SPL" = "purple", 
                "LNG" = "blue", "LLN" = "deepskyblue", 
                "MLN" = "forestgreen", "ILN" = "orange", "Titer"="black")

CC <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))(200))


Fig2C <- ggplot(df_corr, aes(x = p_before, y = p_imp)) +
  geom_abline(slope=1, lty=3, col="black")+
  geom_point(aes(color=value_before, fill=value_imp), shape=21, alpha=0.8, size=3, stroke=1.5)+ theme_bw() +
  geom_hline(yintercept=0.05, lty=2, col="darkgray")+
  geom_vline(xintercept=0.05, lty=2, col="darkgray")+
  scale_color_gradientn(colours = CC, limits = c(-1, 1), name = "Spearman coefficient") +
  scale_fill_gradientn(colours = CC, limits = c(-1, 1), name = "Spearman coefficient") +
  # Allow a second color scale for labels
  ggnewscale::new_scale_color() +
  geom_richtext(data = LL, aes(x = x_rep, y = y_rep, label = name1, vjust=vjust, hjust=hjust, color=Tissue), size = 3,
                fontface = "bold", fill=NA, alpha = 0.6, label.size = NA,  nudge_x = 0.05)+
  scale_x_log10() + xlab("P-values of raw dataset")+
  scale_y_log10() + ylab("Pooled P-values of imputed datasets")+
  scale_color_manual(values=mat_colors, guide = "none")+
  theme(legend.position = "inside",  legend.position.inside = c(0.75, 0.2), legend.background = element_rect(fill = scales::alpha("white", 0.8)))


Fig2C_zoomed <- ggplot(df_corr, aes(x = p_before, y = p_imp)) +
  geom_abline(slope=1, lty=3, col="black")+
  geom_point(aes(color=value_before, fill=value_imp), shape=21, alpha=0.8, size=3, stroke=1.5)+ theme_bw() +
  geom_hline(yintercept=0.05, lty=2, col="darkgray")+
  geom_vline(xintercept=0.05, lty=2, col="darkgray")+
  scale_color_gradientn(colours = CC, limits = c(-1, 1), name = "Value") +
  scale_fill_gradientn(colours = CC, limits = c(-1, 1), name = "Value") +
  # Allow a second color scale for labels
  ggnewscale::new_scale_color() +
  geom_richtext(data = LL, aes(x = x_rep, y = y_rep, label = name1, vjust=vjust, hjust=hjust, color=Tissue), size = 3,
                fontface = "bold", fill=NA, alpha = 0.6, label.size = NA,  nudge_x = 0.05)+
  scale_color_manual(values=mat_colors)+
  #geom_magnify(from = log10(c(1e-6, 1.1, 1e-2, 1.2)), to = c(1e-20, 1e-15, 1e-10, 1e-2))+
  scale_x_log10() + xlab("P-values of raw dataset")+
  scale_y_log10() + ylab("Pooled P-values of imputed datasets")+
  coord_cartesian(xlim = c(1e-6, 1.1), ylim = c(1e-2, 1.2))+
  theme(legend.position="none")

  
#Alternative plot
# Non-significant points (alpha 0.3)
Fig2C_Alt<- ggplot(df_corr, aes(x = value_before, y = value_imp)) +
  geom_point(
    data = subset(df_corr, p_before > 0.05 & p_imp > 0.05),
    shape = 1, color = "grey", alpha = 0.3, show.legend = TRUE
  ) +
  # Significant points (alpha 1)
  geom_point(
    data = subset(df_corr, !(p_before > 0.05 & p_imp > 0.05)),
    aes(
      #shape = factor(ifelse(p_before > 0.05, "sig_imp", "sig_before")),
      color = factor(ifelse(p_imp <= 0.05, "sig_imp", "sig_before"))
    ),
    alpha = 0.8, show.legend = FALSE
  ) +
  scale_shape_manual(values = c("sig_before" = 19, "sig_imp" = 19)) +
  scale_color_manual(values = c("sig_before" = "black", "sig_imp" = "red")) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  scale_x_continuous(limits = c(-1, 1)) + xlab("Spearman coefficients in raw dataset")+
  scale_y_continuous(limits = c(-1, 1)) + ylab("Pooled Spearman coefficients in imputed datasets")


#Overview of imputation quality for main subsets - Figure 2B

#CREATING VARIANCE RATIO HEATMAP of imputed data (only measurements, no categorical covariants)
Ratio_matrix <- matrix(NA, length(data$donor), length(c(Predictors, Predictors_B, Predictors_8, Predictors_4)))
Ratio_matrix <- as.data.frame(Ratio_matrix)
mc<- 1

S<- c(Predictors, Predictors_B, Predictors_8, Predictors_4)

for (i in S){
  #c<- names(imp_data$imp)[i]
  REF_sd <- var(unlist(data[i]), na.rm=TRUE)
  imputed <- as.data.frame(imp_data$imp[i])
  
  names(Ratio_matrix)[mc] <- i
  
  for (j in 1:dim(imputed)[1]){
    r<- row.names(imputed)[j]  
    Cell_sd <- var(as.numeric(imputed[r,1:100]))
    Ratio_matrix[r,i] <- Cell_sd/REF_sd
  }
  
  mc<- mc+1
}

Ratio_matrix_1 <- as.matrix(Ratio_matrix)
Ratio_matrix_1[Ratio_matrix_1 > 1] <- 1 #necessary to get above 1 into the scale  could also be set to 2

colfunc <- colorRampPalette(c("white", "grey75", "grey30"))(500)

#Identify censored data - 
#data into the right order
data_1 <- data[, S]

#Modify data to identify censored and missing data
data_1[, Predictors] <- lapply(data_1[, Predictors], function(x) {
  x[x == 0.001] <- -1
  x
})

data_1[!(is.na(data_1) | data_1 == -1)] <- 0
data_1[data_1 == -1] <- 1

#reorder rows and cols
row_order <- order(rowSums(is.na(data_1[,Predictors])))

# Data frame with column annotations.
mat_col <- data.frame(Tissue = c(rep("BLD",5), rep("BM",5), rep("SPL",5),rep("LNG",5),rep("LLN",5),rep("MLN",5),rep("ILN" ,5)))
rownames(mat_col) <- colnames(data_1[, Predictors])

# List with colors for each annotation.
mat_colors <- list(Tissue =  c("BLD" = "red", "BM" = "darkred", "SPL" = "purple", "LNG" = "blue",  "LLN" = "deepskyblue", "MLN" = "forestgreen", "ILN"="orange"))


#Dimension and order of col and row has to be fixed!!!!
#data_1 comes form 0_Heatmap_Completness.R
Ratio_matrix_1[which(is.na(Ratio_matrix_1))] <- 0
display_mat <-  matrix("", nrow = nrow(Ratio_matrix_1), ncol = ncol(Ratio_matrix_1))
display_mat[which(data_1 == 1)] <- "x" 

colnames(display_mat) <- colnames(Ratio_matrix_1)

P<- pheatmap::pheatmap(as.matrix(Ratio_matrix_1[row_order,Predictors]), 
                       scale = "none", 
                       color = colfunc, 
                       display_numbers = display_mat[row_order,Predictors],
                       cluster_rows = FALSE, 
                       cluster_cols = FALSE, 
                       show_rownames = FALSE, 
                       show_colnames = TRUE, na_col="white",
                       number_color = "black", 
                       legend = FALSE, 
                       labels_col = c(rep(c("mBc", "CD8", "CD4", "TFH", "Treg"),7)),
                       annotation_col = mat_col,
                       annotation_colors = mat_colors,
                       annotation_legend = FALSE,
                       annotation_names_col = FALSE,
                       angle_col = c("90"),
                       #treeheight_row = 0, 
                       #treeheight_col = 0,
                       border_color = "grey90"
)

P[[4]]$grobs[[2]]$gp$col <- c(rep("red",5), rep("darkred",5), rep("purple",5),rep("blue",5),rep("deepskyblue",5),rep("forestgreen",5),rep("orange" ,5))

gt <- P$gtable
layout <- gt$layout
anno_idx <- which(layout$name == "col_annotation")

annotation <- grid.text(names(mat_colors$Tissue), x=seq(0.075,0.925,0.85/6), y=unit(0.5, "npc"), just="center",gp = gpar(fontsize = 11))
gt <- gtable_add_grob(gt,  grobs = annotation,  t = layout[anno_idx, "t"],   l = layout[anno_idx, "l"],   b = layout[anno_idx, "b"],
                      r = layout[anno_idx, "r"],
                      name = "custom_annotation_text"
)


P <- as.ggplot(gt, scale=0.975) +theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Build color bar as list of grobs
legend_boxes <- list(rectGrob(x = unit(0.1, "npc"), y = unit(0.675, "npc"),
                              width = unit(0.025, "npc"), height = unit(0.2, "npc"),
                              gp = gpar(fill = "white", col = "black")))


legend_x <- seq(0.3, 0.9, length.out = length(colfunc))
legend_y <- -0.05  # slightly below bottom edge

# Build color bar as list of grobs
legend_rects <- lapply(seq_along(colfunc), function(i) {
  rectGrob(x = unit(legend_x[i], "npc"), y = unit(0.5, "npc"),
           width = unit(0.012, "npc"), height = unit(0.2, "npc"),
           gp = gpar(fill = colfunc[i], col = NA))
})

# Add labels (optional)
legend_labels <- list(
  textGrob("censored", x = unit(0.1, "npc"), y = unit(0.275, "npc"), gp = gpar(fontsize = 9)),
  textGrob("x", x = unit(0.1, "npc"), y = unit(0.7, "npc"), gp = gpar(fontsize = 9)),
  textGrob("0", x = unit(0.3, "npc"), y = unit(0.1, "npc"), gp = gpar(fontsize = 9)),
  textGrob("0.5", x = unit(0.6, "npc"), y = unit(0.1, "npc"), gp = gpar(fontsize = 9)),
  textGrob(">1", x = unit(0.9, "npc"), y = unit(0.1, "npc"), gp = gpar(fontsize = 9)),
  textGrob("Ratio of imputation variance and population variance", x = 0.6, y = unit(0.85, "npc"), gp = gpar(fontsize = 10, fontface = "bold"))
)

# Combine into one legend grob
legend_grob <- grobTree(gList(do.call(gList, legend_boxes), do.call(gList, legend_rects), do.call(gList, legend_labels)))


Fig2B <- plot_grid(
  P,
  ggdraw(legend_grob),
  ncol = 1,
  rel_heights = c(1, 0.1)  # Adjust legend height as needed
)

##Figure 2

Fig2A<- plot_grid(NULL, labels = c('A'), nrow=1)

Fig2BC <- plot_grid(Fig2B, Fig2C, ncol=2, labels = c('B', 'C'))

Fig2D<- plot_grid(NULL, labels = c('D'), nrow=1)

Figure2 <- plot_grid(Fig2A, Fig2BC, Fig2D, rel_heights = c(1,1.5,1), nrow=3)

save_plot("Manuscript/Figures/Figure2_Final.pdf", Figure2, nrow=3.5,ncol=2)
save_plot("Manuscript/Figures/Figure2C.pdf", Fig2C, base_width = 8, base_height = 8)
save_plot("Manuscript/Figures/Figure2C_zoomed.pdf", Fig2C_zoomed, nrow=1,ncol=1)

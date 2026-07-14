#DimensioReduction - Figure 3

#Packages
source("Code/0_Packages.R")

## load data
#Raw dataset
load("Data/Dataset_DARPA.RData")
data<- as.data.frame(Dataset_DARPA)

#Load imputation dataset
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

#####################################################################
####PCA

PCA_load <- list()
PCA_var <- list()
PCA_val <- list()

for (i in 1:100) {
  m.data<- mice::complete(imp_data, i)
  #Data log transformation
  m.data[,outcome] <- log1p(m.data[,outcome])
  m.data[,c(Predictors, Predictors_B, Predictors_8, Predictors_4)] <- lapply(m.data[,c(Predictors, Predictors_B, Predictors_8, Predictors_4)], function(x) ifelse(x == 0, log10(0.05), log10(x)))
  set.seed(123)
  PP<- prcomp(m.data[, c(outcome, Predictors, Predictors_B, Predictors_8, Predictors_4)], center = TRUE, scale. = TRUE)
  PCA_load[[i]] <- PP$rotation[,1:3]
  PCA_var[[i]] <- summary(PP)$importance[2, 1:3]
  PCA_val[[i]] <- PP$x[,1:3]
}



#Sign adjustment
align_signs <- function(ref, target) {
  aligned <- target
  for (i in 1:ncol(ref)) {
    if (cor(ref[, i], target[, i]) < 0) {
      aligned[, i] <- -target[, i]
    }
  }
  return(aligned)
}

# Apply align_signs
aligned <- list()
for (i in 1:100){
  aligned[[i]] <- align_signs( PCA_load[[1]][, 1:3], PCA_load[[i]][, 1:3])
}


pca_df_load<- data.frame(do.call(rbind, aligned))
pca_df_var<- data.frame(do.call(rbind, PCA_var))
#pca_df_val<- data.frame(do.call(rbind, pca_X_rot))


summary(pca_df_var)

pca_df_load$Names <- rep(c(outcome,Predictors, Predictors_B, Predictors_8, Predictors_4),100)

names(pca_df_load)<- c("PC1", "PC2", "PC3", "Names")


I<- which(pca_df_load$Names %in% c(outcome, Predictors))

my_custom_labels <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron", c(rep(c("mBc", "CD8", "CD4", "TFH", "Treg"),7)))
my_custom_colours <- c(rep("black",5), rep("red",5), rep("darkred",5), rep("purple",5),rep("blue",5),rep("deepskyblue",5),rep("forestgreen",5),rep("orange" ,5))
my_custom_headings <- c("Titers", "Blood", "Bone marrow", "Spleen", "Lung", "Lung lymph nodes", "Mesenteric lymph nodes", "Inguinal lymph nodes")


#Figure 3A
CC<- as.data.frame(do.call(rbind, tapply(pca_df_load$PC1[I]*(-1), factor(pca_df_load$Names[I], levels=unique(pca_df_load$Names[I])), summary)))
CC$colour <- ifelse (CC[,2] >0 & CC[,5] > 0,
                     ifelse( CC[,1] > 0 & CC[,6] > 0,  "red", "lightcoral"), ifelse( CC[,2] < 0  & CC[,5] < 0, ifelse(CC[,1] < 0 & CC[,6] < 0, "blue", "lightskyblue"),"grey70"))

Fig3A <- ggplot(pca_df_load[I, ], aes(x=factor(Names,levels=unique(Names)), y=PC1*(-1)))+ 
  scale_y_continuous(limits=c(-0.32,0.32), breaks=c(-0.3,-0.2, -0.1,0, 0.1, 0.2, 0.3))+
  geom_boxplot(fill=c(CC$colour), col=c(CC$colour), alpha=0.5)+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 90), labels = my_custom_labels)+
  xlab("")+ylab("Factor loading")+ggtitle("PC1 Loading Scores for Main cell populations (~7.12%)")+
  theme(axis.text = element_text(size=9),axis.title=element_text(size=9), axis.text.x = element_text(color = my_custom_colours))+
  geom_hline(yintercept = 0, lty=2)


#Figure 3B

CC<- as.data.frame(do.call(rbind, tapply(pca_df_load$PC2[I]*(1), factor(pca_df_load$Names[I], levels=unique(pca_df_load$Names[I])), summary)))
CC$colour <- ifelse (CC[,2] >0 & CC[,5] > 0,
                     ifelse( CC[,1] > 0 & CC[,6] > 0, "red", "lightcoral"), ifelse( CC[,2] < 0  & CC[,5] < 0, ifelse(CC[,1] < 0 & CC[,6] < 0, "blue", "lightskyblue"),"grey70"))


Fig3B <- ggplot(pca_df_load[I, ], aes(x=factor(Names,levels=unique(Names)), y=PC2*(1)))+ 
  #geom_jitter(width=0.1)+
  scale_y_continuous(limits=c(-0.32,0.32), breaks=c(-0.3,-0.2, -0.1,0, 0.1, 0.2, 0.3))+
  geom_boxplot(fill=c(CC$colour), col=c(CC$colour), alpha=0.5)+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 90), labels = my_custom_labels)+
  xlab("")+ylab("Factor loading")+ggtitle("PC2 Loading Scores for Main cell populations (~5.24%)")+
  theme(axis.text = element_text(size=9),axis.title=element_text(size=9), axis.text.x = element_text(color = my_custom_colours))+
  geom_hline(yintercept = 0, lty=2)



#Figure 3C

I<- which(pca_df_load$Names %in% c(Predictors_B, Predictors_8, Predictors_4))


my_custom_labels <- c(rep(c("CD69", "IgM", "IgG", "IgA"),7), rep(c("Naive", "TEM", "TCM", "TEMRA", "CD49a", "CD103", "CXCR6"),9))
my_custom_colours <- c(rep("red",4), rep("darkred",4), rep("purple",4),rep("blue",4),rep("deepskyblue",4),rep("forestgreen",4),rep("orange" ,4),rep("purple",7),rep("deepskyblue",7) ,rep("red",7), rep("darkred",7), rep("purple",7),rep("blue",7),rep("deepskyblue",7),rep("forestgreen",7),rep("orange" ,7))

CC<- as.data.frame(do.call(rbind, tapply(pca_df_load$PC1[I]*(-1), factor(pca_df_load$Names[I], levels=unique(pca_df_load$Names[I])), summary)))
CC$colour <- ifelse (CC[,2] >0 & CC[,5] > 0,
                     ifelse( CC[,1] > 0 & CC[,6] > 0, "red", "lightcoral"), ifelse( CC[,2] < 0  & CC[,5] < 0, ifelse(CC[,1] < 0 & CC[,6] < 0, "blue", "lightskyblue"),"grey70"))

Fig3C <- ggplot(pca_df_load[I, ], aes(x=factor(Names,levels=unique(Names)), y=PC1*(-1)))+ 
  #geom_jitter(width=0.1)+
  scale_y_continuous(limits=c(-0.32,0.32), breaks=c(-0.3,-0.2, -0.1,0, 0.1, 0.2, 0.3))+
  geom_boxplot(fill=c(CC$colour), col=c(CC$colour), alpha=0.5)+
  theme_bw()+scale_x_discrete(guide = guide_axis(angle = 90), labels = my_custom_labels)+
  xlab("")+ylab("Factor loading")+ggtitle("PC1 Loading Scores for Subset cell populations (~7.12%)")+
  theme(axis.text = element_text(size=9),axis.title=element_text(size=9), axis.text.x = element_text(color = my_custom_colours))+
  geom_hline(yintercept = 0, lty=2)


#######################################
##CCA
########################################

cca_ress_X <- list()
cca_ress_Y <- list()
cca_cor <- list()
cca_X <- list()
cca_Y <- list()
loadings_X <- list()
loadings_Y <- list()
loadings_crossY <- list()


for (i in 1:100) {
  m.data<- mice::complete(imp_data, i)
  m.data[,outcome] <- log1p( m.data[,outcome])
  m.data[,c(Predictors, Predictors_B, Predictors_8, Predictors_4)] <- lapply(m.data[,c(Predictors, Predictors_B, Predictors_8, Predictors_4)], function(x) ifelse(x == 0, log10(0.05), log10(x)))
  
  X_imp <- m.data[,c(Predictors, Predictors_B, Predictors_8, Predictors_4)] 
  Y_imp <-  m.data[,outcome]
  
  X_imp <- scale(X_imp, center = TRUE, scale = TRUE)
  Y_imp <- scale(Y_imp, center = TRUE, scale = TRUE)
  
  # do regularized CCA on the imputed data
  set.seed(1)
  cca_res <- CCA::rcc(X_imp, Y_imp, 1.0, 1.0)
  #cca_res <- rcc(X_imp, Y_imp, method = 'shrinkage')
  
  coords_X <- cca_res$scores$corr.X.xscores
  coords_Y <- cca_res$scores$corr.Y.xscores
  
  cca_ress_X[[i]] <- coords_X
  cca_ress_Y[[i]] <- coords_Y  
  cca_cor[[i]] <- cca_res$cor
  cca_X[[i]] <- cca_res$scores$xscores
  cca_Y[[i]] <- cca_res$scores$yscores
  
  loadings_X[[i]] <- cor(X_imp, cca_res$scores$xscores)
  loadings_Y[[i]] <- cor(Y_imp, cca_res$scores$yscores)
  loadings_crossY[[i]] <- cor(Y_imp, cca_res$scores$xscores)
}

# rotate or mirror the datasets such that they are maximally aligned

X_baseline <- cca_ress_X[[1]][, 1:3]

cca_ress_X_rot <- list()
cca_ress_Y_rot <- list()

for (i in 1:length(cca_ress_X)) {
  X_rot <- cca_ress_X[[i]][, 1:3]
  Y_rot <- cca_ress_Y[[i]][, 1:3]
  # find the rotation matrix that aligns the two datasets
  set.seed(10)
  rot_mat <- procrustes(X_baseline, X_rot)
  
  X <- data.frame(X_rot %*% rot_mat$rotation)
  Y <- data.frame(Y_rot %*% rot_mat$rotation)
  
  # TESTING!! not rotating should not work well
  # X <- data.frame(X_rot)
  # Y <- data.frame(Y_rot)
  
  # add the variable name as a column
  X$var <- as.factor(rownames(X))
  Y$var <- as.factor(rownames(Y))
  
  cca_ress_X_rot[[i]] <- X
  cca_ress_Y_rot[[i]] <- Y
}


# convert the list to a data frame
CCA_cor <- data.frame(do.call(rbind, cca_cor))
summary(CCA_cor)


cca_df_X <- data.frame(do.call(rbind, cca_ress_X_rot))
cca_df_Y <- data.frame(do.call(rbind, cca_ress_Y_rot))

#plot(cca_df_X[, 1], cca_df_X[, 2], xlim=c(-1,1), ylim=c(-1,1), cex=0.1)
#points(cca_df_Y[, 1], cca_df_Y[, 2], col="red", cex=0.1)

# compute the medians of the X and Y data per variable

loc_X <- aggregate(cbind(X1, X2, X3) ~ var, cca_df_X, median)
loc_Y <- aggregate(cbind(X1, X2, X3) ~ var, cca_df_Y, median)

# compute the length of the X and Y means

radius_X <- sqrt(loc_X$X1^2 + loc_X$X2^2)
radius_Y <- sqrt(loc_Y$X1^2 + loc_Y$X2^2)

radius_X1 <- sqrt(loc_X$X1^2 + loc_X$X3^2)
radius_Y1 <- sqrt(loc_Y$X1^2 + loc_Y$X3^2)

radius_X2 <- sqrt(loc_X$X2^2 + loc_X$X3^2)
radius_Y2 <- sqrt(loc_Y$X2^2 + loc_Y$X3^2)


# create data frame for ggplot2 plotting

ellipses_X <- confidence_ellipse(cca_df_X, x=1, y=2, .group_by=var, conf_level=0.5, robust=TRUE)
ellipses_Y <- confidence_ellipse(cca_df_Y, x=1, y=2, .group_by=var, conf_level=0.5, robust=TRUE)

ellipses_X1 <- confidence_ellipse(cca_df_X, x=1, y=3, .group_by=var, conf_level=0.5, robust=TRUE)
ellipses_Y1 <- confidence_ellipse(cca_df_Y, x=1, y=3, .group_by=var, conf_level=0.5, robust=TRUE)

ellipses_X2 <- confidence_ellipse(cca_df_X, x=2, y=3, .group_by=var, conf_level=0.5, robust=TRUE)
ellipses_Y2 <- confidence_ellipse(cca_df_Y, x=2, y=3, .group_by=var, conf_level=0.5, robust=TRUE)

# add radius to the data frame ellipses_X and ellipses_Y

ellipses_X$radius <- rep(0, nrow(ellipses_X))
ellipses_Y$radius <- rep(0, nrow(ellipses_Y))


ellipses_X1$radius <- rep(0, nrow(ellipses_X1))
ellipses_Y1$radius <- rep(0, nrow(ellipses_Y1))


ellipses_X2$radius <- rep(0, nrow(ellipses_X2))
ellipses_Y2$radius <- rep(0, nrow(ellipses_Y2))


for ( var in unique(ellipses_X$var) ) {
  idx <- which(ellipses_X$var == var)
  r <- radius_X[which(loc_X$var == var)]
  ellipses_X[idx,]$radius <- r
}

for ( var in unique(ellipses_Y$var) ) {
  idx <- which(ellipses_Y$var == var)
  r <- radius_Y[which(loc_Y$var == var)]
  ellipses_Y[idx,]$radius <- r
}

for ( var in unique(ellipses_X1$var) ) {
  idx <- which(ellipses_X1$var == var)
  r <- radius_X1[which(loc_X$var == var)]
  ellipses_X1[idx,]$radius <- r
}

for ( var in unique(ellipses_Y1$var) ) {
  idx <- which(ellipses_Y1$var == var)
  r <- radius_Y1[which(loc_Y$var == var)]
  ellipses_Y1[idx,]$radius <- r
}

for ( var in unique(ellipses_X2$var) ) {
  idx <- which(ellipses_X2$var == var)
  r <- radius_X2[which(loc_X$var == var)]
  ellipses_X2[idx,]$radius <- r
}

for ( var in unique(ellipses_Y2$var) ) {
  idx <- which(ellipses_Y2$var == var)
  r <- radius_Y2[which(loc_Y$var == var)]
  ellipses_Y2[idx,]$radius <- r
}


#Colourscheme
mat_colors <- c("BLD" = "red", "BM" = "darkred", "SPL" = "purple", 
                "LNG" = "blue", "LLN" = "deepskyblue", 
                "MLN" = "forestgreen", "ILN" = "orange", "Titer"="black")

cca_df_X$Tissue <- toupper(sub("_.*", "", cca_df_X$var))      # Extract part before "_"

loc_X$Tissue <- toupper(sub("_.*", "", loc_X$var)) 

loc_Y$Tissue <- "Titer"

ellipses_X$Tissue <- toupper(sub("_.*", "", ellipses_X$var))

ellipses_X1$Tissue <- toupper(sub("_.*", "", ellipses_X1$var))

ellipses_X2$Tissue <- toupper(sub("_.*", "", ellipses_X2$var))


#Figure 3D
# plot confidence ellipse for X

E_X<- which(ellipses_X$var %in% c(Predictors, Predictors_B, Predictors_4, Predictors_8))
#E_X<- which(ellipses_X$var %in% c(outcome, Predictors))
E_Y<- which(ellipses_Y$var %in% c(outcome))

# add confidence ellipse for Y

#geom_point(data=loc_Y, aes(x=X1, y=X2, shape=var), size=0.5, alpha=1, col="black")


# add a circle of radius 0.5 to the ggplot

circle1 <- data.frame(x=cos(seq(0, 2*pi, length.out=100)) * 0.5, y=sin(seq(0, 2*pi, length.out=100)) * 0.5)
circle2 <- data.frame(x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)))

# add variable names to the plot using the means of the data
# only select those variables with radius > 0.4

rmin <- 0.4

# add text for titers and cell freqs together with repel to avoid overlapping text

loc_data <- rbind(loc_X[loc_X$var %in% loc_X$var[radius_X > rmin], ], loc_Y[loc_Y$var %in% loc_Y$var[radius_Y > rmin], ])
loc_data$Tissue <- toupper(sub("_.*", "", loc_data$var))      # Extract part before "_"
loc_data$Tissue[(nrow(loc_data)-4):nrow(loc_data)] <- "Titer"

loc_data$vartype <- factor(ifelse(loc_data$var %in% outcome, "Y", "X"))

I_X<- which(cca_df_X$var %in% c(loc_data$var))
I_Y<- which(cca_df_Y$var %in% c(outcome))


#L_X<- which(loc_data$var %in% c(outcome, Predictors))
L_X<- which(loc_data$var %in% c(outcome, Predictors, Predictors_B, Predictors_8, Predictors_4))

LE_X<- which(ellipses_X$var %in% c(loc_data$var))

loc_data$var1 <- c("CD8", "CD69 mBc", "IgG mBc", "IgM mBc", "CD4 CD49a",  "mBc", "Delta", "Omicron", "RBD-IgG", "S-IgG", "WA1")

#########################Name colours accordingly to tissues

Fig3D<- ggplot(loc_X, aes(x=X1, y=X2, color=Tissue)) +
  geom_point(size=1, alpha=0.2, shape=4) +
  geom_path(data=ellipses_X[E_X, ], aes(x=x, y=y, group=var), alpha=0.1) +
  theme_minimal() +
  theme(legend.position="none") +
  xlab("CCA1") + ylab("CCA2") + 
  #ggtitle("CCA on imputed data (all celltypes)") +
  geom_path(data=ellipses_X[LE_X, ], aes(x=x, y=y, group=var), alpha=0.5) +
  geom_point(data=cca_df_X[I_X, ], aes(x=X1, y=X2), size=0.2, alpha=0.3)+
  geom_point(data=loc_data, aes(x=X1, y=X2, color=Tissue), size=1, alpha=0.8, shape=4)+
  xlim(-1, 1) + ylim(-1, 1) + scale_color_manual(values=mat_colors) + coord_fixed()+
  geom_path(data=ellipses_Y[E_Y, ], aes(x=x, y=y, group=var), col="black", alpha=0.5) +
  geom_point(data=cca_df_Y[I_Y, ], aes(x=X1, y=X2), size=0.2, alpha=0.3, col="black")+
  geom_point(data=loc_Y, aes(x=X1, y=X2), size=1, alpha=0.8, col="black", shape=4)+
  geom_path(data=circle1, aes(x=x, y=y), color="darkgrey", linewidth=0.5, lty=2) +
  geom_path(data=circle2, aes(x=x, y=y), color="darkgrey", size=0.5, lty=2)+
  geom_label_repel(data=loc_data[L_X, ], aes(x=X1, y=X2, label=var1, color=Tissue), 
                   size=3,  seed=144,   max.overlaps = 100, fontface = "bold",label.padding=.1, fill = "white", alpha=0.6, label.size = NA)




#################################################

#Figure 3E
# plot confidence ellipse for X

E_X<- which(ellipses_X2$var %in% c(Predictors, Predictors_B, Predictors_8, Predictors_4))
#E_X<- which(ellipses_X$var %in% c(outcome, Predictors))
E_Y<- which(ellipses_Y2$var %in% c(outcome))

# add confidence ellipse for Y

# add a circle of radius 0.5 to the ggplot

circle1 <- data.frame(x=cos(seq(0, 2*pi, length.out=100)) * 0.5, y=sin(seq(0, 2*pi, length.out=100)) * 0.5)
circle2 <- data.frame(x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)))

# add variable names to the plot using the means of the data
# only select those variables with radius > 0.4

rmin <- 0.4

# add text for titers and cell freqs together with repel to avoid overlapping text

loc_data <- rbind(loc_X[loc_X$var %in% loc_X$var[radius_X2 > rmin], ], loc_Y[loc_Y$var %in% loc_Y$var, ])
loc_data$Tissue <- toupper(sub("_.*", "", loc_data$var))      # Extract part before "_"
loc_data$Tissue[(nrow(loc_data)-4):nrow(loc_data)] <- "Titer"

loc_data$vartype <- factor(ifelse(loc_data$var %in% outcome, "Y", "X"))

I_X<- which(cca_df_X$var %in% c(loc_data$var))
I_Y<- which(cca_df_Y$var %in% c(outcome))

#L_X<- which(loc_data$var %in% c(outcome, Predictors))
L_X<- which(loc_data$var %in% c(outcome, Predictors, Predictors_B, Predictors_8, Predictors_4))

LE_X<- which(ellipses_X2$var %in% c(loc_data$var))

loc_data$var1 <- c("CD4", "CD8", "CD69 mBc", "Delta", "Omicron", "RBD-IgG", "S-IgG",  "WA1")

#loc_data$var <- c("CD4", "CD4 TEM", "TFH", "CD4","CD4 CD49a", "CD8", "CD8 TEM" ,"CD69 B-cells", "Delta", "Omicron", "S-IgG", "WA1")

#########################Name colours accordingly to tissues

Fig3E<- ggplot(loc_X, aes(x=X2, y=X3, color=Tissue)) +
  geom_point(size=1, alpha=0.2, shape=4) +
  geom_path(data=ellipses_X2[E_X, ], aes(x=x, y=y, group=var), alpha=0.1) +
  theme_minimal() +
  theme(legend.position="none") +
  xlab("CCA2") + ylab("CCA3") + 
  #ggtitle("CCA on imputed data (all celltypes)") + 
  xlim(-1, 1) + ylim(-1, 1) + scale_color_manual(values=mat_colors) + coord_fixed()+
  geom_path(data=ellipses_X2[LE_X, ], aes(x=x, y=y, group=var), alpha=0.5) +
  geom_point(data=cca_df_X[I_X, ], aes(x=X2, y=X3), size=0.2, alpha=0.3)+
  geom_point(data=loc_data, aes(x=X2, y=X3, color=Tissue), size=1, alpha=0.8, shape=4)+
  geom_path(data=ellipses_Y2[E_Y, ], aes(x=x, y=y, group=var), col="black", alpha=0.5) +
  geom_point(data=cca_df_Y[I_Y, ], aes(x=X2, y=X3), size=0.2, alpha=0.3, col="black")+
  geom_point(data=loc_Y, aes(x=X2, y=X3), size=1, alpha=0.8, col="black", shape=4)+
  geom_path(data=circle1, aes(x=x, y=y), color="darkgrey", linewidth=0.5, lty=2) +
  geom_path(data=circle2, aes(x=x, y=y), color="darkgrey", size=0.5, lty=2)+
  geom_label_repel(data=loc_data[L_X, ], aes(x=X2, y=X3, label=var1, color=Tissue), size=3,  seed=144,   max.overlaps = 100, fontface = "bold", label.padding=.1, fill = "white", alpha=0.6, label.size = NA)


# Figure 3
Fig3AB <- plot_grid(Fig3A, Fig3B,  labels = c('A', 'B'), nrow=1)

Fig3ABC <- plot_grid(Fig3AB, Fig3C,  labels = c('', 'C'), nrow=2)

Fig3DE<- plot_grid(Fig3D, Fig3E,  labels = c('D', 'E'), nrow=1)

Figure3 <-plot_grid(Fig3ABC, Fig3DE,  nrow = 2, rel_heights = c(2,1.5))

save_plot("Manuscript/Figures/Figure3_Final.pdf", Figure3, nrow=3, ncol=2)





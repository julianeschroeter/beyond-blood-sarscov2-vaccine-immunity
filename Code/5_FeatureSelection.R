#Feature Selection including covariates
###Figure 4

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

#Add variables
data$deathdate <- as.Date(paste(data$deathyear,data$deathmonth,1, sep='-'), format='%Y-%B-%d')
data$pandemic <- as.numeric(data$deathdate-as.Date('2019-12-01'))

long <- mice::complete(imp_data, action='long', include=TRUE)
# Generate new variable
long$deathdate <- as.Date(paste(long$deathyear,long$deathmonth,1, sep='-'), format='%Y-%B-%d')
long$pandemic <- as.numeric(long$deathdate-as.Date('2019-12-01'))
# Convert back to Mids
imp <- as.mids(long)

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


var_cat <- c("sex", "age", "infected", "vaccine_brand", "dose_number", "tpv", "pandemic")

#####################################################
###Correlations with covariates - Figure 4A
#####################################################

#Titers
Val_cat<- data.frame()
PVal_cat<- data.frame()

i<-1
j<-1
m<-1
n<-1
for (i in outcome){
  for (j in var_cat){
    expr<- expression(lm(as.formula(paste("log1p(", i, ") ~ ", j, sep = ""))))
    fit_full <- with(imp, expr)
    summary(mice::pool(fit_full))
    S<- D1(fit_full)
    Val_cat[m,n]<- S$result[1]
    PVal_cat[m,n]<- S$result[4]
    m <- m+1
  }
  n <- n+1
  m <-1
}

rownames(Val_cat) <- var_cat
colnames(Val_cat)<-c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")
rownames(PVal_cat) <- var_cat
colnames(PVal_cat)<-c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron")

#Tissue main
TVal_cat<- data.frame()
TPVal_cat<- data.frame()
i<-1
j<-1
m<-1
n<-1
for (i in Predictors){
  for (j in var_cat){
    expr<- expression(lm(as.formula(paste( "log10(ifelse(", i,  "==0, 0.05,",i, ")) ~ ", j, sep = ""))))
    fit_full <- with(imp, expr)
    summary(mice::pool(fit_full))
    S<- D1(fit_full)
    TVal_cat[m,n]<- S$result[1]
    TPVal_cat[m,n]<- S$result[4]
    
    m <- m+1
  }
  n <- n+1
  m <-1
}

rownames(TVal_cat) <- var_cat
colnames(TVal_cat)<-  Predictors
rownames(TPVal_cat) <- var_cat
colnames(TPVal_cat)<-  Predictors

Values_cat <- cbind(Val_cat, TVal_cat)
PValues_cat <- cbind(PVal_cat, TPVal_cat)

save(Values_cat, file = "Data/Fstatistics.RData")
save(PValues_cat, file = "Data/Fstatistics_P.RData")

load("Data/Fstatistics.RData")
load("Data/Fstatistics_P.RData")

#Figure 4A
my_labels <- c("S-IgG", "RBD-IgG", "WA1", "Delta", "Omicron", rep(c("mBc", "CD8", "CD4", "TFH", "Treg"),7))

colnames(Values_cat) <- rep("", ncol(Values_cat))
rownames(Values_cat) <- c("Sex", "Age", "Infected", "Vaccine brand", "Dose number", "Time post vaccination", "Time into pandemic")

my_custom_colours <- c(rep("black",5),rep("red",5), rep("darkred",5), rep("purple",5),rep("blue",5),rep("deepskyblue",5),rep("forestgreen",5),rep("orange" ,5))


CC <- colorRampPalette(brewer.pal(3, "YlGnBu"))(200)

pdf("Manuscript/Figures/Fig4A.pdf", width = 40, height = 13)
# Increase bottom margin to make space for sample size labels
par(mar = c(4, 2, 2, 5), xpd=TRUE, ps=12) 

corrplot(as.matrix(Values_cat), p.mat = as.matrix(PValues_cat), is.corr = FALSE, col.lim = c(0, 25), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "b")
n_col <- ncol(as.matrix(Values_cat))
for(i in 1:n_col){
  text(x=i-0.25, y=nrow(as.matrix(Values_cat))+1, labels=my_labels[i], srt=45, adj=c(0,1), xpd=TRUE, cex=2, col=my_custom_colours[i])
}
mtext("", at=-3, side = 3, line = -5, cex = 3, font = 2, adj=0)
text(20.5, -0.75, "Model Improvement (F-statistic)", cex=2, xpd = TRUE)
dev.off()



### Supplementary Figure S4

#Tissue subset
TVal_cat_sub<- data.frame()
TPVal_cat_sub<- data.frame()
i<-1
j<-1
m<-1
n<-1
for (i in c(Predictors_B, Predictors_8, Predictors_4)){
  for (j in var_cat){
    expr<- expression(lm(as.formula(paste( "log10(ifelse(", i,  "==0, 0.05,",i, ")) ~ ", j, sep = ""))))
    fit_full <- with(imp, expr)
    summary(mice::pool(fit_full))
    S<- D1(fit_full)
    TVal_cat_sub[m,n]<- S$result[1]
    TPVal_cat_sub[m,n]<- S$result[4]
    
    m <- m+1
  }
  n <- n+1
  m <-1
}

rownames(TVal_cat_sub) <- var_cat
colnames(TVal_cat_sub)<-  c(Predictors_B, Predictors_8, Predictors_4)
rownames(TPVal_cat_sub) <- var_cat
colnames(TPVal_cat_sub)<-  c(Predictors_B, Predictors_8, Predictors_4)

save(TVal_cat_sub, file = "Data/Fstatistics_sub.RData")
save(TPVal_cat_sub, file = "Data/Fstatistics_P_sub.RData")

load("Data/Fstatistics_sub.RData")
load("Data/Fstatistics_P_sub.RData")


#B cell
sub_data <- TVal_cat_sub[var_cat, Predictors_B] 
sub_data_P <-  TPVal_cat_sub[var_cat, Predictors_B]


colnames(sub_data) <- rep("", ncol(sub_data))
rownames(sub_data) <- c("sex", "age", "infected", "vaccine brand", "dose number", "time post vaccination", "time into pandemic")

my_labels <- c(rep(c("CD69", "IgM", "IgG", "IgA"),7))
my_custom_colours <- c(rep("red",4), rep("darkred",4), rep("purple",4),rep("blue",4),rep("deepskyblue",4),rep("forestgreen",4),rep("orange" ,4))

CC <- colorRampPalette(brewer.pal(3, "YlGnBu"))(200)

pdf("Manuscript/Figures/FigS4A.pdf", width = 25, height = 11)
# Increase bottom margin to make space for sample size labels
par(mar = c(4, 2, 2, 5), xpd=TRUE, ps=12) 

corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P), is.corr = FALSE, col.lim = c(0, 25), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "n", title="B-cell subsets")
n_col <- ncol(as.matrix(sub_data))
for(i in 1:n_col){
  text(x=i-0.25, y=nrow(as.matrix(sub_data))+1, labels=my_labels[i], srt=45, adj=c(0,1), xpd=TRUE, cex=2, col=my_custom_colours[i])
}
dev.off()


#CD8 T cells 

sub_data <- TVal_cat_sub[var_cat, Predictors_8] 
sub_data_P <-  TPVal_cat_sub[var_cat, Predictors_8]

colnames(sub_data) <- rep("", ncol(sub_data))
rownames(sub_data) <- c("sex", "age", "infected", "vaccine brand", "dose number", "time post vaccination", "time into pandemic")

my_labels <- c(rep(c("Naive", "TEM", "TCM", "TEMRA", "CD49a", "CD103", "CXCR6"),2))
my_custom_colours <- c(rep("purple",7),rep("deepskyblue",7))

CC <- colorRampPalette(brewer.pal(3, "YlGnBu"))(200)

pdf("Manuscript/Figures/FigS4B.pdf", width = 20, height = 12)
# Increase bottom margin to make space for sample size labels
par(mar = c(4, 2, 2, 7), xpd=TRUE, ps=12) 

corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P), is.corr = FALSE, col.lim = c(0, 25), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "n", title="CD8 T-cell subsets")
n_col <- ncol(as.matrix(sub_data))
for(i in 1:n_col){
  text(x=i-0.25, y=nrow(as.matrix(sub_data))+1, labels=my_labels[i], srt=45, adj=c(0,1), xpd=TRUE, cex=2, col=my_custom_colours[i])
}
dev.off()


#CD4 cell
sub_data <- TVal_cat_sub[var_cat, Predictors_4] 
sub_data_P <-  TPVal_cat_sub[var_cat, Predictors_4]

colnames(sub_data) <- rep("", ncol(sub_data))
rownames(sub_data) <- c("sex", "age", "infected", "vaccine brand", "dose number", "time post vaccination", "time into pandemic")

my_labels <- c(rep(c("Naive", "TEM", "TCM", "TEMRA", "CD49a", "CD103", "CXCR6"),7))
my_custom_colours <- c(rep("red",7), rep("darkred",7), rep("purple",7),rep("blue",7),rep("deepskyblue",7),rep("forestgreen",7),rep("orange" ,7))

CC <- colorRampPalette(brewer.pal(3, "YlGnBu"))(200)

pdf("Manuscript/Figures/FigS4C.pdf", width = 45, height = 13)
# Increase bottom margin to make space for sample size labels
par(mar = c(4, 2, 2, 6), xpd=TRUE, ps=12) 

corrplot(as.matrix(sub_data), p.mat = as.matrix(sub_data_P), is.corr = FALSE, col.lim = c(0, 25), tl.col = "black", tl.srt = 45, tl.cex=2, cl.cex=2,
         number.cex = 1, pch.cex= 2, method = 'circle', col=CC, sig.level = c(0.001, 0.01, 0.05),  insig = 'label_sig',
         na.label = "X", na.label.col = "black", cl.pos = "b", title="CD4 T-cell subsets")
n_col <- ncol(as.matrix(sub_data))
for(i in 1:n_col){
  text(x=i-0.25, y=nrow(as.matrix(sub_data))+1, labels=my_labels[i], srt=45, adj=c(0,1), xpd=TRUE, cex=2, col=my_custom_colours[i])
}
mtext("", at=-3, side = 3, line = -5, cex = 3, font = 2, adj=0)
text(20.5, -0.75, "Model Improvement (F-statistic)", cex=2, xpd = TRUE)
dev.off()


##########################################
#Elastic Net - Figure 4B
##########################################
##!!!Long computation time!!!#############
##########################################

e_fit_wa <- list()
e_para_wa <- list()
e_bestpara_wa <- list()

for (i in 1:100){
  m.data<- mice::complete(imp, i)
  data <- m.data[, c("wa1", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)]
  data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)] <- lapply(data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)], function(j)  ifelse(j == 0, log10(0.05), log10(j)))
  data[c("age", "tpv", "pandemic", Predictors, Predictors_B, Predictors_8, Predictors_4)] <- scale(data[c("age", "tpv", "pandemic", Predictors, Predictors_B, Predictors_8, Predictors_4)])
  
  data$sex <- relevel(data$sex, ref = "M")
  data$vaccine_brand <- relevel(data$vaccine_brand, ref = "pfizer")
  
  set.seed(123)
  cv_5 <- trainControl(method = "cv",  number = 5)
  
  # Fit Elastic Net model
  hit_elnet <- train(
    log1p(wa1) ~ .,
    data = data,
    method = "glmnet",
    trControl = cv_5,
    tuneGrid = expand.grid(alpha = seq(0,1,by = 0.025), lambda = seq(0.001,2,by = 0.001))
  )
  
  e_para_wa[[i]] <- hit_elnet$results[which(hit_elnet$results$alpha == hit_elnet$bestTune$alpha),]
  e_bestpara_wa[[i]] <- hit_elnet$bestTune
  e_fit_wa[[i]] <- coef(hit_elnet$finalModel, hit_elnet$bestTune$lambda)
  print(i)
  
}

saveRDS(e_fit_wa, file="Data/e_fit_wa.RData")
saveRDS(e_para_wa, file="Data/e_para_wa.RData")
saveRDS(e_bestpara_wa, file="Data/e_bestpara_wa.RData")


e_fit_delta <- list()
e_para_delta <- list()
e_bestpara_delta <- list()

for (i in 1:100){
  m.data<- mice::complete(imp, i)
  data <- m.data[, c("delta", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)]
  data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)] <- lapply(data[c(Predictors, Predictors_B, Predictors_8, Predictors_4)], function(j)  ifelse(j == 0, log10(0.05), log10(j)))
  data[c("age", "tpv", "pandemic", Predictors, Predictors_B, Predictors_8, Predictors_4)] <- scale(data[c("age", "tpv", "pandemic", Predictors, Predictors_B, Predictors_8, Predictors_4)])
  
  data$sex <- relevel(data$sex, ref = "M")
  data$vaccine_brand <- relevel(data$vaccine_brand, ref = "pfizer")
  
  set.seed(123)
  cv_5 <- trainControl(method = "cv",  number = 5)
  
  # Fit Elastic Net model
  hit_elnet <- train(
    log1p(delta) ~ .,
    data = data,
    method = "glmnet",
    trControl = cv_5,
    tuneGrid = expand.grid(alpha = seq(0,1,by = 0.025), lambda = seq(0.001,2,by = 0.001))
  )
  
  e_para_delta[[i]] <- hit_elnet$results[which(hit_elnet$results$alpha == hit_elnet$bestTune$alpha),]
  e_bestpara_delta[[i]] <- hit_elnet$bestTune
  e_fit_delta[[i]] <- coef(hit_elnet$finalModel, hit_elnet$bestTune$lambda)
  print(i)
}

saveRDS(e_fit_delta, file="Data/e_fit_delta.RData")
saveRDS(e_para_delta, file="Data/e_para_delta.RData")
saveRDS(e_bestpara_delta, file="Data/e_bestpara_delta.RData")




e_fit_omicron <- list()
e_para_omicron <- list()
e_bestpara_omicron <- list()

for (i in 1:100){
  m.data<- mice::complete(imp, i)
  data <- m.data[, c("omicron_ba2.12.1", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)]
  data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)] <- lapply(data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)], function(j)  ifelse(j == 0, log10(0.05), log10(j)))
  data[c("age", "tpv", "pandemic", Predictors, Predictors_B, Predictors_8, Predictors_4)] <- scale(data[c("age", "tpv", "pandemic", Predictors, Predictors_B, Predictors_8, Predictors_4)])
  
  data$sex <- relevel(data$sex, ref = "M")
  data$vaccine_brand <- relevel(data$vaccine_brand, ref = "pfizer")
  
  set.seed(123)
  cv_5 <- trainControl(method = "cv",  number = 5)
  
  # Fit Elastic Net model
  hit_elnet <- train(
    log1p(omicron_ba2.12.1) ~ .,
    data = data,
    method = "glmnet",
    trControl = cv_5,
    tuneGrid = expand.grid(alpha = seq(0,1,by = 0.025), lambda = seq(0.001,2,by = 0.001))
  )
  
  e_para_omicron[[i]] <- hit_elnet$results[which(hit_elnet$results$alpha == hit_elnet$bestTune$alpha),]
  e_bestpara_omicron[[i]] <- hit_elnet$bestTune
  e_fit_omicron[[i]] <- coef(hit_elnet$finalModel, hit_elnet$bestTune$lambda)
  print(i)
}

saveRDS(e_fit_omicron, file="Data/e_fit_omicron.RData")
saveRDS(e_para_omicron, file="Data/e_para_omicron.RData")
saveRDS(e_bestpara_omicron, file="Data/e_bestpara_omicron.RData")




e_fit_s <- list()
e_para_s <- list()
e_bestpara_s <- list()

for (i in 1:100){
  m.data<- mice::complete(imp, i)
  data <- m.data[, c("s_igg", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)]
  data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)] <- lapply(data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)], function(j)  ifelse(j == 0, log10(0.05), log10(j)))
  data[c("age", "tpv", "pandemic", Predictors, Predictors_B, Predictors_8, Predictors_4)] <- scale(data[c("age", "tpv", "pandemic", Predictors, Predictors_B, Predictors_8, Predictors_4)])
  
  data$sex <- relevel(data$sex, ref = "M")
  data$vaccine_brand <- relevel(data$vaccine_brand, ref = "pfizer")
  
  set.seed(123)
  cv_5 <- trainControl(method = "cv",  number = 5)
  
  # Fit Elastic Net model
  hit_elnet <- train(
    log1p(s_igg) ~ .,
    data = data,
    method = "glmnet",
    trControl = cv_5,
    tuneGrid = expand.grid(alpha = seq(0,1,by = 0.025), lambda = seq(0.001,2,by = 0.001))
  )
  
  e_para_s[[i]] <- hit_elnet$results[which(hit_elnet$results$alpha == hit_elnet$bestTune$alpha),]
  e_bestpara_s[[i]] <- hit_elnet$bestTune
  e_fit_s[[i]] <- coef(hit_elnet$finalModel, hit_elnet$bestTune$lambda)
  print(i)
}

saveRDS(e_fit_s, file="Data/e_fit_s.RData")
saveRDS(e_para_s, file="Data/e_para_s.RData")
saveRDS(e_bestpara_s, file="Data/e_bestpara_s.RData")




e_fit_rbd <- list()
e_para_rbd <- list()
e_bestpara_rbd <- list()

for (i in 1:100){
  m.data<- mice::complete(imp, i)
  data <- m.data[, c("rbd_igg", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)]
  data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)] <- lapply(data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)], function(j)  ifelse(j == 0, log10(0.05), log10(j)))
  data[c("age", "tpv", "pandemic", Predictors, Predictors_B, Predictors_8, Predictors_4)] <- scale(data[c("age", "tpv", "pandemic", Predictors, Predictors_B, Predictors_8, Predictors_4)])
  
  data$sex <- relevel(data$sex, ref = "M")
  data$vaccine_brand <- relevel(data$vaccine_brand, ref = "pfizer")
  
  set.seed(123)
  cv_5 <- trainControl(method = "cv",  number = 5)
  
  # Fit Elastic Net model
  hit_elnet <- train(
    log1p(rbd_igg) ~ .,
    data = data,
    method = "glmnet",
    trControl = cv_5,
    tuneGrid = expand.grid(alpha = seq(0,1,by = 0.025), lambda = seq(0.001,2,by = 0.001))
  )
  
  e_para_rbd[[i]] <- hit_elnet$results[which(hit_elnet$results$alpha == hit_elnet$bestTune$alpha),]
  e_bestpara_rbd[[i]] <- hit_elnet$bestTune
  e_fit_rbd[[i]] <- coef(hit_elnet$finalModel, hit_elnet$bestTune$lambda)
  print(i)
}

saveRDS(e_fit_rbd, file="Data/e_fit_rbd.RData")
saveRDS(e_para_rbd, file="Data/e_para_rbd.RData")
saveRDS(e_bestpara_rbd, file="Data/e_bestpara_rbd.RData")





###########################################################
### Boruta Random Forrest - Figure 4C
###########################################################


Feat_var_wa1 <- list()

for (i in 1:100) {
  m.data<- mice::complete(imp, i)
  data <- m.data[, c("wa1", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)]
  
  data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)] <- lapply(data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)], function(j)  ifelse(j == 0, log10(0.05), log10(j)))
  #data[c("age", "tpv", "pandemic", Predictors, Predictors_B, Predictors_8, Predictors_4)] <- scale(data[c("age", "tpv", "pandemic", Predictors, Predictors_B, Predictors_8, Predictors_4)])
  
  set.seed(1)
  B <- Boruta(log1p(wa1) ~., data=data)
  Feat_var_wa1[[i]] <- attStats(B)$decision
  print(i)
}

#Features<- data.frame(do.call(rbind, Feat_var))
#Features %>% summarise_all(~ sum(.== 1))

Features_wa<- data.frame(do.call(rbind, Feat_var_wa1))
colnames(Features_wa)<- c("wa1", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)
Features_wa %>% summarise_all(~ sum(.== 2))

saveRDS(Feat_var_wa1, file="Data/Feat_var_wa1.RData")


Feat_var_delta <- list()

for (i in 1:100) {
  m.data<- mice::complete(imp, i)
  data <- m.data[, c("delta", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)]
  
  data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)] <- lapply(data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)], function(j)  ifelse(j == 0, log10(0.05), log10(j)))
 
  set.seed(1)
  B <- Boruta(log1p(delta) ~., data=data)
  Feat_var_delta[[i]] <- attStats(B)$decision
  print(i)
}

#Features<- data.frame(do.call(rbind, Feat_var))
#Features %>% summarise_all(~ sum(.== 1))

Features_delta<- data.frame(do.call(rbind, Feat_var_delta))
colnames(Features_delta)<- c("delta", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)
Features_delta %>% summarise_all(~ sum(.== 2))

saveRDS(Feat_var_delta, file="Data/Feat_var_delta.RData")

Feat_var_omicron <- list()

for (i in 1:100) {
  m.data<- mice::complete(imp, i)
  data <- m.data[, c("omicron_ba2.12.1", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)]
  
  data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)] <- lapply(data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)], function(j)  ifelse(j == 0, log10(0.05), log10(j)))

  set.seed(1)
  B <- Boruta(log1p(omicron_ba2.12.1) ~., data=data)
  Feat_var_omicron[[i]] <- attStats(B)$decision
  print(i)
}

#Features<- data.frame(do.call(rbind, Feat_var))
#Features %>% summarise_all(~ sum(.== 1))

Features_omicron<- data.frame(do.call(rbind, Feat_var_omicron))
colnames(Features_omicron)<- c("omicron_ba2.12.1", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)
Features_omicron %>% summarise_all(~ sum(.== 2))

saveRDS(Feat_var_omicron, file="Data/Feat_var_omicron.RData")

Feat_var_rbd <- list()

for (i in 1:100) {
  m.data<- mice::complete(imp, i)
  data <- m.data[, c("rbd_igg", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)]
  
  data[c( Predictors, Predictors_B, Predictors_4, Predictors_8)] <- lapply(data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)], function(j)  ifelse(j == 0, log10(0.05), log10(j)))

  set.seed(1)
  B <- Boruta(log1p(rbd_igg) ~., data=data)
  Feat_var_rbd[[i]] <- attStats(B)$decision
  print(i)
}

#Features<- data.frame(do.call(rbind, Feat_var))
#Features %>% summarise_all(~ sum(.== 1))

Features_rbd<- data.frame(do.call(rbind, Feat_var_rbd))
colnames(Features_rbd)<- c("rbd_igg", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)
Features_rbd %>% summarise_all(~ sum(.== 2))

saveRDS(Feat_var_rbd, file="Data/Feat_var_rbd.RData")

Feat_var_s <- list()

for (i in 1:100) {
  m.data<- mice::complete(imp, i)
  data <- m.data[, c("s_igg", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)]
  
  data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)] <- lapply(data[c( Predictors, Predictors_B, Predictors_8, Predictors_4)], function(j)  ifelse(j == 0, log10(0.05), log10(j)))

  set.seed(1)
  B <- Boruta(log1p(s_igg) ~., data=data)
  Feat_var_s[[i]] <- attStats(B)$decision
  print(i)
}

#Features<- data.frame(do.call(rbind, Feat_var))
#Features %>% summarise_all(~ sum(.== 1))

Features_s<- data.frame(do.call(rbind, Feat_var_s))
colnames(Features_s)<- c("s_igg", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)
Features_s %>% summarise_all(~ sum(.== 2))

saveRDS(Feat_var_s, file="Data/Feat_var_s.RData")


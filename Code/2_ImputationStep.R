###IMPUTATION STEP
#Packages
source("Code/0_Packages.R")

#Load dataset
load("Data/Dataset_DARPA.RData")
data<- as.data.frame(Dataset_DARPA)

imp <- mice(data, printFlag = FALSE, seed = 1, method="rf", m=100, pred = quickpred(data, method = "spearman", mincor= .3, minpuc = 0.25), maxit=5, remove.collinear = FALSE)
write.mice.imputation(mi.res=imp, name="imp", mids2spss=F)
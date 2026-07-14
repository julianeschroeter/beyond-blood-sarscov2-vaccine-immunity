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



###Figure 4B
#read in files
e_fit_wa<- readRDS("Data/e_fit_wa.RData")
e_fit_delta<- readRDS("Data/e_fit_delta.RData")
e_fit_omicron<- readRDS("Data/e_fit_omicron.RData")
e_fit_s<- readRDS("Data/e_fit_s.RData")
e_fit_rbd<- readRDS("Data/e_fit_rbd.RData")


EN_plot <- function(x, title){

#x<- e_fit_s
#title <- "S_IgG"
    
output_list <- lapply(x, function(m) as.matrix(m))

# Combine into a data frame
en_output <- as.data.frame(do.call(cbind, output_list))
colnames(en_output)<- c(1:100)

N_predictors <- summary(colSums(en_output != 0)-1)

plus <- rowSums(en_output > 0)
minus <- rowSums(en_output < 0)

# Combine into a data frame
V <- data.frame(
  name = rownames(en_output),
  Positive = plus,
  Negative = minus
)
V<- V[-1, ]

# Reshape for ggplot (long format)
V_long <- reshape2::melt(V, id.vars = "name", variable.name = "Type", value.name = "Count")

# Order by total appearance (descending)
V$Majority <- ifelse(V$Positive >= V$Negative, V$Positive, V$Negative)
V$total <- V$Positive + V$Negative
V <- V[order(V$Majority), ]

last_25_names <- tail(V$name, 25)

# Subset V_long to include only those last 20 factors
V_long_sub <- subset(V_long, name %in% last_25_names)

# Reset factor levels to preserve the order for plotting
V_long_sub$name <- factor(V_long_sub$name, levels = last_25_names)

V_long_sub$Count_signed <- ifelse(V_long_sub$Type == "Negative", -V_long_sub$Count, V_long_sub$Count)

#Effect size 
en_output[en_output == 0] <- NA
en_output_t <- t(en_output)

en_sub <- en_output_t[, last_25_names, drop = FALSE]

WW<- lapply(as.data.frame(en_sub), median, na.rm=T)

WW_long <- melt(WW)
WW_long$L1 <- factor(WW_long$L1, levels = last_25_names)
WW_long$value <- as.numeric(round(WW_long$value,2))
WW_long <- WW_long[order(WW_long$L1), ]
WW_long$label_y <- ifelse(WW_long$value >= 0, 115, 115) 
WW_long$hjust   <- ifelse(WW_long$value >= 0, 1, 1)  



p <- ggplot(V_long_sub, aes(x = name, y = Count, fill = Type)) +
  geom_hline(yintercept=0, lty=1, col="darkgrey")+
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Positive" = "lightcoral", "Negative" = "lightskyblue"))+
  coord_flip() +
  scale_y_continuous(limits = c(0,115), breaks=c(0,20,40,60,80,100))+
  ylab("Appearance Score") +
  xlab("") +
  theme_bw() + ggtitle(title)+
  geom_hline(yintercept=90, lty=2, col="darkgrey")+
  theme(axis.text.y = element_text(size = 8), legend.position = "none", axis.text.x = element_text(size = 8))+
  geom_text(data = WW_long,  aes(x = L1, y =  label_y, label = value, hjust = hjust),  # fixed y-position slightly outside bars
    inherit.aes = FALSE, size = 2)

return(list(plot = p, name_25 = last_25_names, N_predictors = N_predictors))
}

E1 <- EN_plot(e_fit_s, "S-IgG")
E2 <- EN_plot(e_fit_rbd, "RBD-IgG")
E3 <- EN_plot(e_fit_wa, "WA1")
E4 <- EN_plot(e_fit_delta, "Delta")
E5 <- EN_plot(e_fit_omicron, "Omicron")

#annotated by hand
e1 <- E1[[1]]+ scale_x_discrete(labels = c("CD4 TEMRA", "mBC" , "Time post vaccination", "TFH", "Treg", "CD8", "CD4", "Sex", "CD8 CD103", "CD69 mBc", "Age", "Dose number", "mBc", "CD4 naive", "IgG mBc", "IgM mBc", "CD8", "mBc", "IgA mBc", "CD4 TCM", "CD4", "Infected", "CD4 CD49a", "Vaccine brand", "IgG mBc"))+
  theme(axis.text.y = element_text(colour =c("deepskyblue", "deepskyblue", "black", "purple", "blue", "red", "darkred", "black", "deepskyblue", "purple", "black", "black", "blue", "deepskyblue", "darkred", "deepskyblue","purple", "purple", "purple", "blue", "red", "black", "blue", "black", "deepskyblue")))

#e1 <- E1[[1]]+ scale_x_discrete(labels = c("CD4 CXCR6", "IgA mBC", "CD4 TEMRA", "Sex", "CD4 CD49a", "TFH",  "CD4 naive", "CD4", "Age", "IgG mBc", "mBc", "CD4 TCM", "IgM mBC", "CD8 CD103", "Dose number", "mBc", "CD69 mBc", "CD4 CD49a", "CD4", "IgA mBc", "CD4 naive",  "Infected",  "Vaccine brand", "IgG mBc", "CD8"))+
#  theme(axis.text.y = element_text(colour =c("darkred", "blue", "red", "black", "forestgreen","purple", "blue", "deepskyblue", "black", "darkred", "blue", "blue", "deepskyblue", "deepskyblue", "black", "purple","purple", "blue", "red", "purple", "deepskyblue", "black", "black", "deepskyblue",  "purple")))

#e1 <- E1[[1]]+ scale_x_discrete(labels = c("mBc", "CD4", "CD4", "IgA+", "CD4+ CD103+", "Treg", "mBc", "CD8 TCM", "CD4 TCM", "Age", "CD4 naive", "Dose number", "mBc CD69+", "CD4 naive", "TFH", "mBc", "IgA+", "CD8+ CD103+", "CD4 naive", "TFH", "Infected", "CD4+ CD49a+", "Vaccine brand", "CD8", "IgG+"))+
 # theme(axis.text.y = element_text(colour =c("blue", "red", "deepskyblue", "blue", "red", "blue", "deepskyblue", "deepskyblue", "blue", "black", "blue", "black","purple", "red", "red", "purple", "purple", "deepskyblue", "deepskyblue", "purple", "black", "blue", "black", "purple", "deepskyblue")))

e2 <- E2[[1]]+ scale_x_discrete(labels = c("CD8 TEM", "CD8 naive", "CD4 TEM", "CD4 TCM", "mBc", "CD8", "CD4", "IgG mBc", "CD4", "CD69 mBc", "IgM mBc", "CD8 CD49a", "CD4 CD49a", "TFH", "CD4 TEMRA", "Vaccine brand", "IgA mBc", "CD4 naive", "mBc", "CD8", "IgG mBc", "Sex", "mBc", "CD4 CD49a", "Infected"))+
  theme(axis.text.y = element_text(colour =c("purple", "purple", "deepskyblue", "blue", "darkred", "deepskyblue", "red", "darkred","deepskyblue", "purple", "deepskyblue", "purple","purple", "purple", "deepskyblue","black", "purple", "deepskyblue", "blue", "purple", "deepskyblue", "black", "purple", "blue", "black")))

#e2 <- E2[[1]]+ scale_x_discrete(labels = c("IgG mBc", "CD4 naive", "CD4 naive", "CD4 CD49a", "CD8 CD103", "CD4", "Age", "CD4", "CD69 mBc", "IgG mBc", "CD4 CD49a", "Sex", "TFH", "mBc", "CD4 TEMRA", "IgM mBc", "CD4 CD49a", "mBc", "CD4", "IgA mBc", "Vaccine brand", "CD4 naive", "IgG mBc", "Infected", "CD8"))+
 # theme(axis.text.y = element_text(colour =c("red", "blue", "forestgreen", "forestgreen", "deepskyblue", "red", "black", "purple", "purple","darkred", "purple", "black", "purple","purple", "deepskyblue", "deepskyblue","blue", "blue", "deepskyblue", "purple", "black", "deepskyblue", "deepskyblue", "black", "purple" )))

#e2 <- E2[[1]]+ scale_x_discrete(labels = c("IgA+", "CD4+ CXCR6+", "mBc CD69+", "mBc", "mBc", "Age", "IgA+", "CD4 TEMRA", "IgG+", "CD4 naive", "Sex", "CD4+ CXCR6+", "mBc", "TFH", "CD4", "CD8 naive", "CD8+ CD103+", "TFH", "Vaccine brand", "CD4+ CD49a+", "mBc", "CD4 naive", "CD8", "Infected", "IgG+"))+
  #theme(axis.text.y = element_text(colour =c("blue", "deepskyblue", "purple", "deepskyblue", "darkred", "black", "purple", "deepskyblue", "red", "red", "black", "purple","blue", "red", "deepskyblue", "purple", "deepskyblue", "purple", "black", "blue", "purple", "deepskyblue", "purple", "black", "deepskyblue")))

e3 <- E3[[1]]+ scale_x_discrete(labels = c("CD4 CXCR6", "CD4 CD49a", "CD4", "CD8 CD49a", "CD4 TEMRA", "CD4 TCM", "CD4 CD49a", "CD4 TEMRA", "CD69 mBc", "CD4 CXCR6", "CD8 CD103", "CD4", "CD8", "IgA mBc", "Sex", "mBc", "CD4 CD49a" ,"Dose number", "mBc", "Vaccine brand", "IgM mBc", "mBc", "IgG mBc", "Infected" , "IgG mBc"))+
  theme(axis.text.y = element_text(colour =c("purple", "forestgreen", "purple", "purple", "deepskyblue", "red", "darkred", "purple", "purple", "deepskyblue", "deepskyblue","red", "blue", "purple", "black", "deepskyblue", "blue", "black",  "blue",  "black" ,"deepskyblue", "purple","darkred", "black", "deepskyblue")))

#e3 <- E3[[1]]+ scale_x_discrete(labels = c("CD4", "CD8+ CD103+", "CD4+ CD49a+", "CD4+ CD103+", "CD4 TEMRA", "CD4 TEMRA", "CD4 naive", "TFH", "CD8+ CD49a+", "CD8", "CD69+ mBc", "mBc", "CD4+ CD49a+", "Sex", "IgM+ mBc", "Dose number", "CD8+ CD103+", "Vaccine brand", "mBc", "CD4", "IgA+ mBc", "IgG+ mBc", "mBc", "IgG+ mBc", "Infected"))+
 # theme(axis.text.y = element_text(colour =c("purple", "deepskyblue", "darkred", "darkred", "deepskyblue", "purple", "deepskyblue", "deepskyblue","deepskyblue", "purple", "purple", "deepskyblue", "blue", "black", "deepskyblue", "black" ,"purple", "black","blue", "red","purple","darkred","purple", "deepskyblue", "black")))

#e3 <- E3[[1]]+ scale_x_discrete(labels = c("CD4+ CD49a+", "CD4 TCM", "TFH", "CD4+ CD103+", "CD4 TEM", "CD4+ CXCR6+", "IgA+", "IgM+", "CD8", "CD4 naive", "IgA+", "Sex", "CD4+ CXCR6+", "CD8+ CD103+", "mBc", "TFH", "Dose number", "CD4+ CD49a+", "Vaccine brand", "CD4+ naive", "mBc", "CD8+ CD103+", "mBc", "IgG+" , "Infected"))+
 # theme(axis.text.y = element_text(colour =c("darkred", "red", "purple", "red", "forestgreen", "purple", "darkred", "deepskyblue",  "purple", "deepskyblue", "purple", "black", "deepskyblue", "purple", "blue", "deepskyblue", "black", "blue", "black", "red", "deepskyblue","deepskyblue", "purple", "deepskyblue", "black")))

e4 <- E4[[1]]+ scale_x_discrete(labels = c("IgM mBc", "Dose number", "TFH", "CD8", "CD4 CD103", "mBc", "Sex",  "CD4 TCM", "CD4", "Vaccine brand", "CD4 CD49a", "CD4 CD49a", "CD69 mBc", "IgG mBc", "IgA mBc", "Time into pandemic", "mBc", "CD4", "IgM mBc", "CD4 naive", "CD4 TEMRA" , "CD8", "IgG mBc", "CD8", "Infected"))+
  theme(axis.text.y = element_text(colour =c("purple", "black", "purple", "blue", "blue", "deepskyblue","black", "blue", "red", "black", "blue", "purple","deepskyblue", "darkred", "purple", "black",  "purple", "purple", "deepskyblue","deepskyblue","deepskyblue","deepskyblue", "deepskyblue", "purple", "black")))

#e4 <- E4[[1]]+ scale_x_discrete(labels = c("CD4 CD49a", "CD4 TCM", "CD4 CD49a", "TFH", "CD4", "CD4 TEMRA", "CD4 CD49a", "CD4 TCM", "mBc", "CD69 mBc", "CD8 CD49a", "CD4 CD49a", "CD4 naive", "mBc", "CD4", "CD4 naive", "Vaccine brand", "CD8", "IgM mBc", "IgG mBc" , "CD8", "IgA mBc", "Time into pandemic", "IgG mBc", "Infected"))+
  #theme(axis.text.y = element_text(colour =c("forestgreen", "purple", "purple", "purple", "purple", "deepskyblue","red", "blue",  "deepskyblue", "deepskyblue", "deepskyblue", "blue", "forestgreen", "purple","deepskyblue","deepskyblue",  "black", "deepskyblue" ,"deepskyblue", "darkred", "purple", "purple", "black", "deepskyblue", "black")))

#e4 <- E4[[1]]+ scale_x_discrete(labels = c("TFH","IgM+", "CD4+ CD49a+", "CD4 TCM", "CD4 TEM", "mBc CD69+", "IgA+", "CD4 naive", "CD4+ CD49a+", "CD4+ CD49a+", "CD4", "CD4", "CD4+ CXCR6+", "mBc", "mBc CD69+", "CD4 naive", "CD8+ CD103+", "Vaccine brand", "IgG+", "CD8", "mBc", "IgG+" , "CD8", "Time into pandemic", "Infected"))+
 # theme(axis.text.y = element_text(colour =c("purple", "forestgreen", "purple", "blue", "forestgreen", "purple", "purple", "red",  "blue", "deepskyblue", "deepskyblue", "purple", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "black", "darkred", "purple", "purple","deepskyblue", "deepskyblue", "black", "black")))

e5 <- E5[[1]]+ scale_x_discrete(labels = c("CD8", "CD4 TEM", "CD4 CD103",  "mBc", "mBc", "CD69 mBc", "CD4 CD103", "CD4 CD49a", "Vaccine brand" ,"IgA mBc" , "IgG mBc", "CD4 TEMRA", "Sex", "CD4 CD49a", "CD4", "CD69 mBc" , "TFH", "Time into pandemic", "IgM mBc", "IgA mBC",  "mBc", "CD4 CD49a", "IgG mBc", "CD8", "Infected")) +
  theme(axis.text.y = element_text(colour =c("purple", "forestgreen","purple", "blue", "deepskyblue", "deepskyblue", "deepskyblue", "forestgreen", "black", "darkred", "darkred",  "deepskyblue", "black", "purple", "red","purple", "purple", "black","deepskyblue", "purple", "purple", "blue", "deepskyblue","deepskyblue", "black")))

#e5 <- E5[[1]]+ scale_x_discrete(labels = c("CD4+ CXCR6+", "TFH", "CD4+ CD49a+", "Sex", "CD4+ CD49a+", "CD4", "CD69+ mBc", "CD4+ naive", "CD4+ CD49a+", "IgA+ mBc", "Vaccine brand" , "CD8", "CD69+ mBc", "IgG+ mBc", "CD4+ CD49a+", "CD4 TEMRA", "TFH", "IgA+ mBc" , "mBc", "Time into pandemic", "IgM+ mBc", "CD8", "IgG+ mBc", "CD4 TEMRA", "Infected")) +
 # theme(axis.text.y = element_text(colour =c("purple", "deepskyblue","deepskyblue", "black", "purple", "red","deepskyblue", "deepskyblue", "forestgreen","darkred",  "black", "deepskyblue", "purple", "darkred",  "blue", "red", "purple", "purple","purple", "black","deepskyblue", "purple","deepskyblue","deepskyblue", "black")))

#e5 <- E5[[1]]+ scale_x_discrete(labels = c("IgM+", "mBc", "TFH", "mBc CD69+", "CD4+ CD49a+", "CD8+ CD103+", "mBc CD69+", "CD4+ CXCR6+", "CD4+ CD103+", "Sex", "CD8", "IgM+", "CD4 naive", "CD4+ CD49a+", "CD4 TEMRA", "CD4 TEM", "Vaccine brand", "IgA+", "IgA+" , "TFH", "Time into pandemic", "IgG+", "CD8", "mBc", "Infected")) +
 # theme(axis.text.y = element_text(colour =c("forestgreen", "deepskyblue","deepskyblue","deepskyblue", "blue","deepskyblue", "purple", "purple","deepskyblue",  "black", "deepskyblue", "deepskyblue", "deepskyblue",  "purple", "deepskyblue", "forestgreen", "black", "purple", "darkred", "purple", "black","deepskyblue", "purple", "purple", "black")))

#Figure 4B
Fig4B_rename <- plot_grid(E1[[1]], E2[[1]], E3[[1]], E4[[1]], E5[[1]], ncol=5, align = "hv")

Fig4B <- plot_grid(e1, e2, e3, e4, e5, ncol=5, align = "hv")


#Figure S4A
EN_effectsize <- function(x, title, names_25){
  
  output_list <- lapply(x, function(m) as.matrix(m))
  
  # Combine into a data frame
  en_output <- as.data.frame(do.call(cbind, output_list))
  colnames(en_output)<- c(1:100)
  
  en_output[en_output == 0] <- NA
  en_output_t <- t(en_output)
  
  en_sub <- en_output_t[, names_25, drop = FALSE]
  
  en_sub_long <- melt(en_sub)
  colnames(en_sub_long) <- c("Row", "Factor", "Value") 
  
  p <- ggplot(en_sub_long, aes(y = Factor, x = Value)) +
    #
    geom_vline(xintercept=0, lty=2, col="black")+
    geom_boxplot(na.rm = TRUE)+
    scale_x_continuous(limits = c(-0.75, 2))+
    #scale_x_continuous(limits = c(-log10(0.75), log10(2.0)))+
    ylab("") +
    xlab("Effect size") +
    theme_bw() + 
    ggtitle(title)+
    theme(axis.text.y = element_text(size = 6, hjust = 1), legend.position = "none")
  
  p
}


F1 <- EN_effectsize(e_fit_s, "S-IgG", E1[[2]])
F2 <- EN_effectsize(e_fit_rbd,"RBD-IgG",E2[[2]])
F3 <- EN_effectsize(e_fit_wa, "WA1", E3[[2]])
F4 <- EN_effectsize(e_fit_delta, "Delta", E4[[2]])
F5 <- EN_effectsize(e_fit_omicron, "Omicron", E5[[2]])

#annotated by hand
f1 <- F1+ scale_y_discrete(labels = c("CD4 TEMRA", "mBC" , "Time post vaccination", "TFH", "Treg", "CD8", "CD4", "Sex", "CD8 CD103", "CD69 mBc", "Age", "Dose number", "mBc", "CD4 naive", "IgG mBc", "IgM mBc", "CD8", "mBc", "IgA mBc", "CD4 TCM", "CD4", "Infected", "CD4 CD49a", "Vaccine brand", "IgG mBc"))+
  theme(axis.text.y = element_text(colour =c("deepskyblue", "deepskyblue", "black", "purple", "blue", "red", "darkred", "black", "deepskyblue", "purple", "black", "black", "blue", "deepskyblue", "darkred", "deepskyblue","purple", "purple", "purple", "blue", "red", "black", "blue", "black", "deepskyblue")))

f2 <- F2+ scale_y_discrete(labels = c("CD8 TEM", "CD8 naive", "CD4 TEM", "CD4 TCM", "mBc", "CD8", "CD4", "IgG mBc", "CD4", "CD69 mBc", "IgM mBc", "CD8 CD49a", "CD4 CD49a", "TFH", "CD4 TEMRA", "Vaccine brand", "IgA mBc", "CD4 naive", "mBc", "CD8", "IgG mBc", "Sex", "mBc", "CD4 CD49a", "Infected"))+
  theme(axis.text.y = element_text(colour =c("purple", "purple", "deepskyblue", "blue", "darkred", "deepskyblue", "red", "darkred","deepskyblue", "purple", "deepskyblue", "purple","purple", "purple", "deepskyblue","black", "purple", "deepskyblue", "blue", "purple", "deepskyblue", "black", "purple", "blue", "black")))

f3 <- F3+ scale_y_discrete(labels = c("CD4 CXCR6", "CD4 CD49a", "CD4", "CD8 CD49a", "CD4 TEMRA", "CD4 TCM", "CD4 CD49a", "CD4 TEMRA", "CD69 mBc", "CD4 CXCR6", "CD8 CD103", "CD4", "CD8", "IgA mBc", "Sex", "mBc", "CD4 CD49a" ,"Dose number", "mBc", "Vaccine brand", "IgM mBc", "mBc", "IgG mBc", "Infected" , "IgG mBc"))+
  theme(axis.text.y = element_text(colour =c("purple", "forestgreen", "purple", "purple", "deepskyblue", "red", "darkred", "purple", "purple", "deepskyblue", "deepskyblue","red", "blue", "purple", "black", "deepskyblue", "blue", "black",  "blue",  "black" ,"deepskyblue", "purple","darkred", "black", "deepskyblue")))

f4 <- F4+ scale_y_discrete(labels = c("IgM mBc", "Dose number", "TFH", "CD8", "CD4 CD103", "mBc", "Sex",  "CD4 TCM", "CD4", "Vaccine brand", "CD4 CD49a", "CD4 CD49a", "CD69 mBc", "IgG mBc", "IgA mBc", "Time into pandemic", "mBc", "CD4", "IgM mBc", "CD4 naive", "CD4 TEMRA" , "CD8", "IgG mBc", "CD8", "Infected"))+
  theme(axis.text.y = element_text(colour =c("purple", "black", "purple", "blue", "blue", "deepskyblue","black", "blue", "red", "black", "blue", "purple","deepskyblue", "darkred", "purple", "black",  "purple", "purple", "deepskyblue","deepskyblue","deepskyblue","deepskyblue", "deepskyblue", "purple", "black")))

f5 <- F5+ scale_y_discrete(labels = c("CD8", "CD4 TEM", "CD4 CD103",  "mBc", "mBc", "CD69 mBc", "CD4 CD103", "CD4 CD49a", "Vaccine brand" ,"IgA mBc" , "IgG mBc", "CD4 TEMRA", "Sex", "CD4 CD49a", "CD4", "CD69 mBc" , "TFH", "Time into pandemic", "IgM mBc", "IgA mBC",  "mBc", "CD4 CD49a", "IgG mBc", "CD8", "Infected")) +
  theme(axis.text.y = element_text(colour =c("purple", "forestgreen","purple", "blue", "deepskyblue", "deepskyblue", "deepskyblue", "forestgreen", "black", "darkred", "darkred",  "deepskyblue", "black", "purple", "red","purple", "purple", "black","deepskyblue", "purple", "purple", "blue", "deepskyblue","deepskyblue", "black")))


FigS5A <- plot_grid(f1, f2, f3, f4, f5, nrow=1, align = "hv", ncol=5)



###Boruta plots - Figure 4C

Feat_var_s<- readRDS("Data/Feat_var_s.RData")
Feat_var_rbd<- readRDS("Data/Feat_var_rbd.RData")
Feat_var_wa1<- readRDS("Data/Feat_var_wa1.RData")
Feat_var_delta<- readRDS("Data/Feat_var_delta.RData")
Feat_var_omicron<- readRDS("Data/Feat_var_omicron.RData")


Boruta_plot <- function(x, title){

W<-   data.frame(do.call(rbind, x))
colnames(W)<- c("AB", var_cat, Predictors, Predictors_B, Predictors_8, Predictors_4)

Confirmed <- summary(rowSums(W == 2)-1)
N_Predictors <- summary(rowSums(W != 3)-1)

W_sub <- as.data.frame(cbind(as.integer(as.vector(W %>% summarise_all(~ sum(.== 2)))[-1])))
W_sub <- cbind(W_sub, as.integer(as.vector(W %>% summarise_all(~ sum(.== 1)))[-1])) # Tentative
colnames(W_sub) <- c("Confirmed", "Tentative")
W_sub$name <- colnames(W)[-1]

W_sub <- W_sub[order(W_sub$Confirmed), ]

more_50_names <- W_sub$name[which(W_sub$Confirmed+W_sub$Tentative >= 50)]

# Reshape for ggplot (long format)
W_long <- reshape2::melt(W_sub, id.vars = "name", variable.name = "Type", value.name = "Count")

W_long_sub <- subset(W_long, name %in% more_50_names)

# Reset factor levels to preserve the order for plotting
W_long_sub$name <- factor(W_long_sub$name, levels = more_50_names)
W_long_sub$Type <- factor(W_long_sub$Type, levels = c("Tentative", "Confirmed"))

p<- ggplot(W_long_sub, aes(x = name, y = Count, fill = Type)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Confirmed" = "darkgrey", "Tentative" = "lightgrey")) +
  scale_y_continuous(limits = c(0,100), breaks=c(0,20,40,60,80,100))+
  coord_flip() +
  ylab("Appearance Score") +
  xlab("") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8), legend.position = "none")+
  geom_hline(yintercept=0, lty=1, col="darkgrey")+
  geom_hline(yintercept=50, lty=2, col="darkgrey")

legend_plot <- ggplot(W_long_sub, aes(x = name, y = Count, fill = Type)) + theme_bw()+
  geom_col()+
  scale_fill_manual(values = c("Confirmed" = "darkgrey", "Tentative" = "lightgrey"))+
  theme(legend.position = "bottom")+ guides(fill = guide_legend(title = NULL))
  
  return(list(plot = p, legend = legend_plot, Confirmed = Confirmed, N_Predictors=N_Predictors))
}


G1 <- Boruta_plot(Feat_var_s, "S-IgG")
G2 <- Boruta_plot(Feat_var_rbd, "RBD-IgG")
G3 <- Boruta_plot(Feat_var_wa1, "WA1")
G4 <- Boruta_plot(Feat_var_delta, "Delta")
G5 <- Boruta_plot(Feat_var_omicron, "Omicron")


#annotated by hand
g1 <- G1[[1]]+ scale_x_discrete(labels = c("IgM mBc", "CD4", "IgG mBc"))+
  theme(axis.text.y = element_text(colour =c("deepskyblue", "red", "deepskyblue")))

g2 <- G2[[1]]+ scale_x_discrete(labels = c("CD4 naive", "IgM mBc", "IgG mBC", "Infected", "CD8"))+
  theme(axis.text.y = element_text(colour =c("orange", "deepskyblue", "deepskyblue", "black", "deepskyblue" )))

g3 <- G3[[1]]+ scale_x_discrete(labels = c("Time into pandemic", "Infected", "IgG mBc", "mBc", "IgM mBc"))+
  theme(axis.text.y = element_text(colour =c("black",  "black",  "deepskyblue", "purple", "deepskyblue")))

g4 <- G4[[1]]+ scale_x_discrete(labels = c("IgG mBc", "IgM mBc", "CD69 mBc", "CD4 naive", "CD8", "Infected", "Time into pandemic"))+
  theme(axis.text.y = element_text(colour =c("deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "deepskyblue", "black",  "black")))

g5 <- G5[[1]]+ scale_x_discrete(labels = c("mBc", "CD69 mBc", "CD4 CD49a", "IgG mBc", "IgM mBc", "CD8", "Infected", "Time into pandemic")) +
  theme(axis.text.y = element_text(colour =c("purple", "deepskyblue", "blue","deepskyblue", "deepskyblue", "deepskyblue","black","black")))


Fig4C_rename <- plot_grid(G1[[1]], G2[[1]], G3[[1]], G4[[1]], G5[[1]], nrow=1, align = "vh", ncol=5)

Fig4C <- plot_grid(g1, g2, g3, g4, g5, nrow=1, align = "vh", ncol=5)

shared_legend4 <- get_legend(G1[[2]])


#####################################
###Figure 4
#####################################

Fig4C <- plot_grid(Fig4C, shared_legend4, ncol = 1, rel_heights = c(1, 0.1))

Fig4BC <- plot_grid(Fig4B, Fig4C,  labels = c('B', 'C'), nrow=2, rel_heights = c(2,1))

Fig4A <- plot_grid(NULL,  labels = c('A'), nrow=1, rel_heights = c(1), align="hv", rel_weidth = c(3))

Figure4 <- plot_grid(Fig4A, Fig4BC, nrow=2, labels = c('A', ''), rel_heights = c(4,8))

save_plot("Manuscript/Figures/Figure4_Final.pdf", Figure4, nrow=3, ncol=2.5)


#####################################
###Figure S4
#####################################

save_plot("Manuscript/Figures/FigureS5_Final.pdf", FigS5A, nrow=1, ncol=3)


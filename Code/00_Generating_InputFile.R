# GENERATE INPUT FILE
source("Code/0_Packages.R")

# Load Demographics and Antibody titers
Data_Demo <- read_csv("Data/Data_DARPA_1.csv")

#Make typo changes
Data_Demo$race_ethnicity[which(Data_Demo$race_ethnicity == "whie")] <- "white"

#Dates
Data_Demo[Data_Demo == "augus"] <- "august"
Data_Demo[Data_Demo == "sepember"] <- "september"
Data_Demo[Data_Demo == "ocober"] <- "october"

Data_Demo$deathmonth <- capitalize(Data_Demo$deathmonth)

#Ensure that data type is correct
Data_Demo$donor <- as.character(Data_Demo$donor)
Data_Demo$age <- as.numeric(Data_Demo$age)
Data_Demo$tpv <- as.numeric(Data_Demo$tpv)
Data_Demo$sex <- as.factor(Data_Demo$sex)
Data_Demo$race_ethnicity <- as.factor(Data_Demo$race_ethnicity)
Data_Demo$batch <- as.factor(Data_Demo$batch)
Data_Demo$dose_number <- as.factor(Data_Demo$dose_number)
Data_Demo$vaccine_brand <- as.factor(Data_Demo$vaccine_brand)
Data_Demo$vaccinated <- as.factor(Data_Demo$vaccinated)
Data_Demo$infected <- as.factor(Data_Demo$infected)

Dem_var<- c("donor", "age", "sex", "race_ethnicity", "batch", "deathmonth", "deathyear", "tpv", "dose_number", "vaccine_brand", "vaccinated", "infected")
outcome <- c("s_igg", "rbd_igg", "wa1", "delta", "omicron_ba2.12.1")

Data1 <- Data_Demo[,c(Dem_var, outcome)]

#Read in tissue frequencies
Data_CD4 <- read_excel("Data/Julia/MasterRemake_040326.xlsx", sheet = "CD4")
Data_CD8 <- read_excel("Data/Julia/MasterRemake_040326.xlsx", sheet = "CD8")
Data_B <- read_excel("Data/Julia/MasterRemake_040326.xlsx", sheet = "B")

##Tissue variable
Tissue <- c("bld", "bm", "spl", "lng", "lln", "mln", "iln")

#Subsets
Subset_4 <- c("scd4", "scd4_naive", "scd4_tem", "scd4_tcm", "scd4_temra","scd4_cd49a", "scd4_cd103", "scd4_cxcr6", "stfh", "streg")
Predictors_4 <- c()

for (i in Tissue){
  for(j in Subset_4){
    V<- paste(i, j, sep="_")
    Predictors_4 <- c(Predictors_4, V)
  }
}

Subset_8 <- c("scd8_naive", "scd8_tem", "scd8_tcm", "scd8_temra","scd8_cd49a", "scd8_cd103", "scd8_cxcr6")
Predictors_8 <- c()

for (i in Tissue){
  for(j in "scd8"){
    V<- paste(i, j, sep="_")
    Predictors_8 <- c(Predictors_8, V)
  }
}

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


###Merge the right Experiments together

D<- unique(c(Data_CD4$Donor, Data_CD8$Donor))
D_N <- grep("_N$", D, value = TRUE)

for (i in D_N){
  j <- sub("_N$", "", i)
  #Checked for a set of tissue in CD4 and CD8 manually
  #Identify the right rows
  i4 <- which(Data_CD4$Donor ==i)
  j4 <- which(Data_CD4$Donor == j)
  
  i8 <- which(Data_CD8$Donor ==i)
  j8 <- which(Data_CD8$Donor == j)
  
  #Replace the other values with what is already in the dataset
  Data_CD4[i4,] <- ifelse(is.na(Data_CD4[i4,]), Data_CD4[j4,], Data_CD4[i4,])
  Data_CD8[i8,] <- ifelse(is.na(Data_CD8[i8,]), Data_CD8[j8,], Data_CD8[i8,])
  
  #Copy the new vlaues in the right palace
  Data_CD4[j4,] <- Data_CD4[i4,]
  Data_CD8[j8,] <- Data_CD8[i8,]
  
  #Change Donor ID
  Data_CD4[j4,"Donor"] <- j
  Data_CD8[j8,"Donor"] <- j
}

#Remove row with donor ID "Dxxx_N"
Data_CD4 <- Data_CD4 %>%
  filter(!grepl("_N$", Donor))

Data_CD8 <- Data_CD8 %>%
  filter(!grepl("_N$", Donor))

###########################################################
###########################################################

#Merge full dataset for Input 
Data_T <- merge(Data_CD4[,c("Donor",Predictors_4)], Data_CD8[,c("Donor",Predictors_8)], by="Donor", all=TRUE)
Data_cell <- merge(Data_T, Data_B[,c("donor",Predictors_B)], by.x="Donor", by.y="donor", all=TRUE)

#Exclude donor with no neutralization titers: D554
#SARS_COV2 <- SARS_COV2[-23,]

#replace censored values with something to work with
#SARS_COV2[SARS_COV2 == "<5"] <- "0" #Replaced already in the correction step
Data_cell <- Data_cell %>% mutate_all(str_replace_all, "\\S+\\*$", "NA")
##################################################################
###!!!! Quality Check !!!
#Data_cell <- Data_cell %>% mutate_all(str_replace_all, "t", "NA")
##################################################################
Data_cell <- Data_cell %>% mutate_all(str_replace_all, "t", "")

Data_cell[is.na(Data_cell)] <- "NA"
Data_cell[Data_cell == "NA"] <- NA
Data_cell[Data_cell == "<0.001"] <- "0.001"
#New Threshold
#Data_cell[Data_cell == "<0.001"] <- "0.0025"

#Merge Demographics, Antibodies and cell frequencies
Data2 <- merge(Data1, Data_cell, by.x="donor", by.y="Donor", all=TRUE)

#Ensure numeric valus are numeric
var_num <- c(outcome, Predictors_4, Predictors_8, Predictors_B) 
Data2[ , var_num] <- apply(Data2[ , var_num], 2, function(x) as.numeric(as.character(x)))

Data2$donor <- as.factor(Data2$donor)

#Save dataset
Dataset_DARPA <- Data2

save(Dataset_DARPA, file = "Data/Dataset_DARPA.RData")
write.csv(Dataset_DARPA, "Data/Dataset_DARPA.csv", row.names = FALSE)




##########################################################################

#Compare new with old dataset "SARS_COV2"
# Reorder df2 to match df1
df2 <- SARS_COV2[
  rownames(Data2),        # row order
  colnames(Data2)         # column order
]


which(Data2[,13:143] != df2[,13:143], arr.ind = TRUE)
tol <- 1e-13
diff_idx <- which(abs(Data2[,13:143] - df2[,13:143]) > tol, arr.ind = TRUE)





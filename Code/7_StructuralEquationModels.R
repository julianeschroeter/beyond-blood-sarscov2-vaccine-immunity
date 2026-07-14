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

long$s_igg_log <- as.numeric(scale(log1p(long$s_igg)))
long$rbd_igg_log <- as.numeric(scale(log1p(long$rbd_igg)))
long$wa1_log <- as.numeric(scale(log1p(long$wa1)))
long$delta_log <- as.numeric(scale(log1p(long$delta)))
long$omicron_ba2.12.1_log <- as.numeric(scale(log1p(long$omicron_ba2.12.1)))

long$sex_fac <- as.numeric(as.factor(long$sex))-1
long$age_std <-  as.numeric(scale(long$age))
long$pandemic_std <-  as.numeric(scale(long$pandemic))
long$infected_fac <- as.numeric(as.factor(long$infected))-1

long$spl_smbc_log <- as.numeric(scale(ifelse(long$spl_smbc == 0, log10(0.05), log10(long$spl_smbc))))

long$lng_scd4_cd49a_log <- as.numeric(scale(ifelse(long$lng_scd4_cd49a == 0, log10(0.05), log10(long$lng_scd4_cd49a))))

long$lln_smbc_cd69_log <- as.numeric(scale(ifelse(long$lln_smbc_cd69 == 0, log10(0.05), log10(long$lln_smbc_cd69))))
long$lln_smbc_igg_log <- as.numeric(scale(ifelse(long$lln_smbc_igg == 0, log10(0.05), log10(long$lln_smbc_igg))))
long$lln_smbc_igm_log <- as.numeric(scale(ifelse(long$lln_smbc_igm == 0, log10(0.05), log10(long$lln_smbc_igm))))
long$lln_scd8_log <- as.numeric(scale(ifelse(long$lln_scd8 == 0, log10(0.05), log10(long$lln_scd8))))

# Convert back to Mids
imp_l <- as.mids(long)

#Model 1A
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log + wa1_log + delta_log + omicron_ba2.12.1_log

  antibodies ~ pandemic_std

  antibodies ~~ 1*antibodies
'
fit1A <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit1A)
fitMeasures(fit1A)


#Model 1B
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log

  antibodies ~ infected_fac

  antibodies ~~ 1*antibodies
'
fit1B <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit1B)
fitMeasures(fit1B)


#Model 1C
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_scd8_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log

  antibodies ~ cell

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
'
fit1C <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit1C)
fitMeasures(fit1C)


#Model 2
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_scd8_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log

  antibodies ~ pandemic_std + infected_fac + cell

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
'
fit2 <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit2)
fitMeasures(fit2)


#Model 2A
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_scd8_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log

  infected_fac ~ pandemic_std
  antibodies ~ pandemic_std + infected_fac + cell

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
'
fit2A <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit2A)
fitMeasures(fit2A)


#Model 2B
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_scd8_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log

  infected_fac ~ pandemic_std
  cell ~ infected_fac
  antibodies ~  infected_fac + cell

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
'
fit2B <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit2B)
fitMeasures(fit2B)


#Model 2C
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_scd8_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log

  infected_fac ~ pandemic_std
  cell ~ infected_fac
  antibodies ~ pandemic_std + infected_fac + cell

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
'
fit2C <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit2C)
fitMeasures(fit2C)


#Model 2D
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_scd8_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log

  infected_fac ~ pandemic_std
  cell ~ infected_fac + pandemic_std
  antibodies ~ pandemic_std + infected_fac + cell

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
'
fit2D <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit2D)
fitMeasures(fit2D)





#Model 3
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_scd8_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log

  infected_fac ~ pandemic_std
  cell ~ infected_fac
  antibodies ~ cell

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
'
fit3 <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit3)
fitMeasures(fit3)


#Model 3A
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_scd8_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log

  infected_fac ~ pandemic_std
  cell ~ infected_fac + pandemic_std
  antibodies ~ cell

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
'
fit3A <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit3A)
fitMeasures(fit3A)



#Model 3B
model <- '
antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_scd8_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log 

  infected_fac ~ pandemic_std
  cell ~ infected_fac
  antibodies ~ cell
 
  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
'
fit3B <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit3B)
fitMeasures(fit3B)



#Model 3C
model <- '
antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_scd8_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log

  infected_fac ~ pandemic_std
  cell ~ infected_fac
  antibodies ~ cell
 
  s_igg_log ~ pandemic_std
  rbd_igg_log ~ pandemic_std
  wa1_log ~ pandemic_std
  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
'
fit3C <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit3C)
fitMeasures(fit3C)



#Model 4
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  Bcell =~ NA*spl_smbc_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log + lng_scd4_cd49a_log


  antibodies ~ Bcell
  
  Bcell ~~ lln_scd8_log 
  antibodies ~~ lln_scd8_log 

  antibodies ~~ 1*antibodies
  Bcell ~~ 1*Bcell
'
fit4 <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit4)
fitMeasures(fit4)



#Model 4A
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log

  antibodies ~ lln_scd8_log 

  antibodies ~~ 1*antibodies
'
fit4A <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit4A)
fitMeasures(fit4A)


#Model 4B
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log

  antibodies ~~ lln_scd8_log 

  antibodies ~~ 1*antibodies
'
fit4B <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit4B)
fitMeasures(fit4B)

#Model 4C
model <- '
  Bcell =~ NA*spl_smbc_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log + lng_scd4_cd49a_log

  lln_scd8_log ~ infected_fac + pandemic_std + Bcell

  Bcell ~~ 1*Bcell
'
fit4C <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit4C)
fitMeasures(fit4C)


#Model 5
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  Bcell =~ NA*spl_smbc_log + lln_smbc_igg_log + lln_smbc_igm_log + lng_scd4_cd49a_log + lln_smbc_cd69_log
  
  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std 

  infected_fac ~ pandemic_std
  Bcell ~ infected_fac
  lln_scd8_log ~ pandemic_std
  antibodies ~ Bcell
  
  antibodies ~~ 1*antibodies
  Bcell ~~ 1*Bcell
'
fit5 <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit5)
fitMeasures(fit5)



#Model 5.1
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  Bcell =~ NA*spl_smbc_log + lln_smbc_igg_log + lln_smbc_igm_log + lng_scd4_cd49a_log + lln_smbc_cd69_log
  
  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std 

  infected_fac ~ pandemic_std
  Bcell ~ infected_fac
  lln_scd8_log ~ infected_fac + pandemic_std
  antibodies ~ Bcell
  
  antibodies ~~ 1*antibodies
  Bcell ~~ 1*Bcell
'
fit5.1 <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit5.1)
fitMeasures(fit5.1)



#################################################################33


#Model 5A
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  Bcell =~ NA*spl_smbc_log + lln_smbc_igg_log + lln_smbc_igm_log 
  CD8 =~ NA*spl_scd8_log + lln_scd8_log
   
  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std

  infected_fac ~ pandemic_std
  Bcell ~ infected_fac
  antibodies ~ Bcell
  
  CD8 ~ antibodies + Bcell

  antibodies ~~ 1*antibodies
  Bcell ~~ 1*Bcell
  CD8 ~~ 1*CD8
'
fit5A <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit5A)
fitMeasures(fit5A)





#Model 6
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  Bcell =~ NA*spl_smbc_log + lln_smbc_igg_log + lln_smbc_igm_log
  CD8 =~ NA*spl_scd8_log + lln_scd8_log 
  
  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std

  infected_fac ~ pandemic_std
  Bcell ~ infected_fac
  antibodies ~ Bcell + CD8

  antibodies ~~ 1*antibodies
  Bcell ~~ 1*Bcell
  CD8 ~~ 1*CD8
'
fit6 <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit6)
fitMeasures(fit6)








#Model 6
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log + lln_scd4_naive_log + lln_scd4_cd49a_log
  CD8 =~ NA*spl_scd8_log + lln_scd8_log
  
  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std

  infected_fac ~ pandemic_std
  cell ~ infected_fac
  antibodies ~  cell
  
  CD8 ~ cell

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
  CD8 ~~ 1*CD8
'
fit6 <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit6)
fitMeasures(fit6)



#Model 4A
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + spl_scd8_log + lng_scd4_cd49a_log +  lln_scd8_log + lln_smbc_cd69_log + lln_smbc_igg_log +lln_smbc_igm_log + lln_scd4_naive_log + lln_scd4_cd49a_log

  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std

  infected_fac ~ pandemic_std
  cell ~ infected_fac
  antibodies ~ cell

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
'
fit4A <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit4A)
fitMeasures(fit4A)


### Singular dataset
dat1 <- mice::complete(imp_l, 1)
fit_single <- sem(model, data = dat1)
modindices(fit_single, sort = TRUE, maximum.number = 15)



#Model 4B
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log + lln_scd4_naive_log + lln_scd4_cd49a_log
  CD8 =~ NA*spl_scd8_log + lln_scd8_log
  
  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std

  infected_fac ~ pandemic_std
  cell ~ infected_fac
  antibodies ~ cell
  
  CD8 ~~ antibodies

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
  CD8 ~~ 1*CD8
'
fit4B <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit4B)
fitMeasures(fit4B)



#Model 4C
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  Bcell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log + lln_scd4_naive_log + lln_scd4_cd49a_log
  CD8 =~ NA*spl_scd8_log + lln_scd8_log + lng_scd4_cd49a_log + lln_scd4_naive_log + lln_scd4_cd49a_log
  
  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std

  infected_fac ~ pandemic_std
  Bcell ~ infected_fac
  antibodies ~ Bcell
  
  CD8 ~ pandemic_std

  antibodies ~~ 1*antibodies
  Bcell ~~ 1*Bcell
  CD8 ~~ 1*CD8
'
fit4C <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit4C)
fitMeasures(fit4C)




#Model 4
model <- '
  cell =~ NA*spl_smbc_log + spl_scd8_log + lng_scd4_cd49a_log +  lln_scd8_log + lln_smbc_cd69_log + lln_smbc_igg_log +lln_smbc_igm_log + lln_scd4_naive_log + lln_scd4_cd49a_log

  infected_fac ~ pandemic_std
  cell ~ infected_fac
  s_igg_log ~ cell
  rbd_igg_log ~ cell
  wa1_log ~ cell 
  delta_log ~ cell + pandemic_std
  omicron_ba2.12.1_log ~ cell + pandemic_std

  cell ~~ 1*cell
'
fit4 <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit4)
fitMeasures(fit4)


#Model 5
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log + lln_scd4_naive_log + lln_scd4_cd49a_log
  CD8 =~ NA*spl_scd8_log + lln_scd8_log
  
  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std

  infected_fac ~ pandemic_std
  cell ~ infected_fac
  antibodies ~ cell
  CD8 ~ pandemic_std + infected_fac + cell

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
  CD8 ~~ 1*CD8
'
fit5 <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit5)
fitMeasures(fit5)



#Model 5A
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log + lln_scd4_naive_log + lln_scd4_cd49a_log
  CD8 =~ NA*spl_scd8_log + lln_scd8_log
  
  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std

  infected_fac ~ pandemic_std
  cell ~ infected_fac
  antibodies ~ cell
  CD8 ~ pandemic_std

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
  CD8 ~~ 1*CD8
'
fit5A <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit5A)
fitMeasures(fit5A)

### Singular dataset
dat1 <- mice::complete(imp_l, 1)
fit_single <- sem(model, data = dat1)
modindices(fit_single, sort = TRUE, maximum.number = 5)


#Model 5B
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  cell =~ NA*spl_smbc_log + lng_scd4_cd49a_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log + lln_scd4_naive_log + lln_scd4_cd49a_log
  CD8 =~ NA*spl_scd8_log + lln_scd8_log
  
  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std

  infected_fac ~ pandemic_std
  cell ~ infected_fac
  antibodies ~ cell
  CD8 ~ pandemic_std

  antibodies ~~ 1*antibodies
  cell ~~ 1*cell
  CD8 ~~ 1*CD8
'
fit5B <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit5B)
fitMeasures(fit5B)


#Model 6
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  Bcell =~ NA*spl_smbc_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log
  CD4 =~ NA*lng_scd4_cd49a_log + lln_scd4_naive_log + lln_scd4_cd49a_log
  CD8 =~ NA*spl_scd8_log + lln_scd8_log
  
  delta_log ~ pandemic_std
  omicron_ba2.12.1_log ~ pandemic_std

  infected_fac ~ pandemic_std
  CD4 ~ infected_fac
  Bcell ~ CD4
  antibodies ~ Bcell
  CD8 ~ pandemic_std + CD4

  antibodies ~~ 1*antibodies
  Bcell ~~ 1*Bcell
  CD4 ~~ 1*CD4
  CD8 ~~ 1*CD8
'
fit6 <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit6)
fitMeasures(fit6)



#Model 7A
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  CD8 =~ NA*spl_scd8_log + lln_scd8_log
  
  antibodies ~ CD8

  antibodies ~~ 1*antibodies
  CD8 ~~ 1*CD8
'
fit7A <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit7A)
fitMeasures(fit7A)



#Model 7AA
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  CD8 =~ NA*spl_scd8_log + lln_scd8_log
  
  antibodies ~~ CD8

  antibodies ~~ 1*antibodies
  CD8 ~~ 1*CD8
'
fit7AA <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit7AA)
fitMeasures(fit7AA)


#Model 7B
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  CD4 =~ NA*lng_scd4_cd49a_log + lln_scd4_naive_log + lln_scd4_cd49a_log
  
  antibodies ~ CD4

  antibodies ~~ 1*antibodies
  CD4 ~~ 1*CD4
'
fit7B <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit7B)
fitMeasures(fit7B)


#Model 7BB
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  CD4 =~ NA*lng_scd4_cd49a_log + lln_scd4_naive_log + lln_scd4_cd49a_log
  
  antibodies ~~ CD4

  antibodies ~~ 1*antibodies
  CD4 ~~ 1*CD4
'
fit7BB <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit7BB)
fitMeasures(fit7BB)



#Model 7C
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  Bcell =~  NA*spl_smbc_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log
  
  antibodies ~ Bcell

  antibodies ~~ 1*antibodies
  Bcell ~~ 1*Bcell
'
fit7C <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit7C)
fitMeasures(fit7C)



#Model 7CC
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  Bcell =~  NA*spl_smbc_log + lln_smbc_cd69_log + lln_smbc_igg_log + lln_smbc_igm_log + lng_scd4_cd49a_log + lln_scd4_naive_log + lln_scd4_cd49a_log
  CD8 =~ NA*spl_scd8_log + lln_scd8_log
  
  antibodies ~ Bcell
  antibodies ~~ CD8

  antibodies ~~ 1*antibodies
  CD8 ~~ 1*CD8
  Bcell ~~ 1*Bcell
'
fit7CC <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit7CC)
fitMeasures(fit7CC)




#Model 7D
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  LNG =~  NA*lng_scd4_cd49a_log
  
  antibodies ~ LNG

  antibodies ~~ 1*antibodies
  LNG ~~ 1*LNG
'
fit7D <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit7D)
fitMeasures(fit7D)


#Model 7E
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  LLN =~  NA*lln_smbc_igg_log + lln_smbc_igm_log + lln_scd8_log 
  
  antibodies ~ LLN

  antibodies ~~ 1*antibodies
  LLN ~~ 1*LLN
'
fit7E <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit7E)
fitMeasures(fit7E)



#Model 7F
model <- '
  antibodies =~ NA*s_igg_log + rbd_igg_log +  wa1_log + delta_log + omicron_ba2.12.1_log
  SPL =~  NA*spl_scd8_log + spl_smbc_log
  
  antibodies ~ SPL

  antibodies ~~ 1*antibodies
  SPL ~~ 1*SPL
'
fit7F <- lavaan.mi::sem.mi(model, data = imp_l)
summary(fit7F)
fitMeasures(fit7F)









### Singular dataset
dat1 <- mice::complete(imp_l, 1)
fit_single <- sem(model, data = dat1)
semPaths(fit_single, what = "std", layout = "tree", edge.label.cex = 0.8)             


summary(fit, fit.measures = TRUE, standardized = TRUE)

modindices(fit_single, sort = TRUE, maximum.number = 5)


summary(fit_single)

#data <- cbind(data, model.matrix(~ gender - 1, data)

fit1 <- parameterEstimates(fit)
params_sig <- subset(fit1, pvalue < 0.05 & op == "~")
semPlotObj <- semPlotModel(fit1)
semPaths(fit1, whatLabels = "est", layout = "tree", edge.label.cex = 0.9)

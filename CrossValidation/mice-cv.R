##############################################################
# script for cross-validation of the MICE imputation method. #
# allow running on cluster by taking a command line argument #
# equal to the index on the list of known values             #
# This script is for developing and "smoke testing"          #
##############################################################

## NOTE: D631 bm_scd4 is not censored, but eq to 0.001 
## LOD is now set to 0.0025 for parent variables


# check R version

if ( getRversion() != "4.5.2" ) {
    stop("incorrect R version")
}


library(mice)
library(tidyverse)
#devtools::install_github(repo = "amices/mice")

NUM_CORES <- 2
NUM_IMPS <- 2 ## TESTING
NUM_TREES <- 10
LOG_TRANS <- FALSE
BALANCE <- TRUE
CO_MASK_SUBSETS <- TRUE


find_co_mask <- function(variable) {
    # return variable names that should be co-masked
    tissues <- c("spl", "bld", "lln", "mln", "iln", "lng")
    major_cell_types <- c("scd4", "scd8", "stfh", "streg", "smbc")
    tcell_subsets <- c("tcm", "tem", "temra", "naive")
    bcell_subsets <- c("igm", "igg", "iga")
    variable_parts <- strsplit(variable, "_")[[1]]
    if ( length(variable_parts) != 3 ) return(variable)
    tissue <- variable_parts[1]
    if ( !(tissue %in% tissues) ) return(variable)
    cell_type <- variable_parts[2]
    if ( cell_type %in% c("scd4", "scd8") ) {
        subset <- variable_parts[3]
        if ( !(subset %in% tcell_subsets) ) return(variable)
        co_mask_vars <- unlist(lapply(tcell_subsets, function(x) paste(tissue, cell_type, x, sep="_")))
        return(co_mask_vars)
    } else if ( cell_type == "smbc" ) {
        subset <- variable_parts[3]
        if ( !(subset %in% bcell_subsets) ) return(variable)
        co_mask_vars <- unlist(lapply(bcell_subsets, function(x) paste(tissue, cell_type, x, sep="_")))
        return(co_mask_vars)
    } else {
        return(variable)
    }
}



find_balancing_loo <- function(data, masked_donor, masked_variable) {
    # find donor with observed data that balances the left out value for given donor
    
    data_var <- data |> select(all_of(c("donor", masked_variable))) |>
        filter(!is.na(.data[[masked_variable]]))

    # sort data by value of masked_variable 
    
    sorted_var <- data_var |> arrange(.data[[masked_variable]])
    
    # find index of the masked_donor
    
    idx_masked <- which(sorted_var$donor == masked_donor)
    num_obs <- nrow(sorted_var)
    
    balancing_idx <- num_obs - idx_masked + 1
    
    # return balancing donor 
    
    balancing_donor <- sorted_var$donor[balancing_idx]
    
    return(balancing_donor)
}

######## functions ###########

run_mice <- function(data, num_cores, num_imps, masked_donor=NULL, masked_variable=NULL, 
                     ignored_vars=NULL, balance=FALSE, co_mask_subsets=FALSE) {
    # copy data to avoid modifying original
    data <- data %>% as_tibble()
    
    # mask a value for CV purposes
    if (!is.null(masked_variable) && !is.null(masked_donor)) {
        if ( balance ) {
            balancing_donor <- find_balancing_loo(data, masked_donor, masked_variable)
            cat("Balancing donor for", masked_donor, "is", as.character(balancing_donor), "\n")
            data[data$donor==balancing_donor, masked_variable] <- NA
        }
        
        if ( co_mask_subsets ) {
            mask_vars <- find_co_mask(masked_variable)
            cat("Co-masking variables:", paste(mask_vars, collapse=", "), "\n")
            for ( var in mask_vars ) {
                data[data$donor==masked_donor, var] <- NA
            }
        }
        
        data[data$donor==masked_donor, masked_variable] <- NA
    }
    
    ## make a predictor matrix: exclude ignored variables (don't use them for imputations)
    pred <- quickpred(data, method = "spearman", mincor= .3, minpuc = 0.25, exclude = ignored_vars)
    
    run_parallel <- num_cores > 1
    mp <- if ( run_parallel ) {
        futuremice(data, parallelseed = 144169, method="rf",
                   m=num_imps, pred = pred, maxit=5, remove.collinear = FALSE, n.core=num_cores, ntree=NUM_TREES)
    } else {
        mice(data, printFlag=FALSE, method="rf", m=num_imps, seed = 144169,
             pred = pred, maxit=5, remove.collinear = FALSE, ntree=NUM_TREES)
    }
    
    return(mp)
}


loo_cv_mice <- function(data, num_cores, num_imps, masked_donor, masked_variable, ignored_vars=NULL, balance=FALSE, 
                        co_mask_subsets=FALSE) {
    # run mice with masked variable and donor
    mp <- run_mice(data, num_cores, num_imps, 
                   masked_donor=masked_donor, masked_variable=masked_variable, 
                   ignored_vars=ignored_vars, balance=balance, co_mask_subsets=co_mask_subsets)
    
    # extract imputed values for masked variable and donor
    
    imp_vals <- mice::complete(mp, "long") |>
        filter(donor == masked_donor) |>
        select(all_of(masked_variable))
    return(imp_vals)
}


transform_data <- function(data) {
    varnames <- colnames(data)
    
    ## different observations have to be transformed differently
    meta_vars <- c(
        "donor", "age", "sex", "race_ethnicity", "deathyear", "deathmonth", 
        "tip", "tpv", "dose_number", "vaccine_brand", "vaccinated", "infected",
        "batch"
    )
    
    titer_vars <- c("s_igg", "rbd_igg", "wa1", "delta", "omicron_ba2.12.1")
    
    cell_vars <- varnames[!varnames %in% c(meta_vars, titer_vars)]
    
    is_parent <- function(varname) length(strsplit(varname, "_")[[1]]) == 2
    
    parent_vars <- cell_vars[sapply(cell_vars, is_parent)]
    subset_vars <- cell_vars[!sapply(cell_vars, is_parent)]
    
    parent_LOD <- 0.0025
    subset_LOD <- 0.05
    
    parent_trans <- function(x) {
        if (is.na(x)) return(NA)
        if (x == 0) return(log10(parent_LOD))
        return(log10(x))
    }
    
    subset_trans <- function(x) {
        if (is.na(x)) return(NA)
        if (x == 0) return(log10(subset_LOD))
        return(log10(x))
    }
    
    titer_trans <- function(x) {
        if (is.na(x)) return(NA)
        return(log10(1+x))
    }
    
    
    trans_data <- data |> 
        mutate(across(all_of(parent_vars), Vectorize(parent_trans))) |>
        mutate(across(all_of(subset_vars), Vectorize(subset_trans))) |>
        mutate(across(all_of(titer_vars), Vectorize(titer_trans)))
    
    return(trans_data)
}


######### main script #########

# import RData File

#load("data/feb26/Dataset_DARPA.RData")
load("data/apr26/Dataset_DARPA.RData")

data <- as_tibble(Dataset_DARPA)

# add time into pandemic, remove year and month

ref_date <- as.Date("2019-12-01", format="%Y-%m-%d")
death_date <- as.Date(paste(data$deathyear, data$deathmonth,1, sep='-'), format='%Y-%B-%d')
days_post_ref <- as.numeric(death_date - ref_date)
data <- add_column(data, tip=days_post_ref, .before="tpv")

# transform data to log10 scale
if ( LOG_TRANS ) {
    cat("log-transforming data")
    data <- transform_data(data)
} else {
    cat("not log-transforming data")
}

ignored_vars <- c("deathyear", "deathmonth", "donor")

# create a data frame with all known entries: the donor, the variable name and the ground truth value. This will be used for masking and testing.
# only include numerical values. 

numeric_data <- data |> select(-all_of(ignored_vars)) |> select(where(is.numeric))
# add donor ID as an index
numeric_data <- numeric_data %>% add_column(donor = data$donor, .before=1)
# create a long format data frame with donor, variable, and value
long_data <- numeric_data %>%
    pivot_longer(cols = -donor, names_to = "variable", values_to = "ground_truth") %>%
    filter(!is.na(ground_truth))

cat("Total known values for CV:", nrow(long_data), "\n")

# save long_data as a rds file for later use in analysis
if ( LOG_TRANS ) {
    saveRDS(long_data, "data/long_data_trans.rds")
} else {
    saveRDS(long_data, "data/long_data.rds")
}


################ how many predictors do we want per node? #################


pred <- quickpred(data, method = "spearman", mincor= .3, minpuc = 0.25, exclude = ignored_vars)

heatmap(pred)


# sum the rows of the predictor matrix to get the number of predictors used for each variable

num_pred <- rowSums(pred)

hist(num_pred, breaks=20, main="Distribution of number of predictors per variable", xlab="Number of predictors")

which(num_pred == 0)

# check that complete obs have no predictors 

is.complete <- data |> select(-all_of(ignored_vars)) |>
    select(where(is.numeric)) |>
    select(where(~!anyNA(.))) |> names()

if ( all(is.complete %in% colnames(pred)[num_pred == 0]) ) {
    cat("Complete observations have no predictors, as expected\n")
} else {
    cat("Warning: some complete observations have predictors\n")
}



masked_donor <- "D511"
masked_variable <- "spl_scd4_temra"

imp <- run_mice(data, num_cores=NUM_CORES, num_imps=NUM_IMPS, masked_donor, 
                masked_variable, ignored_vars, balance=BALANCE, 
                co_mask_subsets=CO_MASK_SUBSETS)

imp_vals <- loo_cv_mice(data, num_cores=NUM_CORES, num_imps=NUM_IMPS, 
                        masked_donor, masked_variable, ignored_vars, 
                        balance=BALANCE, co_mask_subsets=CO_MASK_SUBSETS)

imp_vals



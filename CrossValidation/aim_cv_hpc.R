##############################################################
# script for cross-validation of the MICE imputation method. #
# allow running on cluster by taking a command line argument #
# equal to the index on the list of known values             #
#               file version for use on HPC                  #
##############################################################

# check R version

if ( getRversion() != "4.5.2" ) {
   stop("incorrect R version")
}


library(mice)
library(tidyverse)
library(ranger)

NUM_CORES <- 20
NUM_IMPS <- 20 ## TESTING
NUM_TREES <- 50 ## TESTING
LOG_TRANS <- FALSE
BALANCE <- TRUE
CO_MASK_SUBSETS <- TRUE


######## functions ###########


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



run_mice <- function(data, num_cores, num_imps, masked_donor=NULL, masked_variable=NULL, 
                     ignored_vars=NULL, balance=FALSE, co_mask_subsets=FALSE) {
    # copy data to avoid modifying original
    data <- data %>% as_tibble()
    
    # mask a value for CV purposes
    if (!is.null(masked_variable) && !is.null(masked_donor)) {
        if ( balance ) {
            balancing_donor <- find_balancing_loo(data, masked_donor, masked_variable)
            cat("Balancing donor for", as.character(masked_donor), "is", as.character(balancing_donor), "\n")
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
                   m=num_imps, pred = pred, maxit=5, remove.collinear = FALSE, 
                   n.core=num_cores, ntree=NUM_TREES)
    } else {
        mice(data, printFlag=FALSE, method="rf", m=num_imps, seed = 144169,
             pred = pred, maxit=5, remove.collinear = FALSE, ntree=NUM_TREES)
    }
    
    return(mp)
}


loo_cv_mice <- function(data, num_cores, num_imps, masked_donor, masked_variable, ignored_vars=NULL,
                        balance=FALSE, co_mask_subsets=FALSE) {
    # run mice with masked variable and donor
    mp <- run_mice(
        data, num_cores, num_imps, 
        masked_donor=masked_donor, masked_variable=masked_variable,
        ignored_vars=ignored_vars, balance=balance, co_mask_subsets=co_mask_subsets
    )
    
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

num_vals <- nrow(long_data)
cat("Total known values for CV:", num_vals, "\n")

task_index <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
num_tasks <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

if ( is.na(task_index) || is.na(num_tasks) ) {
    stop("must provide task index and number of tasks")
}

all_mask_idxs <- 1:num_vals
chunks <- split(all_mask_idxs, cut(seq_along(all_mask_idxs), num_tasks, labels = FALSE))
mask_idxs <- chunks[[task_index]]

cat("mask indices for this task", mask_idxs, "\n")

imp_val_list <- list()

for ( mask_index in mask_idxs ) {
    masked_donor <- long_data$donor[mask_index]
    masked_variable <- long_data$variable[mask_index]
    
    cat("Masked donor:", as.character(masked_donor), "\n")
    cat("Masked variable:", masked_variable, "\n")

    imp_vals <- loo_cv_mice(data, num_cores=NUM_CORES, num_imps=NUM_IMPS, 
                            masked_donor, masked_variable, ignored_vars,
                            balance=BALANCE, co_mask_subsets=CO_MASK_SUBSETS)
    
    # add masked donor and variable to imp_vals for later analysis
    df <- imp_vals |> add_column(variable = masked_variable, .before=1) |>
        add_column(donor = masked_donor, .before=1) |>
        rename(imputed_value = 3) |>
        # add ground truth
        left_join(long_data, by = c("donor" = "donor", "variable" = "variable"))
    
    
    imp_val_list[[length(imp_val_list)+1]] <- df
}

#### combine the list of data frames 

imp_val_df <- bind_rows(imp_val_list)


######## write imputed values to file for later analysis

# WARNING: added "balance" to filename: testing balanced imputation
# WARNING: added "co-mask" to filename: testing co-masking of subsets
output_file <- paste0("results/imputed_values_apr26_balance_co-mask_task_", task_index, ".csv")
write_csv(imp_val_df, output_file)
cat("Imputed values written to:", output_file, "\n")








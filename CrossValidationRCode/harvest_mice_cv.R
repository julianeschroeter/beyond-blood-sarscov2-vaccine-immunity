library(mice)
library(ggplot2)
library(tidyverse)

LOG_TRANS <- FALSE
BALANCED <- TRUE
CO_MASK_SUBSETS <- TRUE

version_id <- "apr26"

# import dataset 

load("Data/Dataset_DARPA.RData")
#load("data/apr26/Dataset_DARPA.RData")
data <- as_tibble(Dataset_DARPA)

# add time into pandemic, remove year and month

ref_date <- as.Date("2019-12-01", format="%Y-%m-%d")
death_date <- as.Date(paste(data$deathyear, data$deathmonth,1, sep='-'), format='%Y-%B-%d')
days_post_ref <- as.numeric(death_date - ref_date)
data <- data |> add_column(tip=days_post_ref, .after="age")


ignored_vars <- c("deathyear", "deathmonth", "donor")


# count the number of donors for which the value of a variable is observed
num_observed_vals <- data |> select(-all_of(ignored_vars)) |> 
    select(where(is.numeric)) |> 
    summarise_all(~sum(!is.na(.))) |> 
    pivot_longer(everything(), names_to="variable", values_to="num_observed")

# import long-format ground-truth values 

if ( LOG_TRANS ) {
    long_data <- readRDS("Data/long_data_trans.rds")
} else {
    long_data <- readRDS("Data/long_data.rds")
}

# import result files WARNING: replace directory name!

result_dir <- "Data/results/"
namebase = paste("imputed_values", version_id, sep="_")

if ( !CO_MASK_SUBSETS ) {
    if ( LOG_TRANS && !BALANCED ) {
        cat("importing log-transformed and unbalanced results\n")
        result_files <- list.files(result_dir, pattern = paste0(namebase, "_trans_task_.*\\.csv"), full.names = TRUE)
    } else if ( LOG_TRANS && BALENCED ) {
        stop("log-transfomed and balanced not implemented yet")
    } else if ( !LOG_TRANS && !BALANCED ) {
        cat("importing untransformed and unbalanced results\n")
        result_files <- list.files(result_dir, pattern = paste0(namebase, "_task_.*\\.csv"), full.names = TRUE)
    } else if ( !LOG_TRANS && BALANCED ) {
        cat("importing untransformed and balanced results\n")
        result_files <- list.files(result_dir, pattern = paste0(namebase, "_balance_task_.*\\.csv"), full.names = TRUE)
    } else {
        stop("invalid combination")
    }
} else {
    if ( LOG_TRANS || !BALANCED ) {
        stop("co-masking subsets not implemented for log-transformed or unbalanced data")
    } else {
        cat("importing co-masked subsets results\n")
        result_files <- list.files(result_dir, pattern = paste0(namebase, "_balance_co-mask_task_.*\\.csv"), full.names = TRUE)
    }
}

# read all result files into a single data frame

sel_result_files <- result_files

imputed_values <- lapply(sel_result_files, read_csv) |> bind_rows()

## make sure the order of variables is the same as in the data set (the columns)
## then sort by donor (again using the order in the donor column of data)

imputed_values <- imputed_values |> 
    mutate(variable = factor(variable, levels = colnames(data))) |> 
    mutate(donor = factor(donor, levels = unique(data$donor))) |>
    arrange(variable, donor)




## check that the number of variables is as expected 

# count the number of unique (donor, variable) pairs in the imputed_values data
num_imputed_pairs <- imputed_values |> select(donor, variable) |> distinct() |> nrow()
print(num_imputed_pairs)



# plot ground truth vs. imputation for the variable "age"
# make boxplots per donor

varname <- "lng_scd4_cd49a"

lod_subsets <- 0.025

if ( LOG_TRANS ) {
    safe_trans <- function(x) x
} else {
    safe_trans <- function(x, limit=1.0) {
        ifelse(x > lod_subsets, log10(x), log10(limit))
    }
}

sel_df <- imputed_values |> 
    filter(variable == varname) |> 
    mutate(ground_truth=safe_trans(ground_truth), imputed_value=safe_trans(imputed_value))



ggplot(sel_df, aes(x=ground_truth, y=imputed_value)) +
    geom_point() +
    geom_boxplot(aes(group=donor)) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    labs(title="CV", x=paste("Ground Truth", varname) , y=paste("Imputed", varname)) +
    theme_minimal()




# compute p-value for correlation between imputed and ground truth values for each variable
corrs <- imputed_values |> 
    group_by(variable) |> 
    summarise(
        correlation = cor(ground_truth, imputed_value, method="spearman"), 
        p_value = cor.test(ground_truth, imputed_value, method="spearman")$p.value
    )

# alternatively, take the mean/median of the imputed values for each (donor, variable) pair and then compute correlation with ground truth
aggr_func <- median
aggr_imputed <- imputed_values |>
    group_by(donor, variable) |>
    summarise(aggr_imputed_value = aggr_func(imputed_value), .groups = "drop") |>
    left_join(long_data, by=c("donor", "variable")) |>
    mutate(variable = factor(variable, levels = colnames(data))) |> 
    mutate(donor = factor(donor, levels = unique(data$donor))) |>
    arrange(variable, donor)

test_vars <- aggr_imputed[aggr_imputed$variable == "lng_scd4_cd49a", c("ground_truth", "aggr_imputed_value")]
test_vars <- aggr_imputed[aggr_imputed$variable == "lng_scd4_cd49a", c("donor", "ground_truth", "aggr_imputed_value")]

plot(test_vars)

print(test_vars, n=200)


cor(test_vars$ground_truth, test_vars$aggr_imputed_value, method="spearman")


cor_method <- "spearman"

corrs_aggr <- aggr_imputed |>
    group_by(variable) |>
    summarise(
        correlation = cor(ground_truth, aggr_imputed_value, method=cor_method), 
        p_value = cor.test(ground_truth, aggr_imputed_value, method=cor_method, continuity=TRUE)$p.value
    )

# replace NA with 0 for correlations 
corrs_aggr <- corrs_aggr |> mutate(correlation = ifelse(is.na(correlation), 0, correlation))


## compute the expected LOOCV bias: for each variable and donor, 
## compute the mean of the remaining values after leaving out a single donor

aggr_func <- median
loo_bias <- long_data |>
    group_by(variable, donor) |>
    summarise(naive_imputation = 
                  aggr_func(long_data$ground_truth[long_data$donor != cur_group()$donor & long_data$variable == cur_group()$variable]), .groups = "drop") |>
    left_join(long_data, by=c("donor", "variable")) |>
    mutate(variable = factor(variable, levels = colnames(data))) |>
    mutate(donor = factor(donor, levels = unique(data$donor))) |>
    arrange(variable, donor)


corrs_bias <- loo_bias |> 
    group_by(variable) |>
    summarise(
        correlation = cor(ground_truth, naive_imputation, method=cor_method), 
        p_value = cor.test(ground_truth, naive_imputation, method=cor_method, continuity=TRUE)$p.value
    )


signif_marks <- function(p) {
    if ( is.na(p) || is.nan(p) ) {
        return("ns")
    }
    if (p < 0.001) {
        return("***")
    } else if (p < 0.01) {
        return("**")
    } else if (p < 0.05) {
        return("*")
    } else {
        return("ns")
    }
}

get_pretty_varname <- function(varname) {
    varname <- as.character(varname)
    # map some values to more readable names for plotting
    pv <- switch(
        varname,
        age = "Age",
        tip = "Time into pandemic",
        tpv = "Time post vaccination",
        s_igg = "S-IgG",
        rbd_igg = "RBD-IgG",
        wa1 = "WA1",
        delta = "Delta",
        omicron_ba2.12.1 = "Omicron"
    )
    
    if ( !is.null(pv) ) return(pv)
    
    parts <- strsplit(varname, "_")[[1]]
    
    tissue <- toupper(parts[1])
    major_cell <- switch(
        parts[2],
        smbc = "mBc",
        scd8 = "CD8",
        scd4 = "CD4",
        stfh = "TFH",
        streg = "Treg",
    )
    if ( length(parts) == 2 ) {
        return(paste(tissue, major_cell))
    }

    subset_cell <- switch(
        parts[3],
        igm = "IgM",
        igg = "IgG",
        iga = "IgA",
        cd49a = "CD49a+",
        cd69 = "CD69+",
        cd103 = "CD103+",
        cxcr6 = "CXCR6+",
        tem = "TEM",
        tcm = "TCM",
        temra = "TEMRA",
        naive = "Naive",
    )
    return(paste(tissue, major_cell, subset_cell))
}


corrs <- corrs |> mutate(significance_class = sapply(p_value, signif_marks))
corrs_aggr <- corrs_aggr |> mutate(significance_class = sapply(p_value, signif_marks))
corrs_aggr <- corrs_aggr |> add_column(num_observed = num_observed_vals$num_observed, .before=1)

corrs_aggr <- corrs_aggr |> add_column(varname = sapply(corrs_aggr$variable, get_pretty_varname))
corrs_aggr <- corrs_aggr |> mutate(varname = factor(varname, levels=varname))

# add a plot group to the rows so that we can plot 3 different panels

corrs_aggr$idx <- 1:nrow(corrs_aggr)
num_vars_per_grp <- ceiling(nrow(corrs_aggr) / 3)
corrs_aggr <- corrs_aggr |> mutate(plot_group = idx %/% num_vars_per_grp)

corrs_bias <- corrs_bias |> mutate(significance_class = sapply(p_value, signif_marks))

palette <- c('lightcoral', 'red', '#67001F', 'lightgray')

# plot corr values and color by p-value
# add the number of observed values to the plot as text labels on the bars
p1 <- ggplot(corrs_aggr, aes(y=fct_rev(varname), x=correlation, fill=significance_class)) +
    geom_bar(stat="identity") +
    scale_x_continuous(breaks = c(-0.25, 0, 0.25, 0.5, 0.75), lim=c(-0.4, 0.85)) +
    geom_text(aes(y=fct_rev(varname), x=correlation, label=num_observed), inherit.aes=FALSE, size=3, hjust=0) +
    scale_fill_discrete(name="p-value", palette=palette) +
    labs(title="Correlation between Imputed and Ground Truth Values", y="Variable", x="Correlation") +
    theme_minimal() +
    facet_wrap(~plot_group, scales="free_y") +
    theme(strip.text.x = element_blank())


p1


ggsave(paste("cv_corrs", version_id, "balanced_co-mask.pdf", sep='_'), plot=p1, width=15, height=8)





# show example of mean imputation vs. ground truth

varname <- "lng_scd4_cd49a"

lod_subsets <- 0.025

safe_log <- function(x, limit=1.0) {
    ifelse(x > lod_subsets, log10(x), log10(limit))
}


sel_df <- aggr_imputed |> 
    filter(variable == varname) |> 
    mutate(ground_truth=safe_log(ground_truth, limit=0.001), imputed_value=safe_log(aggr_imputed_value, limit=0.001))



ggplot(sel_df, aes(x=ground_truth, y=imputed_value)) +
    geom_point() +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    labs(title="CV", x=paste("Ground Truth", varname) , y=paste("Imputed", varname)) +
    theme_minimal()



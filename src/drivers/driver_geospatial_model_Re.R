# Script Name: driver_geospatial_model_Re.R
# Purpose: Fits a model to effective reproduction estimates (Re) using weighted regression
# Author: Michael McPhail
# Date: 2024-07-09
# Dependencies: raster, Matrix, dplyr, rgdal, sp, sf, optparse, jsonlite 
# Usage: Rscript driver_geospatial_model_Re.R --config='PATH_TO_CONFIG'

set.seed(42)  #Currently removing randomness: good for validation, 
library(raster)
library(Matrix)
library(dplyr)
library(rgdal)
library(sp)
library(sf)
library(MASS)
library(glmmTMB)

library(optparse)
library(jsonlite)

options(warn=-1)  

#' Title: Load config
#' 
#' Description: This function loads the .json configuration file
#' 
#' @param config_path the path to the configuration file
#' @return A list of file paths.
load_config <- function(config_path) {
  config <- fromJSON(config_path)
  return(config)
}

option_list <- list(
  make_option(c("--config"), type="character", default=NULL, 
              help="Path to the configuration file", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$config)) {
  print_help(opt_parser)
  stop("The --config argument is required", call.=FALSE)
}

#---------------------------------------------------------------------------------------------------------------#
#             Start: Configuration
#---------------------------------------------------------------------------------------------------------------#
config <- load_config(opt$config)
#Path to input data
data_in_fn <- config$case_data_withReCovs_fn
#Output file paths
data_out_suff = config$Re_prediction_suffix
data_out_pf_prefix = config$Re_prediction_pf_prefix
data_out_pv_prefix = config$Re_prediction_pv_prefix
if(!dir.exists(dirname(data_out_pf_prefix))){
  dir.create(dirname(data_out_pf_prefix))
}
#Year range for prediction
Re_prediction_min_year = config$Re_prediction_min_year
Re_prediction_max_year = config$Re_prediction_max_year

#Country information for subsetting 
country_extent <- extent(config$country_extent)
country_iso = config$country_iso

#source utility functions
source('./src/utils/utilities_geospatial_modelling.R')
source('./config/data_lookup.R')

#---------------------------------------------------------------------------------------------------------------#
#             Start: Modelling 
#---------------------------------------------------------------------------------------------------------------#
parasite_list <- c('P.falciparum', 'P.vivax')
#Covariates to include with monthly lags
dyn_cov_name_list = c('CHIRPS', 'EVI','LST_day', 'LST_delta', 'LST_night','TCB', 'TCW', 'TSI_Pf', 'TSI_Pv')
#Covariates that never change
static_cov_name_list = c('ELEVATION', 'ACCESSIBILITY') #Easier to add forrest factor in the driver script # treat population as static (use 2020)
#Covariates that change annually (must include population for population weighting)
annual_cov_name_list = c('CHIRPS_synoptic', 'POPULATION')
use_forest = TRUE
full_base_cov_name_list <- c(dyn_cov_name_list, static_cov_name_list,annual_cov_name_list)
cov_names_full_list = get_cov_names(dyn_cov_name_list, n_lagged=1)

#Loop model fit and prediction over parasite_list
for(parasite in parasite_list){
  message(paste0('Starting ', parasite, ' iteration.'))
  #Stage 1: load covariates (load, wrangle into appropriate form)
  input_db <- read.csv(data_in_fn)
  # input_db$Re =  input_db$Re + 0.01
  cov_names = cov_names_full_list
  if(parasite == 'P.falciparum'){
    cov_names = cov_names[!grepl('TSI_Pv', cov_names)]
  }else{
    cov_names = cov_names[!grepl('TSI_Pf', cov_names)]
  }
  input_db = wrangle_input_data(input_db, cov_names, parasite)

  cov_names_to_normalise <-cov_names
  reduced_cov_names <- c(cov_names_to_normalise, 'FOREST')
  retained_cols <- c('Re', 'coordinate_index','LAT', 'LONG','inv_probability', reduced_cov_names)

  #Normalise covariates
  normalisation_data = normalise_input_covariates(input_db, cov_names_to_normalise)
  input_db = normalisation_data$input_db
  cov_scales_df = normalisation_data$cov_scales_df

  #Currently scales and truncates weights (to reduce outlier effect)
  input_db = get_regression_weights(input_db)

  #Feature selection loop
  retained_covariates <- reduced_cov_names
  retained_covariates_store <- retained_covariates[retained_covariates!='FOREST']
  retained_covariates_store<- retained_covariates_store[unname(sapply(retained_covariates_store, get_abs_corr, input_db=input_db))>0.05]

  message(paste0('Starting with ', toString(length(c(retained_covariates_store, 'FOREST'))), ' covariates'))
  retained_covariates = feature_select_greedy_wrapper(input_db, c(retained_covariates_store, 'FOREST'))
  message(paste0('Ending with ', toString(length(retained_covariates)), ' covariates'))

  model <- as.formula(paste("Re ~ ",paste(retained_covariates, collapse = " + ")))

  #Fitting final model
  fitted_model <- glmmTMB(model, 
                            ziformula = ~ ., 
                            family = ziGamma(link='log'), 
                            data = input_db,
                            weights=inv_probability)
  final_summary = summary(fitted_model)
  
  if(parasite=='P.falciparum'){
    write.csv(as.data.frame(final_summary$coefficients$cond), paste0(data_out_pf_prefix,'slopes_cond' , data_out_suff))
    write.csv(as.data.frame(final_summary$coefficients$zi), paste0(data_out_pf_prefix,'slopes_zi' , data_out_suff))
  }else{
    write.csv(as.data.frame(final_summary$coefficients$cond), paste0(data_out_pv_prefix, 'slopes_cond', data_out_suff))
    write.csv(as.data.frame(final_summary$coefficients$zi), paste0(data_out_pv_prefix, 'slopes_zi', data_out_suff))
  }
  # fitted_model <- glm.convert(fitted_model)

  #Make predictions for each year
  for(year_int in Re_prediction_min_year:Re_prediction_max_year){
    message(paste0('Predicting for year: ', year_int))
    valid_coord_df = predict_single_year(fitted_model, year_int, retained_covariates, cov_scales_df, country_extent, country_iso, max_cov_year = 2022)
    # valid_coord_df = add_prediction_intervals(fitted_model, valid_coord_df)
    if(parasite=='P.falciparum'){
      write.csv(valid_coord_df, paste0(data_out_pf_prefix, toString(year_int), data_out_suff))
    }else{
      write.csv(valid_coord_df, paste0(data_out_pv_prefix, toString(year_int), data_out_suff))
    }

  }
}

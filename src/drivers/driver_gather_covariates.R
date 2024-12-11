# Script Name: driver_gather_covariates.R
# Purpose: Appends covariate values to each row of data for subsequent regression
# Author: Michael McPhail 
# Date: 2024-07-09
# Dependencies: dplyr, ggplot2
# Usage: Rscript my_script.R input_file.csv output_file.csv

library(raster)
library(dplyr)
library(tidyr)
# library(rgdal)
library(sp)
library(terra)
library(sf)
library(optparse)
library(jsonlite)

#The following code takes a structured data base, builds a covariate data base for each admin unit, for each month, 
#for each year. Then it attatches the selected covariates (and for some, lagged covariates) to the data base.

#Generally, you only need to change the configuration part. Eventually this should enter as a config file.

# Input data frame format (data_base_in) should be provided by the precurser 'gather' script: 
# c('source_id','hfname','ISO','ADMIN0','ADMIN1','ADMIN2','ADMIN3','year','month_int','admin_level','LAT','LONG','demographic'
# 'tested','confirmed','test_type','source_file')
# But, some of these colums are only needed downstream. For this code in particular, we need:
# c('source_id','ISO','ADMIN0','ADMIN1','ADMIN2','ADMIN3','year','month_int','admin_level','LAT','LONG')
# LAT and LONG don't need values, they will be replaced anyway.  

#Run in the terminal as follows: 
#Rscript driver_gather_covariates.R --config='path_to_config' --skip_covariate_build
#the config path is necessary, skip_covariate_build is FALSE by default, but can be set to TRUE 
# if desired (to save run time if already built)


#Note, the covariate filepaths assume the code is being run on linux. If running on a windows machine, 
# make sure to change the start of the file paths. 

load_config <- function(config_path) {
  config <- fromJSON(config_path)
  return(config)
}

option_list <- list(
  make_option(c("--config"), type="character", default=NULL, 
              help="Path to the configuration file", metavar="character"),
  make_option(c("--skip_covariate_build"), type="logical", default=FALSE, 
              help="Skip (re-)building covariate databases")
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

#Build the covariate data base from scratch (TRUE or FALSE)?
skip_covariate_build <- opt$skip_covariate_build

# Data base to load (i.e. input data)
case_data_withRe_fn <- config$case_data_withRe_fn


#File name to save covariate data base. format will be 'covariate_file_prefix'_ADMINX_'covariate_file_suffix', where 
# X is the admin level.
covariate_file_prefix <- config$covariate_file_prefix 
covariate_file_suffix<- config$covariate_file_suffix

#File name to save values with covariates to. 
case_data_withReCovs_fn <- config$case_data_withReCovs_fn


# Data base to load (i.e. input data)
# data_base_in <- '/mnt/Z/Michael/Data/hf_data/eth_hf_data_for_seasonality_study_170423.csv'

#Load functions for processing
source('./src/utils/utilities_covariate.R')

#Load file path information
source('./config/data_lookup.R')

# nlagged <- 3
#Number of months before and after for which we collect covariates
nlagged <- 1 

# metric_of_interest <- 'Entomological'
# metric_of_interest <- 'Incidence'
#---------------------------------------------------------------------------------------------------------------#
#             End: Configuration ----- Start: preliminary wrangling and covariate path specification 
#---------------------------------------------------------------------------------------------------------------#

limited_df <- read.csv(case_data_withRe_fn)
limited_df <- limited_df[!((limited_df$year==max(limited_df$year))&(limited_df$month==12)),]


# pop_root_str <- '/mnt/Z/mastergrids/Other_Global_Covariates/Population/WorldPop/v3/PopulationCounts_DRC_fixed/5km/WorldPop_UNAdj_v3_DRC_fix.'
# pop_suffix_str <- '.Annual.Data.5km.sum.tif'




root_str <- dyn_root_str_list[annual_cov_name_list[1]]
suff_str <- dyn_suff_str_list[annual_cov_name_list[1]]
month_name_list <- c('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC')
month_idx_list <- c('01', '02', '03', '04', '05', '06','07', '08', '09', '10', '11', '12')

# shape_file_paths = list('prefix'='/mnt/Z/master_geometries/Admin_Units/Global/MAP/2022/admin2022_', 
#                         'suffix'='.shp')

#------------------------------------------------------------------------------------------#
#       Step 1: build point data frame
#------------------------------------------------------------------------------------------#
# retain original_row_num, if point use LAT LONG, if higher admin level sample points.

# Sub step 1: - make empty points_df: columns = 'site_id', 'lat', 'long'
#             - make column data frame id_map: columns = 'original_row_num', 'site_id', 'iso', 'ADMIN_LEVEL', 'ADMIN_UNIT'
# Sub step 2: - for i in original_row_num: if ADMIN_LEVEL==POINT, add lat long + 'site_id' to points_df
#             - add 'site_id' in id_map
#             - if ADMIN_LEVEL!=POINT: check for match in id_map (add data to id_map), 
#             - else load appropriate shape file, sample points, populate points_df and id_map.
#             - return list(id_map=id_map, points_df = points_df)

limited_df['admin_level'] <- 'POINT'
admin_levels <- limited_df$admin_level %>% unique()
limited_df[,'source_id'] <- NA
source_id <- 1
for (admin_level in admin_levels){
     if(admin_level =='POINT'){
        limited_df[limited_df$admin_level=='POINT',]$source_id <-  source_id:nrow(limited_df[limited_df$admin_level=='POINT',])
         source_id <- source_id + nrow(limited_df[limited_df$admin_level=='POINT',])
     }else{
         admin_units <- limited_df[limited_df$admin_level==admin_level,admin_level] %>% unique()
         for (admin_unit in admin_units){
             limited_df[limited_df[,admin_level]==admin_unit, 'source_id'] <-source_id
             source_id<- source_id + 1
         }
    }
}

limited_df <- limited_df[!is.na(limited_df$LAT)&(!is.na(limited_df$LONG)),]


#------------------------------------------------------------------------------------------#
#       Step 2: Build yymm dataframe for each row (original_row_num)
#------------------------------------------------------------------------------------------#
#Columns are: c(original_row_num, yymm_nolag, lag1, ...) yymm_nolag is used for static covariates
#pivoted_df has all the same columns as limited_df, plus 4 yymm columns for dates of covariates
message('Adding yymm columns')
yymm_df <- build_yymm_df(limited_df, nlagged)

#------------------------------------------------------------------------------------------#
#       Step 3: pivot original dataframe and add covariates column by column
#------------------------------------------------------------------------------------------#
#Here is the time-consuming part of the code
#Substep 1: build a 'covariate data base'. This gets the population-weighted mean covariate values
# for each admin unit, for each admin level for all years and months. 

#Save intermediate file with all points, years etc. Then extract required covariates from this file
covariate_meta <- list()
covariate_meta$min_year <- min(limited_df$year) -1
covariate_meta$max_year <- max(limited_df$year)

covariate_meta$dyn_cov_name_list <- dyn_cov_name_list
covariate_meta$annual_cov_name_list <- annual_cov_name_list 
covariate_meta$static_cov_name_list <- static_cov_name_list 

# covariate_meta$month_name_list <- month_name_list


covariate_meta$out_fn_prefix  <- covariate_file_prefix
covariate_meta$out_fn_suffix  <- covariate_file_suffix

#Build covariate data base and save (file name accounts for admin level automatically)
if(!skip_covariate_build){
    build_full_covariate_database(limited_df, covariate_meta) 
}

for (i in -nlagged:nlagged){
  if(i==0){
    limited_df['covariate_no_lag'] <- yymm_df['covariate_no_lag']
  }else if(i > 0){
    limited_df[paste0('covariate_lag',toString(i))] <- yymm_df[paste0('covariate_lag',toString(i))]
  }else{
    limited_df[paste0('covariate_lead',toString(-i))] <- yymm_df[paste0('covariate_lead',toString(-i))]   
  }
}

covariate_df <-read.csv(paste0(covariate_meta$out_fn_prefix,admin_level,covariate_meta$out_fn_suffix)) 
na_rows <- covariate_df[rowSums(is.na(covariate_df[dyn_cov_name_list])) > 0, ]
limited_df <- limited_df[!(paste(limited_df$LAT, limited_df$LONG) %in% 
                          paste(na_rows$LAT, na_rows$LONG)), ]

# admin_levels <- limited_df$admin_level %>% unique() 
# for(admin_level in admin_levels){
#     covariate_df <-read.csv(paste0(covariate_meta$out_fn_prefix,admin_level,covariate_meta$out_fn_suffix)) 
#     missing_isos <- covariate_df[is.na(covariate_df[dyn_cov_name_list[1]]),]$ISO %>% unique()
#     for (iso in missing_isos){
#         missing_admin_units <- covariate_df[is.na(covariate_df[dyn_cov_name_list[1]]),admin_level] %>% unique()
#         limited_df <- limited_df[((limited_df$admin_level==admin_level) &!(limited_df[[admin_level]]%in%missing_admin_units))|(limited_df$admin_level!=admin_level),]
#     }
# }


#Write checkpoint (just in case of a crash)
write.csv(limited_df, case_data_withReCovs_fn,  row.names = FALSE)
#Append covariates from covariate database
limited_df <- append_covariate_columns(limited_df, covariate_meta, nlagged=nlagged)
write.csv(limited_df, case_data_withReCovs_fn,  row.names = FALSE)

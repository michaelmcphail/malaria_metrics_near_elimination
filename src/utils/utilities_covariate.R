# Script Name: utilities_covariate.R
# Purpose: Utility functions for covariate assembly.
# Author: Michael McPhail
# Date: 2024-07-09
# Dependencies: dplyr, rgdal, sf, sp, spatstat, terra, tidyr
# Usage: source(PATH_TO_FILE)

library(dplyr)
library(rgdal)
library(raster)
library(sf)
library(sp)
library(spatstat)
library(terra)
library(tidyr)

source('./config/data_lookup.R')

#' For a data frame with a column called 'admin_level' indicating data admin level (e.g. 'ADMIN0', ...,'ADMIN3'), find matching admin unit in MAP
#' shape file. Data frame/MAP names may be messy (with '-', '/', spaces, capitalisations). So matching involves stripping all of these characters
#' from the string. Also match the ISO becasue admin level names are non-unique internationally. Warning, assumes that there are no special 
#' characters (e.g. accented characters). It won't be able to match these. This function also replaces the ADMINx name with the MAP name.
#'
#' @param weighted_mean_df A data frame with columns ['admin_level', 'ISO', 'ADMINx' (for each admin level required), 'LAT', 'LONG']
#' @param shape_file_paths A list with entries shape_file_paths$prefix (the first part of the shape file path string)
#' and shape_file_paths$suffix where paste0(shape_file_paths$prefix,admin_level_number,shape_file_paths$suffix) is the path to a shape file for 
#' admin_level_number one of {0, 1, 2, 3}.

# set_admin_unit_coords_to_centroids <- function(weighted_mean_df, shape_file_paths){
# #TODO: still many missing locations, could be string matching, but I think the issue is deeper (completely different names)
#     missing_admin_units = c()
#     for(admin_level_number in 0:3){
#         # sf_use_s2(FALSE) #Previously used because of edge intercepts, not ideal because of wrong coordinate system, st_make_valid() is a better fix
#         message(paste('Admin',toString(admin_level_number)))
#         admin_level_str <- paste0('ADMIN', toString(admin_level_number))
#         admin_rows <- weighted_mean_df[,'admin_level']==admin_level_str 
#         iso_list <- weighted_mean_df[admin_rows,'ISO'] %>% unique() #Get list of all unique ISOs in data frame for this admin_level.
#         #Read full shape file (could restrict to africa here, but not really necessary)
#         full_sf <- st_read(paste0(shape_file_paths$prefix,admin_level_number,shape_file_paths$suffix))

#         for(iso in iso_list){
#             message(paste0('ISO: ', iso,' for admin_level', admin_level_str))
#             admin_rows <- (weighted_mean_df[, 'admin_level']==admin_level_str) &(weighted_mean_df[, 'ISO']==iso) 
#             admin_unit_list <- weighted_mean_df[admin_rows,admin_level_str] %>% unique() #Get all admin units in ISO at admin level. 
#             #TODO: Ugly repeated code, compress to single taking admin_level_number as input.
#             for(admin_unit in admin_unit_list){
#                 if(admin_level_number==0){
#                     matching_admin_level <- tolower(gsub('-','',gsub(' ', '', gsub('/', '',full_sf$Name_0))))==tolower(gsub('-','',gsub(' ', '', gsub('/', '',admin_unit))))
#                     matching_iso <- full_sf$ISO == iso
#                     centroid_coords <- full_sf[matching_iso&matching_admin_level,'geometry'] %>% st_make_valid() %>% st_centroid() %>% st_coordinates()
#                     shp_file_name <- full_sf[matching_admin_level,]$Name_0 
#                 }else if(admin_level_number==1){
#                     matching_admin_level <- tolower(gsub('-','',gsub(' ', '', gsub('/', '',full_sf$Name_1))))==tolower(gsub('-','',gsub(' ', '', gsub('/', '',admin_unit))))
#                     matching_iso <- full_sf$ISO == iso
#                     centroid_coords <- full_sf[matching_iso&matching_admin_level,'geometry'] %>% st_make_valid() %>% st_centroid() %>% st_coordinates()
#                     shp_file_name <- full_sf[matching_admin_level,]$Name_1 
#                 }else if(admin_level_number==2){
#                     matching_admin_level <- tolower(gsub('-','',gsub(' ', '', gsub('/', '',full_sf$Name_2))))==tolower(gsub('-','',gsub(' ', '', gsub('/', '',admin_unit))))
#                     matching_iso <- full_sf$ISO == iso
#                     centroid_coords <- full_sf[matching_iso&matching_admin_level,'geometry'] %>% st_make_valid() %>% st_centroid() %>% st_coordinates()
#                     shp_file_name <- full_sf[matching_admin_level,]$Name_2 
#                 }else if(admin_level_number==3){
#                     matching_admin_level <- tolower(gsub('-','',gsub(' ', '', gsub('/', '',full_sf$Name_3))))==tolower(gsub('-','',gsub(' ', '', gsub('/', '',admin_unit))))
#                     matching_iso <- full_sf$ISO == iso
#                     centroid_coords <- full_sf[matching_iso&matching_admin_level,'geometry'] %>% st_make_valid() %>% st_centroid() %>% st_coordinates()
#                     shp_file_name <- full_sf[matching_admin_level,]$Name_3 

#                 }
#                 if(nrow(centroid_coords)>0){
#                     weighted_mean_df[admin_rows&(weighted_mean_df[,admin_level_str]==admin_unit), 'LAT']<-centroid_coords[1,2] 
#                     weighted_mean_df[admin_rows&(weighted_mean_df[,admin_level_str]==admin_unit), 'LONG']<-centroid_coords[1,1]
#                     weighted_mean_df[admin_rows&(weighted_mean_df[,admin_level_str]==admin_unit),  admin_level_str]<-shp_file_name[1]
#                 }else{
#                     #No coordinates => no match
#                     message(paste(unique(weighted_mean_df[admin_rows&(weighted_mean_df[,admin_level_str]==admin_unit), 'ISO']), admin_unit))
#                     missing_admin_units <- c(missing_admin_units, admin_unit) #Could return this as well if you want to study the failed matches.
#                 }

#             }
#         }
#     }

#   #Return data frame with Latitude, Longitude, admin unit names replaced for all matching admin units.
#   return(weighted_mean_df)
# }

#' Build yymm data frame 
#'
#' Builds a data frame indicating which years (YYYY) and months (MM) to pull each (vanilla, lagged, lead) covariate from
#' for each row of data
#' 
#' @param long_df The data frame in long form (1 row for each observation).
#' @param n_lagged Number of months to lag and lead by.
#' @return the input data frame with additional columns appended
#' @export
build_yymm_df <- function(long_df, nlagged){
# add yyyy.mm columns associated with the year (yyyy) and month (mm) of the covariates (with lagged dates for dynamic covariates)
    month_idx_list <- c('01', '02', '03', '04', '05', '06','07', '08', '09', '10', '11', '12')
    get_yymm_str <- function(year, month_int, lag) {
        month <-month_int 
        year <- year
        if ((month_int - lag) <1){
            year <- year -1
        }else if ((month_int - lag) >12 ){
            year <- year + 1
        }
        month_lag <- (month -1 - lag + 12)%%12 +1
        return(paste0(toString(year), '.', month_idx_list[month_lag]))
    }
    for (i in (-nlagged):nlagged){
      if(i==0){
        long_df['covariate_no_lag'] <- mapply(get_yymm_str, long_df$year, long_df$month_int, i)
      }else{
        if(i <0){
          long_df[paste0('covariate_lead',toString(-i))] <- mapply(get_yymm_str, long_df$year, long_df$month_int, i)
        }else{
          long_df[paste0('covariate_lag',toString(i))] <- mapply(get_yymm_str, long_df$year, long_df$month_int, i)
        }
      }
    }
    return(long_df)
}

#' Build full covariate data base
#'
#' Build the covariate data base for the indicated year range, which will later be used to extract the covariates 
#' The admin form hasn't been updated from the seasonality code, so at the moment only the point extraction will work (I think)
#'
#' @param limited_df the data frame containing the geometric units of interest (coordinates or admin unit names)
#' @param covariate_meta A list of covariate meta data (min_year, max_year, dyn_cov_name_list, annual_cov_name_list, static_cov_name_list)
#' the file prefixes/suffixes for saving data out.
#' @return nothing, but saves covariate csv's at each stage.
#' @export
build_full_covariate_database <- function(limited_df, covariate_meta){
  admin_levels <- limited_df$admin_level %>% unique()
  for (admin_level in admin_levels){
      if (admin_level =='POINT'){
          append_covariates_points(limited_df, covariate_meta)
      }
      else{
          append_covariates_admin_units(limited_df, covariate_meta, admin_level)
      }
  }
}

#' Append covariate points
#'
#' Build the covariate data set for points
#'
#' @param limited_df The input data set
#' @param covariate_meta The covariate meta data
#' @return no return but saves csv of assembled data.
#' @export
append_covariates_points <- function(limited_df, covariate_meta){

    min_year <- covariate_meta$min_year 
    max_year <- covariate_meta$max_year 

    dyn_cov_name_list <- covariate_meta$dyn_cov_name_list 
    annual_cov_name_list <- covariate_meta$annual_cov_name_list 
    static_cov_name_list <- covariate_meta$static_cov_name_list

    out_fn <- covariate_meta$out_fn 
    out_fn_prefix <- covariate_meta$out_fn_prefix 
    out_fn_suffix <- covariate_meta$out_fn_suffix 
    out_fn <- paste0(out_fn_prefix, 'POINT', out_fn_suffix) 

    years <- min_year:max_year
    month_idx_list<-c('01', '02', '03', '04', '05', '06','07', '08', '09', '10', '11', '12')

    new_col_names = c()
    for (year in years){
        for (month in month_idx_list){
        new_col_names = c(new_col_names, paste0(toString(year), '_',month, sep=''))
        }
    }

    source_ids <- limited_df[limited_df$admin_level=='POINT',]$source_id %>% unique()
    admin_unit_df = unique(limited_df[(limited_df$admin_level=='POINT'),c('source_id', 'ISO','ADMIN0', 'LAT', 'LONG')])
    # admin_unit_df = unique(limited_df[(limited_df$admin_level=='POINT'),c('source_id','Name_3', 'Name_2', 'Name_1', 'LAT', 'LONG')])

    admin_unit_df[new_col_names] = NA
    admin_unit_df<- admin_unit_df %>% 
        pivot_longer(
        cols=all_of(new_col_names),
        names_to = c("year", "month"),
        names_sep="_",
        values_to='covariate')

    admin_unit_df<-admin_unit_df[c('source_id', 'ISO','ADMIN0', 'LAT', 'LONG','year','month')]


    for (cov_name in dyn_cov_name_list){
        message(paste('Extracting dynamic covariate', cov_name))
        admin_unit_df[cov_name] <-NA

        for (year in years){
            for (month in month_idx_list){
                rows_to_extract<- (admin_unit_df$year==year)&(admin_unit_df$month==month)
                if((cov_name %in% c('LST_day', 'LST_delta', 'LST_night')) & (year>2021)) {# a temporary fix because 2022 LST values are in a different directory
                    cov_fn = paste0(dyn_root_post2021_str_list[cov_name], toString(year),'.', month , dyn_suff_str_list[cov_name]) 
                }else{
                    cov_fn = paste0(dyn_root_str_list[cov_name], toString(year),'.', month , dyn_suff_str_list[cov_name]) 
                    # cov_fn <- paste0(dynamic_file_name_df[dynamic_file_name_df$cov_name_list==cov_name,'root_str_list'], toString(year),'.', month,dynamic_file_name_df[dynamic_file_name_df$cov_name_list==cov_name,'suff_str_list']) 
                }
                cov_rast <- terra::rast(cov_fn)
                covariate_sub_values <- terra::extract(cov_rast, 
                                                        admin_unit_df[rows_to_extract,c('LONG', 'LAT')])
                admin_unit_df[rows_to_extract,cov_name]<- covariate_sub_values[,2]

            }
            message(paste('Completed',year))

        }
        write.csv(admin_unit_df, out_fn)
    }

    #Add annual covariates such as population (remove month loop)
    for (cov_name in annual_cov_name_list){
        message(paste('Extracting static covariate', cov_name))
        admin_unit_df[cov_name] <-NA
        for (year in years){
                rows_to_extract<- (admin_unit_df$year==year)
                cov_fn <- paste0(annual_root_str_list[cov_name], toString(year), annual_suff_str_list[cov_name]) 
                if((year >2020)&(cov_name=='POPULATION')){
                    cov_fn <- paste0(annual_root_str_list[cov_name], '2020', annual_suff_str_list[cov_name])  #Note pop_roo/suffix is defined in main text, correct this.
                }
                cov_rast <- terra::rast(cov_fn)
                covariate_sub_values <- terra::extract(cov_rast, 
                                                        admin_unit_df[rows_to_extract,c('LONG', 'LAT')])
                admin_unit_df[rows_to_extract,cov_name]<- covariate_sub_values[,2]
            message(paste('Completed',year))

        }
        write.csv(admin_unit_df, out_fn)
    }

    #Add static covariates such as elevation Remove month and year loop
    for (cov_name in static_cov_name_list){
        message(paste('Extracting static covariate', cov_name))
        admin_unit_df[cov_name] <-NA
        cov_fn <- static_fn_str_list[cov_name]
        cov_rast <- terra::rast(cov_fn)
        covariate_sub_values <- terra::extract(cov_rast, 
                                                admin_unit_df[,c('LONG', 'LAT')])
        admin_unit_df[,cov_name]<- covariate_sub_values[,2]

        write.csv(admin_unit_df, out_fn)
    }
    message('Finished adding covariates')
    return(admin_unit_df)
}

#' Append covariate ADMINX
#'
#' Build covariate data set for admin levels (not relevant for Vietnam study so not tested)
#'
#' @param limited_df The input data set
#' @param covariate_meta The covariate meta data
#' @param admin_level the specific admin level to collect
#' @return no return but saves csv of assembled data.
#' @export
append_covariates_admin_units <- function(limited_df, covariate_meta, admin_level){
    min_year <- covariate_meta$min_year 
    max_year <- covariate_meta$max_year 

    dyn_cov_name_list <- covariate_meta$dyn_cov_name_list 
    annual_cov_name_list <- covariate_meta$annual_cov_name_list 
    static_cov_name_list <- covariate_meta$static_cov_name_list

    out_fn_prefix <- covariate_meta$out_fn_prefix 
    out_fn_suffix <- covariate_meta$out_fn_suffix 

    shp_prefix <- shape_file_paths$prefix
    shp_suffix <- shape_file_paths$suffix

    out_fn <- paste0(out_fn_prefix, admin_level, out_fn_suffix) 

    years <- min_year:max_year
    month_idx_list<-c('01', '02', '03', '04', '05', '06','07', '08', '09', '10', '11', '12')
    new_col_names = c()
    for (year in years){
        for (month in month_idx_list){
        new_col_names = c(new_col_names, paste0(toString(year), '_',month, sep=''))
        }
    }

    #NOTE, this part will break with current specifications of pipeline, as ADMIN1...3 columns are dropped at an earlier stage
    admin_unit_df = unique(limited_df[(limited_df$admin_level==admin_level),c('source_id', 'ISO','ADMIN0', 'ADMIN1', 'ADMIN2', 'ADMIN3', 'LAT', 'LONG')])
    admin_unit_df[new_col_names] = NA
    admin_unit_df<- admin_unit_df %>% 
        pivot_longer(
        cols=all_of(new_col_names),
        names_to = c("year", "month"),
        names_sep="_",
        values_to='covariate')


    admin_unit_df<-admin_unit_df[c('source_id', 'ISO','ADMIN0', 'ADMIN1', 'ADMIN2', 'ADMIN3', 'LAT', 'LONG','year','month')]

    polygon_names<- c()
    admin_idx <- substr(admin_level,nchar(admin_level), nchar(admin_level))
    shp_count <-1

    admin_shp <- st_read(paste0(shp_prefix,admin_idx,shp_suffix))
    admin_shp <- admin_shp[admin_shp[[paste0('Name_', admin_idx)]] %in% unique(unlist(admin_unit_df[,admin_level])),]
    polygon_names <- as.data.frame(admin_shp[,paste0('Name_', admin_idx)]) %>% select(paste0('Name_', admin_idx))

    message(unique(unlist(admin_unit_df[,admin_level])))
    for (cov_name in dyn_cov_name_list){
        message(paste('Extracting dynamic covariate', cov_name))
        admin_unit_df[cov_name] <-NA

        for (year in years){
                if(year <=2020){
                    pop_fn <- paste0(annual_root_str_list['POPULATION'], toString(year), annual_suff_str_list['POPULATION']) 
                }else{
                    pop_fn <- paste0(annual_root_str_list['POPULATION'], '2020', annual_suff_str_list['POPULATION'])
                }
                pop_rast <- terra::rast(pop_fn)
                pop_rast[is.na(pop_rast)]<-0
                for (month in month_idx_list){
                   if((cov_name %in% c('LST_day', 'LST_delta', 'LST_night')) & (year > 2021)){
                      cov_fn <- paste0(dyn_root_post2021_str_list[cov_name], toString(year),'.', month, dyn_suff_str_list[cov_name]) 
                   }else{
                      cov_fn <- paste0(dyn_root_str_list[cov_name], toString(year),'.', month, dyn_suff_str_list[cov_name]) 
                   }
                   cov_rast <- terra::rast(cov_fn)
                   cov_out<- exactextractr::exact_extract(cov_rast, admin_shp, 
                                                        function(values, coverage_frac, weights) {
                                                        weighted.mean(values, coverage_frac * weights, na.rm=TRUE)
                                                        },  
                                                        weights=pop_rast,
                                                       stack_apply=TRUE)
                    #This could probably be done faster (vector assign), but I need to check if polygon_names is in the same order
                    for (i in 1:length(cov_out)){
                      admin_unit_df[(admin_unit_df[,admin_level]==polygon_names[i,])&(admin_unit_df$year==year)&(admin_unit_df$month==month),cov_name]<- cov_out[i]
                    }
                }
                message(paste('Completed',year))

            }
        write.csv(admin_unit_df, out_fn)
    }

    #Add annual covariates such as population (remove month loop)
    for (cov_name in annual_cov_name_list){
        message(paste('Extracting static covariate', cov_name))
        admin_unit_df[cov_name] <-NA


        for (year in years){
            if(year <=2020){
                pop_fn <- paste0(annual_root_str_list['POPULATION'], toString(year), annual_suff_str_list['POPULATION']) 
            }else{
                pop_fn <- paste0(annual_root_str_list['POPULATION'], '2020', annual_suff_str_list['POPULATION'])
            }
            pop_rast <- terra::rast(pop_fn)
            pop_rast[is.na(pop_rast)]<-0
            if (cov_name != 'POPULATION'){
              cov_fn <- paste0(annual_root_str_list[cov_name], toString(year), annual_suff_str_list[cov_name]) 
              cov_rast <- terra::rast(cov_fn)
              cov_out<- exactextractr::exact_extract(cov_rast, admin_shp, 
                                                function(values, coverage_frac, weights) {
                                                weighted.mean(values, coverage_frac * weights, na.rm=TRUE)
                                                },  
                                                weights=pop_rast,
                                                stack_apply=TRUE)
            }else{
                cov_out<- exactextractr::exact_extract(pop_rast, admin_shp,'mean', stack_apply=TRUE)                
            }
            #This could probably be done faster (vector assign), but I need to check if polygon_names is in the same order
            for (i in 1:length(cov_out)){
                admin_unit_df[(admin_unit_df[,admin_level]==polygon_names[i,])&(admin_unit_df$year==year),cov_name]<- cov_out[i]
            }
            message(paste('Completed',year))
        }
        write.csv(admin_unit_df, out_fn)

    }

  #Add annual covariates such as population (remove month loop)
    for (cov_name in static_cov_name_list){
        message(paste('Extracting static covariate', cov_name))
        admin_unit_df[cov_name] <-NA


        for (year in years){ #demographic changes means that 'static' population weighted covariates change in time.
            if(year <=2020){
                pop_fn <- paste0(annual_root_str_list['POPULATION'], toString(year), annual_suff_str_list['POPULATION']) 
            }else{
                pop_fn <- paste0(annual_root_str_list['POPULATION'], '2020', annual_suff_str_list['POPULATION'])
            }
            pop_rast <- terra::rast(pop_fn)
            pop_rast[is.na(pop_rast)]<-0
            cov_fn <- paste0(static_fn_str_list[cov_name])

            cov_rast <- terra::rast(cov_fn)
            if (cov_name != 'POPULATION'){
               cov_out<- exactextractr::exact_extract(cov_rast, admin_shp, 
                                                function(values, coverage_frac, weights) {
                                                weighted.mean(values, coverage_frac * weights, na.rm=TRUE)
                                                },  
                                                weights=pop_rast,
                                                stack_apply=TRUE)
            }else{
                cov_out<- exactextractr::exact_extract(cov_rast, admin_shp,'mean', stack_apply=TRUE)                
            }
            #This could probably be done faster (vector assign), but I need to check if polygon_names is in the same order
            for (i in 1:length(cov_out)){
                admin_unit_df[(admin_unit_df[,admin_level]==polygon_names[i,])&(admin_unit_df$year==year),cov_name]<- cov_out[i]
            }
            message(paste('Completed',year))
        }
        write.csv(admin_unit_df, out_fn)

    }


    message('Finished adding covariates')
    return(admin_unit_df)
}


#' Append covariate column
#'
#' Append the covariates to the data for the regression
#'
#' @param limited_df The input data set
#' @param covariate_meta The covariate meta data
#' @param nlagged The number of months to include that are laged/leading
#' @return an augmentation of the input data with additional covariate columns.
#' @export
append_covariate_columns<- function(limited_df, covariate_meta, nlagged=1){
    dyn_cov_name_list <- covariate_meta$dyn_cov_name_list 
    annual_cov_name_list <- covariate_meta$annual_cov_name_list 
    static_cov_name_list <- covariate_meta$static_cov_name_list
    cov_fn_prefix <- covariate_meta$out_fn_prefix 
    cov_fn_suffix <- covariate_meta$out_fn_suffix     

    # This should only ever return the vector c('POINT') unless the code base is extended to include admin-units
    admin_levels <- limited_df$admin_level %>% unique()

    base_and_lagged_names <- c()
    for (i in (-nlagged):nlagged){
      if(i==0){
        base_and_lagged_names<- c(base_and_lagged_names, 'covariate_no_lag')
      }else if(i>0){
        base_and_lagged_names <- c(base_and_lagged_names, paste0('covariate_lag',toString(i)))
      }else{
         base_and_lagged_names <- c(base_and_lagged_names, paste0('covariate_lead',toString(-i)))
      }
    } 

    for(cov_name in dyn_cov_name_list){
        for(i in 1:length(base_and_lagged_names)){
            if(grepl('_no_', base_and_lagged_names[i])){
                column_name <- cov_name
            }else if(grepl('_lag', base_and_lagged_names[i])){
                column_name <- paste0(cov_name,substring(base_and_lagged_names[i],  nchar(base_and_lagged_names[i]) - 4, nchar(base_and_lagged_names[i])))
            }else{
                 column_name <- paste0(cov_name,substring(base_and_lagged_names[i],  nchar(base_and_lagged_names[i]) - 5, nchar(base_and_lagged_names[i])))
            }
            limited_df[column_name]<- NA
        }
    }

    for(cov_name in c(annual_cov_name_list,static_cov_name_list)){
        limited_df[cov_name]<- NA
    }

    for(admin_level in admin_levels){
        cov_fn <- paste0(cov_fn_prefix, admin_level, cov_fn_suffix)
        covariate_db <- read.csv(cov_fn)

        if(admin_level=='POINT'){
            #match by source_id
            admin_rows <- limited_df$admin_level=='POINT'
            source_ids <- limited_df[admin_rows,]$source_id %>% unique()
            pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(source_ids), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to build progress bar
            finished_source_ids <-0

            for(source_id in source_ids){
                source_id_rows <- limited_df$source_id == source_id
                cov_source_id_rows <- covariate_db$source_id == source_id
                for(cov_name in dyn_cov_name_list){
                    for(i in 1:length(base_and_lagged_names)){
                        ref_col_name <- base_and_lagged_names[i]
                        if(grepl('_no_', base_and_lagged_names[i])){
                            column_name <- cov_name
                            limited_df[source_id_rows, column_name] <- add_covariate_column_subset(limited_df[source_id_rows, ref_col_name], cov_name, covariate_db[cov_source_id_rows,])

                        }else if(grepl('_lag', base_and_lagged_names[i])){
                            column_name <- paste0(cov_name,substring(ref_col_name,  nchar(ref_col_name) - 4, nchar(ref_col_name)))
                            limited_df[source_id_rows, column_name] <- add_covariate_column_subset(limited_df[source_id_rows, ref_col_name], cov_name, covariate_db[cov_source_id_rows,])
                        }else{
                            column_name <- paste0(cov_name,substring(base_and_lagged_names[i],  nchar(base_and_lagged_names[i]) - 5, nchar(base_and_lagged_names[i])))
                            limited_df[source_id_rows, column_name] <- add_covariate_column_subset(limited_df[source_id_rows, ref_col_name], cov_name, covariate_db[cov_source_id_rows,])
                        }
                     }#End lag loop
                }#end dynamic covariate loop
                for(cov_name in c(annual_cov_name_list,static_cov_name_list)){
                        # ref_col_name <- base_and_lagged_names[1]
                        ref_col_name <-'covariate_no_lag'
                        column_name <- cov_name
                        limited_df[source_id_rows, column_name] <- add_covariate_column_subset(limited_df[source_id_rows, ref_col_name], cov_name, covariate_db[cov_source_id_rows,])
                }#end static/annual covariate loop
            setTxtProgressBar(pb, finished_source_ids)            
            finished_source_ids<-finished_source_ids+1
            }

        }else{
            na_cov_admin_units <- c()
            for (cov in c(dyn_cov_name_list,annual_cov_name_list,static_cov_name_list)){
                na_admin_units <- covariate_db[is.na(cov), admin_level] %>% unique()
                na_cov_admin_units <-c(na_cov_admin_units, na_admin_units)
            }
            na_cov_admin_units <- na_cov_admin_units %>% unique()
            rows_to_drop <- (limited_df$admin_level==admin_level)&(limited_df[[admin_level]]%in%na_cov_admin_units)
            limited_df<- limited_df[!rows_to_drop,]
            admin_rows <- limited_df$admin_level==admin_level
            for (iso in limited_df[admin_rows, 'ISO'] %>% unique()){
                admin_units <- limited_df[admin_rows&(limited_df$ISO==iso), admin_level] %>% unique()
                for(admin_unit in admin_units){
                    admin_unit_rows <- (limited_df$ISO==iso) & admin_rows & limited_df[[admin_level]]==admin_unit
                    #If there is more than 1 admin unit with that name in the covariate_db, check higher admin level for match
                    #TODO: what if there is more than one in limited_df?
                    cov_admin_unit_rows <- (covariate_db[[admin_level]] == admin_unit) &(covariate_db$ISO==iso)
                    if(admin_level=='ADMIN2'){
                        higher_admin_level <- 'ADMIN1'
                    }else{
                        higher_admin_level <- 'ADMIN2' #obviously doesn't cover all cases, but iso~ADMIN0 is already considered
                    }
                    if(length(covariate_db[cov_admin_unit_rows,higher_admin_level] %>% unique())>1){
                        #more than one adminX unit with the same name
                        higher_admin_unit <- limited_df[admin_unit_rows, higher_admin_level] %>% unique() 
                        cov_admin_unit_rows <- cov_admin_unit_rows&(covariate_db[[higher_admin_level]] == higher_admin_unit)
                    }
                    
                    message(paste0('adding dynamic covariates for ', admin_unit))

                    for(cov_name in dyn_cov_name_list){
                        for(i in 1:length(base_and_lagged_names)){
                            ref_col_name <- base_and_lagged_names[i]
                            if(grepl('_no_', base_and_lagged_names[i])){
                                column_name <- cov_name
                                limited_df[admin_unit_rows, column_name] <- add_covariate_column_subset(limited_df[admin_unit_rows, ref_col_name], cov_name, covariate_db[cov_admin_unit_rows,])
                                
                            }else if(grepl('_lag', base_and_lagged_names[i])){
                                column_name <- paste0(cov_name,substring(ref_col_name,  nchar(ref_col_name) - 4, nchar(ref_col_name)))
                                limited_df[admin_unit_rows, column_name] <- add_covariate_column_subset(limited_df[admin_unit_rows, ref_col_name], cov_name, covariate_db[cov_admin_unit_rows,])
                            }else{
                                column_name <- paste0(cov_name,substring(base_and_lagged_names[i],  nchar(base_and_lagged_names[i]) - 5, nchar(base_and_lagged_names[i])))
                                limited_df[admin_unit_rows, column_name] <- add_covariate_column_subset(limited_df[admin_unit_rows, ref_col_name], cov_name, covariate_db[cov_admin_unit_rows,])
                            } 
                        }#End lag loop
                    
                    }#end dynamic covariate loop
                    message(paste0('finished adding dynamic covariates for ', admin_level))
                    for(cov_name in c(annual_cov_name_list,static_cov_name_list)){
                        # ref_col_name <- base_and_lagged_names[1]
                        ref_col_name <-'covariate_no_lag' 
                        column_name <- cov_name
                        limited_df[admin_unit_rows, column_name] <- add_covariate_column_subset(limited_df[admin_unit_rows, ref_col_name], cov_name, covariate_db[cov_admin_unit_rows,])
                    }#end static covariate loop
                }#end admin_unit loop
          } #END ISO loop (needed because multiple isos have equal adminx units)
        }
        message(paste('Finshed with', admin_level))
    }#end admin_level loop
    return(limited_df)

}

#'  Add covariate column subset
#'
#' Append the specified covariate for the subset (based on admin unit) of rows
#'
#' @param yy_mm_column_subset The subset of limited_df (based on admin_unit name, in this case a specific point), and only the yy.mm column of interest
#' @param cov_name Single covariate name.
#' @param covariate_db_subset The subset of the full covariate data set corrresponding to a specific admin level (in this case POINT).
#' @return an augmentation of the input data with additional covariate columns.
#' @export
add_covariate_column_subset <- function(yy_mm_column_subset, cov_name,covariate_db_subset){
    month_idx_to_int_map <- list()
    month_idx_list <- c('01', '02', '03', '04', '05', '06','07', '08', '09', '10', '11', '12')
    month_int_list <- 1:12
    for (i in 1:length(month_idx_list)){
        month_idx_to_int_map[month_idx_list[i]] = month_int_list[i]
    }
    covariate_column_subset <- c()
    for(j in 1:length(yy_mm_column_subset)){
        year <- substring(yy_mm_column_subset[j], 1, 4)
        month_str <- substring(yy_mm_column_subset[j], 6, 7)
        month_int <- month_idx_to_int_map[[month_str]]
        if (nrow(covariate_db_subset[(covariate_db_subset$year==year)&(covariate_db_subset$month==month_int), ])==0){
            covariate_column_subset<-c(covariate_column_subset, NA)
        }else{
        covariate_column_subset <- c(covariate_column_subset, covariate_db_subset[(covariate_db_subset$year==year)&(covariate_db_subset$month==month_int), cov_name])}
    }
    return(covariate_column_subset)
}


# aggregate_to_admin_level <- function(limited_df, iso_list){
#     aggregated_df <- limited_df[0,]
    
#     for(iso in iso_list){
#       iso_rows <-  limited_df$ISO==iso
#       admin_levels <- limited_df[iso_rows, ]$admin_level %>% unique()
#         for (admin_level in admin_levels){
#           admin_level_rows<- iso_rows&(limited_df$admin_level==admin_level)
#           admin_units <- limited_df[admin_level_rows, admin_level] %>% unique()
#           for(admin_unit in admin_units){
#             admin_unit_rows<- admin_level_rows&(limited_df[[admin_level]]==admin_unit)
#             demographics <- limited_df[admin_unit_rows, ]$demographic %>% unique()
#               for(demographic in demographics){
#                 demographic_rows<- admin_unit_rows&(limited_df$demographic==demographic)
#                 test_types <- limited_df[demographic_rows, ]$test_type %>% unique()
#                   for(test_type in test_types){
#                     test_type_rows<- demographic_rows&(limited_df$test_type==test_type)
#                     years <- limited_df[test_type_rows, ]$year %>% unique()
#                     for(year in years){
#                        year_rows_rows<- demographic_rows&(limited_df$year==year)
#                        months <- limited_df[year_rows_rows, ]$month_int %>% unique()
#                         for (month in months){
#                             month_rows <- year_rows_rows &(limited_df$month_int==month)
#                             new_row <- limited_df[month_rows, ][1,]
#                             new_row['hfname']<-""
#                             new_row['tested']<-limited_df[month_rows&(!is.na(limited_df$tested)), ]$tested %>% sum()
#                             new_row['confirmed']<-limited_df[month_rows&(!is.na(limited_df$confirmed)), ]$confirmed %>% sum()
                            
#                             aggregated_df <- rbind(aggregated_df, new_row)
#                         }    

#                     }#end year loop
                      
#                   }#end test_type loop
                  
#               } #end demographic loop
#             } # end admin unit loop
#         } #end admin_level loop
#     } #end iso loop
    
#     return(aggregated_df)
# }
#consider subsetting df each time, might be much faster



# Script Name: utilities_geospatial_modelling.R
# Purpose: Contains a number of functions useful for geospatial modelling 
# Author: Michael McPhail
# Date: 2024-07-09
# Dependencies: raster, Matrix, dplyr, rgdal, sp, sf
# Usage: souruce(PATH_TO_SCRIPT) in driver script

library(raster)
library(Matrix)
library(dplyr)
library(rgdal)
library(sp)
library(sf) # For loading shape file from Z drive
library(MASS) #Remove if glm is no longer used
library(glmmTMB)

#-----------------------------------------------------#
#Start: List data (file paths, covariate names ...)
#-----------------------------------------------------#
month_idx_list = c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12')
#-----------------------------------------------------#
#Start: Data processing
#-----------------------------------------------------#
#' Build vector of covariate names with lags
#'
#' Takes in the vector of base names and returns the augmented list
#'
#' @param dyn_cov_name_list is a vector of base covariate names
#' @return a vector of covariate names
#' @export
get_cov_names <- function(dyn_cov_name_list, n_lagged=1){
  cov_names_full_list = c(dyn_cov_name_list)
  if(n_lagged > 0){
    for(lag in 1:n_lagged){
      for(displacement in c(paste0('_lag', toString(lag)), paste0('_lead', toString(lag)))){
        for(cov_name in dyn_cov_name_list){
          cov_names_full_list = c(cov_names_full_list, paste0(cov_name, displacement))
        }
      }
    }
  }  
  return(cov_names_full_list)
}

#' Wrangle input data
#'
#' Takes the input data and covariate names, adds Forest indicator and performs some data processing
#'
#' @param input_db Must have columns c('Re', 'LAT', 'LONG', cov_names, 'month_int', 'year' (int), 'mp_species')
#' @param cov_names A vector of covariate names
#' @return A wrangled/subsetted data frame.
#' @export
wrangle_input_data  <- function(input_db, cov_names, parasite){
  input_db<- input_db[!is.na(input_db$LONG)|(!is.na(input_db$LAT)), ]

  uncropped_raster = raster(static_fn_str_list['FOREST'])
  input_db['FOREST'] <- factor(terra::extract(uncropped_raster, cbind(input_db$LONG, input_db$LAT)))
  if(nrow(input_db[input_db$FOREST==3,])>0){
    input_db[input_db$FOREST==3,]$FOREST <- 2
  }

  #Retain only data points with complete covariate information
  for(cov in cov_names){
      input_db <- input_db[!is.na(input_db[,cov]),]
  }
  input_db <- input_db[!is.na(input_db$Re),]

  #Drop first month (often Re is too high because of missing sources for the first month)
  input_db <- input_db[!((input_db$month_int==1)&(input_db$year==2019)),]
  input_db <- input_db[(input_db$year!=2018),]
  # input_db$Re = input_db$Re
  input_db['logR_value'] <- log(input_db$Re+0.000001) 

  #Drop outlier values (there may be none).
  input_db <- input_db[input_db$Re < 10,] 
  # input_db <- input_db[input_db$Re > 0,] 

  
  input_db <- input_db[(input_db$mp_species==parasite),]

  return(input_db)
}

#' Normalise the covariates
#'
#' Normalise the specified list of covariates and replace input_db
#'
#' @param input_db Un-normalised input data.
#' @param cov_names The list of covariates to normalise, rainfall data is transformed to log.
#' @return A normalised data frame.
#' @export
normalise_input_covariates <- function(input_db, cov_names_to_normalise){
  cov_scales_df = data.frame(matrix(nrow=length(cov_names_to_normalise), ncol=3))
  names(cov_scales_df) = c('cov_name', 'mean', 'sd')
  cov_scales_df$cov_name = cov_names_to_normalise
  for(col_name in cov_names_to_normalise){
    if(grepl('CHIRPS', col_name)){
      cov_scales_df[cov_scales_df$cov_name==col_name, 'mean'] = mean(log(input_db[,col_name]))
      cov_scales_df[cov_scales_df$cov_name==col_name, 'sd']  = sd(log(input_db[,col_name]))
      input_db[,col_name] <-(log(input_db[,col_name]) - mean(log(input_db[,col_name])))/sd(log(input_db[,col_name]))
    }else{
      cov_scales_df[cov_scales_df$cov_name==col_name, 'mean'] = mean(input_db[,col_name])
      cov_scales_df[cov_scales_df$cov_name==col_name, 'sd']  = sd(input_db[,col_name])
      input_db[,col_name] <-(input_db[,col_name] - mean(input_db[,col_name]))/sd(input_db[,col_name])
    }
  }
  return(list('input_db' = input_db, 'cov_scales_df' = cov_scales_df))
}

#' Get regression weights
#'
#' Extract the inverse probabilities from input_db and convert to regression weights
#'
#' @param input_db Must have columns c('inv_probability')
#' @param scale_fac A scale factor for transforming weights, has no significant effect (no real effect, only that normalised values can be small for lots of obs)
#' @param truncate_quant At which quantile to truncate weights. Default is c(5%, 95%)
#' @param no_weight Setting to TRUE makes all weights equal (1); effectively no weight.
#' @return data frame with weights in place of inv_probability column
#' @export
get_regression_weights <- function(input_db, scale_fac = 1000, truncate_quant = 0.95, no_weight=FALSE){
  if(no_weight){
    input_db$inv_probability = 1
    return(input_db)
  }
  input_db$inv_probability <- input_db$inv_probability/sum(input_db$inv_probability)
  upper_lim = quantile(input_db$inv_probability, truncate_quant)
  lower_lim = quantile(input_db$inv_probability, 1-truncate_quant)

  input_db[input_db$inv_probability >=upper_lim,]$inv_probability = upper_lim
  input_db[input_db$inv_probability <=lower_lim,]$inv_probability = lower_lim
  #Scale values for better conditioning in optimisation
  input_db$inv_probability <- input_db$inv_probability/mean(input_db$inv_probability)
  
  
  return(input_db)
}


#-----------------------------------------------------#
#End: Data processing

#Start: Vulnerability/LGCM
#-----------------------------------------------------#
#' Get the dual mesh 
#'
#' Takes in a mesh, generates the dual mesh
#'
#' @param mesh an INLA mesh object
#' @return  a dual mesh (SpatialPolygon vector)
#' @examples
#' add_numbers(3, 5)
#' add_numbers(10, 20)
#' @export
book.mesh.dual <- function(mesh) {
  if (mesh$manifold=='R2') {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0) 
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      }
      else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}


#-----------------------------------------------------#
#End: Vulnerability/LGCM 

#Start: Regression modelling functions
#-----------------------------------------------------#
#' Add Two Numbers
#'
#' This function takes two numeric arguments and returns their sum.
#'
#' @param x A numeric value.
#' @param y A numeric value.
#' @return A numeric value representing the sum of \code{x} and \code{y}.
#' @examples
#' add_numbers(3, 5)
#' add_numbers(10, 20)
#' @export
# model_predict <- function(mesh, coordinate_covariate_data, nsamples = 100, model){
# #mesh is an inla mesh object
# #Columns in coordinate_covariate_data are LAT and LONG of prediction points, and associated covariate values
#   A.val <- inla.spde.make.A(mesh,
#                             loc=cbind(coordinate_covariate_data$x, coordinate_covariate_data$y))

#   set.seed(1)
#   samp <- inla.posterior.sample(nsamples, model)

#   pred <- matrix(NA, nrow = dim(A.val)[1], ncol = nsamples)
#   covar_pred <- model$names.fixed[-1] # Remove intercept.
#   covar_col <- as.matrix(coordinate_covariate_data[, covar_pred])
#   k <- length(covar_pred)

#   for (i in 1:nsamples){
#     field <- samp[[i]]$latent[grep('field', rownames(samp[[i]]$latent))]
#     intercept <- samp[[i]]$latent[grep('Intercept', rownames(samp[[i]]$latent))]
#     beta <- rep(NA, k)
#     for (j in 1:k){
#       #NOTE: MM changed here to take substring -2 characters from row names because it seems to add ':1' to the end of covariates
#       beta[j] <- samp[[i]]$latent[substr(rownames(samp[[i]]$latent),1,nchar(rownames(samp[[i]]$latent))-2) == covar_pred[j], ] # fixed = TRUE to avoid problems with lags and without.
#     }
#     pred[, i] <-  intercept + c(covar_col%*%beta) + drop(A.val%*%field)
#   }

#   return(pred)

# }

#' Get absolute value of correlation
#'
#' Get absolute value of correlation
#'
#' @param cov_name covariate name  
#' @param input_db data frame of covariate values (must have columns c('cov_name','Re'))
#' @return the absolute value of the correlation between Re and covariate
#' @export
get_abs_corr <- function(cov_name,input_db){
  return(abs(cor(input_db[,cov_name],input_db$Re )))
}

#' get AIC
#'
#' Gets the Aikake Information  (AIC) value of a proposed model.
#'
#' @param retained_covariates a vector of covariate names (possibly some subset of all available covariates).
#' @param input_db data frame of covariate values (must have columns c('Re', retained_covariates)).
#' @return An AIC value
#' @export
get_aic<- function(retained_covariates, input_db){
    model <- as.formula(paste("Re ~ ",paste(retained_covariates, collapse = " + ")))

    # fitted_model <- glm(model, data = input_db, family=Gamma(link='log'), weights=input_db$inv_probability)
    
    fitted_model <- glmmTMB(model, 
                            ziformula = ~ ., 
                            family = ziGamma(link='log'), 
                            data = input_db,
                            weights=inv_probability)
    summary = summary(fitted_model)
    return(unname(summary$AICtab['AIC']))
    # return(fitted_model$aic)
}

#' Feature selection
#'
#' Use a greedy-method-based wrapper to perform feature selection. Remove based on AIC improvement.
#'
#' @param input_db data frame of covariate values (must have columns c('Re', retained_covariates))
#' @param retained_covariates a vector of covariate names (possibly some subset of all available covariates).
#' @return A subset of the original covariate list corresponding to the approximate best model.
#' @export
feature_select_greedy_wrapper <- function(input_db, retained_covariates){
  best_aic <- get_aic(retained_covariates, input_db)
  try_drop_each_cov <- TRUE
  message(paste0('Starting with ', toString(length(retained_covariates)), ' Covariates'))
  while(try_drop_each_cov){
      cov_to_drop <- 'None'
      AIC_was_improved <- FALSE
      for(covariate in retained_covariates){
         candidate_covariates <- retained_covariates[retained_covariates!=covariate]    
         new_aic <- get_aic(candidate_covariates, input_db)
         if(new_aic<best_aic){
             cov_to_drop <- covariate
             best_aic = new_aic
             AIC_was_improved<-TRUE
         }
         message(new_aic)
          
      }
      if(AIC_was_improved){
          retained_covariates <- retained_covariates[retained_covariates!=cov_to_drop] 
      }else{
          try_drop_each_cov <- FALSE
      }
  }
  return(retained_covariates)
}


#-----------------------------------------------------#
#End:  Regression modelling functions

#Start: Prediction functions
#-----------------------------------------------------#

simulate_Re  <- function(mu, alpha = 1/0.5, nsims = 100){
  if(is.na(mu)){
    uci=NA
    lci=NA
  }else{
    beta <- alpha / mu
    simulated_vals = rgamma(nsims, shape = alpha, rate = beta)
    uci = quantile(simulated_vals, 0.975) %>% unname()
    lci = quantile(simulated_vals, 0.025) %>% unname()
  }
  return(c(uci, lci))
}

#' Predict single year
#'
#' Makes a country-wide prediction for a specified year (all 12 months).
#'
#' @param fitted_model A fitted glm() model.
#' @param year_int the specified year (integer).
#' @param cov_names the vector of covariate values.
#' @param cov_scales_df a data frame of covariate mean/sd values for transformation.
#' @param country_extent the approximate limits of the country in question for raster cropping.
#' @param country_iso the iso3 code of the target country (e.g. VNM for Vietnam) for shape file subsetting.
#' @param max_cov_year the maximum year for which we have covariate values.
#' @return A data frame with the predictions, and corresponding covariates.
#' @export
predict_single_year <- function(fitted_model, year_int, cov_names, cov_scales_df, country_extent = extent(c(100, 110, 5, 25)), country_iso = 'VNM',  max_cov_year = 9999){
  
  phi = sigma(fitted_model)^2
  
  #Load and subset shape file, get prediction points.
  full_sf <- st_read(paste0(shape_file_paths$prefix, '0', shape_file_paths$suffix)) 
  full_shp <- as(full_sf, Class = 'Spatial')
  sub_shp <-  full_shp[full_shp$ISO %in% c(country_iso),]
  sub_sf <- sub_shp
  poly.shapefile <- st_as_sf(sub_sf) 
  
  #Get coordinates for prediction
  if((gsub('_lead.', '', gsub('_lag.', '', 'CHIRPS')) %in% c('LST_day', 'LST_delta', 'LST_night')) & (year_int>2021)) {
    cropped_raster <- raster(paste0(dyn_root_post2021_str_list['CHIRPS'], toString(year_int),'.01',dyn_suff_str_list['CHIRPS']))%>% crop(country_extent)
  }else{
    cropped_raster <- raster(paste0(dyn_root_str_list['CHIRPS'], toString(year_int),'.01',dyn_suff_str_list['CHIRPS']))%>% crop(country_extent)
  }
  buffer.image <- cropped_raster %>% mask(poly.shapefile %>% st_zm())
  buffer.values   <- getValues(buffer.image)
  valid.pix     <- which(!is.na(buffer.values) & buffer.values >= 0)
  valid.coords     <- coordinates(buffer.image)[valid.pix,]
  valid_coord_df <- as.data.frame(valid.coords)

  #Get predictions for all months
  for(month_int in 1:12){
    for(cov_name in cov_names[!grepl('FOREST', cov_names)]){
        base_cov_name <-  gsub('_lead.', '', gsub('_lag.', '', cov_name))
        month_of_interest <- month_int
        year_of_interest <- year_int

        if(grepl('lag', cov_name)){
            offset<- substr(cov_name,nchar(cov_name), nchar(cov_name))%>% strtoi()
            month_of_interest <- month_of_interest - offset
        }else if(grepl('lead', cov_name)){
            offset<- substr(cov_name,nchar(cov_name), nchar(cov_name))%>% strtoi()
            month_of_interest <- month_of_interest + offset
        }
        # Connect Dec to Jan
        if(month_of_interest > 12){
            year_of_interest <- year_of_interest + 1
            if(year_of_interest > max_cov_year){
                year_of_interest <- max_cov_year
            }
            month_of_interest<- (month_of_interest-1)%%12 + 1
        }else if (month_of_interest < 1){
            year_of_interest <- year_of_interest - 1
            month_of_interest<- (month_of_interest -1 + 12)%%12
        }

        month_str <- month_idx_list[month_of_interest]
        year_str<- toString(year_of_interest)

        if(base_cov_name %in% dyn_cov_name_list){
            if((base_cov_name %in% c('LST_day', 'LST_delta', 'LST_night')) & (year_of_interest>2021)) {# a temporary fix because 2022 LST values are in a different directory
              cov_fn = paste0(dyn_root_post2021_str_list[base_cov_name], year_str, '.', month_str,dyn_suff_str_list[base_cov_name])
            }else{
              cov_fn = paste0(dyn_root_str_list[base_cov_name], year_str, '.', month_str, dyn_suff_str_list[base_cov_name])
            }
        }else if(base_cov_name %in% annual_root_str_list){
            cov_fn = paste0(annual_root_str_list[base_cov_name], year_str, annual_suff_str_list[base_cov_name])
            if((year > 2020)&(cov_name=='POPULATION')){
                cov_fn <- paste0(annual_root_str_list['POPULATION'], '2020', annual_suff_str_list['POPULATION']) 
            }
        }else if(base_cov_name %in% static_fn_str_list){
            cov_fn <- static_fn_str_list[base_cov_name]
        }else{
            message('not a valid covariate')
            break
        }

        cropped_raster <-cov_fn %>% raster() %>% crop(country_extent)

        valid_coord_df[cov_name] <- terra::extract(cropped_raster, valid.coords)    

    }

    if('FOREST' %in% cov_names){
      cropped_raster <- raster(static_fn_str_list['FOREST'])%>% crop(country_extent)
      valid_coord_df['FOREST'] <- factor(terra::extract(cropped_raster, valid.coords))
      valid_coord_df[(valid_coord_df$FOREST==3)&(!is.na(valid_coord_df$FOREST)),]$FOREST <-2
    }
    #Transform covariates appropriately (consistent with training data)
    for(col_name in cov_names[!grepl('FOREST', cov_names)]){
      mean_val = cov_scales_df[cov_scales_df$cov_name==col_name, 'mean']
      sd_val = cov_scales_df[cov_scales_df$cov_name==col_name, 'sd']
      if(grepl('CHIRPS', col_name)){
        valid_coord_df[,col_name] <-(log(valid_coord_df[,col_name]) - mean_val)/sd_val
      }else{
        valid_coord_df[,col_name] <-(valid_coord_df[,col_name] - mean_val)/sd_val
      }
      
    }
    if(!('inv_probability' %in% (valid_coord_df %>% names()))){
      valid_coord_df$inv_probability = 1
    }
    
    rows_with_sefit = rowSums(is.na(valid_coord_df[,cov_names]))==0
    #Get predictions
    predictions <- predict(fitted_model, valid_coord_df, type='response')
    zi_predictions <- predict(fitted_model, valid_coord_df, type='zprob')
    cont_predictions = predict(fitted_model,newdata = valid_coord_df, type = 'conditional')
    
    valid_coord_df[paste0('prediction_', month_idx_list[month_int])] <- predictions 
    valid_coord_df[paste0('prediction_zi_', month_idx_list[month_int])] <- zi_predictions
    
    valid_coord_df[paste0('prediction_se_', month_idx_list[month_int])]  = NA
    valid_coord_df[paste0('prediction_zi_se_', month_idx_list[month_int])] <- NA
    predictions <- predict(fitted_model, valid_coord_df[rows_with_sefit,], type='response', se.fit=TRUE)
    zi_predictions <- predict(fitted_model, valid_coord_df[rows_with_sefit,], type='zprob', se.fit=TRUE)
    valid_coord_df[rows_with_sefit, paste0('prediction_se_', month_idx_list[month_int])] <- predictions$se.fit %>% unname()
    valid_coord_df[rows_with_sefit, paste0('prediction_zi_se_', month_idx_list[month_int])] <- zi_predictions$se.fit %>% unname()
    
    simulated_cis = sapply(cont_predictions, simulate_Re, alpha = 1/phi, nsims=1000) %>% t()  
    valid_coord_df[paste0('prediction_pi_uci_', month_idx_list[month_int])] = simulated_cis[,1]*(1-valid_coord_df[paste0('prediction_zi_', month_idx_list[month_int])])
    valid_coord_df[paste0('prediction_pi_lci_', month_idx_list[month_int])] = simulated_cis[,2]*(1-valid_coord_df[paste0('prediction_zi_', month_idx_list[month_int])])
  }
  return(valid_coord_df[,-which(names(valid_coord_df) %in% c('inv_probability'))])
}



add_prediction_intervals <- function(fitted_model, valid_coord_df){
  valid_coord_df$inv_probability = 1
  phi = sigma(fitted_model)^2
  for(month_int in 1:12){
    message(month_int)
    #TODO change to conditional prediction (rather than condition*zi)
    # cont_predictions = predict(fitted_model,newdata = valid_coord_df, type = 'conditional')
    cont_predictions = valid_coord_df[[paste0('prediction_', month_idx_list[month_int])]]/(1-valid_coord_df[[paste0('prediction_zi_', month_idx_list[month_int])]]) 
    simulated_cis = sapply(cont_predictions, simulate_Re, alpha = 1/phi, nsims=1000) %>% t()  
    valid_coord_df[paste0('prediction_pi_uci_', month_idx_list[month_int])] = simulated_cis[,1]*(1-valid_coord_df[paste0('prediction_zi_', month_idx_list[month_int])])
    valid_coord_df[paste0('prediction_pi_lci_', month_idx_list[month_int])] = simulated_cis[,2]*(1-valid_coord_df[paste0('prediction_zi_', month_idx_list[month_int])])
    
  }
  return(valid_coord_df[,-which(names(valid_coord_df) %in% c('inv_probability'))])
}
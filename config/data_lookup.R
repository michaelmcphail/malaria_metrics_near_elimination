#Covariates to include with monthly lags
dyn_cov_name_list = c('CHIRPS')
# dyn_cov_name_list = c('LST_day', 'LST_delta', 'LST_night')
#Covariates that never change
static_cov_name_list = c('ELEVATION') #Easier to add forrest factor in the driver script # treat population as static (use 2020)
#Covariates that change annually (must include population for population weighting)
annual_cov_name_list = c('CHIRPS_synoptic')

# look up table data structures, will depend on specific machines access to data.
shape_file_paths = list('prefix'='SHAPE_NAME_PREFIX', 
                        'suffix'='SHAPE_NAME_SUFFIX')

# dynamic covariates

dyn_root_str_list = list('CHIRPS' = 'PATH_TO_COVARIATE_PREFIX')

dyn_suff_str_list = c('CHIRPS' = 'COVARIATE_SUFFIX')

annual_root_str_list = c('CHIRPS_synoptic' = 'PATH_TO_COVARIATE_PREFIX')

annual_suff_str_list = c('CHIRPS_synoptic' = 'COVARIATE_SUFFIX')

static_fn_str_list = c('ELEVATION' = 'PATH_TO_COVARIATE') 


# Script Name: my_script.R
# Purpose: Brief description of what the script does.
# Author: Your Name
# Date: YYYY-MM-DD
# Dependencies: dplyr, ggplot2
# Usage: Rscript my_script.R input_file.csv output_file.csv

library(sf)
library(spocc)
library(rnaturalearth)
library(sf)
library(terra)
library(INLA)
library(rgeos)
library(raster)
library(fields) 
library(optparse)
library(jsonlite)



#TODO: reference template script

load_config <- function(config_path) {
  config <- fromJSON(config_path)
  return(config)
}

option_list <- list(
  make_option(c("--config"), type="character", default=NULL, 
              help="Path to the configuration file", metavar="character"),
  make_option(c("--all_index_cases"), type="logical", default=FALSE, 
              help="Set to true to include orphaned cases from network diffusion model")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$config)) {
  print_help(opt_parser)
  stop("The --config argument is required", call.=FALSE)
}


config <- load_config(opt$config)
all_index_cases <- opt$all_index_cases

data_in_fn = config$case_data_withVulnCovs_fn

if(all_index_cases){
  pf_vulnerability_raster_fn = config$pf_vulnerability_index_raster_fn
  pv_vulnerability_raster_fn = config$pv_vulnerability_index_raster_fn
}else{
  pf_vulnerability_raster_fn = config$pf_vulnerability_imported_raster_fn
  pv_vulnerability_raster_fn = config$pv_vulnerability_imported_raster_fn
}

country_name = config$country_name


parasite_list <- c('P.falciparum', 'P.vivax')

source('./src/utils/utilities_geospatial_modelling.R')
source('./config/data_lookup.R')

projUTM <-'+proj=tmerc +lat_0=0 +lon_0=104 +k=0.9999 +x_0=500000 +y_0=0 +ellps=WGS84 +towgs84=-191.90441429,-39.30318279,-111.45032835,0.00928836,-0.01975479,0.00427372,0.252906278 +units=km +no_defs'

for(parasite in parasite_list){
  message(paste0('Estimating vulnerability for: ', parasite))

  #------------------------ 
  input_db <- read.csv(data_in_fn)
  input_db<- input_db[!is.na(input_db$LONG)|(!is.na(input_db$LAT)), ]
      
  if(all_index_cases){    
      point_conditions <- (input_db$is_index_case ==1)&(input_db$mp_species == parasite)
  }else{
      point_conditions <- (input_db$imported =='Y')&(input_db$mp_species == parasite) 
  }
  input_db<- input_db[point_conditions,]
  
  unique_coords_obs <- input_db[c('LAT', 'LONG')] %>% unique()

  cov_fn = static_fn_str_list['POPULATION']
  cov_rast <- terra::rast(cov_fn)
  
  obs_points_sf <- st_as_sf(unique_coords_obs, coords = c("LONG", "LAT"),
                crs = st_crs(cov_rast)) %>% st_transform(crs=projUTM) %>% st_transform()
  
  map <- ne_countries(type = "countries", country = country_name,
                      scale = "medium", returnclass = "sf") %>% st_transform(crs=projUTM)
  # map <- st_transform(map, crs = projUTM)
  
  
  # raster grid covering map
  grid <- terra::rast(map, nrows = 250, ncols = 250)
  # coordinates of all cells
  xy <- terra::xyFromCell(grid, 1:ncell(grid))
  
  # transform points to a sf object
  grid_points_sf <- st_as_sf(as.data.frame(xy), coords = c("x", "y"),
                  crs = st_crs(map))
  
  # indices points within the map
  points_within_idx <- which(st_intersects(grid_points_sf, map,
                                              sparse = FALSE))
  
  # points within the map
  grid_points_sf <- st_filter(grid_points_sf, map)
  grid_coords_sf <- st_coordinates(grid_points_sf)
  
  domain_points <- cbind(st_coordinates(map)[, 1], st_coordinates(map)[, 2])
  mesh <- inla.mesh.2d(loc.domain = domain_points, max.edge = c(50, 100),
                        offset = c(50, 100), cutoff = 1)
  
  mesh_n <- mesh$n
  points_n <- nrow(obs_points_sf)
  
  spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
  
  dual_mesh <- book.mesh.dual(mesh)
  
  # Domain polygon is converted into a SpatialPolygons
  domain.polys <- Polygons(list(Polygon(domain_points)), '0')
  domain_sp <- SpatialPolygons(list(domain.polys))
  domain_sp2 <- gBuffer(domain_sp, byid=TRUE, width=0)
  domain_sp2 <- gSimplify(domain_sp2, tol = 0.00001)
  
  # Because the mesh is larger than the study area, we need to
  # compute the intersection between each polygon
  # in the dual mesh and the study area
  #TODO: speed up by parallelising
  w <- sapply(1:length(dual_mesh), function(i) {
    if (gIntersects(dual_mesh[i, ], domain_sp2))
      return(gArea(gIntersection(dual_mesh[i, ], domain_sp2)))
    else return(0)
  })
  
  y.pp <- rep(0:1, c(mesh_n, points_n))
  e.pp <- c(w, rep(0, points_n))
  
  # Projection matrix for the integration points (mesh vertices)
  A.int <- Diagonal(mesh_n, rep(1, mesh_n))
  # Projection matrix for observed points (event locations)
  A.y <- inla.spde.make.A(mesh = mesh, loc = obs_points_sf)
  # Projection matrix for mesh vertices and event locations
  A.pp <- rbind(A.int, A.y)
  Ap.pp <- inla.spde.make.A(mesh = mesh, loc = grid_points_sf)
  
  
  mesh_points <- cbind(mesh$loc[,1],mesh$loc[,2])
  mesh_coords <- as.data.frame(mesh_points)
  names(mesh_coords)<-c('x', 'y')
  mesh_coords_sf <- st_as_sf(mesh_coords, coords = c("x", "y"),
                              crs = st_crs(map))
  
  # stack for estimation
  stk.e.pp <- inla.stack(tag = "est.pp",
                          data = list(y = y.pp, e = e.pp),
                          A = list(1, A.pp),
                          effects = list(list(b0 = rep(1, mesh_n + points_n)),
                                        list(s = 1:mesh_n)))
  
  # stack for prediction stk.p
  stk.p.pp <- inla.stack(tag = "pred.pp",
                          data = list(y = rep(NA, nrow(grid_coords_sf)), e = rep(0, nrow(grid_coords_sf))),
                          A = list(1, Ap.pp),
                          effects = list(data.frame(b0 = rep(1, nrow(grid_coords_sf))),
                                        list(s = 1:mesh_n)))
  formula <- y ~ 0 + b0  + f(s, model = spde)
  
  # stk.full has stk.e and stk.p
  stk.full.pp <- inla.stack(stk.e.pp, stk.p.pp)
  
  
  
  res <- inla(formula,  family = 'poisson',
              data = inla.stack.data(stk.full.pp),
              control.predictor = list(compute = TRUE, link = 1,
                                        A = inla.stack.A(stk.full.pp)),
              E = inla.stack.data(stk.full.pp)$e)
  
  index <- inla.stack.index(stk.full.pp, tag = "pred.pp")$data
  pred_mean <- res$summary.fitted.values[index, "mean"]
  pred_ll <- res$summary.fitted.values[index, "0.025quant"]
  pred_ul <- res$summary.fitted.values[index, "0.975quant"]
  
  grid$mean <- NA
  grid$ll <- NA
  grid$ul <- NA
  
  grid$mean[points_within_idx] <- pred_mean
  grid$ll[points_within_idx] <- pred_ll
  grid$ul[points_within_idx] <- pred_ul
  
  #Load and subset shape file
  full_sf <- st_read(paste0(shape_file_paths$prefix, '0', shape_file_paths$suffix)) 
  full_shp <- as(full_sf, Class = 'Spatial')
  sub_sf <-  full_shp[full_shp$ISO %in% c('VNM'),]
  poly.shapefile <- st_as_sf(sub_sf) 
  
  grid_as_df = as.data.frame(grid, xy = TRUE)
  points_to_plot <- st_as_sf(grid_as_df, coords = c("x", "y"),crs = st_crs(map)) %>% st_transform(crs=st_crs(sub_sf))
  
  points_to_plot$x <- st_coordinates(points_to_plot)[,'X']
  points_to_plot$y <- st_coordinates(points_to_plot)[,'Y']
  
  start.time <- Sys.time()
  message('Downscaling and saving')
  raster_template <- raster(extent(points_to_plot), res = 0.01)
  crs(raster_template) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  message('Building interpolator')
  message(Sys.time() - start.time)
  spline <- Tps(cbind(points_to_plot$x, points_to_plot$y), points_to_plot$mean)
  message('Interpolating')
  message(Sys.time() - start.time)
  splined = interpolate(raster_template, spline)
  
  message('Saving')
  message(Sys.time() - start.time)
  
  if(parasite=='P.falciparum'){
    writeRaster(splined, pf_vulnerability_raster_fn, overwrite=TRUE)
  }else{
    writeRaster(splined, pv_vulnerability_raster_fn, overwrite=TRUE)
  }
}

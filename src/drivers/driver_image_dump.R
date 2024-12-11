# Script Name: driver_driver_image_dump.R
# Purpose: Dumps a bunch of output images for the paper 
# Author: Michael McPhail 
# Date: 2024-07-16
# Dependencies: ggplot2, raster, dplyr, rgdal, sp, optparse, jsonlite
# Usage: Rscript driver_image_dump.R --config='PATH_TO_CONFIG_FILE'
rm(list = ls())
library(ggplot2)
library(raster)
library(dplyr)
library(rgdal)
# library(sp)
library(sf)
library(scales)
library(ggrepel)
library(cowplot)
library(ggpubr)#For saving legend
library(optparse)
library(jsonlite)
library(tmaptools)

# options(warning=1)
#---------------------------------------------------------------------------------------------------------------#
#             Start: Configuration setup
#---------------------------------------------------------------------------------------------------------------#
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


midline_decimal <- function(x) { gsub("\\.", "·", format(x)) }

# "#FFB562"
custom_gradient <- scale_fill_gradientn(
  colours = c("#F9F2ED",  "#3AB0FF","#0C359E", "#011f4b",  "#7CFC00","red"),
  values = scales::rescale(c(0, 0.4,0.8, 1, 1.001,2)),  # specify the values where colors change
  limits = c(0, 2),
  name=NULL,
  # name=expression(R[e]),
  labels=midline_decimal
)

custom_gradient_truncated <- scale_fill_gradientn(
  colours = c("#F9F2ED",  "#3AB0FF","#0C359E", "#011f4b"),
  values = scales::rescale(c(0, 0.4,0.8, 1)),  # specify the values where colors change
  limits = c(0, 1),
  name=NULL,
  # name=expression(R[e]),
  labels=midline_decimal
)

custom_gradient_zi <- scale_fill_gradientn(
  colours = c("red", "#7CFC00"),
  values = scales::rescale(c(0, 1)),  # specify the values where colors change
  limits = c(0, 1),
  name=NULL,
  # name=expression(paste("P[R"[e], "=0]")),
  labels=midline_decimal
)

#---------------------------------------------------------------------------------------------------------------#
#             End: Configuration setup
#             Start: Configuration
#---------------------------------------------------------------------------------------------------------------#
config <- load_config(opt$config)
Re_prediction_min_year = config$Re_prediction_min_year
Re_prediction_max_year = config$Re_prediction_max_year


if(!dir.exists(config$image_out_root)){
  dir.create(config$image_out_root)
}

parasite_list <- c('P.falciparum', 'P.vivax')
study_year_list <- seq(Re_prediction_min_year, Re_prediction_max_year)
year_str_list = sapply(study_year_list, toString)

country_extent <- extent(config$country_extent)
country_iso = config$country_iso
surrounding_country_isos = config$surrounding_country_isos

vietnam_subregions <- list(
  "Red River Delta" = c("Ha Noi", "Hai Phong", "Quang Ninh", "Hai Duong", "Hung Yen", "Hai Phong", "Bac Ninh", "Ha Nam", "Vinh Phuc", "Thai Binh", "Nam Dinh", "Ninh Binh", 'TP. Ha Noi', 'TP. Hai Phon'),
  "Northeast" = c("Cao Bang", "Bac Kan", "Lang Son", "Ha Giang", "Tuyen Quang", "Quang Ninh", "Dien Bien", "Lai Chau", "Lao Cai", "Yen Bai", "Phu Tho", "Bac Can", 'Bac Giang', 'Thai Nguyen'),
  "North Central Coast" = c("Thanh Hoa", "Nghe An", "Ha Tinh", "Quang Binh", "Quang Tri", "Thua Thien H"),
  "Northwest" = c("Dien Bien", "Lai Chau", "Lao Cai", "Yen Bai", "Hoang Lien Son", 'Hoa Binh','Son La'),  
  "South Central Coast" = c("TP. Da Nang", "Quang Nam", "Quang Ngai", "Binh Dinh", "Phu Yen", "Khanh Hoa", "Ninh Thuan", "Binh Thuan"),
  "Central Highlands" = c("Kon Tum", "Gia Lai", "Dak Lak", "Dak Nong", "Lam Dong"),
  "Southeast" = c("TP. Ho Chi M", "Ba Ria-Vung", "Binh Duong", "Binh Phuoc", "Dong Nai", "Tay Ninh"),
  "Mekong Delta" = c("TP. Can Tho", "Long An", "Tien Giang", "Ben Tre", "Dong Thap", "Vinh Long", "Tra Vinh", "An Giang", "Kien Giang", "Hau Giang", "Soc Trang", "Bac Lieu", "Ca Mau")
)

#Load and subset shape files
full_sf <- st_read(paste0('/mnt/Z/master_geometries/Admin_Units/Global/MAP/2022/admin2022_0.shp')) 
full_shp <- as(full_sf, Class = 'Spatial')
sub_sf <-  full_shp[full_shp$ISO %in% c(country_iso),]
poly.shapefile <- st_as_sf(sub_sf) 

subplot_sf <- full_sf[full_sf$ISO %in% c(country_iso,surrounding_country_isos),] 
asp_ratio = get_asp_ratio(subplot_sf)

full_sf1 <- st_read(paste0('/mnt/Z/master_geometries/Admin_Units/Global/MAP/2022/admin2022_1.shp')) 
admin1_sf <- st_as_sf(full_sf1[full_sf1$ISO %in% c(country_iso),])%>% st_make_valid()

full_sf2 <- st_read(paste0('/mnt/Z/master_geometries/Admin_Units/Global/MAP/2022/admin2022_2.shp')) 
admin2_sf <- st_as_sf(full_sf2[full_sf2$ISO %in% c(country_iso),])%>% st_make_valid()

admin1_copy = full_sf1[full_sf1$ISO %in% c(country_iso),] %>% st_make_valid()
admin1_copy$Region = 'NA'
for(name in vietnam_subregions%>% names()){
  admin1_copy[admin1_copy$Name_1 %in% vietnam_subregions[[name]],'Region'] <- name
}


nyears = 4 

yyyy_mm_float_list = c()
for(year in Re_prediction_min_year:Re_prediction_max_year){
  for(month_int in 1:12){
    yyyy_mm_float_list = c(yyyy_mm_float_list, year + month_int/100)
  }
}

# get_legend <- function(myplot) {
#   tmp <- ggplotGrob(myplot)
#   legend <- tmp$grobs[which(sapply(tmp$grobs, function(x) x$name) == "guide-box")]
#   legend
# }

get_legend_custom <- function(myplot) {
  tmp <- ggplotGrob(myplot)
  legend <- gtable::gtable_filter(tmp, "guide-box")
  return(legend)
}


parasite = 'P.falciparum'
# parasite_list = c('P.vivax')
for(parasite in parasite_list){ 
  #-------------------------------------------
  # Load data
  case_data_withRe_fn = read.csv(config$case_data_withRe_fn)
  case_data_withRe_fn = case_data_withRe_fn[(case_data_withRe_fn$mp_species==parasite)&(!grepl('2018', case_data_withRe_fn$date_onset)),]
  case_data_withRe_fn$yyyy_mm_float = case_data_withRe_fn$year + case_data_withRe_fn$month_int/100
  
  case_data_base = read.csv(config$case_data_base_fn)
  case_data_base = case_data_base[(case_data_base$mp_species == parasite)&(!grepl('2018',case_data_base$date_onset)),]
  case_data_base$yyyy_mm_float = 2019.01 #Initialise
  yyyy_mm_str_list = c()
  value_counts = c()
  mm_str_list = c('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12')
  for(year in Re_prediction_min_year:Re_prediction_max_year){
    for(month_int in 1:12){
      new_yyyy_mm_str = paste0(year, '-', mm_str_list[month_int])
      yyyy_mm_str_list = c(yyyy_mm_str_list, new_yyyy_mm_str)
      value_counts = c(value_counts, nrow(case_data_base[grepl(new_yyyy_mm_str,case_data_base$date_onset),]))
      value_counts = c(value_counts, nrow(case_data_base[grepl(new_yyyy_mm_str,case_data_base$date_onset),]))
      case_data_base[grepl(new_yyyy_mm_str,case_data_base$date_onset),]$yyyy_mm_float=year + month_int/100
    }
  }
  
  Re_data_with_covs = read.csv(config$case_data_withReCovs_fn)
  Re_data_with_covs = Re_data_with_covs[(Re_data_with_covs$mp_species==parasite)&(!grepl('2018', Re_data_with_covs$date_onset)),]
  Re_data_with_covs$relative_date = Re_data_with_covs$relative_date - 7
  
  
  #Load an example data frame to build on
  if(parasite=='P.falciparum'){
    seasonal_average_df <- read.csv(paste0(config$Re_prediction_pf_prefix, year_str_list[1], config$Re_prediction_suffix))
    abs_max<- 4 #zlim for receptivity
  }else{
    seasonal_average_df <- read.csv(paste0(config$Re_prediction_pv_prefix, year_str_list[1], config$Re_prediction_suffix))
    abs_max<- 1.1 #zlim for receptivity
  }
  prediction_names <- seasonal_average_df %>% names()
  zi_names <- prediction_names[grepl('prediction_zi', prediction_names)]
  zi_names <- zi_names[!grepl('prediction_zi_se_', zi_names)]
  
  prediction_names <- prediction_names[grepl('prediction', prediction_names)]
  prediction_names <- prediction_names[!grepl('zi', prediction_names)]
  prediction_names <- prediction_names[!grepl('prediction_se_', prediction_names)]
  prediction_names <- prediction_names[!grepl('prediction_pi', prediction_names)]
  
  seasonal_average_preds_df <- seasonal_average_df[, c('x', 'y', prediction_names)]
  seasonal_average_zi_df <- seasonal_average_df[, c('x', 'y', zi_names)]
  count <- 1
  for(year_str in year_str_list[2:length(year_str_list)]){
    if(parasite=='P.falciparum'){
      new_df <- read.csv(paste0(config$Re_prediction_pf_prefix, year_str, config$Re_prediction_suffix))
    }else{
      new_df <- read.csv(paste0(config$Re_prediction_pv_prefix, year_str, config$Re_prediction_suffix))
    }
    seasonal_average_preds_df[,prediction_names] = (count*seasonal_average_preds_df[,prediction_names] + new_df[, prediction_names])/(count + 1)
    seasonal_average_zi_df[,zi_names] = (count*seasonal_average_zi_df[,zi_names] + new_df[, zi_names])/(count + 1)
    count <- count + 1
  }
  
  for(prediction_name in prediction_names){
    if (nrow(seasonal_average_preds_df[(!is.na(seasonal_average_preds_df[,prediction_name]))&seasonal_average_preds_df[,prediction_name]>abs_max,])>0){
      seasonal_average_preds_df[(!is.na(seasonal_average_preds_df[,prediction_name]))&seasonal_average_preds_df[,prediction_name]>abs_max,prediction_name] = abs_max
    }
  }
  
  for(zi_name in zi_names){
    if (nrow(seasonal_average_zi_df[(!is.na(seasonal_average_zi_df[,zi_name]))&seasonal_average_zi_df[,zi_name]>abs_max,])>0){
      seasonal_average_zi_df[(!is.na(seasonal_average_zi_df[,zi_name]))&seasonal_average_zi_df[,zi_name]>abs_max,zi_name] = abs_max
    }
  }
  
  seasonal_average_preds_df$annual_mean <- seasonal_average_preds_df[,prediction_names] %>% rowMeans()
  seasonal_average_preds_df$annual_max <- apply(seasonal_average_preds_df[, prediction_names], 1, function(x) max(c(x[!is.na(x)],0)))
  seasonal_average_zi_df$annual_mean <- seasonal_average_zi_df[,zi_names] %>% rowMeans()
  
  seasonal_average_preds_df = seasonal_average_preds_df[!is.na(seasonal_average_preds_df$annual_mean), ]
  if (length(seasonal_average_preds_df[seasonal_average_preds_df$annual_mean>abs_max,]$annual_mean)>0){
    seasonal_average_preds_df[seasonal_average_preds_df$annual_mean>abs_max,]$annual_mean = abs_max
  }
  

  # Load vulnerability raster    
  if(parasite=='P.falciparum'){
    vuln_raster = raster(config$pf_vulnerability_imported_raster_fn)
    index_raster = raster(config$pf_vulnerability_index_raster_fn)
  }else{
    vuln_raster = raster(config$pv_vulnerability_imported_raster_fn)
    index_raster = raster(config$pv_vulnerability_index_raster_fn)
  }
  seasonal_average_preds_df$vulnerability <- terra::extract(vuln_raster, 
                                                            seasonal_average_preds_df[,c('x', 'y')])/nyears
  seasonal_average_preds_df$vulnerability_index <- terra::extract(index_raster, 
                                                            seasonal_average_preds_df[,c('x', 'y')])/nyears
  
  # Calculate malariogenic potential
  seasonal_average_preds_df$malariogenic <- seasonal_average_preds_df$vulnerability*seasonal_average_preds_df$annual_mean
  seasonal_average_preds_df$malariogenic_index <- seasonal_average_preds_df$vulnerability_index*seasonal_average_preds_df$annual_mean
  
  
  pop_rast = paste0('/mnt/Z/mastergrids/Other_Global_Covariates/Population/WorldPop/v3/PopulationCounts_DRC_fixed/5km/WorldPop_UNAdj_v3_DRC_fix.', '2020', '.Annual.Data.5km.sum.tif') %>% raster()
  coords_sf = st_as_sf(seasonal_average_preds_df[, c('x', 'y')], coords = c("x", "y"), crs = st_crs(admin1_sf))
  seasonal_average_preds_df$population = terra::extract(pop_rast, coords_sf)
  breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6)
  cumulative_mean_pop = c()
  for(i in 1:(length(breaks)-2)){
    conds = (seasonal_average_preds_df$annual_mean>=breaks[i])&(seasonal_average_preds_df$annual_mean<breaks[i+1])&(!is.na(seasonal_average_preds_df$annual_mean))&(!is.na(seasonal_average_preds_df$population))
    cumulative_mean_pop = c(cumulative_mean_pop, sum(seasonal_average_preds_df[conds, ]$population))
  }
  conds = (seasonal_average_preds_df$annual_mean>=breaks[length(breaks)-1])&(!is.na(seasonal_average_preds_df$annual_mean))&(!is.na(seasonal_average_preds_df$population))
  cumulative_mean_pop = c(cumulative_mean_pop, sum(seasonal_average_preds_df[conds, ]$population))
  
  
  cumulative_max_pop = c()
  for(i in 1:(length(breaks)-2)){
    conds = (seasonal_average_preds_df$annual_max>=breaks[i])&(seasonal_average_preds_df$annual_max<breaks[i+1])&(!is.na(seasonal_average_preds_df$annual_max))&(!is.na(seasonal_average_preds_df$population))
    cumulative_max_pop = c(cumulative_max_pop, sum(seasonal_average_preds_df[conds, ]$population))
  }
  conds = (seasonal_average_preds_df$annual_max>=breaks[length(breaks)-1])&(!is.na(seasonal_average_preds_df$annual_max))&(!is.na(seasonal_average_preds_df$population))
  cumulative_max_pop = c(cumulative_max_pop, sum(seasonal_average_preds_df[conds, ]$population))
  
  pop_total_df = rbind(data.frame(breaks=breaks[2:length(breaks)], pop = cumulative_mean_pop, metric=rep('Mean', length(breaks)-1)),
                       data.frame(breaks=breaks[2:length(breaks)], pop = cumulative_max_pop, metric=rep('Max', length(breaks)-1)))
  
  # ggplot(pop_total_df, aes(fill=metric, y=pop, x=breaks)) + 
  #   geom_bar(position="dodge", stat="identity")
  
  #-------------------------------------------
  # Plots -- Region map
  ggplot(data = admin1_copy) + 
    geom_sf(aes(fill=Region))
  
  ggsave(paste0(config$image_out_root, 'region_map', config$image_out_suff),  height=5*asp_ratio, width=3.8)
  #-------------------------------------------
  # Plots -- Bar charts
  plt = ggplot() + geom_bar(data=case_data_base, 
                     aes(factor(yyyy_mm_float), 
                         fill=factor(imported, levels=c('N', 'Y'), labels=c('Indigenous', 'Imported') ))) + 
    scale_x_discrete(breaks = c(2019.01,2019.07,2020.01,2020.07, 2021.01,2021.07,2022.01,2022.07), 
                     # labels = c("2019-01","2019-07","2020-01","2020-07", "2021-01","2021-07","2022-01","2022-07"))+
                     labels = c("","","","", "","","",""))+
  # ylim(0, 400, expand=c(0, 0)) + 
    scale_y_continuous(limits=c(0, 450), expand=c(0, 0))+
    ylab("") + 
    xlab("") + 
    # labs(fill='Imported', title=(parasite)) + 
    labs(fill='', title='') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.7, hjust=0.7), 
            axis.line = element_line(colour = "black"),
            axis.text.y=element_blank(),
            # legend.position=c(0.85, 0.8),
            legend.position="none",
            plot.title = element_text(vjust = -7.5, hjust = 0.95),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # panel.border = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA))
  
    my_legend = get_legend(plt)  
    my_legend = as_ggplot(my_legend) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) 
    plt = plt + theme(legend.position='none')

    if(parasite=='P.falciparum'){
      # ggsave(paste0(config$image_out_root, 'Pf_Re_hist_legend', config$image_out_suff), plot=my_legend, height=1.5, width=2)
      ggsave(paste0(config$image_out_root, 'Pf_Re_hist', config$image_out_suff), plot = plt, height=5, width=5)
    }else{
      ggsave(paste0(config$image_out_root, 'Pv_Re_hist_legend', config$image_out_suff), plot=my_legend, height=1.5, width=2)
      ggsave(paste0(config$image_out_root, 'Pv_Re_hist', config$image_out_suff),plot=plt, height=5, width=5)
    }
  
  #-------------------------------------------
  # Plots -- Re predictions 
  plt = ggplot(Re_data_with_covs) + 
    geom_point(aes(x=relative_date, y=Re, colour=factor(imported, levels=c('N','Y'), labels=c('Indigenous', 'Imported'))), shape=1) +
    theme_bw() +
    # labs(title=(parasite))+
    scale_y_continuous(limits=c(0, 8), expand=c(0.01, 0))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5),
          plot.title = element_text(vjust = -7.5, hjust=0.5),
          panel.border = element_rect(colour = "black", fill=NA)
          )+
    labs(y='', x='', colour='Imported')
  
  
  if(parasite=='P.falciparum'){
    plt = plt +theme(legend.position="none") + 
    # + 
            # theme(legend.position=c(0.9, 0.85)) + 
            scale_x_continuous(breaks = c(0, 181, 365, 547, 731,  921, 1096, 1277, 1461),
                                   labels = c("", "", "", "", "", "", "", "", ""),
                                   expand=c(0.01, 0))
    ggsave(paste0(config$image_out_root, 'Pf_Re_timeSeries', config$image_out_suff), height=3, width=9)
  }else{
    outliers = Re_data_with_covs[Re_data_with_covs$Re > 8,]
    outliers$store_Re = outliers$Re
    outliers$Re_labs = outliers$store_Re %>% round(1)
    outliers$Re_labs = midline_decimal(outliers$Re_labs)
    outliers$Re = 7.7 + (outliers$store_Re - 7.7)/20
    plt = plt + geom_point(data = outliers,aes(x=relative_date, y=Re,
                                         colour=factor(imported, levels=c('N','Y'), labels=c('Indigenous', 'Imported'))),
                     shape=1)+
          # geom_text_repel(data=outliers,aes(relative_date, y=Re),label=outliers$Re_labs) + 
          # theme(legend.position="none")+
          scale_x_continuous(breaks = c(0, 181, 365, 547, 731,  921, 1096, 1277, 1461), 
                                   # labels = c("2019-01-01","2019-07-01","2020-01-01","2020-07-01", "2021-01-01","2021-07-01","2022-01-01","2022-07-01", "2023-01-01"),
                                   labels = c("","","","", "","","","", ""),
                                   expand=c(0.01, 0))
    my_legend = get_legend(plt)  
    my_legend = as_ggplot(my_legend) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) 
    plt = plt + theme(legend.position='none')
    
    ggsave(paste0(config$image_out_root, 'Pv_Re_timeSeries_legend', config$image_out_suff), plot=my_legend, height=1.5, width=2)
    ggsave(paste0(config$image_out_root, 'Pv_Re_timeSeries', config$image_out_suff), plot = plt,height=3, width=9)
  }
  
  #TO add time series with covariate
  
  #-------------------------------------------
  # Re zif motivation
  
  plt = ggplot() + geom_histogram(data=case_data_withRe_fn[case_data_withRe_fn$Re < 2.1,], 
                                  aes(Re, fill=imported), 
                                  binwidth=0.05)+
    scale_x_continuous(expand=c(0, 0))+
    scale_y_continuous(limits=c(0, 330), expand=c(0, 0))+
    ylab("") + 
    xlab("") + 
    labs(fill='') + 
    theme(axis.line = element_line(colour = "black"),
          legend.position=c(0.85, 0.85),
          # legend.title="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # panel.border = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA))
  if(parasite=='P.falciparum'){
    plt = plt + theme(legend.position="none")
    ggsave(paste0(config$image_out_root, 'Pf_Re_count_hist', config$image_out_suff), height=5, width=5)
  }else{
    my_legend = get_legend(plt)  
    my_legend = as_ggplot(my_legend) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) 
    plt = plt + theme(legend.position='none')
    ggsave(paste0(config$image_out_root, 'Pv_Re_count_hist_legend', config$image_out_suff), plot=my_legend, height=1.5, width=2)
    
    ggsave(paste0(config$image_out_root, 'Pv_Re_count_hist', config$image_out_suff), plot = plt, height=5, width=5)
  }
  #-------------------------------------------
  # Receptivity 
  prediction_month_name_map = c('January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December')
  month_count = 1
  for(prediction_name in prediction_names){
    plt<-  ggplot(data = subplot_sf) + 
      geom_sf() + 
      geom_tile(data = seasonal_average_preds_df, mapping=aes(x = x, y = y, fill = get(prediction_name)))+ 
      labs(x="", y="", fill='') + 
      coord_sf(xlim = c(unname(st_bbox(poly.shapefile)$xmin), unname(st_bbox(poly.shapefile)$xmax)),
               ylim = c(unname(st_bbox(poly.shapefile)$ymin), unname(st_bbox(poly.shapefile)$ymax))) +
      # geom_label(aes(x = 102.8+ nchar(prediction_month_name_map[month_count])/10, y = 23.65, label = prediction_month_name_map[month_count]), fill = "white")+
      custom_gradient +
      geom_sf(data = admin1_sf, fill=NA) +
      coord_sf(xlim = c(102, 110), ylim = c(8, 24), expand = FALSE) + theme_void()
    if(parasite=='P.falciparum'){
      plt = ggdraw(plt + theme(legend.position='none', plot.margin=margin(c(0, 0, 0, 0))))
      ggsave(paste0(config$image_out_root, 'Pf_', prediction_name, config$image_out_suff), height=5*asp_ratio, width=2.9)
    }else{
      if(month_count==1){
        my_legend = get_legend(plt)  
        my_legend = as_ggplot(my_legend) + theme(plot.margin=grid::unit(c(0,0,4,0), "pt"))
        ggsave(paste0(config$image_out_root, 'Pv_', prediction_name, '_legend', config$image_out_suff), plot=my_legend, height=1.5, width=2)
      }
      plt = ggdraw(plt + theme(legend.position='none', plot.margin=margin(c(0, 0, 0, 0))))
      
      ggsave(paste0(config$image_out_root, 'Pv_', prediction_name, config$image_out_suff),plot = plt, height=5*asp_ratio, width=2.9)
    }
    month_count = month_count + 1
  }
  
  plt = ggplot(data = subplot_sf) + 
    geom_sf() + 
    geom_tile(data = seasonal_average_preds_df, mapping=aes(x = x, y = y, fill = annual_mean), linewidth=0)+ 
    labs(x="", y="", fill='') + 
    # geom_point(data = to_plot, mapping=aes(x = LONG, y = LAT, colour = pred))+ 
    coord_sf(xlim = c(unname(st_bbox(poly.shapefile)$xmin), unname(st_bbox(poly.shapefile)$xmax)),
             ylim = c(unname(st_bbox(poly.shapefile)$ymin), unname(st_bbox(poly.shapefile)$ymax))) +
    # scale_fill_viridis_c(option="turbo" , name=NULL, limits = c(0, abs_max))+ # limits = c(0, 4.5)
    custom_gradient_truncated +
    # labs(title=(parasite)) + 
    geom_sf(data = admin1_sf, fill=NA) +
    coord_sf(xlim = c(102, 110), ylim = c(8, 24), expand = FALSE) + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position=c(0.5, 0.5),
          plot.title = element_text(vjust = -5.9)) + scale_x_continuous(expand=c(0, 0))+
    scale_y_continuous(expand=c(0, 0))
    # theme_void()  +
    # theme(plot.title = element_text(vjust = -5.4))  
  # annotate("text", x=107, y=23, label= "Sep")
  if(parasite=='P.falciparum'){
    plt = ggdraw(plt + theme(legend.position='none', plot.margin=margin(c(0, 0, 0, 0))))
    ggsave(paste0(config$image_out_root, 'Pf_Mean', config$image_out_suff), height=5*asp_ratio, width=2.9)
  }else{
    my_legend = get_legend(plt)  
    my_legend = as_ggplot(my_legend) + theme(plot.margin=grid::unit(c(0,0,4,0), "pt"))
    ggsave(paste0(config$image_out_root, 'Pv_', prediction_name, '_legend', config$image_out_suff), plot=my_legend, height=1.5, width=2)
    plt = ggdraw(plt + theme(legend.position='none', plot.margin=margin(c(0, 0, 0, 0))))
    
    ggsave(paste0(config$image_out_root, 'Pv_Mean_legend', config$image_out_suff), plot=my_legend, height=1.5, width=2)
    ggsave(paste0(config$image_out_root, 'Pv_Mean', config$image_out_suff), plot=plt, height=5*asp_ratio, width=2.9)
  }
  
  #-------------------------------------------
  # Zero inflated probability
  
  month_count = 1
  for(zi_name in zi_names){
    plt<-  ggplot(data = subplot_sf) + 
      geom_sf() + 
      geom_tile(data = seasonal_average_zi_df, mapping=aes(x = x, y = y, fill = get(zi_name)))+ 
      # geom_point(data = to_plot, mapping=aes(x = LONG, y = LAT, colour = pred))+ 
      coord_sf(xlim = c(unname(st_bbox(poly.shapefile)$xmin), unname(st_bbox(poly.shapefile)$xmax)),
               ylim = c(unname(st_bbox(poly.shapefile)$ymin), unname(st_bbox(poly.shapefile)$ymax))) +
      # scale_fill_viridis_c(option="turbo" , limits = c(0, abs_max), name=NULL) +
      # geom_label(aes(x = 102.8+ nchar(prediction_month_name_map[month_count])/10, y = 23.65, label = prediction_month_name_map[month_count]), fill = "white")+
      custom_gradient_zi +
      geom_sf(data = admin1_sf, fill=NA) +
      coord_sf(xlim = c(102, 110), ylim = c(8, 24), expand = FALSE) + theme_void()
    if(parasite=='P.falciparum'){
      plt = plt + theme(legend.position='none')
      ggsave(paste0(config$image_out_root, 'Pf_zi_', zi_name, config$image_out_suff),height=5*asp_ratio, width=2.9)
    }else{
      if(month_count==1){
        my_legend = get_legend(plt)  
        my_legend = as_ggplot(my_legend) + theme(plot.margin=grid::unit(c(0,0,4,0), "pt"))
        ggsave(paste0(config$image_out_root, 'Pv_zi_', zi_name, '_legend', config$image_out_suff), plot=my_legend, height=1.5, width=2)
      }
      plt = plt + theme(legend.position='none')
      
      ggsave(paste0(config$image_out_root, 'Pv_zi_', zi_name, config$image_out_suff), plot=plt, height=5*asp_ratio, width=2.9)
    }
    month_count = month_count + 1
  }
  
  ggplot(data = subplot_sf) + 
    geom_sf() + 
    geom_tile(data = seasonal_average_zi_df, mapping=aes(x = x, y = y, fill = annual_mean))+ 
    # geom_point(data = to_plot, mapping=aes(x = LONG, y = LAT, colour = pred))+ 
    coord_sf(xlim = c(unname(st_bbox(poly.shapefile)$xmin), unname(st_bbox(poly.shapefile)$xmax)),
             ylim = c(unname(st_bbox(poly.shapefile)$ymin), unname(st_bbox(poly.shapefile)$ymax))) +
    # scale_fill_viridis_c(option="turbo" , name=NULL, limits = c(0, abs_max))+ # limits = c(0, 4.5)
    custom_gradient_zi +
    geom_sf(data = admin1_sf, fill=NA) +
    coord_sf(xlim = c(102, 110), ylim = c(8, 24), expand = FALSE) + theme_void()
  if(parasite=='P.falciparum'){
    plt = plt + theme(legend.position='none')
    ggsave(paste0(config$image_out_root, 'Pf_zi_Mean', config$image_out_suff), height=5*asp_ratio, width=2.9)
  }else{
    my_legend = get_legend(plt)  
    my_legend = as_ggplot(my_legend) + theme(plot.margin=grid::unit(c(0,0,4,0), "pt"))
    ggsave(paste0(config$image_out_root, 'Pv_zi_Mean_legend', config$image_out_suff), plot=my_legend, height=1.5, width=2)
    plt = plt + theme(legend.position='none')
    ggsave(paste0(config$image_out_root, 'Pv_zi_Mean', config$image_out_suff), plot = plt, height=5*asp_ratio, width=2.9)
  }

  #-------------------------------------------
  # Vulnerability 
  if(parasite=='P.falciparum'){
    image_out_fn <- paste0(config$image_out_root, 'vuln_Pf', config$image_out_suff)
    vuln_lim<- 0.02/4
    vuln_gradient <- scale_fill_gradientn(
      colours = c("#f8f8f8", "#f87c7c", "#c80815"),
      values = scales::rescale(c(0, 0.4*vuln_lim, vuln_lim)),  # specify the values where colors change
      limits = c(0, vuln_lim),
      name=NULL,
      labels=midline_decimal
    )  
  }else{
    image_out_fn <- paste0(config$image_out_root, 'vuln_Pv', config$image_out_suff)
    # vuln_lim<- 0.005/4
    vuln_lim<- 0.02/4
    vuln_gradient <- scale_fill_gradientn(
      colours = c("#f8f8f8", "#f87c7c", "#c80815"),
      values = scales::rescale(c(0,0.4*vuln_lim,  vuln_lim)),  # specify the values where colors change
      limits = c(0, vuln_lim),
      # name='Importation\n     Rate',
      name=NULL,
      labels=midline_decimal
    )  
  } 
  

  plt <- ggplot(data = subplot_sf) + 
    geom_sf() + 
    geom_tile(data = seasonal_average_preds_df, aes(x=x, y=y, fill = vulnerability))+ 
    labs(x="", y="", fill='') + 
    coord_sf(xlim = c(unname(st_bbox(poly.shapefile)$xmin), unname(st_bbox(poly.shapefile)$xmax)),
             ylim = c(unname(st_bbox(poly.shapefile)$ymin), unname(st_bbox(poly.shapefile)$ymax))) +
    vuln_gradient+
    geom_sf(data = admin1_sf, fill=NA) +
    # labs(title=(parasite)) + 
    coord_sf(xlim = c(102, 110), ylim = c(8, 24), expand = FALSE) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position=c(0.5, 0.5),
          plot.title = element_text(vjust = -5.9)) + scale_x_continuous(expand=c(0, 0))+
          scale_y_continuous(expand=c(0, 0))
  if(parasite=='P.falciparum'){
    plt <- ggdraw(plt + theme(legend.position="none", plot.margin=margin(c(0, 0, 0, 0))))
    ggsave(image_out_fn, height=5*asp_ratio, width=2.9)
  }else{
    my_legend = get_legend(plt)  
    my_legend = as_ggplot(my_legend) + theme(plot.margin=grid::unit(c(0,0,4,0), "pt"))
    ggsave(paste0(config$image_out_root, 'vuln_Pv_legend', config$image_out_suff), plot=my_legend, height=1.5, width=2)
    plt = ggdraw(plt + theme(legend.position='none', plot.margin=margin(c(0, 0, 0, 0))))
    # plt <- plt + theme(legend.position=c(1.2, 0.45), plot.margin=unit(c(0, 2.5, 0, 0), "cm"))
    ggsave(image_out_fn,height=5*asp_ratio,plot=plt, width=2.9)
  } 
    

  #-------------------------------------------
  # Vulnerability 
  if(parasite=='P.falciparum'){
    image_out_fn <- paste0(config$image_out_root, 'vuln_index_Pf', config$image_out_suff)
    vuln_lim<- 0.02/4
    vuln_gradient <- scale_fill_gradientn(
      colours = c("#f8f8f8", "#f87c7c", "#c80815"),
      values = scales::rescale(c(0, 0.4*vuln_lim, vuln_lim)),  # specify the values where colors change
      limits = c(0, vuln_lim),
      name=NULL,
      labels=midline_decimal
    )  
  }else{
    image_out_fn <- paste0(config$image_out_root, 'vuln_index_Pv', config$image_out_suff)
    # vuln_lim<- 0.005/4
    vuln_lim<- 0.02/4
    vuln_gradient <- scale_fill_gradientn(
      colours = c("#f8f8f8", "#f87c7c", "#c80815"),
      values = scales::rescale(c(0,0.4*vuln_lim,  vuln_lim)),  # specify the values where colors change
      limits = c(0, vuln_lim),
      # name='Importation\n     Rate',
      name=NULL,
      labels=midline_decimal
    )  
  } 
  
  
  plt <- ggplot(data = subplot_sf) + 
    geom_sf() + 
    geom_tile(data = seasonal_average_preds_df, aes(x=x, y=y, fill = vulnerability_index))+ 
    labs(x="", y="", fill='Importation rate') + 
    coord_sf(xlim = c(unname(st_bbox(poly.shapefile)$xmin), unname(st_bbox(poly.shapefile)$xmax)),
             ylim = c(unname(st_bbox(poly.shapefile)$ymin), unname(st_bbox(poly.shapefile)$ymax))) +
    vuln_gradient+
    geom_sf(data = admin1_sf, fill=NA) +
    coord_sf(xlim = c(102, 110), ylim = c(8, 24), expand = FALSE) +
    # labs(title=(parasite)) + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(vjust = -5.9))
  if(parasite=='P.falciparum'){
    plt <- plt + theme(legend.position="none")
    ggsave(image_out_fn, height=5*asp_ratio, width=2.9)
  }else{
    my_legend = get_legend(plt)  
    my_legend = as_ggplot(my_legend) + theme(plot.margin=grid::unit(c(0,0,4,0), "pt"))
    ggsave(paste0(config$image_out_root, 'vuln_index_Pv_legend', config$image_out_suff), plot=my_legend, height=1.5, width=2)
    plt = plt + theme(legend.position='none')
    # plt <- plt + theme(legend.position=c(1.2, 0.45), plot.margin=unit(c(0, 2.5, 0, 0), "cm"))
    ggsave(image_out_fn, height=5*asp_ratio, width=2.9)
  } 
  
  #-------------------------------------------
  # Malariogenic potential 
  mgc_lim<- 0.0022
  legend_out_fn = paste0(config$image_out_root, 'mgp_legend', config$image_out_suff)
  if(parasite=='P.falciparum'){
    image_out_fn <- paste0(config$image_out_root, 'mgp_Pf', config$image_out_suff)
    mgc_gradient <- scale_fill_gradientn(
      colours = c("#f8f8f8", "#f87c7c", "#c80815"),
      values = scales::rescale(c(0, 0.4*mgc_lim, mgc_lim)),  # specify the values where colors change
      limits = c(0, mgc_lim),
      name=NULL,
      labels=midline_decimal
    )  
  }else{
    image_out_fn <- paste0(config$image_out_root, 'mgp_Pv', config$image_out_suff)
    mgc_gradient <- scale_fill_gradientn(
      colours = c("#f8f8f8", "#f87c7c", "#c80815"),
      values = scales::rescale(c(0, 0.4*mgc_lim, mgc_lim)),  # specify the values where colors change
      limits = c(0, mgc_lim),
      labels=midline_decimal,
      name=NULL)
      # name='Malariogenic\n   Potential')
  }
  plt<- ggplot(data = subplot_sf) + 
    geom_sf() + 
    geom_tile(data = seasonal_average_preds_df, aes(x=x, y=y, fill = malariogenic))+ 
    labs(x="", y="", fill='') + 
    coord_sf(xlim = c(unname(st_bbox(poly.shapefile)$xmin), unname(st_bbox(poly.shapefile)$xmax)),
             ylim = c(unname(st_bbox(poly.shapefile)$ymin), unname(st_bbox(poly.shapefile)$ymax))) +
    mgc_gradient+
    # scale_fill_viridis_c(option="magma" , name=NULL, aesthetics = 'fill', limits = c(0, mgc_lim))+ # limits = c(0, 4.5)
    # scale_fill_viridis_c(option="magma" , name=NULL, aesthetics = 'colour', limits = c(0, mgc_lim))+ # limits = c(0, 4.5)
    geom_sf(data = admin1_sf, fill=NA) +
    coord_sf(xlim = c(102, 110), ylim = c(8, 24), expand = FALSE) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = c(0.5, 0.5),
          plot.title = element_text(vjust = -5.9))+ scale_x_continuous(expand=c(0, 0))+
    scale_y_continuous(expand=c(0, 0))
  if(parasite=='P.falciparum'){
    plt <- ggdraw(plt + theme(legend.position="none", plot.margin=margin(c(0, 0, 0, 0))))
    ggsave(image_out_fn, height=5*asp_ratio, width=2.9)
  }else{
    my_legend = get_legend(plt)  
    my_legend = as_ggplot(my_legend) + theme(plot.margin=grid::unit(c(0,0,4,0), "pt"))
    ggsave(paste0(config$image_out_root, 'mgp_Pv_legend', config$image_out_suff), plot=my_legend, height=1.5, width=2)
    plt = ggdraw(plt + theme(legend.position='none', plot.margin=margin(c(0, 0, 0, 0))))
    ggsave(image_out_fn,plot=plt, height=5*asp_ratio, width=2.9)
  } 
  
  #-------------------------------------------
  # Malariogenic potential (based on "index cases") 
  mgc_lim<- 0.0022
  # legend_out_fn = paste0(config$image_out_root, 'mgp_legend', config$image_out_suff)
  if(parasite=='P.falciparum'){
    image_out_fn <- paste0(config$image_out_root, 'mgp_index_Pf', config$image_out_suff)
    mgc_gradient <- scale_fill_gradientn(
      colours = c("#f8f8f8", "#f87c7c", "#c80815"),
      values = scales::rescale(c(0, mgc_lim)),  # specify the values where colors change
      limits = c(0, mgc_lim),
      name=NULL,
      labels=midline_decimal
    )  
  }else{
    image_out_fn <- paste0(config$image_out_root, 'mgp_index_Pv', config$image_out_suff)
    mgc_gradient <- scale_fill_gradientn(
      colours = c("#f8f8f8", "#f87c7c", "#c80815"),
      values = scales::rescale(c(0, 0.4*mgc_lim, mgc_lim)),  # specify the values where colors change
      limits = c(0, mgc_lim),
      labels=midline_decimal,
      name=NULL)
      # name='Malariogenic\n   Potential')
  }
  plt<- ggplot(data = subplot_sf) + 
    geom_sf() + 
    geom_tile(data = seasonal_average_preds_df, aes(x=x, y=y, fill = malariogenic_index))+ 
    # labs(x="", y="", title=(parasite)) + 
    labs(x="", y="",fill='') + 
    coord_sf(xlim = c(unname(st_bbox(poly.shapefile)$xmin), unname(st_bbox(poly.shapefile)$xmax)),
             ylim = c(unname(st_bbox(poly.shapefile)$ymin), unname(st_bbox(poly.shapefile)$ymax))) +
    mgc_gradient+
    geom_sf(data = admin1_sf, fill=NA) +
    coord_sf(xlim = c(102, 110), ylim = c(8, 24), expand = FALSE) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(vjust = -5.9))
  if(parasite=='P.falciparum'){
    plt <- ggdraw(plt + theme(legend.position="none"))
    ggsave(image_out_fn, height=5*asp_ratio, width=2.9)
  }else{
    my_legend = get_legend(plt)  
    my_legend = as_ggplot(my_legend) + theme(plot.margin=grid::unit(c(0,0,4,0), "pt"))
    ggsave(paste0(config$image_out_root, 'mgp_index_Pv_legend', config$image_out_suff), plot=my_legend, height=1.5, width=2)
    plt = ggdraw(plt + theme(legend.position='none'))
    # plt <- plt + theme(legend.position=c(1.22, 0.45), plot.margin=unit(c(0, 3, 0, 0), "cm"))
    ggsave(image_out_fn,plot=plt, height=5*asp_ratio, width=2.9)
  } 

  
}

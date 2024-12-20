
image_out_suff <- paste0(date_out_str, '.png')
image_out_root <- paste0('/mnt/Z/Michael/Vietnam/Modelling/Re_and_descriptive/Re_images_',date_out_str, '/')
if(!dir.exists(image_out_root)){
  dir.create(image_out_root)
}







  #----------------------------------------------------------------------------#
  #            Geospatial Rc plots
  #----------------------------------------------------------------------------#
#TODO: parasite mush be chosen
if(parasite=='P.falciparum'){
    seasonal_average_df <- read.csv(paste0(data_out_root, 'Re_gamma_Pf_2019', data_out_suff))
}else{
    seasonal_average_df <- read.csv(paste0(data_out_root, 'Re_gamma_Pv_2019', data_out_suff))
}

prediction_names <- seasonal_average_df %>% names()
prediction_names <- prediction_names[grepl('prediction', prediction_names)]
seasonal_average_df <- seasonal_average_df[, c('x', 'y', prediction_names)]

count <- 1
for(year_str in c('2020', '2021', '2022')){
    if(parasite=='P.falciparum'){
    new_df <- read.csv(paste0(data_out_root, 'Re_gamma_Pf_', year_str, data_out_suff))
    }else{
    new_df <- read.csv(paste0(data_out_root, 'Re_gamma_Pv_', year_str, data_out_suff))
    }
    seasonal_average_df[,prediction_names] = (count*seasonal_average_df[,prediction_names] + new_df[, prediction_names])/(count + 1)
    count <- count + 1
}

row_max <- apply(seasonal_average_df[,prediction_names]  , 1, max)
abs_max <- row_max[!is.na(row_max)] %>% max()
if(parasite=='P.falciparum'){
    new_df <- read.csv(paste0(data_out_root, 'Re_gamma_Pf_', year_str, data_out_suff))
    abs_max<- 4
}else{
    new_df <- read.csv(paste0(data_out_root, 'Re_gamma_Pv_', year_str, data_out_suff))
    abs_max<- 1.1
}
seasonal_average_df$annual_mean <- seasonal_average_df[,prediction_names] %>% rowMeans()

for(prediction_name in prediction_names){
    if (nrow(seasonal_average_df[(!is.na(seasonal_average_df[,prediction_name]))&seasonal_average_df[,prediction_name]>abs_max,])>0){
    seasonal_average_df[(!is.na(seasonal_average_df[,prediction_name]))&seasonal_average_df[,prediction_name]>abs_max,prediction_name] = abs_max
    }
    
}

for(prediction_name in prediction_names){
    plt<-  ggplot(data = subplot_sf) + 
    geom_sf() + 
    geom_tile(data = seasonal_average_df, mapping=aes(x = x, y = y, fill = get(prediction_name)))+ 
    # geom_point(data = to_plot, mapping=aes(x = LONG, y = LAT, colour = pred))+ 
    coord_sf(xlim = c(unname(st_bbox(poly.shapefile)$xmin), unname(st_bbox(poly.shapefile)$xmax)),
                ylim = c(unname(st_bbox(poly.shapefile)$ymin), unname(st_bbox(poly.shapefile)$ymax))) +
    scale_fill_viridis_c(option="turbo" , limits = c(0, abs_max), name=NULL) +
    geom_sf(data = admin2_sf, fill=NA) +
    coord_sf(xlim = c(102, 110), ylim = c(8, 24), expand = FALSE) + theme_void()
    # print(plt) #TODO: save
    if(parasite=='P.falciparum'){
    ggsave(paste0(image_out_root, 'Pf_', prediction_name,image_out_suff), height=6, width=3)
    }else{
    ggsave(paste0(image_out_root, 'Pv_', prediction_name, image_out_suff), height=6, width=3)
    }
}

seasonal_average_df = seasonal_average_df[!is.na(seasonal_average_df$annual_mean), ]
if (length(seasonal_average_df[seasonal_average_df$annual_mean>abs_max,]$annual_mean)>0){
    seasonal_average_df[seasonal_average_df$annual_mean>abs_max,]$annual_mean = abs_max
}


ggplot(data = subplot_sf) + 
    geom_sf() + 
    geom_tile(data = seasonal_average_df, mapping=aes(x = x, y = y, fill = annual_mean))+ 
    # geom_point(data = to_plot, mapping=aes(x = LONG, y = LAT, colour = pred))+ 
    coord_sf(xlim = c(unname(st_bbox(poly.shapefile)$xmin), unname(st_bbox(poly.shapefile)$xmax)),
            ylim = c(unname(st_bbox(poly.shapefile)$ymin), unname(st_bbox(poly.shapefile)$ymax))) +
    scale_fill_viridis_c(option="turbo" , name=NULL, limits = c(0, abs_max))+ # limits = c(0, 4.5)
    geom_sf(data = admin2_sf, fill=NA) +
    coord_sf(xlim = c(102, 110), ylim = c(8, 24), expand = FALSE) + theme_void()
if(parasite=='P.falciparum'){
    ggsave(paste0(image_out_root, 'Pf_Mean', image_out_suff), height=6, width=3)
}else{
    ggsave(paste0(image_out_root, 'Pv_Mean', image_out_suff), height=6, width=3)
}
    
  
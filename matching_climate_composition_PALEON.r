#get historic climate data paleon.
rm(list = ls())

#using x-y csv files provided by Ann Raiho of historic and contemporary climate.
data.dir <- '/fs/data3/caverill/paleo_myco_2018/'
    historic <- read.csv(paste0(data.dir,'biomass_prediction_v0.9-10_bam.csv'))
contemporary <- read.csv(paste0(data.dir,'FIA_biomass_predictions_v0.1.csv'))
#load historic climate
h.clim <- read.csv(paste0(data.dir,'paleon_models_environment_master.csv'))


#convert x-y albers coordiantes to WGS 84 for historic and contemporary products. Code by Ann Raiho.
#historic
albers <- cbind(as.numeric(as.character(historic$x)),as.numeric(as.character(historic$y)))
albers.df = data.frame(albers)
colnames(albers.df) = c('y', 'x')
coordinates(albers.df) <- ~ x + y
proj4string(albers.df) <- CRS('+init=epsg:3175')
lat.long <- spTransform(albers.df, CRS('+proj=longlat +ellps=WGS84'))
lat.long <- (data.frame(lat.long))
colnames(lat.long) <- c('longitude','latitude')
historic <- cbind(lat.long,historic)
#repeat for contemporary
albers <- cbind(as.numeric(as.character(contemporary$x)),as.numeric(as.character(contemporary$y)))
albers.df = data.frame(albers)
colnames(albers.df) = c('y', 'x')
coordinates(albers.df) <- ~ x + y
proj4string(albers.df) <- CRS('+init=epsg:3175')
lat.long <- spTransform(albers.df, CRS('+proj=longlat +ellps=WGS84'))
lat.long <- (data.frame(lat.long))
colnames(lat.long) <- c('longitude','latitude')
contemporary <- cbind(lat.long,contemporary)


###How do I assign historic climate to historic or contemporary composition products?
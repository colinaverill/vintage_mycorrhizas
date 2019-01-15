#calculate historic and contemporary composition
library(raster)
library(sp)
library(rgdal)
library(data.table)
library(boot)
library(betareg)
library(wesanderson)
source('/home/caverill/NEFI_microbe/NEFI_functions/extract_ndep.r')
source('/home/caverill/NEFI_microbe/NEFI_functions/crib_fun.r')
source('/home/caverill/vintage_mycorrhizas/project_functions/prism_query.r')

#source paths
source('paths.r')

#read historic data for a species
    historic <- raster(    historic.composition.path,varname = 'Hemlock')
contemporary <- raster(contemporary.composition.path,varname = 'Hemlock')
myc <- readRDS('/fs/data3/caverill/myc_traits/myc_assignments.rds')


#using x-y csv files provided by Ann Raiho
    historic <- read.csv('/fs/data3/caverill/paleo_myco_2018/biomass_prediction_v0.9-10_bam.csv')
contemporary <- read.csv('/fs/data3/caverill/paleo_myco_2018/fia_contemporary_biomass/FIA_biomass_predictions_v0.1.csv')

#convert x-y albers coordiantes to WGS 84
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

#data.table
    historic <- data.table(historic)
contemporary <- data.table(contemporary)

#See if you can match up by cell
    historic <- historic[cell %in% contemporary$cell,]
contemporary <- contemporary[cell %in% historic$cell,]

#assign mycorrhizal associations
spp <- colnames(historic)[6:27]
#still unsure on Salix (willow) for North America. Need species level list. 
h.em <- c(0,0,1,1,1,0,1,1,1,1,1,1,1,0,1,NA,1,0,1,1,0,NA)
h.am <- 1 - h.em
historic.basal <- historic[,6:27]
historic.em <- t(t(historic.basal)*h.em)
historic.am <- t(t(historic.basal)*h.am)
historic.em <- rowSums(historic.em, na.rm = T)
historic.am <- rowSums(historic.am, na.rm = T)
historic.total <- rowSums(historic.basal, na.rm=T)
historic[,AM := historic.am]
historic[,EM := historic.em]
historic[,Total := historic.total]
historic[,AM.EM := AM + EM]
contemporary.basal <- contemporary[,6:24]
c.em <- c(0,1,1,1,0,1,1,1,1,1,0,1,NA,1,NA,1,1,0,0)
c.am <- 1-c.em
contemporary.em <- t(t(contemporary.basal)*c.em)
contemporary.am <- t(t(contemporary.basal)*c.am)
contemporary.em <- rowSums(contemporary.em, na.rm=T)
contemporary.am <- rowSums(contemporary.am, na.rm=T)
contemporary.total <- rowSums(contemporary.basal)
contemporary[,AM := contemporary.am]
contemporary[,EM := contemporary.em]
contemporary[,Total := contemporary.total]
contemporary[,AM.EM := AM + EM]

#subset to cells where AM and EM trees make up 90% of composition in both datasets.
#More limited by historic composition than contemporary, which makes sense.
contemporary <- contemporary[AM.EM/Total > 0.9,]
    historic <-     historic[AM.EM/Total > 0.9,]
contemporary <- contemporary[cell %in% historic$cell,]
    historic <-     historic[cell %in% contemporary$cell,]
    
#Get a dataframe that is historic+contemporary AM, EM and total abundance, lat and lon.
c.out <- contemporary[,.(latitude,longitude,cell,AM,EM,AM.EM,Total)]
names(c.out)[4:7] <- c('c.AM','c.EM','c.AM.EM','c.Total')
h.out <-     historic[,.(cell,AM,EM,AM.EM,Total)]
names(h.out) <- c('cell','h.AM','h.EM','h.AM.EM','h.Total')
out <- merge(c.out,h.out)

#Get N deposition.
Ndep <- extract_ndep(longitude = out$longitude, latitude = out$latitude)
out <- cbind(out,Ndep)
#Get PRISM climate
climate <- prism_query(longitude = out$longitude,latitude = out$latitude)
out <- cbind(out,climate)

#Get historic and contemporary compsoition.
out[, c.relEM := c.EM / c.Total]
out[, h.relEM := h.EM / h.Total]


#get historic climate and soils data.
h.clim <- read.csv(historic.climate.path)

#function to get historic climate data for given lat/lon
get_ind_coors <- function(h.clim,data){
  centers = h.clim[,c('lon','lat')]#cbind(historic$x,historic$y)
  idx_cores = vector(length=nrow(data))
  for(i in seq_along(idx_cores)){
    core_site = data[i,c('longitude','latitude')]
    d = fields::rdist(matrix(core_site, ncol=2), as.matrix(centers))
    #find smallest distance to climate coordinate
    idx_cores[i] = which.min(d)
  }
  return(idx_cores)
}

#get row numbers of historic climate product closest to plot centers from PALEON data.
hist_idx <- get_ind_coors(h.clim,out)
out$key <- hist_idx
h.clim$key <- seq(1,nrow(h.clim))
out <- data.table(out)
h.clim <- data.table(h.clim)
setkey(out, key)
setkey(h.clim, key)
out <- merge(out, h.clim)
#drop unncessary columns
dropcols <- names(h.clim)[1:6]
out[,c(dropcols) := NULL]

#convert from degrees kelvin to degrees celsius
out$tair.yr.set <- out$tair.yr.set - 273.15
out$tair.jja.set <- out$tair.jja.set - 273.15
out$tair.yr.cru <- out$tair.yr.cru - 273.15
out$tair.jja.cru <- out$tair.jja.cru - 273.15

#convert from kg / m2 / s rainfall to mm / yr rainfall
out$precip.jja.cru <- out$precip.jja.cru *  60*60*24*365.25
out$precip.jja.set <- out$precip.jja.set *  60*60*24*365.25
out$precip.yr.cru  <- out$precip.yr.cru  *  60*60*24*365.25
out$precip.yr.set  <- out$precip.yr.set  *  60*60*24*365.25

#save output
saveRDS(out, historic_contemporary_merge.path)
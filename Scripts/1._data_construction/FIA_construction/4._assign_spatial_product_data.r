#Assign mat, map and n-deposition to sites.
rm(list=ls())
library(data.table)
library(raster)
library(rgdal)
library(sp)
detach("package:runjags", unload=TRUE) #make sure runjags is not loaded.
source('required_products_utilities/extract_ndep.r')
source('required_products_utilities/prism_query.r')
source('required_products_utilities/arid_extract.r')
source('required_products_utilities/worldclim2_grab.r')
source('paths.r')


p1 <- data.table(readRDS(Product_1.path))
p2 <- data.table(readRDS(Product_2.path))
p3 <- readRDS(Product_3.path); p3 <- data.table(p3)
#soil paths.
p1.soil <- readRDS(Product_1.soil.path)
p2.soil <- readRDS(Product_2.soil.path)
p3.soil <- readRDS(Product_3.soil.path); p3.soil <- data.table(p3.soil)

#remove spatial columns if already present.
to.drop <- c('n.dep','wet.dep','dry.dep','mat','map','mat_CV','map_CV','mat_sd','map_sd','mdr','aridity')
p1 <- p1[,c(to.drop) := NULL]
p2 <- p2[,c(to.drop) := NULL]
p3 <- p3[,c(to.drop) := NULL]
p1.soil <- p1.soil[,c(to.drop) := NULL]
p2.soil <- p2.soil[,c(to.drop) := NULL]
p3.soil <- p3.soil[,c(to.drop) := NULL]

#too much RAM required to do huge products 1 by 1, but all PLT_CN values are in p1.
#pull N deposition data. ignore warnings.
p1 <- cbind(p1,extract_ndep(p1$longitude, p1$latitude))
p2 <- merge(p2,p1[,.(PLT_CN,n.dep,dry.dep,wet.dep)])
p3 <- merge(p3,p1[,.(PLT_CN,n.dep,dry.dep,wet.dep)])
p1.soil <- cbind(p1.soil, extract_ndep(p1.soil$longitude, p1.soil$latitude))
p2.soil <- cbind(p2.soil, extract_ndep(p2.soil$LON      , p2.soil$LAT     ))
p3.soil <- cbind(p3.soil, extract_ndep(p3.soil$longitude, p3.soil$latitude))

#pull worldclim2 climate data and aridity
p1 <- cbind(p1,worldclim2_grab(p1$latitude,p1$longitude))
p1$aridity <- arid_extract(p1$latitude,p1$longitude)
p2 <- merge(p2, p1[,.(PLT_CN,mat,map,mat_CV,map_CV,mdr,aridity)])
p3 <- merge(p3, p1[,.(PLT_CN,mat,map,mat_CV,map_CV,mdr,aridity)])
p1.soil <- cbind(p1.soil, worldclim2_grab(p1.soil$latitude, p1.soil$longitude))
p2.soil <- cbind(p2.soil, worldclim2_grab(p2.soil$LAT     , p2.soil$LON      ))
p3.soil <- cbind(p3.soil, worldclim2_grab(p3.soil$latitude, p3.soil$longitude))
p1.soil$aridity <- arid_extract(p1.soil$latitude, p1.soil$longitude)
p2.soil$aridity <- arid_extract(p2.soil$LAT     , p2.soil$LON      )
p3.soil$aridity <- arid_extract(p3.soil$latitude, p3.soil$longitude)

#save output.
saveRDS(p1, Product_1.path)
saveRDS(p2, Product_2.path)
saveRDS(p3, Product_3.path)
saveRDS(p1.soil, Product_1.soil.path)
saveRDS(p2.soil, Product_2.soil.path)
saveRDS(p3.soil, Product_3.soil.path)

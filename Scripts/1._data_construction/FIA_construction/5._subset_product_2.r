#Subset p1 and p2 to get a random stratified sampling for fitting models quicky.
rm(list=ls())
source('paths.r')
library(ggplot2)
library(ggalt)
library(data.table)

#set output path.----
output.path <- Product_2.subset.path

#load data.----
p2 <- readRDS(Product_2.path)
states <- read.csv('required_products_utilities/FIA_state_codes_regions.csv')

#Grab a plot table with PLT_CN, lat-lon and STATECD.----
d <- p2[,.(PLT_CN,LAT,LON,STATECD,relEM,REMPER)]
setkey(d, 'PLT_CN')
d <- unique(d)
d <- d[d$REMPER >=4.9 & d$REMPER <= 5.1,]


#Subset.----
set.seed(42069)
n.plots <- 4000 #set number of plots to subsample.
d <- d[sample(nrow(d), n.plots),]
p2 <- p2[PLT_CN %in% d$PLT_CN,]

#plot subset to make sure its representative.----
plot = T
if(plot == T){
  lon <- d$LON
  lat <- d$LAT
  world <- map_data('usa')
  map <- ggplot() + geom_cartogram(data = world, map = world, 
                                   aes(x=long, y = lat, group = group, map_id=region))
  map <- map + coord_proj("+proj=wintri", ylim = c(20,50))
  map <- map + geom_point(aes(x = lon, y = lat), color = "yellow"    , size = .2)
  map <- map + theme(axis.line=element_blank(),
                     axis.text.x=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks=element_blank(),
                     axis.title.x=element_blank(),
                     axis.title.y=element_blank(),
                     legend.position="none",
                     panel.border=element_blank()
  )
  map
}

#Save output.----
saveRDS(p2, output.path)

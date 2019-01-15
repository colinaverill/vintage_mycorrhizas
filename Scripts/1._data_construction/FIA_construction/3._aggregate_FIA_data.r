#This script uses the output of "1. FIA Extract and Filter.r". This output represents all data at the tree-level that meet our filtering criteria, and are plots that are represented in the soils data set. 
#This script uses the output of "2. FIA soil data extraction.r". This output is the profile scale aggregated soil data for each plot that met all filtering criteria. 
#This script builds 3 products
#Product_1: The total basal area of all trees in every plot at the time of soil measurement, and then subsets this basal area by mycorrhizal type and PFT, pairs with soils.
#This is for the soil C storage analysis, and the relative abundance of AM and EM trees analysis.
#Product_2: This is the tree level data for mortality analysis, paired with soils data. Uses the 'future' FIA plot data, generates a death vector, and identifies the mycorrhizal status of each tree.
#Product_3: Plot level growth and recruitment data- plot level basal area increment of all trees that survived the remeasurement interval. I do not subtract death here. Also all new trees that cross 5in DIA threshold counted as new recruits.

#s.FIA is soils
#a.FIA is all FIA (not subsetted to soils).
#s.FIA.1 or a.FIA.1 is the 1st measurement.
#s.FIA.2 or a.FIA.1 is the 2nd measurement.

rm(list=ls())
require(data.table)
source('paths.r')

#load data from Tree and Soil queries
s.FIA.1 <- readRDS(FIA_extraction_out.path)
s.FIA.2 <- readRDS(FIA_extraction_out_FUTURE.path)
a.FIA.1 <- readRDS(all.past1.path)
a.FIA.2 <- readRDS(all.present.path)
Soils      <- read.csv(soil_data.processed.path)
FIA.states <- data.table(read.csv('required_products_utilities/FIA_state_codes_regions.csv'))
#put data frames in a list. This will make your life easier.
data.list <- list(s.FIA.1,s.FIA.2,a.FIA.1,a.FIA.2)
names(data.list) <- c('s.FIA.1','s.FIA.2','a.FIA.1','a.FIA.2')

#remove quotes from CN values for both FIA data sets because they mess everything up.
for(i in 1:length(data.list)){
  data.list[[i]]$PLT_CN      <- as.numeric(gsub('"', "",data.list[[i]]$PLT_CN     ))
  data.list[[i]]$PREV_PLT_CN <- as.numeric(gsub('"', "",data.list[[i]]$PREV_PLT_CN))
  data.list[[i]]$TRE_CN      <- as.numeric(gsub('"', "",data.list[[i]]$TRE_CN     ))
  data.list[[i]]$PREV_TRE_CN <- as.numeric(gsub('"', "",data.list[[i]]$PREV_TRE_CN))
}

#remove any plots that have "clear evidence of artificial regeneration." STDORGCD == 1. 
for(i in 1:length(data.list)){
  to.remove <- unique(data.list[[i]][STDORGCD == 1,]$PLT_CN)
  data.list[[i]] <- data.list[[i]][!(PLT_CN %in% to.remove),]
}

#remove any plots that have a tree with STATUSCD = 3. These are plots where humans cut down a tree.
#125 of 4155 intiial, and 285/2918 future plots.
for(i in 1:length(data.list)){
  to.remove <- unique(data.list[[i]][STATUSCD == 3,]$PLT_CN)
  data.list[[i]] <- data.list[[i]][!(PLT_CN %in% to.remove),]
}

####Remove all saplings (DIA < 5inches) based on microplot samplings.
for(i in 1:length(data.list)){
  data.list[[i]] <- data.list[[i]][!(TPA_UNADJ == 74.965282),]
}

####Remove one random site that has very strange growth/recruitment numbers.
####Fairly confident this is recovering from a recent clearcut, but was not indicated in other filters.
for(i in 1:length(data.list)){
  data.list[[i]] <- data.list[[i]][!(PLT_CN == 65355954010538),]
}

###############################################################
#####Product 1. Basal area of each plot paired with soils######
###############################################################

#Calculate number of species in each plot.
for(i in 1:length(data.list)){
  data.list[[i]][, spp.count := uniqueN(SPCD), by = PLT_CN]
}

#generate lists of myc types and PFTs.
 em.list <- levels(a.FIA.2$MYCO_ASSO)
pft.list <- levels(a.FIA.2$PFT)
all.list <- c(em.list,pft.list)

#convert some things to numeric that should be.
for(i in 1:length(data.list)){
  data.list[[i]][,PREVDIA := as.numeric(PREVDIA)]
  data.list[[i]][, REMPER := as.numeric(REMPER) ]
}

#if you are currently dead, your current basal area is assigned NA
for(i in 1:length(data.list)){
  data.list[[i]]$BASAL <- ifelse(!is.na(data.list[[i]]$AGENTCD), NA, data.list[[i]]$BASAL)
}

#current and previous basal area in cm2
for(i in 1:length(data.list)){
  data.list[[i]][,BASAL    := pi*((2.54*DIA    )/2)^2]
  data.list[[i]][,PREVBASAL:= pi*((2.54*PREVDIA)/2)^2]
}

#calculate current basal area of all trees by mycorrhizal type and PFT
#myc type
for(i in 1:length(data.list)){
  for(k in 1:length(em.list)){
    name <- paste0('BASAL.',em.list[k])
    data.list[[i]][MYCO_ASSO == em.list[k],new := BASAL]
    setnames(data.list[[i]], 'new', name)
  }
}
#pft
for(i in 1:length(data.list)){
  for(k in 1:length(pft.list)){
    name <- paste0('BASAL.',pft.list[k])
    data.list[[i]][MYCO_ASSO == pft.list[k],new := BASAL]
    setnames(data.list[[i]], 'new', name)
  }
}

#begin aggregation of tree-level data to plot level. 
scaled.list <- list()
for(i in 1:length(data.list)){
  scaled <- aggregate(data.list[[i]]$BASAL ~ data.list[[i]]$PLT_CN, FUN = 'sum', na.rm = T, na.action = na.pass)
  names(scaled) <- c('PLT_CN','plot.BASAL')
  scaled.list[[i]] <- data.table(scaled)
}
names(scaled.list) <- names(data.list)

#aggregate basal area per plot by myctype and PFT
for(i in 1:length(scaled.list)){
  for(k in 1:length(all.list)){
    name <- paste0('BASAL.',all.list[k])
    scaled.list[[i]]$new <- aggregate(data.list[[i]][[name]] ~ data.list[[i]]$PLT_CN, FUN='sum',na.rm=T, na.action=na.pass)[,2]
    setnames(scaled.list[[i]],'new',name)
  }
}

#pop in relevant data from plot table by taking medians
for(i in 1:length(scaled.list)){
  scaled.list[[i]]$latitude    <- aggregate(data.list[[i]]$LAT     ~ data.list[[i]]$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
  scaled.list[[i]]$longitude   <- aggregate(data.list[[i]]$LON     ~ data.list[[i]]$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
  scaled.list[[i]]$elevation   <- aggregate(data.list[[i]]$ELEV    ~ data.list[[i]]$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
  scaled.list[[i]]$INVYR       <- aggregate(data.list[[i]]$INVYR   ~ data.list[[i]]$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
  scaled.list[[i]]$STATECD     <- aggregate(data.list[[i]]$STATECD ~ data.list[[i]]$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
  scaled.list[[i]]$STDAGE      <- aggregate(data.list[[i]]$STDAGE  ~ data.list[[i]]$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
  scaled.list[[i]]$REMPER      <- aggregate(data.list[[i]]$REMPER  ~ data.list[[i]]$PLT_CN, FUN='median',na.rm=T,na.action=na.pass)[,2]
  scaled.list[[i]]$n.trees     <- data.list[[i]][, .N, by = PLT_CN][,2]
  scaled.list[[i]]$PREV_PLT_CN <- aggregate(data.list[[i]]$PREV_PLT_CN ~ data.list[[i]]$PLT_CN, FUN='unique', na.action = na.pass)[,2]
  #calculate relative abundance of EM trees at time of soil sampling and total EM + AM
  scaled.list[[i]]$relEM    <- scaled.list[[i]]$BASAL.ECM / scaled.list[[i]]$plot.BASAL
  scaled.list[[i]]$relEM.AM <- (scaled.list[[i]]$BASAL.ECM + scaled.list[[i]]$BASAL.AM) / scaled.list[[i]]$plot.BASAL
}

#remove sites that are less than 90% AM+ECM trees.
for(i in 1:length(scaled.list)){
  scaled.list[[i]] <- scaled.list[[i]][!(relEM.AM < 0.9),]
}

#Product_1 is a list of scaled forest data for soil and all 1 and 2 samplings.
#Product_1.soil is the first soil sampling paired with soils data.
Product_1      <- scaled.list[[4]]
Product_1.soil <- scaled.list[[1]]
Product_1.soil <- merge(Product_1.soil, Soils, by = 'PLT_CN')

#remove sites with CN > 90
Product_1.soil$cn <- Product_1.soil$C.storage / Product_1.soil$N.storage
Product_1.soil    <- Product_1.soil[cn < 90,]

#save output.
saveRDS(Product_1     , file=Product_1.path     )
saveRDS(Product_1.soil, file=Product_1.soil.path)


##################################################################
##### Product 2. Individual-level Growth and Mortality data ######
##################################################################

#Mortality is modeled at the tree level. Take most recent data.
mort.list <- list(data.list[[2]], data.list[[4]])

#flag whether a tree died or not, for any reason. If it did, there will be a value associated with "AGENTCD" grater than 0.
for( i in 1:length(mort.list)){
  mort.list[[i]][,death := ifelse(AGENTCD > 0, 1,  0)]
}

#Trees that are new (cross 5in threshold) during the remeasurement period don't count as surviving the measurement period. give them NA values.
for(i in 1:length(mort.list)){
  mort.list[[i]][,death := ifelse(PREV_TRE_CN > 0, death, NA)]
}


#get diameter in centimeters
for(i in 1:length(mort.list)){
  mort.list[[i]][,    DIA.cm :=     DIA*2.54]
  mort.list[[i]][,PREVDIA.cm := PREVDIA*2.54]
  mort.list[[i]][,inc.cm2.yr := (DIA.cm - PREVDIA.cm)/REMPER]
}

#must remove plots that don't have a previous plot number or a remper value of zero.
for(i in 1:length(mort.list)){
  mort.list[[i]] <- mort.list[[i]][!(is.na(PREV_PLT_CN)),]
  mort.list[[i]] <- mort.list[[i]][!(REMPER == 0),]
}

#Break out all and soils subset.
Product_2.soil <- mort.list[[1]]
Product_2      <- mort.list[[2]]

#Need to use previous PLT_CN values. Must also merge in relative abundance EM for downstream filtering.
#This will remove sites that didn't make it through Product_1 filtering, which is great.
Product_2      <- merge(Product_2     ,scaled.list[[4]][,.(relEM,relEM.AM,plot.BASAL,PLT_CN)], by = 'PLT_CN') 
Product_2.soil <- merge(Product_2.soil,  Product_1.soil[,.(relEM,relEM.AM,plot.BASAL,PLT_CN)], by.x = 'PREV_PLT_CN', by.y = 'PLT_CN')

#save the output! Tree-level mortality data paired with soils!
saveRDS(Product_2     , Product_2.path     )
saveRDS(Product_2.soil, Product_2.soil.path)

############################################################
#####      Product 3. Plot-level Recruitment Data     ######
############################################################
##****NOTE: Not sure if this needs to be its own product anymore. 
##Don't see why it can't just be tacked onto site level data in Product_1.

#In this section I calculate:
#A. Plot area from TPA_UNADJ column (expansion factors)
#B. Calculate total recruitment and by mycorrhizal type. 
#C. Aggregate individual tree data to the plot scale. Add site level data. 

###########################################################
#A. Plot area from TPA_UNADJ column (expansion factors)
###########################################################
#Note- all microplot observations of saplings have already been removed. 
#assuming they sample everything within the plot at the biggest resolution, the correct TPA_UNADJ is the smallest non-zero, non-NA TPA_UNADJ number.
recr.list <- list(data.list[[2]],data.list[[4]])
names(recr.list) <- c(names(data.list)[2], names(data.list)[4])

for(i in 1:length(recr.list)){
  recr.list[[i]][, TPA.fixed := ifelse(TPA_UNADJ == 0, NA, TPA_UNADJ)]
  recr.list[[i]][, TPA.fixed := ifelse(TPA_UNADJ == "", NA, TPA.fixed)]
  recr.list[[i]][, TPA.fixed := min(TPA.fixed, na.rm=T), by = PLT_CN]
}

#17 sites do not have a TPA_UNADJ values that is not NA. 239/23697 sites in all trees subset.
length(unique(recr.list[[2]]$PLT_CN))
length(unique(recr.list[[2]][is.na(TPA.fixed)]$PLT_CN))
for(i in 1:length(recr.list)){
  recr.list[[i]] <- recr.list[[i]][!(is.na(TPA.fixed)),]
  recr.list[[i]][,TPA.fixed := as.numeric(TPA.fixed)]
}

##calculate plot area as 1/TPA_UNADJ, which returns the area in acres. Convert to m2 by multiplying by 4046.86
for(i in 1:length(recr.list)){
  recr.list[[i]][,area.m2 := 1/TPA.fixed * 4046.86]
} 


###########################################################
#######   B. calculate recruitment, by MYC   ##############
###########################################################

#get a recruitment vector. Anything that does not have a PREV_TRE_CN value.
for (i in 1:length(recr.list)){
  recr.list[[i]][,   recruit := ifelse(is.na(PREV_TRE_CN), 1, 0)]
  recr.list[[i]][,recruit.em := ifelse(MYCO_ASSO == 'ECM',recruit, 0)]
  recr.list[[i]][,recruit.am := ifelse(MYCO_ASSO ==  'AM',recruit, 0)]
}

########################################################################
##### C. Get recruitment at the plot-scale, pair with other data. ######
########################################################################
r.scaled.list <- list()
for(i in 1:length(recr.list)){
  r.scaled.list[[i]] <- aggregate(recr.list[[i]]$recruit ~ recr.list[[i]]$PLT_CN, FUN = 'sum', na.rm = T, na.action = na.pass)
  colnames(r.scaled.list[[i]]) <- c('PLT_CN','recruit')
}
names(r.scaled.list) <- names(recr.list)

#pop in AM and EM recruit numbers, plot areas.
for(i in 1:length(r.scaled.list)){
  r.scaled.list[[i]]$recruit.am <- aggregate(recr.list[[i]]$recruit.am ~ recr.list[[i]]$PLT_CN, FUN = 'sum'   , na.rm=T, na.action = na.pass)[,2]
  r.scaled.list[[i]]$recruit.em <- aggregate(recr.list[[i]]$recruit.em ~ recr.list[[i]]$PLT_CN, FUN = 'sum'   , na.rm=T, na.action = na.pass)[,2]
  r.scaled.list[[i]]$area.m2    <- aggregate(recr.list[[i]]$area.m2    ~ recr.list[[i]]$PLT_CN, FUN = 'median', na.rm=T, na.action = na.pass)[,2]
}

#break out all and soil product for merging in final plot-scale stuff from Product_1
Product_3      <- r.scaled.list[[2]]
Product_3.soil <- r.scaled.list[[1]]

#merge in plot scale data and soils.
Product_3      <- merge(Product_3     ,scaled.list[[4]][,.(latitude,longitude,elevation,STDAGE,relEM,relEM.AM,plot.BASAL,n.trees,REMPER,INVYR,STATECD,PREV_PLT_CN,PLT_CN)], by = 'PLT_CN') 
Product_3.soil <- merge(Product_3.soil,scaled.list[[2]][,.(latitude,longitude,elevation,STDAGE,relEM,relEM.AM,plot.BASAL,n.trees,REMPER,INVYR,STATECD,PREV_PLT_CN,PLT_CN)], by = 'PLT_CN')
Product_3.soil <- merge(Product_3.soil,Soils, by.x = 'PREV_PLT_CN', by.y = "PLT_CN")

#save this output file. Gross basal increment of surviving trees at the plot scale, paired with soils!
saveRDS(Product_3     , Product_3.path     )
saveRDS(Product_3.soil, Product_3.soil.path)

##end script.
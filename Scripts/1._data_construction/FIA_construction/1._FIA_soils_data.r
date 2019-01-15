# CA working withi FIA Phase 3 soils data
# soils have two different layering structures:
# 1. FF_TOTAL / MIN_1 / MIN_2 <- organic horizon and mineral soil 0-10cm and 10-20cm depth. 
# 2. L_ORG / ORG_1 / ORG_2 <- organic horizon litter layer and then the soil organic horizon 0-10cm and 10-20cm depth. 
#clear R environment, load pacakges.
rm(list=ls())
library(data.table)
source('paths.r')

#where to save soil output file
out.file <- soil_data.processed.path

#kill scientific notation- this avoids problems with long PLT_CN values. 
options(scipen=999)

#load in soil chemical and location data
soil.chem <- data.table(read.csv(FIA.soil.chem.path))       #4997 unique PLT_CN values here.
soil.loc  <- data.table(read.csv(FIA.soil.loc.path))        #8195 unique PLT_CN values here. 


#averaging values within a soil layer by plot
#CA double checked how many rows should be there before/after. There are a few duplicate mineral and organic soils in addition to the expected multiple observations of FF_TOTAL. 
soil.chem<- aggregate(.~ PLT_CN  + LAYER_TYPE,data=soil.chem,FUN=mean,na.action=na.pass)
soil.chem<- soil.chem[order(soil.chem$PLT_CN),] #5403 unique plots

#calcualte missing BD values for MIN1, MIN2, ORG1, ORG2. This doesn't work for FF_TOTAL or L_ORG
test <- soil.chem[soil.chem$LAYER_TYPE %in% c('MIN_1','MIN_2','ORG_1','ORG_2'),]
test2<- soil.chem[soil.chem$LAYER_TYPE %in% c('FF_TOTAL','L_ORG'),]
#calculate BD for relevant soil layers
test$BULK_DENSITY<- ifelse(is.na(test$BULK_DENSITY)==T,test$OVEN_DRY_SOIL_WT/181,test$BULK_DENSITY)
#merge all all layers back together as soil.chem
soil.chem<-rbind(test,test2)
soil.chem<- soil.chem[order(soil.chem$PLT_CN),] #order data set by PLT_CN

#kill all observations missing soil C- 5283 sites remaining. 
to.remove <- soil.chem[is.na(soil.chem$C_ORG_PCT),] 
soil.chem<- soil.chem[!soil.chem$PLT_CN %in% to.remove$PLT_CN,]

#also kill observations that are missing bulk density for mineral or org1/2 horizons (this is 12 plots)
#1001 unique sites missing bulk density for L_ORG or FF_TOTAL
#keep so long as they have oven dry soil weight (only 14 sites excluded in this case)
to.remove <- soil.chem[is.na(soil.chem$BULK_DENSITY),] 
to.remove <- to.remove[to.remove$LAYER_TYPE %in% c('MIN_1','MIN_2','ORG_1','ORG_2'),] 
soil.chem <- soil.chem[!soil.chem$PLT_CN %in% to.remove$PLT_CN,] #down to 5,271 sites

#kill all observations with missing bulk density information that are FF_TOTAL or L_ORG horizons with no oven dry weight
#keep so long as they have oven dry soil weight (only 14 sites excluded in this case)
to.remove <- soil.chem[is.na(soil.chem$BULK_DENSITY),] 
to.remove <- to.remove[to.remove$LAYER_TYPE %in% c('L_ORG','FF_TOTAL'),] 
to.remove <- to.remove[is.na(to.remove$OVEN_DRY_SOIL_WT),]
soil.chem <- soil.chem[!soil.chem$PLT_CN %in% to.remove$PLT_CN,]


#calculate soil C and N storage in kg C / m2 by horizon
#forget depth reporting, try adding in FF and L_ORG horizons at end, after mineral calculation. 
# (%C/100) * bd (g/cm3) * depth (10 cm) * 10,000 cm2 / m2 * (1kg / 1,000g)
soil.chem$C.storage<- soil.chem$BULK_DENSITY * (soil.chem$C_ORG_PCT  /100) * 10 * 10000 * (1/1000)
soil.chem$N.storage<- soil.chem$BULK_DENSITY * (soil.chem$N_TOTAL_PCT/100) * 10 * 10000 * (1/1000)
#sub in the fact that if you are a FF_TOTAL or a L_ORG observation we will calculate this differently (aerial basis w/o depth or bulk density)
#(%C/100) * mass (g) * 1 / (SF cm2 (pi*15.24^2)) * 10,000 cm2 / m2 * (1kg / 1,000g)
soil.chem$C.storage <- ifelse(soil.chem$LAYER_TYPE %in% c('FF_TOTAL','L_ORG'),(soil.chem$C_ORG_PCT  /100)*soil.chem$OVEN_DRY_SOIL_WT*(1/(pi*15.24^2))*10000*(1/1000),soil.chem$C.storage)
soil.chem$N.storage <- ifelse(soil.chem$LAYER_TYPE %in% c('FF_TOTAL','L_ORG'),(soil.chem$N_TOTAL_PCT/100)*soil.chem$OVEN_DRY_SOIL_WT*(1/(pi*15.24^2))*10000*(1/1000),soil.chem$N.storage)

#sum across horizons to get a single soil C storage number per site. 
output           <-aggregate(soil.chem$C.storage~soil.chem$PLT_CN,FUN='sum',na.action=na.pass)
colnames(output) <- c('PLT_CN','C.storage')
output$N.storage <- aggregate(soil.chem$N.storage ~ soil.chem$PLT_CN, FUN='sum'   , na.action=na.pass)[,2]
output$VSTNBR    <- aggregate(soil.chem$VSTNBR    ~ soil.chem$PLT_CN, FUN='median', na.action=na.pass)[,2]
output$VSTNBR <- as.numeric(output$VSTNBR)

#sanity check plot
plot(C.storage~N.storage,data=output, pch=16, cex=0.3)
abline(lm(C.storage~N.storage,data=output),lty=2,lwd=3,col='purple')

#calculate mineral horizon micronutrient storage- mg nutrient / m2
mineral.chem            <- soil.chem[soil.chem$LAYER_TYPE %in% c('MIN_1','MIN_2'),]
mineral.chem$Na.storage <- mineral.chem$EXCHNG_NA * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$K.storage  <- mineral.chem$EXCHNG_K  * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$Mg.storage <- mineral.chem$EXCHNG_MG * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$Ca.storage <- mineral.chem$EXCHNG_CA * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$Al.storage <- mineral.chem$EXCHNG_AL * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$Mn.storage <- mineral.chem$EXCHNG_MN * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$Fe.storage <- mineral.chem$EXCHNG_FE * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$Ni.storage <- mineral.chem$EXCHNG_NI * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$Cu.storage <- mineral.chem$EXCHNG_CU * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$Zn.storage <- mineral.chem$EXCHNG_ZN * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$Cd.storage <- mineral.chem$EXCHNG_CD * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$Pb.storage <- mineral.chem$EXCHNG_PB * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$S.storage  <- mineral.chem$EXCHNG_S  * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$Bray.P.storage  <- mineral.chem$BRAY1_P * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000
mineral.chem$Olsen.P.storage <- mineral.chem$OLSEN_P * mineral.chem$BULK_DENSITY * 10 * (1/1000) * 10000

#aggregate mineral soil C, N and micronutrient storage, as well as depth and average pH. 
m.output              <- aggregate(mineral.chem$C.storage~mineral.chem$PLT_CN,FUN='sum')
colnames(m.output)    <- c('PLT_CN','m.C.storage')
m.output$m.N.storage  <- aggregate(mineral.chem$N.storage~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$pH_H2O       <- aggregate(mineral.chem$PH_H2O   ~mineral.chem$PLT_CN,FUN='mean',na.action=na.pass)[,2]
m.output$pH_CaCl2     <- aggregate(mineral.chem$PH_CACL2 ~mineral.chem$PLT_CN,FUN='mean',na.action=na.pass)[,2]
m.output$m.Na.storage <- aggregate(mineral.chem$Na.storage~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.K.storage  <- aggregate(mineral.chem$K.storage ~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.Mg.storage <- aggregate(mineral.chem$Mg.storage~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.Ca.storage <- aggregate(mineral.chem$Ca.storage~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.Al.storage <- aggregate(mineral.chem$Al.storage~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.Mn.storage <- aggregate(mineral.chem$Mn.storage~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.Fe.storage <- aggregate(mineral.chem$Fe.storage~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.Ni.storage <- aggregate(mineral.chem$Ni.storage~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.Cu.storage <- aggregate(mineral.chem$Cu.storage~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.Zn.storage <- aggregate(mineral.chem$Zn.storage~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.Cd.storage <- aggregate(mineral.chem$Cd.storage~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.Pb.storage <- aggregate(mineral.chem$Pb.storage~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.S.storage  <- aggregate(mineral.chem$S.storage ~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.BrayP.storage   <- aggregate(mineral.chem$Bray.P.storage  ~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]
m.output$m.OlsenP.storage  <- aggregate(mineral.chem$Olsen.P.storage ~mineral.chem$PLT_CN,FUN='sum',na.action=na.pass)[,2]

#merge mineral soil output table (m.output) with total soil data (output), inserting NAs if there is no data for the mineral soil. 
output <- merge(output,m.output,by='PLT_CN',all=T)

#save output as csv file.
write.csv(output,out.file)
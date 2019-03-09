#paths.r
#Main file dump directory.
host <- system('hostname', intern = T)
storage.dir <- '/projectnb/talbot-lab-data/caverill/paleo_myco_2018/'
if(host == 'pecan2'){storage.dir <- '/fs/data3/caverill/paleo_myco_2018/'}
if(host == 'Colins-MacBook-Pro-2.local'){storage.dir <- '/Users/colin/Desktop/paleo_myco_2018/'}
cmd <- paste0('mkdir -p ',storage.dir)
system(cmd)


#data product paths
    historic.composition.path <- paste0(storage.dir,'SetTreeComp_Level2_v1.0.nc'     )
contemporary.composition.path <- paste0(storage.dir,'fia_composition_midwest_v0.1.nc')
        historic.climate.path <- paste0(storage.dir,'paleon_models_environment_master.csv')
historic_contemporary_merge.path <- paste0(storage.dir,'paleon_historic_contemporary_merge.rds')

#model output paths
model.output.path <- paste0(storage.dir,'model_output.rds')

#JAGS model output paths.
historic_contemporary_EM_spatial.path <- paste0(storage.dir,'historic_contemporary_EM_spatial.path')
        
#FIA input paths
FIA7.dir.path <- '/fs/data3/caverill/FIA7/'
        FIAdb.path <- paste0(FIA7.dir.path,'FIA7.sqlite')
FIA.soil.chem.path <- paste0(FIA7.dir.path,'soils/SOILS_LAB.csv')
 FIA.soil.loc.path <- paste0(FIA7.dir.path,'soils/SOILS_SAMPLE_LOC.csv')
        
#FIA output paths
fia.dir <- paste0(storage.dir,'FIA_output/')
#soils
soil_data.processed.path <- paste0(fia.dir,'FIA7soil_output.csv')
#FIA with soils data, visit 1 and visit 2.
FIA_extraction_out.path        <- paste0(fia.dir,'soilC.FIA.out.rds')
FIA_extraction_out_FUTURE.path <- paste0(fia.dir,'soilC.FIA.out.FUTURE.rds')
#All FIA data, regardless of whether it has soil.
all.present.path <- paste0(fia.dir,'FIA.all.present.rds')
all.past1.path <- paste0(fia.dir,'FIA.all.past1.rds')
all.past2.path <- paste0(fia.dir,'FIA.all.past2.rds')
all.past3.path <- paste0(fia.dir,'FIA.all.past3.rds')

#products for downstream analysis
Product_1.path         <- paste0(fia.dir,"Product_1.rds")
Product_1.soil.path    <- paste0(fia.dir,"Product_1.soil.rds")
Product_2.path         <- paste0(fia.dir,"Product_2.rds")
Product_2.soil.path    <- paste0(fia.dir,"Product_2.soil.rds")
Product_3.path         <- paste0(fia.dir,"Product_3.rds")
Product_3.soil.path    <- paste0(fia.dir,"Product_3.soil.rds")
Product_1.all.path     <- paste0(fia.dir,"Product_1.all.rds")
time_series_dat.path   <- paste0(fia.dir,'time_series_dat.rds')

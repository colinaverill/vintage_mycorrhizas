#Querying the FIA database.
#Colin downloaded the most recent PLOT, COND, TREE_GRM_EST, and TREE tables from the FIA on May 12, 2017
#This uses RSQlite to query the respective tables, rather than pSQL as done previously by R. Kelly.
#This script takes a while to run (>20 minutes).
#clear environment, load packages.
rm(list=ls())
library(data.table)
library(RSQLite)
source('paths.r')

#Some functions to make life easier
sumNA  = function(x)  sum(x,na.rm=T)
meanNA = function(x) mean(x,na.rm=T)
maxNA  = function(x)  max(x,na.rm=T)
tic = function() assign("timer", Sys.time(), envir=.GlobalEnv)
toc = function() print(Sys.time()-timer)

#Connect to FIA7 database.
con <- dbConnect(SQLite(), dbname = FIAdb.path)

#File paths to data
file.pft = "required_products_utilities/gcbPFT.csv"
file.myc = "required_products_utilities/mycorrhizal_SPCD_data.csv"
file.soil = read.csv(soil_data.processed.path)
file.soil$plt_cn <- file.soil$PLT_CN

#Where to save filtered outputs. 
file.out        = FIA_extraction_out.path
file.out.future = FIA_extraction_out_FUTURE.path

#state codes
states <- read.csv('required_products_utilities/FIA_state_codes_regions.csv')
states <- states[states$paleon == 1,]$STATECD
states <- paste(states, collapse = ', ')
states <- paste0('(',states,')')

###---Query PLOT table
cat("Query PLOT...\n")
tic()
query <- paste("select 
                CN, STATECD, PREV_PLT_CN, REMPER, LAT, LON, ELEV 
                from 
                PLOT 
                where 
                STATECD IN",states)
PLOT = dbGetQuery(con, query)

PLOT = data.table(PLOT)
toc()
setnames(PLOT,"CN","PLT_CN")
states  = sort(unique(PLOT$STATECD))
n.state = length(states)

# Remove this one miscellaneous plot, per Trevor Andrews
PLOT = PLOT[ PLT_CN!= 134680578010854 ]


###---Query COND table
cat("Query COND...\n")
tic()
query <- paste("select 
                        PLT_CN, CONDID, STDORGCD, STDAGE, CONDPROP_UNADJ, AFFORESTATION_CD
                        from 
                        COND")
COND = dbGetQuery(con, query)
COND = data.table(COND)
toc()

###---Query SUBP_COND table
cat("Query SUBP_COND...\n")
tic()
SUBP_COND = dbGetQuery(con, "select 
                       CN, PLT_CN, CONDID, SUBP, SUBPCOND_PROP
                       from 
                       SUBP_COND")
SUBP_COND = data.table(SUBP_COND)
toc()

#Calculate forested condition proportion.
forest_prop <- SUBP_COND[CONDID == 1,]
forest_prop <- data.table(aggregate(SUBPCOND_PROP ~ PLT_CN, data = forest_prop, FUN = 'sum'))
forest_prop[,SUBPCOND_PROP := SUBPCOND_PROP / 4]
colnames(forest_prop)[2] <- 'forest_proportion'
#merge into COND table
COND <- merge(COND,forest_prop, all.x = T)

#need to have at least one plot with condition =1 (forested)
COND <- COND[CONDID == 1,]

#only keep plots that are 100% forested.
COND <- COND[COND$forest_proportion == 1,]

###---merge PLOT and COND tables
PC = merge(COND, PLOT, by="PLT_CN")
PC <- data.table(PC)

#remove PLOT, COND and SUBP_COND tables from memory.
rm(PLOT,COND, SUBP_COND)

#---Query TREE table
#Only query states in the PalEON domain.
#Colin did this because the full query started crashing pecan2. boring.
# 1. colin has removed p2a_grm_flg!=\'N\' from the query. This killed the west coast.
# 2. Onlty take trees that are alive, or were alive in previous census.
cat("Query TREE...\n")

#first, grab all PLT_CN values where soil was measured, as well as remeasurements
     PC$PLT_CN_filter <- as.numeric(gsub('"', "", PC$PLT_CN))
PC$PREV_PLT_CN_filter <- as.numeric(gsub('"', "", PC$PREV_PLT_CN))
a <- PC[     PLT_CN_filter %in% file.soil$PLT_CN,     PLT_CN] #PC for soil subset.
a.soil <- PC[     PLT_CN_filter %in% file.soil$PLT_CN,     PLT_CN] 
b.soil <- PC[PREV_PLT_CN_filter %in% file.soil$PLT_CN,     PLT_CN]
#of_interest <- data.frame(c(a,b))
#colnames(of_interest)<- c('test')
initial.soil <- data.frame(a.soil)
revisit.soil <- data.frame(b.soil)

#build an empty data.table to store output. 
variables_to_extract <- c('cn','prev_tre_cn','plt_cn','invyr','condid','dia','tpa_unadj','carbon_ag','carbon_bg',
                          'spcd','stocking','statuscd','prevdia','prev_status_cd','p2a_grm_flg','reconcilecd',
                          'agentcd','tpamort_unadj','diahtcd','ht','htcd','actualht','cclcd')
out <- data.frame(matrix(NA, nrow = 0, ncol = length(variables_to_extract)))
colnames(out) <- variables_to_extract
out <- data.table(out)
setnames(out, toupper(names(out)))

#RI is 44, CT is 9

#write a for loop, this should be dropped in parallel.
#really I should just query ones that match the file.soil PLT_CN vector. But. SQL queries hate me. So I'm doing this.
#querying based on the 'of_interest' PLT_CN values would probably speed this up a ton.
#could also query in parallel using the doParallel package and a foreach loop.
tic()
for(i in 1:length(states)){
  query = paste('select CN, PREV_TRE_CN, PLT_CN, INVYR, CONDID, DIA, TPA_UNADJ, CARBON_AG, CARBON_BG,
              SPCD, STOCKING, STATUSCD, PREVDIA, PREV_STATUS_CD, P2A_GRM_FLG, RECONCILECD, AGENTCD, TPAMORT_UNADJ,
              DIAHTCD, HT, HTCD, ACTUALHT, CCLCD
                from TREE 
                WHERE (PREVDIA>5 OR DIA>5) AND (STATUSCD=1 OR PREV_STATUS_CD=1) AND 
                STATECD IN (', paste(states[i],collapse=','), ')')
  pre.tree = as.data.table(dbGetQuery(con, query))
  out <- rbind(out,pre.tree)
  cat(paste0(i,' of ',length(states),' states queried.\n'));toc()
}
TREE <-out

#how many plots and trees?
length(unique(TREE$PLT_CN))
nrow(TREE)


# --- Filter TREE
cat("Filter TREE ...\n")
# By plot/cond criteria
TREE = TREE[ PLT_CN %in% PC$PLT_CN ]

# CONDID ("Remove edge effects" --TA)
TREE[, CONmax := maxNA(CONDID), by=PLT_CN]

# STATUSCD
# *** RK: Next line looks wrong. It's a sum, not max, despite the name. I did rewrite the line but this is equivalent to what Travis had so keeping for now.
TREE[, STATUSCDmax := sumNA(3*as.integer(STATUSCD==3)), by=PLT_CN]

# RECONCILECD. This is just signaling that a tree isn't a new tree to the plot.
TREE[is.na(RECONCILECD), RECONCILECD :=0] # Set NA values to 0 (unused)

# Filter
#TREE = TREE[ CONmax==1 & STATUSCDmax!=3 & STATUSCD!=0 & RECONCILECD<=4 ]
#STATUSCD=0 are remeasured trees that shouldn't be there or something.
#STATUSCD=3 means a tree was cut down by humans.
#RECONCILECD<=4 means these are acceptable reasons to have missed counting trees in the previous inventory.
TREE = TREE[STATUSCDmax!=3 & STATUSCD!=0 & RECONCILECD<=4 ]

#CALCULATE number of trees in a plot. 
TREE[,n.trees  := length(TPA_UNADJ), by=PLT_CN]

# --- Merge in PFTs and mycorrhizal associations
cat("Merge in PFTs and mycorrhizal associations...\n")
MCDPFT     = as.data.table(read.csv(file.pft, header = T)) 
CA_myctype = as.data.table(read.csv(file.myc, header = T)) #colin loads in mycorrhizal associations
CA_myctype = CA_myctype[,c("SPCD","MYCO_ASSO"),with=F]     #colin loads in mycorrhizal associations
TREE = merge(TREE, MCDPFT    , all.x=T, by = "SPCD")
TREE = merge(TREE, CA_myctype, all.x=T, by = "SPCD")


###subset to only include TREE sites with soil profiles, merge with complementary PC keys. 
setnames(TREE, 'CN','TRE_CN')
PC.soil   <-   PC[PLT_CN %in% initial.soil$a,]
TREE.soil <- TREE[PLT_CN %in% initial.soil$a,]

#Link together tempora sequences of trees and PC plots.
#n.past <- sum(!is.na(PC.present$PREV_PLT_CN))
#PC_time_series <- list()
#PC_time_series[[1]] <- PC[!(PLT_CN %in% PREV_PLT_CN),] #newest observations are not a PREV_PLT_CN of anything.
#n.past <- sum(!is.na(PC_time_series[[i]]$PREV_PLT_CN))
#i = 1
#Chain together past measurements.
#while(n.past > 0){
#  PC_time_series[[i+1]] <- PC[PLT_CN %in% PC_time_series[[i]]$PREV_PLT_CN]
#  n.past <- sum(!is.na(PC_time_series[[i+1]]$PLT_CN))
#  i = i + 1
#}
PC.present <- PC[!(PLT_CN %in% PREV_PLT_CN),] #newest observations are not a PREV_PLT_CN of anything.
PC.past1   <- PC[PLT_CN %in% PC.present$PREV_PLT_CN,]
PC.past2   <- PC[PLT_CN %in% PC.past1$PREV_PLT_CN,]
PC.past3   <- PC[PLT_CN %in% PC.past2$PREV_PLT_CN,]
TREE.present <- TREE[PLT_CN %in% PC.present$PLT_CN]
TREE.past1   <- TREE[PLT_CN %in% PC.past1$PLT_CN]
TREE.past2   <- TREE[PLT_CN %in% PC.past2$PLT_CN]
TREE.past3   <- TREE[PLT_CN %in% PC.past3$PLT_CN]
#merge together
all.present <- merge(TREE.present, PC.present, by = 'PLT_CN')
all.past1   <- merge(TREE.past1, PC.past1, by = 'PLT_CN')
all.past2   <- merge(TREE.past2, PC.past2, by = 'PLT_CN')
all.past3   <- merge(TREE.past3, PC.past3, by = 'PLT_CN')

#grab the set of future remeasurements of all soil FIA plots. 
PC.soil.future <-     PC[PLT_CN %in% revisit.soil$b,]
TREE.soil.future <- TREE[PLT_CN %in% revisit.soil$b,]


#merge current soils-trees-plot data
cat("Final merge...\n")
ALL = merge(TREE.soil, PC.soil, by='PLT_CN')
#kill observations that are actually remeasurements. 
#this happens because some sites have had soils measured multiple times. ~186 sites. 
ALL = ALL[!PLT_CN %in% PREV_PLT_CN,]

#merge future soils-trees-plot data
ALL.future = merge(TREE.soil.future, PC.soil.future, by='PLT_CN')

#save outputs
cat("Save.../n")
tic()
saveRDS(ALL       , file = file.out        )
saveRDS(ALL.future, file = file.out.future )
saveRDS(all.present,file = all.present.path)
saveRDS(all.past1  ,file = all.past1.path  )
saveRDS(all.past2  ,file = all.past2.path  )
saveRDS(all.past3  ,file = all.past3.path  )
toc()
###end script.
################################################################################
### Eve Eisemann -- 2022 
### Functions & data preparation for flux calculations for Pea Island, NC
### based on ~1400 elevation profiles from oregon inlet to Salvo NC
### profiles spaced 20 m apart
### 10 lidar datasets used beginning in 2005 ending in 2019
################################################################################

#install.packages('hrbrthemes')

library(tidyverse)
library(reshape)
library(viridis)
options(digits = 10)

## ## ## ---------------------------------------------- ## ## ##
   ## ## ------------------ FUNCTIONS ----------------- ## ##
## ## ## ---------------------------------------------- ## ## ##

## --- euclidean distance --- ##
euclid_dist <- function(a, b) sqrt(sum((a - b)^2)) #defining euc dist function (a b are (x y)vectors)

## --- adding distance column to data frame --- ##
dist_vector <- function(tab) {  
  tdist = rep(NULL, length(tab$OID))  ## empty tdist vector to place distance values in 
  
  for(j in min(tab$Trans_ID):max(tab$Trans_ID)){   ##cycling transects
    
    ind = (tab$Trans_ID == j)  # index locs
    sub = tab[ind,]  #subset only that trns based on index locs
    objs = (min(sub$OID):max(sub$OID))  ## min to max object id
    td = rep(NA, length(objs))  ## empty vector for tdist
    
    for (k in 1:length(objs)){  ## cycling points in transect
      ed0 = 0
      if(k == 1){ ed = 0 } else { ## start of transect = 0 m
        
        ed = euclid_dist(c(sub$POINT_X[sub$OID == objs[k]-1], sub$POINT_Y[sub$OID == objs[k]-1]), 
                         c(sub$POINT_X[sub$OID == objs[k]], sub$POINT_Y[sub$OID == objs[k]])) ##euclidean dist
        
        ed0 = td[sub$OID == objs[k-1]]
      }
      td[sub$OID == objs[k]] = ed+ed0 ##binding ed values together to form vector
    }
    ##placing new distance vector in to index locations for specific year and transect
    tdist[ind] = td
  }
  return(tdist)
}
## ----------------------------------------------------------------------------- ##


# 
# ## --- calculating area change between time a and time b (b - a) per transect --- ##
# transect_change <- function(df, d1, d2)  {
#   ## transect elev data times d1 (first e.g. tab$2005) and d2 (second e.g. tab$2012)
#   ## df indicates the shared transect IDs and points along lines (e.g. tab$Trans_ID)
#   ## Output will be a vector containing the sum of elev changes along each transect
#   ## can be appended as a new column in a dataframe with ntransects rows
#   
# #  dd1 = dist_vector(d1)
# #  dd2 = dist_vector(d2)  ## using previous function to calc dist
#   
# #  dd1 = d1
# #  dd2 = d2
#   
#   #trans = c(min(dd1$Trans_ID):max(dd1$Trans_ID))
#   trans = c(min(df):max(df))
#   totchvec  = vector(length = length(trans))  ## empty vector to fill
#   
#   for(t in 1:length(trans)){
#     
#     tt1 = d1[df == trans[t]]
#     tt2 = d2[df == trans[t]]
#     
#     ed = tt2 - tt1    # elevation differences
#     
#     ted = sum(ed, na.rm = TRUE)   #total elevation difference (sum along TRANSECT)
#     
#     totchvec[t] = ted
#   }
#   return(totchvec)
# }
# ## ----------------------------------------------------------------------------- ##

#d1 = "12"
#d2 = "14"



## --- calculating area change between time a and time b (b - a) per transect --- ##
transect_change <-  function(Pr, d1, d2){
  ## Pr - profile point data frame - including x and y points (POINT_X...), distances (tdist)
  ## d1 and d2 - abbreviated dates eg d1 = '12' or d1 = '19a'
  ## Output will be a vector containing the sum of elev changes along each transect
  ## can be appended as a new column in a dataframe with ntransects rows
  
  df = Pr$Trans_ID
  trans = c(min(df):max(df))
  totchvec  = vector(length = length(trans))  ## empty vector to fill
  
  for(t in 1:length(trans)){
    
    cn = colnames(Pr)
    tt1 = Pr[Pr$Trans_ID == trans[t],cn == paste("elev",d1, sep = "")]
    tt2 = Pr[Pr$Trans_ID == trans[t],cn == paste("elev",d2, sep = "")]
    
    ed = tt2 - tt1    # elevation differences
    
    ted = sum(ed, na.rm = TRUE)   #total elevation difference (sum along transect)
    
    totchvec[t] = ted
    
  }
  return(totchvec)
}
## ----------------------------------------------------------------------------- ##




## --- calculating x-sec area change behind the most landward dune of the two --- ##
transect_change_bd <- function(Pr, Fe, d1, d2)  {
  ## Fe - the big features data frame (includes Year - nnnn_ or nnnna/b, Type - including "FDC", Trans_ID, )
  ## Pr - profile point data frame - including x and y points (POINT_X...), distances (tdist)
  ## d1 and d2 - abbreviated dates eg d1 = '12' or d1 = '19a'
  ## Output will be a vector containing the sum of elev changes along each transect
  ## can be appended as a new column in a dataframe with ntransects rows
  
  #trans = c(min(dd1$Trans_ID):max(dd1$Trans_ID))
  df = Pr$Trans_ID
  trans = c(min(df):max(df))
  totchvec  = vector(length = length(trans))  ## empty vector to fill
  
  if(str_length(d1) == 3){fd1 = Fe[Fe$Year == paste("20",d1, sep="") & Fe$Type == "FDC",]} else {fd1 = Fe[Fe$Year == paste("20",d1,"_",sep="")& Fe$Type == "FDC",]}
  if(str_length(d2) == 3){fd2 = Fe[Fe$Year == paste("20",d2, sep="") & Fe$Type == "FDC",]} else {fd2 = Fe[Fe$Year == paste("20",d2,"_",sep="")& Fe$Type == "FDC",]}
  
  
  for(t in 1:length(trans)){
    Pr1x = Pr$POINT_X[Pr$Trans_ID == trans[t] & Pr$tdist == 0]
    Pr1y = Pr$POINT_Y[Pr$Trans_ID == trans[t] & Pr$tdist == 0]  ## first point in transect t 
    
    f1 = fd1[fd1$Trans_ID == trans[t],]    ## foredunes from time 1 
    f2 = fd2[fd2$Trans_ID == trans[t],]    ## foredunes from time 2 
   
    fds = rbind(f1, f2)  ## foredunes for time 1 and 2
    
    if(length(fds$tdist) == 1){fd = fds[fds$tdist != 0] ## IF one is empty, choose the other, else, pick the biggest one. 
      cn = colnames(Pr)
      tt1 = Pr[Pr$tdist>fd$tdist & Pr$Trans_ID == trans[t],cn == paste("elev",d1, sep = "")]
      tt2 = Pr[Pr$tdist>fd$tdist & Pr$Trans_ID == trans[t],cn == paste("elev",d2, sep = "")]
      
      ed = tt2 - tt1    # elevation differences
      
      ted = sum(ed, na.rm = TRUE)   #total elevation difference (sum along transect)
      
      totchvec[t] = ted
      
      }else if (length(fds$tdist) == 0) {
        ted = NA
        
      }else if (length(fds$tdist) > 2) {
        ted = NA
        
      }else {
        fds$dist = c(euclid_dist(c(f1$POINT_X, f1$POINT_Y), c(Pr1x, Pr1y)), 
                     euclid_dist(c(f2$POINT_X, f2$POINT_Y), c(Pr1x, Pr1y))) ## distances to the two foredunes for time 1 and time 2 
        
        
        ## whichever has the largest distance to the first transect point will be used as the seaward limit for volume calculation 
        fd = fds[fds$dist == max(fds$dist),]
        cn = colnames(Pr)
        tt1 = Pr[Pr$tdist>fd$tdist & Pr$Trans_ID == trans[t],cn == paste("elev",d1, sep = "")]
        tt2 = Pr[Pr$tdist>fd$tdist & Pr$Trans_ID == trans[t],cn == paste("elev",d2, sep = "")]
        
        ed = tt2 - tt1    # elevation differences
        
        ted = sum(ed, na.rm = TRUE)   #total elevation difference (sum along transect)
        
        totchvec[t] = ted
      }
  }
  return(totchvec)
}
## ----------------------------------------------------------------------------- ##



## --- importing features function  --- ##
get_features = function(years){  ## -- list of text years or identifier in files
  results = vector("list", length(years)) ## -- make list
  i = 1
  for (yr in years){
    ## -- import
    yr_features = read.table(paste("Features_",yr,"_FDC.txt", sep = ""), header = TRUE, sep = ',')
    yr_features$Year = yr
    
    ## -- store in list
    results[[i]] = yr_features
    i = i+1
  }
  do.call(rbind, results)
}
## ----------------------------------------------------------------------------- ##




## ------------- adding distance from transect start to features --------------- ##
dist_features <- function(tab,feats) {  ## -- tab: all transects tab (T) with tdist completed, feats: allFeatures tab 
  tdist = rep(NULL, length(feats$OID_))  ## -- empty tdist vector to place distance values in
  
  for(j in min(feats$Trans_ID):max(feats$Trans_ID)){   ##cycling transects
    
    ind = (feats$Trans_ID == j)  # index locs
    sub = feats[ind,]  #subset only that trns based on index locs
#    objs = (min(sub$OID):max(sub$OID))  ## min to max object id
    #td = rep(NA, length(sub$Type))  ## empty vector for tdist
    td = rep(NA, length(sub$OID_))# For the development setback
    
    first = tab[tab$Trans_ID == j & tab$tdist == 0,]  ## -- first point from transects to measure from
    
    for (k in 1:length(td)){  ## cycling points in transect
        ed = euclid_dist(c(first$POINT_X, first$POINT_Y), 
                         c(sub$POINT_X[k], sub$POINT_Y[k])) ##euclidean dist
        
        #ed0 = td[sub$OID == objs[k-1]]
      
      td[k] = ed ##binding ed values together to form vector
    }
    ##placing new distance vector in to index locations for specific year and transect
    tdist[ind] = td
  }
  return(tdist)
}
## ----------------------------------------------------------------------------- ##



## -------------shorelines from elevation profiles------------------------------- ##
## first create shorelines for each year (MHHW = 0.457 m NAVD88 - DUCK station)
## calc dist change, retreat rate between each 
## calc overall retreat rate 
## INPUTS - pr (PROFILES big elev point data frame), yr (e.g. '05' year) 
## Add to Fe (features) with new "type"

#yr = "05"   ## INPUT TO FUNC
#dfSUB = Pr[(Pr$Trans_ID <= 300 & Pr$Trans_ID >= 200),] ### Subsample for testing and better plotting 
#DF = dfSUB

shoreline_points <- function(yr, pr){

d1 = yr
DF = pr
trans = c(min(DF$Trans_ID):max(DF$Trans_ID))
print("Extracting seaward transect points")
print(paste("year: ", "20",yr))
print(paste("transects: ", min(trans), " to ", max(trans)))
SLptDF = data.frame(matrix(ncol = 6, nrow = 0)) ## empty data frame to fill 
colnames(SLptDF) = c("yr","yrelev","POINT_X", "POINT_Y", "Trans_ID", "tdist")

yrCol = data.frame(yr = rep(d1, length(DF$OID)),
                   yrelev = DF[,grep(paste("elev",d1, sep = ""),colnames(DF))], 
                   POINT_X = DF$POINT_X,
                   POINT_Y = DF$POINT_Y,
                   Trans_ID = DF$Trans_ID,
                   tdist = DF$tdist)

## choosing sl point first (smallest tdist) non NA yrelev  - most shoreward
for(t in 1:length(trans)){
  SLpt  = yrCol %>%
    filter(Trans_ID == trans[t]) %>%
    filter(!is.na(yrelev)) %>%
    filter(tdist == min(tdist))
  
  if(length(SLpt$yrelev) == 0){
    SLpt[1,] = c(d1,0, NA, NA, trans[t], NA)
  }
  
  SLptDF[t,] = SLpt  ## shoreline point data frame
}

SLptDF$yrelev = as.numeric(SLptDF$yrelev)
SLptDF$POINT_X = as.numeric(SLptDF$POINT_X)
SLptDF$POINT_Y = as.numeric(SLptDF$POINT_Y)
SLptDF$tdist = as.numeric(SLptDF$tdist)
return(SLptDF)
}

## ----------------------------------------------------------------------------- ##


## -----------------------using elevation to choose SL-------------------------- ##
## ----------------------- Not a function --
# SLptsA = SLptDF
# 
# ## choosing sl pt based on elevation
# for(t in 1:length(trans)){
#   SLpt = yrCol %>%
#     filter(Trans_ID == trans[t]) %>%
#     filter(tdist < 500) %>%
#     filter(yrelev <= 0.6)  ## 0.457 is mhhw elev used
#     
#   if (length(SLpt$yrelev) == 0){
#     print(paste("transect ",trans[t],": NO shoreline point"))
#     SLpt[1,] = c(d1,0, NA, NA, trans[t], NA)
#     #SLpt$yrelev = 0
#     #SLpt$POINT_X = NA
#     #SLpt$POINT_Y = NA
#     #SLpt$Trans_ID = t
#     #SLpt$tdist = NA
#   }else if (length(SLpt$yrelev)>1){
#     print(paste("transect ",trans[t],": MULTIPLE shoreline points, choosing shoreward"))
#     #SLpt = SLpt[SLpt$yrelev == min(SLpt$yrelev),]
#     SLpt = SLpt[SLpt$tdist == min(SLpt$tdist),]
#   }else{
#     print(paste("transect ", trans[t],": ONE shoreline point"))
#     SLpt = SLpt
#   }
#   #filter(yrelev == min(yrelev, na.rm = TRUE))
#   SLptDF[t,] = SLpt  ## shoreline point data frame
# }
# SLptsB = SLptDF

## ----------------------------------------------------------------------------- ##



## ## ## ---------------------------------------------- ## ## ##
## ## ## ---------------------------------------------- ## ## ##
## ## ## --------------- END FUNCTIONS ---------------- ## ## ##
## ## ## ---------------------------------------------- ## ## ##
## ## ## ---------------------------------------------- ## ## ##


setwd("~/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/R_analysis/FluxCalcs") ### mac
setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/R_analysis/FluxCalcs") ### pc




tab05 = read.table("points_2005__20_2000_1m.txt", header = TRUE, sep = ',')
#s05 = tab05[(tab05$Trans_ID >= 134 & tab05$Trans_ID <= 200), ]  #subset of 2005 transects (500-600)

tab09 = read.table("points_2009__20_2000_1m.txt", header = TRUE, sep = ',')

tab11 = read.table("points_2011__20_2000_1m.txt", header = TRUE, sep = ',')

tab12 = read.table("points_2012__20_2000.txt", header = TRUE, sep = ',')

tab14 = read.table("points_2014__20_2000.txt", header = TRUE, sep = ',')

tab16 = read.table("points_2016__20_2000.txt", header = TRUE, sep = ',')

tab17 = read.table("points_2017__20_2000.txt", header = TRUE, sep = ',')

tab18a = read.table("points_2018a_20_2000.txt", header = TRUE, sep = ',')

tab18b = read.table("points_2018b_20_2000.txt", header = TRUE, sep = ',')

tab19a = read.table("points_2019a_20_2000.txt", header = TRUE, sep = ',')

tab19b = read.table("points_2019b_20_2000.txt", header = TRUE, sep = ',')
#s19b = tab19b[(tab19b$Trans_ID >= 134 & tab19b$Trans_ID <= 200), ]  #subset of 2005 transects


## ## ## --- putting all together in one big dataframe (same points)
allTransects = data.frame(OID = tab05$OID_, Trans_ID = tab05$Trans_ID,POINT_X = tab05$POINT_X, POINT_Y = tab05$POINT_Y,
                     elev05 = tab05$RASTERVALU, elev09 = tab09$RASTERVALU, 
                     elev11 = tab11$RASTERVALU,
                     elev12 = tab12$RASTERVALU, 
                     elev14 = tab14$RASTERVALU, elev16 = tab16$RASTERVALU, elev17 = tab17$RASTERVALU, 
                     elev18a = tab18a$RASTERVALU, elev18b = tab18b$RASTERVALU, elev19a = tab19a$RASTERVALU, 
                     elev19b = tab19b$RASTERVALU)


## ## ## --- Remove all the big tables 
rm(tab05,tab09,tab11,tab12, tab14,tab16,tab17,tab18a,tab18b,tab19a,tab19b)


## ## ## ------- creating distance vector & adding distance to features
Pr = allTransects   ### all profiles
Pr$tdist = dist_vector(allTransects)  ## using transect distance (dist_vector) function 
#rm(allTransects)

## ## ## -----------------IMPORT FEATURE POINTS ------------------- ## ## ##
yrs = c("2005_","2009_","2011_","2012_","2014_","2016_","2017_","2018a","2018b","2019a","2019b")
#yrs = c("2005_","2011_","2012_","2014_","2016_","2017_","2018a","2018b","2019a","2019b")
allFeatures = get_features(yrs)  ## -- using new function 
## ## ## ------- ADDING DISTANCE TO FEATURES
Fe = allFeatures
Fe$tdist = dist_features(Pr,allFeatures)
#rm(allFeatures)


## ## ## ------- SUBSET option
# sampleProf = Tr[(Tr$Trans_ID >= 900 & Tr$Trans_ID <= 1000), ] 
# sampleFeat = Fe[(Fe$Trans_ID >= 900 & Fe$Trans_ID <= 1000), ] 

## ## ## -------- DEFINE TIMESLICE NAMES
CN = c("d05_09","d09_11","d11_12","d12_14","d14_16","d16_17","d17_18a","d18a_18b","d18b_19a","d19a_19b","d05_19b") 
#CN = c("d05_11","d11_12","d12_14","d14_16","d16_17","d17_18a","d18a_18b","d18b_19a","d19a_19b","d05_19b")

### --- survey dates
sdates = as.Date(c('24-8-2005','10-08-2009','28-8-2011','5-11-2012','6-01-2014','21-11-2016','24-6-2017',
                   '24-8-2018','2-10-2018','18-6-2019','26-9-2019'), '%d-%m-%Y')

## removing problem years - classification issues (2012 and 2018a)
#CN = c("d05_09","d09_11","d11_14","d14_16","d16_17","d17_18b","d18b_19a","d19a_19b","d05_19b") 

### ---- survey dates problem years removed
#sdates = as.Date(c('24-8-2005','10-08-2009','28-8-2011','6-01-2014','21-11-2016','24-6-2017',
#                  '2-10-2018','18-6-2019','26-9-2019'), '%d-%m-%Y')


i = (1:length(sdates))

dd = difftime(sdates[i+1],sdates[i],units='days')
diffdates = as.numeric(c(dd[1:max(i)-1], difftime(sdates[max(i)], sdates[1], units = 'days'))/365)

t1 = str_remove(substr(CN, 2,4), "_")
t2 = str_remove(substr(CN, 5,8), "_")

#diffdates_no09 = c(diffdates[1] + diffdates[2], diffdates[3:11])
## ## ## ---------------------------------------------------- ## ## ##




## ## ## ---------------------------------------------------- ## ## ##
## ## ## -------------- CALCULATION OF FLUXES -------- ------ ## ## ##
## ## ## ---------------------------------------------------- ## ## ##


## ## ## ---------------    total flux      ----------------- ## ## ##
## summed elevation change 
# TotFluxTab = data.frame(d05_09 = transect_change(Pr$Trans_ID, Pr$elev05, Pr$elev09),
#                            d09_11 = transect_change(Pr$Trans_ID, Pr$elev09, Pr$elev11),
#                            d11_12 = transect_change(Pr$Trans_ID, Pr$elev11, Pr$elev12),
#                            d12_14 = transect_change(Pr$Trans_ID, Pr$elev12, Pr$elev14),
#                            d14_16 = transect_change(Pr$Trans_ID, Pr$elev14, Pr$elev16),
#                            d16_17 = transect_change(Pr$Trans_ID, Pr$elev16, Pr$elev17),
#                            d17_18a = transect_change(Pr$Trans_ID, Pr$elev17, Pr$elev18a),
#                            d18a_18b = transect_change(Pr$Trans_ID, Pr$elev18a, Pr$elev18b),
#                            d18b_19a = transect_change(Pr$Trans_ID, Pr$elev18b, Pr$elev19a),
#                            d19a_19b = transect_change(Pr$Trans_ID, Pr$elev19a, Pr$elev19b))
# 
# TotFluxTab$d05_19b = transect_change(Pr$Trans_ID, Pr$elev05, Pr$elev19b)




## ## ## ----------   Total flux - full transect calculation ----------- ## ## ##
TotFluxTab = data.frame(matrix(ncol = length(CN),nrow = length((min(Pr$Trans_ID) : max(Pr$Trans_ID)))))
colnames(TotFluxTab) = CN
t1 = str_remove(substr(CN, 2,4), "_")
t2 = str_remove(substr(CN, 5,8), "_")


## ## ## ----- Calculating change & flux ------ ## ## ##
for(i in 1:length(CN)){
  print(paste("starting", t1[i], "to", t2[i], sep = " "))
  TotFluxTab[,i] = transect_change(Pr, t1[i], t2[i])
}

TotFluxTab$Transect = (min(Pr$Trans_ID) : max(Pr$Trans_ID))
#slices = colnames(TotFluxTab[1:11])
## ## ## ---------------------------------------------------- ## ## ##


## ## ## ----------   shoreface included - 0-500 m only ----------- ## ## ##
## ## ## ----------   max washover penetration from sandy is 500 m. this will remove some of the inclusion of error from veg. issues

TotFluxTab500 = data.frame(matrix(ncol = length(CN),nrow = length((min(Pr$Trans_ID) : max(Pr$Trans_ID)))))
colnames(TotFluxTab500) = CN
t1 = str_remove(substr(CN, 2,4), "_")
t2 = str_remove(substr(CN, 5,8), "_")

## ## ## ----- Calculating change & flux within 500 m of shoreline ------ ## ## ##
for(i in 1:length(CN)){
  print(paste("starting", t1[i], "to", t2[i], sep = " "))
  TotFluxTab500[,i] = transect_change(Pr[Pr$tdist<500,], t1[i], t2[i])
}

TotFluxTab500$Transect = (min(Pr$Trans_ID) : max(Pr$Trans_ID))
## ## ## ---------------------------------------------------- ## ## ##


## ## ## ---------- ------ behind dune      ---   ----------- ## ## ##
#BDtotFlux = transect_change_bd(Pr,Fe,"05","19b")
BDFluxTab = data.frame(matrix(ncol = length(CN),nrow = length((min(Pr$Trans_ID) : max(Pr$Trans_ID)))))
colnames(BDFluxTab) = CN
t1 = str_remove(substr(CN, 2,4), "_")
t2 = str_remove(substr(CN, 5,8), "_")

for(i in 1:length(CN)){
  paste("starting", t1[i], "to", t2[i], sep = " ")
  BDFluxTab[,i] = transect_change_bd(Pr, Fe, t1[i], t2[i])
}

BDFluxTab$Transect = (min(Pr$Trans_ID) : max(Pr$Trans_ID))

## ## ## -------------------- per year ------------------ ## ## ##

mat = as.matrix(BDFluxTab[,1:length(BDFluxTab-1)])  ## dataframe of fluxes as matrix for sweep function

BDFluxPY = as.data.frame(sweep(mat,2,diffdates,'/'))      ## BEHIND DUNE fluxes PER YEAR
BDFluxPY$Transect = BDFluxTab$Transect
## ## ## ---------------------------------------------------- ## ## ##





Zone = c(rep('Developed', 402), rep('NWR', 996))
SubZone = c(rep('Salvo', 136), rep('Waves', 144), rep('Rodanthe', 122), rep("S-Curves", 135), 
            rep("Stable Zone", 226), rep("New Inlet", 95), rep("Old Sandbag", 160), 
            rep("North Pond", 129), rep("Canal Zone", 133), rep("Inlet Zone", 118))






## ## ## ----------   behind dune - 0-500 m only ----------- ## ## ##
 

## ## ## ----- Calculating change & flux within 500 m of shoreline ------ ## ## ##
for(i in 1:length(CN)){
  print(paste("starting", t1[i], "to", t2[i], sep = " "))
  BDFluxTab500[,i] = transect_change_bd(Pr[Pr$tdist<500,], Fe, t1[i], t2[i])
}

BDFluxTab500$Transect = (min(Pr$Trans_ID) : max(Pr$Trans_ID))

## ## ## --------- Per Year ----------------------------------------------## ## ##
mat = as.matrix(BDFluxTab500[,c(1:ncol(BDFluxTab500)-1)])  ## dataframe of fluxes as matrix for sweep function

BDFPY500 = as.data.frame(sweep(mat,2,diffdates,'/'))      ## BEHIND DUNE fluxes PER YEAR
BDFPY500$Transect = BDFTab500$Transect
## ## ## ---------------------------------------------------- ## ## ##

BDFluxPY$Zone = Zone
BDFluxPY$SubZone = SubZone



## ## ## --------- -------------------------------- ------------------ ----------- ## ## ##  ADDING ZONES: aggregate FLUXES for zones and overall ------------- ## ## ## 

mFluxs = aggregate(BDFluxPY$d05_19b, list(BDFluxPY$SubZone), FUN = mean)
ORD = c(8,10,6,7,9,3,5,4,1,2)   ## ORDER of zones based on the auto alphabetical used in aggregate 
mFluxs = mFluxs[ORD,]



## ## ## ----------------------------------------------------------------------   ## ## ## CUMULATIVE FLUX ## ## ##
CN = c("d05_09","d09_11","d11_12","d12_14","d14_16","d16_17","d17_18a","d18a_18b","d18b_19a","d19a_19b","d05_19b") 
cn = CN[1:length(CN) - 1]
CumulativeVol500= data.frame(matrix(ncol = length(cn),nrow = length((min(BDFTab500$Transect) : max(BDFTab500$Transect)))))
colnames(CumulativeVol500) = cn
CumulativeVol500$Trans_ID = BDFTab500$Transect

for (i in 1:length(cn)){
  CumulativeVol500[,i] = apply(BDFTab500[1:i], 1, sum)
}


## ## ## ----------------------------------------------------------------------- ## ## ##  CUMULATIVE FLUX ALL ## ## ##

CumulativeVol= data.frame(matrix(ncol = length(cn),nrow = length((min(BDFTab$Transect) : max(BDFTab$Transect)))))
colnames(CumulativeVol) = cn
CumulativeVol$Trans_ID = BDFTab500$Transect

for (i in 1:length(cn)){
  CumulativeVol[,i] = apply(BDFTab[1:i], 1, sum)
}

CumulativeVol$Zone = Zone
CumulativeVol$SubZone = SubZone


## ## ## --------- -------------------------------- ------------------ ----------- ## ## ## Shorelines & Change -------------------- ## ## ##
yrss = c("05","09","11","12","14","16","17","18a","18b","19a","19b")

Shorelines = data.frame(matrix(ncol = 6,nrow = length(yrss)*length(min(Pr$Trans_ID) : max(Pr$Trans_ID))))
colnames(Shorelines) = c("yr","yrelev","POINT_X","POINT_Y","Trans_ID", "tdist")

Shorelines$yr[Shorelines$yr == "5"] = "05"   ## fixing the bad labeling for 05 and 09 ## CAUTION ##
Shorelines$yr[Shorelines$yr == "9"] = "09"

for (i in 1:length(yrss)){
  sl = shoreline_points(yrss[i], Pr)
  Shorelines = rbind(Shorelines, sl)
}



## ## ## --------- -------------------------------- ------------------ ----------- ## ## ## SHORELINE CHANGE ---  ## based on Tdist now (dist along transect)
ShorelineChange = data.frame(matrix(ncol = length(CN), nrow = max(Pr$Trans_ID)))
colnames(ShorelineChange) = c("Trans_ID", CN[1:length(CN) -1])
#colnames(ShorelineChange) = c(yrss)
ShorelineChange$Trans_ID = 1: max(Pr$Trans_ID)
ShorelineChange$POINT_X = Pr$POINT_X[Pr$tdist == 0] ## adding xy coords (start of transect)
ShorelineChange$POINT_Y = Pr$POINT_Y[Pr$tdist == 0]

for (i in 1:(length(CN) - 1)){
  YRS = unique(Shorelines$yr)
  yr1 = Shorelines %>% filter(yr == YRS[i]) 
  yr2 = Shorelines %>% filter(yr == YRS[i+1]) 

  for(ii in 1:(length(yr1$POINT_X))){
    ShorelineChange[ii,i+1] = yr1$tdist[ii] - yr2$tdist[ii]    ## negative is retreat, positive is growth 
  }
}

#ShorelineChange$total = rowSums(ShorelineChange[,c(2:11)], na.rm=TRUE)  ## adding total shoreline change

# 
# ## ## TESTING WITH ONE TRANSECT 
# SLCtest = data.frame(matrix(ncol = length(CN), nrow = 1))
# colnames(SLCtest) = c("Trans_ID", CN[1:length(CN) -1])
# #colnames(ShorelineChange) = c(yrss)
# SLCtest$Trans_ID = 1000
# SLCtest$POINT_X =  tst$POINT_X[1]## adding xy coords (start of transect)
# SLCtest$POINT_Y =  tst$POINT_Y[1]
#                      
# for (i in 1:(length(CN) - 1)){
#   YRS = unique(tst$yr)
#   yr1 = tst %>% filter(yr == YRS[i]) 
#   yr2 = tst %>% filter(yr == YRS[i+1])
#   
#   for(ii in 1:(length(yr1$POINT_X))){
#     SLCtest[ii,i+1] = yr1$tdist - yr2$tdist
#   }
# }
## ## ## ---------------------------------------------------- ## ## ##



## ## ## --------- -------------------------------- ------------------ ----------- ## ## ##        adding zones to allTransects ---------- ## ## ##
Pr$Zone = rep("NA", length(Pr$Trans_ID))
Pr$Zone[Pr$Trans_ID <= 402] = "Developed"
Pr$Zone[Pr$Trans_ID > 402] = "NWR"
Pr$SubZone = rep("NA", length(Pr$Trans_ID))
Pr$SubZone[Pr$Trans_ID <= 1398] = "Inlet Zone"
Pr$SubZone[Pr$Trans_ID <= 1280] = "Canal Zone"
Pr$SubZone[Pr$Trans_ID <= 1147] = "North Pond"
Pr$SubZone[Pr$Trans_ID <= 1018] = "Old Sandbag"
Pr$SubZone[Pr$Trans_ID <= 858] = "New Inlet"
Pr$SubZone[Pr$Trans_ID <= 763] = "Stable Zone"
Pr$SubZone[Pr$Trans_ID <= 537] = "S-Curves"
Pr$SubZone[Pr$Trans_ID <= 402] = "Rodanthe"
Pr$SubZone[Pr$Trans_ID <= 280] = "Waves"
Pr$SubZone[Pr$Trans_ID <= 136] = "Salvo"

# 
# BDFPY500$Zone = rep("NA", length(BDFPY500$Transect))
# BDFPY500$Zone[BDFPY500$Transect <= 402] = "Developed"
# BDFPY500$Zone[BDFPY500$Transect > 402] = "NWR"
# BDFPY500$SubZone = rep("NA", length(BDFPY500$Transect))
# BDFPY500$SubZone[BDFPY500$Transect <= 1398] = "Inlet Zone"
# BDFPY500$SubZone[BDFPY500$Transect <= 1280] = "Canal Zone"
# BDFPY500$SubZone[BDFPY500$Transect <= 1147] = "North Pond"
# BDFPY500$SubZone[BDFPY500$Transect <= 1018] = "Old Sandbag"
# BDFPY500$SubZone[BDFPY500$Transect <= 858] = "New Inlet"
# BDFPY500$SubZone[BDFPY500$Transect <= 763] = "Stable Zone"
# BDFPY500$SubZone[BDFPY500$Transect <= 537] = "S-Curves"
# BDFPY500$SubZone[BDFPY500$Transect <= 402] = "Rodanthe"
# BDFPY500$SubZone[BDFPY500$Transect <= 280] = "Waves"
# BDFPY500$SubZone[BDFPY500$Transect <= 136] = "Salvo"



## ## ## --------- -------------------------------- ------------------ ----------- ## ## ##  aggregating elevation statistics for zones and transects  ---------- ## ## ##
ZoneElevs = melt(Pr[,c(2,5:15, 18)], id = c("SubZone", "Trans_ID"))
ORD = c(8,10,6,7,9,3,5,4,1,2)                                                            ## ORDER of zones based on the auto alphabetical used in aggregate 
mElevs = aggregate(ZoneElevs$value, list(ZoneElevs$SubZone), FUN = mean, na.rm = TRUE)   ## Mean elevations of transects in ZONEs 
mElevs = mElevs[ORD, ] 
sds = aggregate(ZoneElevs$value, list(ZoneElevs$SubZone), FUN = sd, na.rm = TRUE)        ## Standard deviations 
sds = sds[ORD, ]  ## orders should be reverse of list above

ZoneElevStats = data.frame(Zone = mElevs$Group.1, MEAN = mElevs$x, SD = sds$x)           ## ZONE ELEV STATS in one data frame to save


mTElevs = aggregate(ZoneElevs$value, list(ZoneElevs$Trans_ID), FUN = mean, na.rm = TRUE) ## mean elevation for Each Transect all years (vs zones as above)
Tsds = aggregate(ZoneElevs$value, list(ZoneElevs$Trans_ID), FUN = sd, na.rm = TRUE)      ## sd for same
moTEelvs =  aggregate(ZoneElevs$value, list(ZoneElevs$Trans_ID), Mode)                   ## mode - provides a bunch for each transect, not really the peak

TransElevStats = data.frame(Trans_ID = c(min(ZoneElevs$Trans_ID):max(ZoneElevs$Trans_ID)), 
                            MEAN = mTElevs$x, SD = Tsds$x)                               ## TRANSECT STATS in one dataframe to save
rm(mTElevs,Tsds,moTEelvs)   ## removing
rm(ZoneElevs,mElevs,sds)    ## removing



## ## ## ------------------------------------------------------------------------ ## ## ## Creating elevation change cumulative data frame 
## ## ##           this includes cumulative elevation change at each point along the transect for looking at overwash penetration 


tst = Pr
tstsub = data.frame(Trans_ID = tst$Trans_ID, POINT_X = tst$POINT_X, POINT_Y = tst$POINT_Y, tdist = tst$tdist)
ECC = data.frame(Trans_ID = tst$Trans_ID, POINT_X = tst$POINT_X, POINT_Y = tst$POINT_Y, tdist = tst$tdist)
ECC$cumch = (tst$elev19b - tst$elev19a) + (tst$elev19a - tst$elev18b) + (tst$elev18b - tst$elev18a)+(tst$elev18a - tst$elev17)+(tst$elev17 - tst$elev16)+(tst$elev16 - tst$elev14)+(tst$elev14 - tst$elev12)+(tst$elev12 - tst$elev11)+(tst$elev11 - tst$elev09)+(tst$elev09 - tst$elev05)
ECC$sldist = rep(NA, length(ECC$Trans_ID), type = "num")
ECC$Zone = rep(NA, length(ECC$Trans_ID))
ECC$Zone[ECC$Trans_ID<=402] = "Developed"
ECC$Zone[ECC$Trans_ID>402] = "NWR"
ECC$SubZone = rep(NA, length(ECC$Trans_ID))
ECC$SubZone[ECC$Trans_ID <= 1398] = "Inlet Zone"
ECC$SubZone[ECC$Trans_ID <= 1280] = "Canal Zone"
ECC$SubZone[ECC$Trans_ID <= 1147] = "North Pond"
ECC$SubZone[ECC$Trans_ID <= 1018] = "Old Sandbag"
ECC$SubZone[ECC$Trans_ID <= 858] = "New Inlet"
ECC$SubZone[ECC$Trans_ID <= 763] = "Stable Zone"
ECC$SubZone[ECC$Trans_ID <= 537] = "S-Curves"
ECC$SubZone[ECC$Trans_ID <= 402] = "Rodanthe"
ECC$SubZone[ECC$Trans_ID <= 280] = "Waves"
ECC$SubZone[ECC$Trans_ID <= 136] = "Salvo"

ECC$mzbl = rep(NA, length(ECC$Trans_ID))
ECC$mzbl[ECC$Trans_ID <= 1398] = "inoi"
ECC$mzbl[ECC$Trans_ID <= 1272] = "nman"
ECC$mzbl[ECC$Trans_ID <= 1009] = "nmas"
ECC$mzbl[ECC$Trans_ID <= 851] = "inni"
ECC$mzbl[ECC$Trans_ID <= 819] = "nnat"
ECC$mzbl[ECC$Trans_ID <= 492] = "insc"
ECC$mzbl[ECC$Trans_ID <= 403] = "tman"
ECC$mzbl[ECC$Trans_ID <= 250] = "tnat"


mzbl.boundaries = c(250, 403, 492, 819, 851, 1009, 1272, 1398)    ## m zone boundaries lumped: more morpho behavior designations. lumping f and g together hummocky/small/no veg dunes in NWR
mzbl.bboundaries = c(0, 250, 403, 492, 819, 851, 1009, 1272) 

names = c("Waves & Salvo", "Rodanthe", "S-Curves Inlet", "NWR C", "New Inlet", "NWR B", "NWR A", "Oregon Inlet")

ECC$mzbl.name =  rep(NA, length(ECC$Trans_ID))
ECC$mzbl.order =  rep(NA, length(ECC$Trans_ID))
for (i  in c(1:length(mzbl.boundaries))) {
  ECC$mzbl.name[(ECC$Trans_ID > mzbl.bboundaries[i])&(ECC$Trans_ID <= mzbl.boundaries[i])] = names[i]
  ECC$mzbl.order[(ECC$Trans_ID > mzbl.bboundaries[i])&(ECC$Trans_ID <= mzbl.boundaries[i])] = i
}



sl19 = Shorelines[Shorelines$yr == "19b", ]    ### DISTANCE FROM SHORELINE
for (i in 1:length(Shorelines$Trans_ID)){
  
  ECC$sldist[ECC$Trans_ID == i] = ECC$tdist[ECC$Trans_ID == i] - sl19$tdist[sl19$Trans_ID == i]
  
}


## ## ## ------------------------------------------------------------------------ ## ## ## Adding SL distance to dune features (Fe)

Fe$sldist = rep(NA, length(Fe$OID_))

sl19 = Shorelines[Shorelines$yr == "19b", ]    ### DISTANCE FROM SHORELINE
for (i in min(Shorelines$Trans_ID):max(Shorelines$Trans_ID)){    ## min to max shoreline Trans ID, one point per transect
  
  Fe$sldist[Fe$Trans_ID == i] = Fe$tdist[Fe$Trans_ID == i] - sl19$tdist[sl19$Trans_ID == i]   
  ## all feature dist from beginning of transect for one transect minus the dist of shoreline in 2019
  
}

Fe$sldist[Fe$Trans_ID <6] = NA      ### a few of the first transects had erroneous features (in the ocean)
Fe$Distance[Fe$Trans_ID <6] = NA      ### a few of the first transects had erroneous features (in the ocean)
Fe$tdist[Fe$Trans_ID <6] = NA      ### a few of the first transects had erroneous features (in the ocean)

Fe$sldist[Fe$Type == "FD" & Fe$sldist >200] = NA      ### some way too far in foredunes

## ## ## ------------------------------------------------------------------------ ## ## ## Adding mzbl zones to dune features (Fe)

mzbl.boundaries = c(250, 403, 492, 819, 851, 1009, 1272, 1398)    ## m zone boundaries lumped: more morpho behavior designations. lumping f and g together hummocky/small/no veg dunes in NWR
mzbl.bboundaries = c(0, 250, 403, 492, 819, 851, 1009, 1272) 

names = c("Waves & Salvo", "Rodanthe", "S-Curves Inlet", "NWR C", "New Inlet", "NWR B", "NWR A", "Oregon Inlet")

Fe$mzbl.name =  rep(NA, length(Fe$Trans_ID))
Fe$mzbl.order =  rep(NA, length(Fe$Trans_ID))

for (i  in c(1:length(mzbl.boundaries))) {
  Fe$mzbl.name[(Fe$Trans_ID > mzbl.bboundaries[i])&(Fe$Trans_ID <= mzbl.boundaries[i])] = names[i]
  Fe$mzbl.order[(Fe$Trans_ID > mzbl.bboundaries[i])&(Fe$Trans_ID <= mzbl.boundaries[i])] = i
}


Fe$Distance[Fe$Type == "FD" & Fe$mzbl.name == "Rodanthe" & Fe$Distance > 125] = NA      ### some way too far in foredunes - scattery points <5 each 
Fe$Distance[Fe$Type == "FD" & Fe$mzbl.name == "S-Curves Inlet" & Fe$Distance > 125] = NA      ### some way too far in foredunes - scattery points <5 each 
Fe$Distance[Fe$Type == "FD" & Fe$mzbl.name == "New Inlet" & Fe$Distance > 250] = NA      ### some way too far in foredunes - scattery points <5 each 
Fe$Distance[Fe$Type == "FD" & Fe$mzbl.name == "NWR B" & Fe$Distance > 250] = NA      ### some way too far in foredunes - scattery points <5 each 
Fe$Distance[Fe$Type == "FD" & Fe$mzbl.name == "NWR A" & Fe$Distance > 250] = NA      ### some way too far in foredunes - scattery points <5 each 

### SUMMARY MEAN OF FOREDUNE Setback - respective years ###
Fe%>%
  filter(Type == "FD")%>%
  group_by(mzbl.order)%>%
  summarise(M = mean(Distance, na.rm = TRUE))%>%
  mutate(across(M,~num(.x,digits = 5)))


## ## ## ------------------------------------------------------------------------ ## ## ## Development intersects - 20 m transects version!

DevInt = read.table("developmentIntersects_20mTransects.txt", header = TRUE, sep = ',')

DevInt$tdist = dist_features(Pr,DevInt)   ### Distance from transect start to each development intersect

DevInt$sldist = rep(NA, length(DevInt$OID_))
sl19 = Shorelines[Shorelines$yr == "19b", ]    ### DISTANCE FROM 2019b SHORELINE
for (i in min(Shorelines$Trans_ID):max(Shorelines$Trans_ID)){    ## min to max shoreline Trans ID, one point per transect
  DevInt$sldist[DevInt$Trans_ID == i] = DevInt$tdist[DevInt$Trans_ID == i] - sl19$tdist[sl19$Trans_ID == i]   
  ## all feature dist from beginning of transect for one transect minus the dist of shoreline in 2019
}

DevInt$mzbl.name =  rep(NA, length(DevInt$Trans_ID))
DevInt$mzbl.order =  rep(NA, length(DevInt$Trans_ID))

for (i  in c(1:length(mzbl.boundaries))) {
  DevInt$mzbl.name[(DevInt$Trans_ID > mzbl.bboundaries[i])&(DevInt$Trans_ID <= mzbl.boundaries[i])] = names[i]
  DevInt$mzbl.order[(DevInt$Trans_ID > mzbl.bboundaries[i])&(DevInt$Trans_ID <= mzbl.boundaries[i])] = i
}

DevInt$DevIntCount =  rep(NA, length(DevInt$Trans_ID))
for(i in c(min(DevInt$Trans_ID):max(DevInt$Trans_ID))){
  DevInt$DevIntCount[DevInt$Trans_ID == i] = length(DevInt$tdist[DevInt$Trans_ID == i])
}


### SUMMARY MEAN OF FIRST DEV INTERSECT ###
DevInt%>%
  group_by(Trans_ID,mzbl.order)%>%summarise(msld = min(sldist))%>%
  group_by(mzbl.order)%>%
  summarise(M = mean(msld, na.rm = TRUE))%>%
  mutate(across(M,~num(.x,digits = 5)))

## ## ## ---------------------------------------------------- ## ## ##
## ## ## -------------------SAVE FILES----------------------- ## ## ##
write.table(Pr, file = "Profiles.txt", sep = ",", col.names = TRUE)                   ## all profile points and elevs, includes zones and sub zones now (9/23)
write.table(ZoneElevStats, file = "ZoneElevStats.txt", sep = ",", col.names = TRUE)   ## elevation statistics by zone, all years (9/23)
write.table(TransElevStats, file = "TransElevStats.txt", sep = ",", col.names = TRUE) ## elevation statistics by transect, all years (9/23)
write.table(Fe, file = "Features.txt", sep = ",", col.names = TRUE)                   ## all features locations and elevs
write.table(BDFluxPY, file="BDFluxPY.txt", sep = ",", col.names = TRUE)               ## m^3/m/yr behind dune 
write.table(BDFluxTab, file="BDFluxTab.txt", sep = ",", col.names = TRUE)             ## m^3/m behind dune
write.table(BDFluxTab500, file="BDFluxTab500.txt", sep = ",", col.names = TRUE)       ## m^3/m dune to 500 m along transect
write.table(BDFPY500, file="BDFluxPY500.txt", sep = ",", col.names = TRUE)            ## m^3/m/yr dune to 500 m along transect
write.table(CumulativeVol500, file = "CumulativeVol500.txt", sep = ",", col.names = TRUE)     ## m^3/m dune to 500 m cumulative vols
write.table(CumulativeVol, file = "CumulativeVol.txt", sep = ",", col.names = TRUE)     ## m^3/m dune to 500 m cumulative vols
write.table(TotFluxTab, file="TotFluxTab.txt", sep = ",", col.names = TRUE)           ## m^3/m whole transect
write.table(TotFluxTab500, file="TotFluxTab500.txt", sep = ",", col.names = TRUE)     ## m^3/m shoreline to 500 m along transect 
write.table(Shorelines, file = "Shorelines.csv", sep = ",", col.names = TRUE)         ## Shorelines 
write.table(ShorelineChange, file = "ShorelineChange.csv", sep = ",", col.names = TRUE)  ## Shoreline Change!
write.table(ECC, file = "ElevChangeCumulative.csv", sep = ",", col.names = TRUE)   ## Elevation change at each point cumulatively 
write.table(DevInt, file = "DevelopmentIntersects.csv", sep = ",", col.names = TRUE)   ## Development Intersects with Tdist etc

## ## ## ---------------------------------------------------- ## ## ##


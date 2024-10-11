################################################################################
### Eve Eisemann -- 08/2023
### calculation of long term fluxes and equilibrium overwash flux
### files to load were created from DataPrepForFlux.R
################################################################################



library(tidyverse)
library(reshape)
library(viridis)
library(dplyr)

# Install
install.packages("wesanderson")
# Load
library(wesanderson)

options(digits = 10)

## ## ## --- MODE FUNC ----- ## ## ##
Mode <- function(x) {
  ux <- na.omit(unique(x) )
  tab <- tabulate(match(x, ux)); ux[tab == max(tab) ]
}

## ## ## --------- -------------------------------- ------------------ ----------- ## ## ##             IMPORTING FILES -------------------- ## ## ##
setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/R_analysis/FluxCalcs") ### pc
setwd("~/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/R_analysis/FluxCalcs") ### mac

Pr = read.table("Profiles.txt", header = TRUE, sep = ',')            ## raw profile elevation data for each lidar dataset (bias corrected)
Fe = read.table("Features.txt", header = TRUE, sep = ',')            ## dune features 
BDFPY = read.table("BDFluxPY.txt", header = TRUE, sep = ',')         ## PER YEAR = true flux : behind dune to extent of data, with zones
BDFPY500 = read.table("BDFluxPY500.txt", header = TRUE, sep = ",")   ## PER YEAR = true flux : behind dune to 500 m inland
BDFTab = read.table("BDFluxTab.txt", header = TRUE, sep = ',')       ## VOLUME CHANGE per meter shoreline (transect change) behind dune to extent of data
BDFTab500 = read.table("BDFluxTab500.txt", header = TRUE, sep = ',') ## VOLUME CHANGE per meter shoreline (transect change) behind dune to 500 m inland 
CumulativeVol500 = read.table("CumulativeVol500.txt", header = TRUE, sep = ',')
CumulativeVol = read.table("CumulativeVol.txt", header = TRUE, sep = ',')
TFTab = read.table("TotFluxTab.txt", header = TRUE, sep = ',')       ## Total VOLUME CHANGE per meter shoreline (transect change) shoreline to extent of data
TFTab500 = read.table("TotFluxTab500.txt", header = TRUE, sep = ',') ## Total VOLUME CHANGE per meter shoreline (transect change) shoreline to 500 m inland
Shorelines = read.table("Shorelines.csv", header = TRUE, sep = ',')  ## shoreline location points from first point above MHW on transect 
ShorelineChange = read.table("ShorelineChange.csv", header = TRUE, sep = ",") ## distance between points from each timestep ## negative is retreat, positive is growth 
ZoneElevStats = read.table("ZoneElevStats.txt", header = TRUE, sep = ",")     ## elevation statistics by zone, all years (9/23)
TransElevStats = read.table("TransElevStats.txt", header = TRUE, sep = ",")   ## elevation statistics by transect, all years (9/23)
ECC = read.table("ElevChangeCumulative.csv", header = TRUE, sep = ',')        ## elevation changes cumulative for each point on transects (m)
## ## ## ----------------------------------------------------- ## ## ##

sdates = as.Date(c('24-8-2005','10-08-2009','28-8-2011','5-11-2012','6-01-2014','21-11-2016','24-6-2017',
                   '24-8-2018','2-10-2018','18-6-2019','26-9-2019'), '%d-%m-%Y')
i = (1:length(sdates))
dd = difftime(sdates[i+1],sdates[i],units='days')
diffdates = as.numeric(c(dd[1:max(i)-1], difftime(sdates[max(i)], sdates[1], units = 'days'))/365)
textdiffdates = colnames(BDFPY)[1:11]

# ## ## ## --------- -------------------------------- ------------------ ----------- ## ## ##   aggregate FLUXES for zones and overall ------------- ## ## ##  all this stuff commented got moved to data prep for flux.r
Zone = c(rep('Developed', 402), rep('NWR', 996))
SubZone = c(rep('Salvo', 136), rep('Waves', 144), rep('Rodanthe', 122), rep("S-Curves", 135), 
            rep("Stable Zone", 226), rep("New Inlet", 95), rep("Old Sandbag", 160), 
            rep("North Pond", 129), rep("Canal Zone", 133), rep("Inlet Zone", 118))
SubZoneCode = c(rep('A', 136), rep('B', 144), rep('C', 122), rep("D", 135), 
               rep("E", 226), rep("F", 95), rep("G", 160), 
               rep("H", 129), rep("I", 133), rep("J", 118))

zoneboundaries = c(136,280,402,537,763,858,1018,1147,1208,1398 )
zonenames = c("Salvo",  "Waves","Rodanthe","S-Curves", "Stable Zone", "New Inlet",   "Old Sandbag",
               "North Pond",  "Canal Zone"  ,"Inlet Zone")

mzoneboundaries = c(0,219, 403, 492, 819, 916, 1009, 1300, 1398) 
mzonenames = c("DD","DND","NND","ND","NND","ND","NND","inlet")  
mzonenames2 =  c("a_DD","b_DND","c_NND","d_ND","e_NND","f_ND","g_NND","h_inlet")  
## Developed + dune
## Developed + no dune
## NWR + dune
## NWR + no dune
## inlet zone (north of bridge connection to jetty)

i = (1:length(mzonenames))
mzone = c(rep(mzonenames[i],(mzoneboundaries[i+1]-mzoneboundaries[i])))
mzone2 = c(rep(mzonenames2[i],(mzoneboundaries[i+1]-mzoneboundaries[i])))

## ## ## ------------ ADDING MZONES --------------------- ## ## ## -------- ## ## ## ------------ ADDING MZONES --------10/2023------------- ## ## ## 

CumulativeVol$mZone = mzone
CumulativeVol500$mZone = mzone
CumulativeVol$mZone2 = mzone2
CumulativeVol500$mZone2 = mzone2

#BDFPY$Zone = Zone
#BDFPY$SubZone = SubZone



## ## ## ------- ---- ----  ---  ---ADDING ZONES TO Fe ---- --- ---- --------------------------- ## ## ## 

Fe$mZone = rep(NA, length(Fe$Trans_ID))
Fe$mZone2 = rep(NA, length(Fe$Trans_ID))


for( i in (1:length(mzonenames2))){
  
  Fe$mZone[Fe$Trans_ID > mzoneboundaries[i] & Fe$Trans_ID <= mzoneboundaries[i+1]] = 
    rep(mzonenames[i], length(Fe$mZone[Fe$Trans_ID > mzoneboundaries[i] & Fe$Trans_ID <= mzoneboundaries[i+1]]))
  
  Fe$mZone2[Fe$Trans_ID > mzoneboundaries[i] & Fe$Trans_ID <= mzoneboundaries[i+1]] = 
    rep(mzonenames2[i], length(Fe$mZone[Fe$Trans_ID > mzoneboundaries[i] & Fe$Trans_ID <= mzoneboundaries[i+1]]))
  
}

## ## # #plotting zone feature distribs. 
#pal = wes_palette("FantasticFox1", n = 5)

x11()
ggplot(Fe[Fe$Type == "HD" & Fe$mZone != "inlet",], aes(mZone, Elevation, col = mZone))+ ### ------------------ Plotting dune heights in zones 
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_boxplot(outlier.shape = NA)+
  #geom_jitter(alpha = 0.3)+
  #facet_wrap(~mZone2, scales = 'free_y')+
  scale_fill_manual(values= wes_palette("GrandBudapest1", n = 4))+
  #scale_color_manual(values=pal)+
 # scale_y_continuous(limits = c(-20,20))+
  ggtitle("Morpho Zone Foredune Elevations")+
  labs(y = "ELevation (NAVD 88)", x =" " )+
  guides(color = guide_legend(title = "Morpho Zone"))




# allTransects = Pr
# 
# ## ## ## --------- -------------------------------- ------------------ ----------- ## ## ##        adding zones to allTransects ---------- ## ## ##
# allTransects$Zone = rep("NA", length(allTransects$Trans_ID))
# allTransects$Zone[allTransects$Trans_ID <= 402] = "Developed"
# allTransects$Zone[allTransects$Trans_ID > 402] = "NWR"
# allTransects$SubZone = rep("NA", length(allTransects$Trans_ID))
# allTransects$SubZone[allTransects$Trans_ID <= 1398] = "Inlet Zone"
# allTransects$SubZone[allTransects$Trans_ID <= 1280] = "Canal Zone"
# allTransects$SubZone[allTransects$Trans_ID <= 1147] = "North Pond"
# allTransects$SubZone[allTransects$Trans_ID <= 1018] = "Old Sandbag"
# allTransects$SubZone[allTransects$Trans_ID <= 858] = "New Inlet"
# allTransects$SubZone[allTransects$Trans_ID <= 763] = "Stable Zone"
# allTransects$SubZone[allTransects$Trans_ID <= 537] = "S-Curves"
# allTransects$SubZone[allTransects$Trans_ID <= 402] = "Rodanthe"
# allTransects$SubZone[allTransects$Trans_ID <= 280] = "Waves"
# allTransects$SubZone[allTransects$Trans_ID <= 136] = "Salvo"
# 
# 
# ## ## ## --------- -------------------------------- ------------------ ----------- ## ## ##  aggregating elevation statistics for zones and transects  ---------- ## ## ##
# ZoneElevs = melt(allTransects[,c(2,5:15, 18)], id = c("SubZone", "Trans_ID"))
# ORD = c(8,10,6,7,9,3,5,4,1,2)                                                            ## ORDER of zones based on the auto alphabetical used in aggregate 
# mElevs = aggregate(ZoneElevs$value, list(ZoneElevs$SubZone), FUN = mean, na.rm = TRUE)   ## Mean elevations of transects in ZONEs 
# mElevs = mElevs[ORD, ] 
# sds = aggregate(ZoneElevs$value, list(ZoneElevs$SubZone), FUN = sd, na.rm = TRUE)        ## Standard deviations 
# sds = sds[ORD, ]  ## orders should be reverse of list above
# 
# ZoneElevStats = data.frame(Zone = mElevs$Group.1, MEAN = mElevs$x, SD = sds$x)
# 
# mTElevs = aggregate(ZoneElevs$value, list(ZoneElevs$Trans_ID), FUN = mean, na.rm = TRUE) ## mean elevation for Each Transect all years (vs zones as above)
# Tsds = aggregate(ZoneElevs$value, list(ZoneElevs$Trans_ID), FUN = sd, na.rm = TRUE)      ## sd for same
# moTEelvs =  aggregate(ZoneElevs$value, list(ZoneElevs$Trans_ID), Mode)                   ## mode - provides a bunch for each transect, not really the peak
# 
# TransElevStats = data.frame(Trans_ID = c(min(ZoneElevs$Trans_ID):max(ZoneElevs$Trans_ID)), MEAN = mTElevs$x, SD = Tsds$x)


mFlux =  (aggregate(CumulativeVol$d19a_19b, list(CumulativeVol$mZone2), FUN = mean, na.rm = TRUE)) ## ## ## Mean fluxes per zones - no truncation 
mFlux$Flux = mFlux$x/diffdates[11]

CumulativeVol500$Zone = Zone
CumulativeVol500$SubZone = SubZone
mFlux500 = (aggregate(CumulativeVol500$d19a_19b, list(CumulativeVol500$mZone), FUN = mean, na.rm = TRUE))## ## ## mean fluxes per zones - 500 m truncation 
mFlux500$Flux = mFlux500$x/diffdates[11]


## ## ## --------- -------------------------------- ------------------ ----------- ## ## ##     PLOT ELEVATIONS IN ZONE FOR ALL YEARS --- ## ## ##
ZoneElevs = melt(Pr[,c(2,5:15, 18)], id = c("SubZone", "Trans_ID"))

x11()
p = ggplot(ZoneElevs[ZoneElevs$SubZone != "Inlet Zone",], aes(x = reorder(SubZone, value, median, na.rm=TRUE), y = value)) +
  #geom_jitter(width = 100, height = 0.5, size = 0.9) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = c(0,3), name = 'Island Elevations') +
  labs(y= "elevation (m NAVD88)", x = "regions")+
  # scale_x_date(date_breaks = "1 year", date_labels = "%Y")
  ggtitle("Island Elevations (all elevs along each transect) 2005 - 2019")
p


## ## ## --------- ---------------- ---------------- ------------------ ----------- ## ## ##  CALCULATE Qow based on TRANSECT RETREAT RATES --- ## ## ##
## ## ## --------- ---------------- ---------------- ------------------ ----------- ## ## ##     1. SHORELINE CHANGE MEASURED FROM LIDAR ## ## ## 
## ## ## --------- ---- QOW = (Dp + Hb)*Rs 


## ## ## -- SHORELINE CHANGE melt-- ## ## ## 
ShorelineChange$Zone = Zone
ShorelineChange$SubZone = SubZone
ShorelineChange$SubZoneCode = SubZoneCode
ShorelineChange$mZone = mzone
ShorelineChange$mZone2 = mzone2
SLCmelt = melt(ShorelineChange, id = c('Trans_ID', 'POINT_X', 'POINT_Y', 'Zone', 'SubZone', 'SubZoneCode', 'mZone', 'mZone2'))  ## ShoreLineChange 
SLCmelt$date = rep(NA,length(SLCmelt$Trans_ID))
SLCmelt$Rate = rep(NA,length(SLCmelt$Trans_ID))  ## dividing change by decimal years from a to b (rate)
SLCmelt$cumulativeRetreat = rep(NA,length(SLCmelt$Trans_ID))
SLCmelt$linearRate = rep(NA,length(SLCmelt$Trans_ID))
SLCmelt$rSquared = rep(NA,length(SLCmelt$Trans_ID))
Labs = (unique(SLCmelt$variable))

#define function to extract overall p-value of model
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

for(i in c(1:10)){   ## adding decimal years from first survey to date
  SLCmelt$date[SLCmelt$variable == Labs[i]] = cumsum(c(diffdates))[i]
  SLCmelt$Rate[SLCmelt$variable == Labs[i]] = SLCmelt$value[SLCmelt$variable == Labs[i]]/diffdates[i]
}

for(i in 1:length(SLCmelt$Trans_ID)){  ## adding cumulative slc                                      ## and linear regression sl change rate
  SLCmelt$cumulativeRetreat[SLCmelt$Trans_ID == i] = cumsum(SLCmelt$value[SLCmelt$Trans_ID == i])
  lengthna = sum(is.na(SLCmelt$value[SLCmelt$Trans_ID == i]))  ## if there are fewer than three retreat values na
  if((length(SLCmelt$value[SLCmelt$Trans_ID == i])-lengthna) < 3){
    SLCmelt$linearRate[SLCmelt$Trans_ID == i] = NA
    SLCmelt$rSquared[SLCmelt$Trans_ID == i] = NA
  }else{  ## if more than three slc values, create linear regression and save slope and r sq
    mod  = lm(SLCmelt$cumulativeRetreat[SLCmelt$Trans_ID == i] ~ 0 + SLCmelt$date[SLCmelt$Trans_ID == i])    ### the 0 locks the trend line through 0,0
    coeff = coef(mod)
    SLCmelt$linearRate[SLCmelt$Trans_ID == i] = coeff[1]   ## meters per year  -- now coeff 1 b/c no longer intercept needed with locked 0.0
    SLCmelt$rSquared[SLCmelt$Trans_ID == i] = summary(mod)$adj.r.squared
    SLCmelt$Pval[SLCmelt$Trans_ID == i] = overall_p(mod)
  }
}

SLCmelt$linearRate[SLCmelt$Trans_ID == 1] = NA   ## ## ## being lazy to get rid of some bad points 
SLCmelt$linearRate[SLCmelt$Trans_ID == 2] = NA
SLCmelt$linearRate[SLCmelt$Trans_ID == 3] = NA
SLCmelt$linearRate[SLCmelt$Trans_ID == 4] = NA
SLCmelt$linearRate[SLCmelt$Trans_ID == 1380] = NA
SLCmelt$linearRate[SLCmelt$Trans_ID == 1381] = NA
SLCmelt$linearRate[SLCmelt$Trans_ID == 1382] = NA
SLCmelt$linearRate[SLCmelt$Trans_ID == 1383] = NA

SLCmelt$rSquared[SLCmelt$Trans_ID == 1] = NA
SLCmelt$rSquared[SLCmelt$Trans_ID == 2] = NA
SLCmelt$rSquared[SLCmelt$Trans_ID == 3] = NA
SLCmelt$rSquared[SLCmelt$Trans_ID == 1380] = NA
SLCmelt$rSquared[SLCmelt$Trans_ID == 1381] = NA
SLCmelt$rSquared[SLCmelt$Trans_ID == 1382] = NA
SLCmelt$rSquared[SLCmelt$Trans_ID == 1383] = NA


aggSLCmelt = SLCmelt %>% group_by(SubZone,date,variable) %>% summarise(mean_SLC = mean(cumulativeRetreat, na.rm = TRUE))  ##mean slc values for zones
aggSLC_meanval = aggSLCmelt %>% as.data.frame()


aggSLCmeltrate = SLCmelt %>% group_by(SubZone,Zone) %>% summarise(mean_SLC = mean(linearRate, na.rm = TRUE))  ##mean slc values for zones
aggSLC_meanvalrate = aggSLCmeltrate %>% as.data.frame()


aggSLCmeltM = SLCmelt %>% group_by(mZone2,date,variable) %>% summarise(mean_SLC = mean(cumulativeRetreat, na.rm = TRUE))  ##mean slc values for mzones
aggSLC_meanvalM = aggSLCmeltM %>% as.data.frame()


aggSLCmeltrateM = SLCmelt %>% group_by(mZone2) %>% summarise(mean_SLC = mean(linearRate, na.rm = TRUE))  ##mean slc values for mzones
aggSLC_meanvalrateM = aggSLCmeltrateM %>% as.data.frame()


## ## ## --------- ---------------- ---------------- ------------------ ----------- ## ## ##PLOTTING SHORELINE CHANGE RATES ##

## ## ## CREATING LABEL DATA FOR ZONES ## ## ##
zonelabs = data.frame(trans_id = zoneboundaries, y = rep(0,length(zoneboundaries)), names = zonenames)
mzonelabs = data.frame(trans_id = mzoneboundaries[2:9], y = rep(0,length(mzoneboundaries)-1), names = mzonenames)

##### ------- -------- --------- ------- ######## #
## Plotting Shorelines map view ##
x11()
PlotData = Shorelines[,]
plt = ggplot(PlotData[PlotData$Trans_ID<250 & PlotData$Trans_ID>150,], aes(x = POINT_X, y = POINT_Y))+
  geom_point(aes(color = yr), size = 1)+
  coord_fixed(ratio = 1)

## Plotting Shoreline position by transect - not map ##
x11()
PlotData = Shorelines[,]
plt = ggplot(PlotData[PlotData$Trans_ID<1380 & PlotData$Trans_ID>4,], aes(x = Trans_ID, y = tdist))+
  geom_line(aes(color = yr), size = 1)+
  scale_y_continuous(limits = c(0,200))

cc = scales::seq_gradient_pal("blue","red", "Lab")(seq(0,1,length.out = 14))
plt + scale_color_manual(values = cc)

## HEATMAP ##
x11()
PlotData = SLCmelt
PlotData$group = cut(PlotData$value, breaks = c(-900,-20,-10,0,5,10,20,900))
ggplot(PlotData, aes(Trans_ID,variable,fill=group))+
   geom_tile()+
   scale_fill_viridis(discrete = TRUE,option = "turbo")+
   ggtitle("shoreline change, negative is retreat")+
   scale_x_continuous(limits = c(0,500))


## ------ LINEAR Shoreline Change rate - ALONGSHORE by transect ------- ##
x11()
Plt = ggplot (SLCmelt, aes(x = Trans_ID, y = linearRate)) +
  geom_vline(xintercept = mzoneboundaries, color = 'gray') +
  geom_hline(yintercept =  0, linetype = "dotted") +
  geom_text(data=mzonelabs, aes( x=trans_id-35, y=y-10, label=names), angle = 90) +    ### LABELING THE ZONES! ###
  geom_point(aes(color = rSquared))

Plt + scale_y_continuous(limits = c(-15,10)) + ggtitle("linear shoreline change rate (m/yr)") + scale_x_continuous(breaks = pretty(SLCmelt$Trans_ID, n = 10))


### ------- small sample area to illustrate smoothing----------- ###
tst = SLCmelt[SLCmelt$variable == "d12_14",]
ggplot(tst, aes(x = Trans_ID, y = value))+geom_point()+geom_smooth(span = 0.005)+
  ggtitle("12 to 14 shoreline change transects")+
  labs(x = "transect number", y = "shoreline change (m)")  ## small sample area to illustrate smoothing


## Shoreline Change Rate - ALONGSHORE by transect ##
Plt = ggplot (SLCmelt, aes(x = Trans_ID, y = Rate))+
  geom_point(aes(color = variable))+
  geom_vline(xintercept = c(zoneboundaries))

Plt + scale_y_continuous(limits = c(-300,300)) + ggtitle("Shoreline Change Rate m/yr")


## Shoreline Cumulative Change over TIME ##
plt = ggplot(SLCmelt[SLCmelt$SubZone == "Canal Zone",], aes (x = date, y = cumulativeRetreat))+
  geom_point(aes(color = Trans_ID))+
  geom_line(aes(color = Trans_ID))

plt + ggtitle("cumulative shoreline change all points Canal Zone")

## Shoreline Cumulative Change over time - means for zones. 
plt = ggplot(aggSLC_meanCUmRet, aes(x = date, y = mean_SLC))+
  geom_point(aes(color = SubZone))+
  geom_line(aes(color = SubZone))
plt + ggtitle("cumulative shoreline change - mean for zones")

## Shoreline Cumulative Change ALONGSHORE by transect ##
Plt = ggplot (SLCmelt, aes(x = Trans_ID, y = cumulativeRetreat))+
  geom_point(aes(color = variable))+
  geom_vline(xintercept = c(zoneboundaries))

Plt + scale_y_continuous(limits = c(-10,800)) + ggtitle("Cumulative Shoreline Change") 

### BOX AND WHISKER FOR ZONES ###  LINEAR SHORELINE CHANGE VALUES ####
x11()
#p = ggplot(SLCmelt, aes(x = reorder(SubZone, linearRate, na.rm=TRUE), y = linearRate)) +
p = ggplot(SLCmelt, aes(x = SubZoneCode, y = linearRate)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-10,10), name = 'shoreline change (m)') +
  labs(y= "sediment flux", x = "regions")+
  # scale_x_date(date_breaks = "1 year", date_labels = "%Y")
  ggtitle("shoreline change rates (m/yr)")
p


## ## ## ----------------------------------------------------------------------- ### BOX AND WHISKER FOR ZONES ### TOTAL CUMULATIVE VOLUME CHANGE ### 
CumulativeVol500$Zone = ShorelineChange$Zone
CumulativeVol500$SubZone = ShorelineChange$SubZone
CumulativeVol500$SubZoneCode = ShorelineChange$SubZoneCode

CVmelt = melt(CumulativeVol500, id = c('Trans_ID', 'Zone', 'SubZone', 'SubZoneCode', 'mZone'))

x11()
#p = ggplot(SLCmelt, aes(x = reorder(SubZone, linearRate, na.rm=TRUE), y = linearRate)) +
p = ggplot(CumulativeVol500, aes(x = mZone, y = d19a_19b/diffdates[11])) +
  geom_boxplot() +
  scale_y_continuous(limits = c(-25,25), name = 'Cumulative Volume Change (m^3/m)') +
  labs(y= 'Cumulative Volume Change (m^3/m)', x = "morpho-zones")+
  # scale_x_date(date_breaks = "1 year", date_labels = "%Y")
  ggtitle("Cumulative Flux Per Transect 2005-2019")
p

## ## ## --------- -------------------------------- ------------------ ----------- ## ## ## --- ### CALCULATING EQUILIBRIUM FLUX ### -------- ## ## ## 
## ## ## ------------- QOW = (Dp + Hb)*Rs 
Dp = rep(1, length(CumulativeVol500$Trans_ID))  ## Shallow Back Bay Depth (WO platform ~ 1m across the island)
Hb = TransElevStats$MEAN - 0.168   ## ## ## MHW at Oregon inlet noaa tidal station. 
HbMin = TransElevStats$MEAN -0.168 - 2*TransElevStats$SD
HbMax = TransElevStats$MEAN -0.168 + 2*TransElevStats$SD    ## barrier height 
Rs = -SLCmelt$linearRate[SLCmelt$variable == 'd19a_19b']

CumulativeVol500$QOW = (Dp + Hb)*Rs
CumulativeVol500$QOW[CumulativeVol500$QOW <= 0] = 0    ### ### ### ---- setting negative Qow to zero here ------  ## ## ##



### ----------------------------------------------------------------------------### PLOT EQUILIBRIUM FLUX FROM SLC rates ###
x11()
#p = ggplot(SLCmelt, aes(x = reorder(SubZone, linearRate, na.rm=TRUE), y = linearRate)) +
p = ggplot(CumulativeVol500, aes(x = SubZoneCode, y = QOW)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.3, color = "tomato")+
  scale_y_continuous(limits = c(-25,25), name = 'Cumulative Volume Change (m^3/m)') +
  labs(y= '"Equilibrium Overwash Flux"', x = "regions")+
  # scale_x_date(date_breaks = "1 year", date_labels = "%Y")
  ggtitle("Equilibrium Overwash Flux")
p

### ----------------------------------------------------------------------------### PLOT EQUILIBRIUM FLUX AND MEASURED FLUX TOGETHER ###

CumulativeVol500$MOW = CumulativeVol500$d19a_19b/diffdates[11] ###-------------- measured overwash total flux (dividing by the time elapsed)
CumulativeVol500$OWDIFF = CumulativeVol500$MOW - CumulativeVol500$QOW ## ------ difference between equilibrium flux and measured flux (Had this mixed up before! fixed now)



## ## SUM IT ## ##
sum(CumulativeVol500$QOW, na.rm = TRUE)*20

sum(CumulativeVol500$MOW, na.rm = TRUE)*20
## ## ----- ## ##


fluxesMelt = melt(CumulativeVol500[,c(11:16)], id = c('Trans_ID', 'OWDIFF', 'mZone', 'mZone2'))

pal = wes_palette("IsleofDogs1", n = 5)

x11()
ggplot(fluxesMelt, aes(variable, value, col = variable))+ ### ------------------ qow and mow box plots and jitter plots 
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_boxplot(outlier.shape = NA)+
 # geom_jitter(alpha = 0.3)+
  facet_wrap(~mZone2, scales = 'free_y')+
 # scale_fill_manual(values= wes_palette("GrandBudapest1", n = 4))+
  scale_color_manual(values=pal[c(1,4)])+
  scale_y_continuous(limits = c(-10,15))+
  ggtitle("Lidar-Measured and Equilibrium Overwash Fluxes")+
  labs(y = "Overwash Volume Flux (m^3/m/yr)", x =" " )+
  guides(color = guide_legend(title = "Flux Dataset"))


x11()
ggplot(fluxesMelt, aes(Trans_ID, value, color = variable))+   ### --------------- plotting q ow and m ow together along transects
  geom_vline(xintercept = mzoneboundaries)+
  scale_x_continuous(breaks = pretty(SLCmelt$Trans_ID, n = 10), name = "Transect ID")+
  scale_y_continuous(limits = c(-10,30), name = 'Sediment Flux')+
  theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  geom_vline(xintercept = mzoneboundaries, color = 'darkgray')+
  scale_color_manual(values=pal[c(3,4)])+
  geom_text(data=mzonelabs, aes( x=trans_id-30, y=y-10, label=names), angle = 90, color = 'black', size = 5)+     ### LABELING THE ZONES! ###
  geom_hline(yintercept = 0, linetype = 'dotted', size = 1.1)+
  geom_point()+
  geom_smooth(method = 'loess',span = 0.05)

  
x11()
p = ggplot(fluxesMelt[fluxesMelt$variable == "MOW",], aes(Trans_ID, value, color = variable))+  ### --------------- plotting  m ow alone and attempting to smooth
  geom_vline(xintercept = mzoneboundaries)+
  scale_x_continuous(breaks = pretty(SLCmelt$Trans_ID, n = 10), name = "Transect ID")+
  scale_y_continuous(limits = c(-25,30), name = 'Sediment Flux')+
  theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  geom_vline(xintercept = mzoneboundaries, color = 'darkgray')+
  geom_text(data=mzonelabs, aes( x=trans_id-30, y=y-20, label=names), angle = 90, color = 'black', size = 5)+     ### LABELING THE ZONES! ###
  geom_hline(yintercept = 0, linetype = 'dotted', size = 1.1)+
  geom_point()+
  geom_smooth(method = 'loess',span = 0.05)+
  scale_color_manual(values=pal[c(4)])

p + geom_segment(aes(x = 1152,y = 4.24, xend = 1182, yend = 4.24), size = 1.5, color = "#F98400")+  ## horiz ##adding shara's data
  geom_segment(aes(x = 1167, y = 3.96, xend = 1167, yend = 4.52), size = 1.5, color = "#F98400")+   ## vert
  geom_segment(aes(x = 406, y = 4.29, xend = 428, yend = 4.29), size = 1.5, color ="#F98400")+    ## horiz
  geom_segment(aes(x = 417, y = 4.23, xend = 417, yend = 4.35), size = 1.5, color = "#F98400")    ## vert


pal = wes_palette(( "Zissou1" ), n = 5)
x11()
ggplot(fluxesMelt, aes(Trans_ID, OWDIFF))+   ### ------------------------------- plotting the difference between equilibrium flux and measured flux ## ## ## 
  #geom_point()+
  scale_y_continuous(limits = c(-40,40), name = 'Measured Flux - Equilibrium Flux') +
  geom_vline(xintercept = mzoneboundaries, color = 'gray')+
  geom_hline(yintercept = 0)+
  geom_text(data=mzonelabs, aes( x=trans_id-30, y=y+30, label=names), angle = 90)+     ### LABELING THE ZONES! ###
  geom_bar(aes(fill = OWDIFF < 0), stat = "identity") + scale_fill_manual(guide = FALSE, breaks = c(TRUE, FALSE), values=c(pal[5], pal[1]))+
  ggtitle("Overwash deficit")+
  scale_x_continuous(breaks = pretty(SLCmelt$Trans_ID, n = 10))

### Compare Shoreline Change to volume /flux values ###
#tst = data.frame(Trans_ID = ShorelineChange$Trans_ID, Retreatd09_11 = ShorelineChange$d09_11, 
#SubZone = ShorelineChange$SubZone, Vol09_11 = BDFTab$d09_11, Zone = ShorelineChange$Zone)

BDFTab$mZone2 = mzone2
BDFTab500$mZone2 = mzone2
BDFmelt = melt(BDFTab[,c(1:10,12)], id = c("Transect"))
tst = SLCmelt
tst$BDVol = BDFmelt$value
tst$BDvolYRS = BDFmelt$variable

#tstPLTDAT = melt(tst, id = c('Trans_ID', 'SubZone'))

p = ggplot(tst, aes(x = linearRate, y = BDVol))+
  geom_point(aes(color = mZone2))+ 
  scale_x_continuous(limits = c(-100,100))+
  scale_y_continuous(limits = c(-200,200))+
  labs(y = "volume change behind dune (m^3/m shoreline)", x = "shoreline change (m), negative is retreat")
p

### LINEAR RATE SLC VS TOTAL vol ###
tst = data.frame(Trans_ID = SLCmelt$Trans_ID, CumVol = CumulativeVol500$d19a_19b, 
                 SLCLinRate = SLCmelt$linearRate[SLCmelt$variable == "d19a_19b"], 
                 Zone = SLCmelt$Zone[SLCmelt$variable == "d19a_19b"],
                 SubZone = SLCmelt$SubZone[SLCmelt$variable == "d19a_19b"])

p = ggplot(tst, aes(x = SLCLinRate, y = CumVol))+
  geom_point(aes(color = SubZone))+ 
  #scale_x_continuous(limits = c(-10,15))+
  #scale_y_continuous(limits = c(0,300))+
  labs(y = "volume change behind dune (m^3/m shoreline)", x = "shoreline change rate (m/yr)")
p


p = ggplot(SLCmelt, aes(x = linearRate, y = MOW))+
  geom_point(aes(color = mZone))+ 
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(-10,25))+
  facet_wrap(~mZone, scales = 'free_y')+
  geom_smooth(method = 'lm')+
  labs(y = "Flux behind dune (m^3/m/yr))", x = "shoreline change rate (m/yr)")
p





## ## ## --------- -------------------------------- ------------------ ----------- ## ## ## 2. SHORELINE RETREAT RATE ZONE ESTIMATES FROM USGS -------- ## ## ## 
## ## ## ----- per transect calc from new shoreline change (08/23) above
## ## ## ------------- CALCULATE Qow 
## ## ## ------------- QOW = (Dp + Hb)*Rs 
mFluxs = aggregate(BDFPY$d05_19b, list(BDFPY$SubZone), FUN = mean)
ORD = c(8,10,6,7,9,3,5,4,1,2)   ## ORDER of zones based on the auto alphabetical used in aggregate 
mFluxs = mFluxs[ORD,]

Rs = c(0.62,   ## Rs = retreat rate 
       -1.31,
       -3.17,
       -2.28,
       -0.75,
       -3.06,
       -1.85,
       -2.58,
       -1.02,
       -6.53)  ## shoreline retreat rate
Dpdeep = c(4.6, 4.7, 4.6, 4.3, 4.1, 4.1,4.5,4.5, 4.5, 4.5)  ## back bay depth (depth pond = Dp)
Dp = rep(1, length(Rs))  ## Shallow Back Bay Depth (WO platform ~ 1m across the island)

HbMin = mElevs$x - 2*sds$x
HbMax = mElevs$x + 2*sds$x

Qow = data.frame(Zone = mElevs$Group.1)

Qow$mFlux = mFluxs$x    ## measured flux
Qow$Dp = Dp
Qow$Dpdeep = Dpdeep
Qow$HbMin = HbMin
Qow$HbMax = HbMax
Qow$HbMean = mElevs$x
Qow$Rs = Rs

Qow$Min = (Dp + HbMin)*-Rs ## Equilibrium OW flux min
Qow$Max = (Dp + HbMax)*-Rs ## EQUILIBRIUM OW FLUX MAX 
Qow$Mean = (Dp + mElevs$x)*-Rs ## EQ OW flux mean (based on min max mean island elvs)
Qow$DeepFlux = (Dpdeep + mElevs$x)*-Rs
Qow$mFlux = mFluxs$x    ## measured flux


#Qow = Qow[ORD,]
## ## ## ---------------------------------------------------- ## ## ##



## ## ## --------- -------------------------------- ------------------ ----------- ## ## ##  PLOT SCRIPTS---------------------- ## ## ##

## ## ## ----------------   cumulative fluxes --------------- ## ## ##
# CN = c("d05_09","d09_11","d11_12","d12_14","d14_16","d16_17","d17_18a","d18a_18b","d18b_19a","d19a_19b","d05_19b") 
# cn = CN[1:length(CN) - 1]
# CumulativeFlux = data.frame(matrix(ncol = length(cn),nrow = length((min(BDFTab$Transect) : max(BDFTab$Transect)))))
# colnames(CumulativeFlux) = cn
# CumulativeFlux$Trans_ID = BDFTab$Transect
# 
# for (i in 1:length(cn)){
#   CumulativeFlux[,i] = apply(BDFTab[1:i], 1, sum)
# }

### prep data
cf2 = data.frame(t(CumulativeFlux))
colnames(cf2) = BDFTab500$Transect
cf2$timeslice = sdates[2:length(sdates)]


PlotData = melt(cf2,id = c("timeslice"))
PlotData$variable = as.numeric(PlotData$variable)
PlotData$Zone = rep("NA", length(PlotData$timeslice))
PlotData$Zone[PlotData$variable <= 402] = "Developed"
PlotData$Zone[PlotData$variable > 402] = "NWR"
PlotData$SubZone = rep("NA", length(PlotData$timeslice))
PlotData$SubZone[PlotData$variable <= 1398] = "Inlet Zone"
PlotData$SubZone[PlotData$variable <= 1280] = "Canal Zone"
PlotData$SubZone[PlotData$variable <= 1147] = "North Pond"
PlotData$SubZone[PlotData$variable <= 1018] = "Old Sandbag"
PlotData$SubZone[PlotData$variable <= 858] = "New Inlet"
PlotData$SubZone[PlotData$variable <= 763] = "Stable Zone"
PlotData$SubZone[PlotData$variable <= 537] = "S-Curves"
PlotData$SubZone[PlotData$variable <= 402] = "Rodanthe"
PlotData$SubZone[PlotData$variable <= 280] = "Waves"
PlotData$SubZone[PlotData$variable <= 136] = "Salvo"


### basic plot 
#CumulativeFlux$Transect = BDFluxTab500$Transect
x11()
plot(sdates[2:length(sdates)], CumulativeFlux[1,], type = 'n', ylim = c(-400,400))
for(i in 1:length(BDFTab500$Transect)){
  points(sdates[2:length(sdates)], CumulativeFlux[i,])
}


### Violin plots
x11()
p = ggplot(PlotData, aes(x = timeslice, y = value)) + 
  geom_point() +
  #geom_smooth(method = lm, color = 'red', fill ="#69b3a2", se=TRUE) 
  geom_smooth(method = lm, color = 'red', se=TRUE, level = 0.95) +
  scale_y_continuous(limits = c(-150,300))
p


#PlotData$timeslice = as.factor(PlotData$timeslice)
x11()
p = ggplot(PlotData, aes(x = timeslice, y = value)) +
  geom_violin(aes(group = cut_width(timeslice,1)), width = 300, fill = NA, size = 1.2) +
  geom_smooth(color = 'red', se=TRUE, level = 0.95, method = 'lm') +
  scale_y_continuous(limits = c(-100,200), name = 'cumulative deposition per m shoreline') +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")
p


#jitter and box and whisker
x11()
p = ggplot(PlotData, aes(x = timeslice, y = value, color = SubZone)) +
  geom_jitter(width = 100, height = 0.5, size = 0.9) +
  geom_boxplot(aes(group = cut_width(timeslice, 1)), outlier.shape = NA) +
  scale_y_continuous(limits = c(-100,200), name = 'cumulative flux') +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")
#ggtitle("2009, 2011, 2012, 2017 removed")
p


#jitter and box and whisker - just two zones 
PDsub = PlotData[PlotData$SubZone == "Old Sandbag" | PlotData$SubZone == "Salvo",]
x11()
p = ggplot(PDsub, aes(x = timeslice, y = value, fill = SubZone, color = SubZone)) +
  geom_jitter(width = 200, height = 0.5, size = 0.9) +
  # geom_boxplot(outlier.shape = NA) +
  geom_smooth(method = lm, color = 'red', se=TRUE, level = 0.95) +
  scale_y_continuous(limits = c(-100,200), name = 'cubic meters per meter shoreline') +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")+
  ggtitle("cumulative sediment deposition")
p


## ## ##  box and whisker (sep. by zones) - FLUX not just DEPOSITION/m   ## ## ##
BDFPY500$Zone = rep("NA", length(BDFPY500$Transect))
BDFPY500$Zone[BDFPY500$Transect <= 402] = "Developed"
BDFPY500$Zone[BDFPY500$Transect > 402] = "NWR"
BDFPY500$SubZone = rep("NA", length(BDFPY500$Transect))
BDFPY500$SubZone[BDFPY500$Transect <= 1398] = "Inlet Zone"
BDFPY500$SubZone[BDFPY500$Transect <= 1280] = "Canal Zone"
BDFPY500$SubZone[BDFPY500$Transect <= 1147] = "North Pond"
BDFPY500$SubZone[BDFPY500$Transect <= 1018] = "Old Sandbag"
BDFPY500$SubZone[BDFPY500$Transect <= 858] = "New Inlet"
BDFPY500$SubZone[BDFPY500$Transect <= 763] = "Stable Zone"
BDFPY500$SubZone[BDFPY500$Transect <= 537] = "S-Curves"
BDFPY500$SubZone[BDFPY500$Transect <= 402] = "Rodanthe"
BDFPY500$SubZone[BDFPY500$Transect <= 280] = "Waves"
BDFPY500$SubZone[BDFPY500$Transect <= 136] = "Salvo"

ow = Qow[order(Qow$mFlux),]

x11()
p = ggplot(BDFPY500[BDFPY500$SubZone != "Inlet Zone",], aes(x = reorder(SubZone, d05_19b, median), y = d05_19b)) +
  #geom_jitter(width = 100, height = 0.5, size = 0.9) +
  geom_boxplot(aes(x = reorder(SubZone, d05_19b, median), y = d05_19b),outlier.shape = NA) +
  #geom_point(data = ow[ow$Zone != "Inlet Zone",], aes(x = Zone, y = DeepFlux), color = 'darkgreen', size = 2)+
  scale_y_continuous(limits = c(-5,20), name = 'sediment flux') +
  labs(y= "sediment flux", x = "regions")+
  # scale_x_date(date_breaks = "1 year", date_labels = "%Y")
  ggtitle("Washover Flux 2005 - 2019 - deep bay depth")
p


## ## ## ------------------ HEAT MAP ------------------------ ## ## ##

sub = subset(BDFTab500[,c(1:10,12)], Transect >= 1158 &  Transect <= 1176)


PlotData = melt(sub, id = c("Transect"))
#PlotData = melt(BDFluxPY, id = c("Transect"))

PlotData1 = PlotData
PlotData1$group = cut(PlotData1$value, breaks = c(150,50,20,-20,-50,-150))


x11()
ggplot(PlotData1, aes(Transect,variable, fill=group))+
  geom_tile()+
  #scale_fill_viridis(discrete=FALSE)+ 
  scale_fill_manual(breaks = levels(PlotData1$group),
                    values = c("#3794bf", "#56B4E9","#FFFFFF",
                               "#E69F00","#df8640"))+
  ggtitle("sediment flux behind the foredune in cubic meters per meter shoreline, shoreline to end of transect")+
  scale_x_continuous(name = "distance along transect",limits = c(0,max(PlotData$Transect)), breaks = seq(0,max(PlotData$Transect),100))


### Continuous
x11()
ggplot(PlotData1, aes(Transect,variable, fill=value))+
  geom_tile()+
  scale_fill_viridis(discrete=FALSE)+ 
  #scale_fill_manual(breaks = levels(PlotData1$group),
  #                  values = c("#3794bf", "#56B4E9","#FFFFFF",
  #                             "#E69F00","#df8640"))+
  ggtitle("sediment flux behind the foredune in cubic meters per meter shoreline, shoreline to end of transect, Fan 3")
#scale_x_continuous(name = "distance along transect",limits = c(0,max(PlotData$Transect)), breaks = seq(0,max(PlotData$Transect),100))


## ## ## ------- line plot: individual timestep - all transects -------- ## ## ##
PlotData = melt(BDFluxTab[,c(5,12)], id = c("Transect"))
PlotData1 = melt(BDFluxTab[,c(6,12)], id = c("Transect"))

x11()
p = ggplot()+
  geom_line(data = PlotData, aes(x = Transect, y = value))+
  geom_line(data = PlotData1, aes(x = Transect, y = value, color = 'darkorange'))+
  ggtitle(paste("sediment flux behind the foredune in cubic meters per meter shoreline",PlotData$variable[1], PlotData1$variable[1]))+
  theme(legend.position="topleft")
p
## ## ## ------------------line plot: individual transect--------------- ## ## ##


## ## ## melting ahead of time
tt = melt(Pr, id = c("OID","Trans_ID","POINT_X","POINT_Y", "tdist"))

pltr = 915 ## --- plot transect
#y1 = "05"
#y2 ="11"
#y3 = "12"

## -- 838 right in new inlet
## -- 832 between two channels - more subaerial 

## ## ## Subsetting to look at specific years ## ## #
#ttt = subset(tt, subset = Trans_ID == pltr & (variable == 'elev05'|variable == "elev11" | variable == "elev12"))
#fff = subset(Fe, subset = Trans_ID == pltr & Type == "FDC" & (Year == "2005_"|Year == "2011_" | Year =="2012_"))

## ## ## All years for transect pltr ## ## ##
ttt = tt[tt$Trans_ID == pltr,]   ### profiles
fff = Fe[Fe$Trans_ID == pltr,]   ### features

x11()
p = ggplot()+
  geom_line(data = ttt, aes(x = tdist, y = value, color = variable))+
  geom_point(data = fff[fff$Type == "FDC",], aes(x = tdist, y = Elevation)) +
  #geom_line(data = plotData[plotData$variable == 'elev18b',], aes(x = tdist, y = value, color = variable))+
  #geom_point(data = fff[fff$Type == "FDC"& fff$Year == "2018b",], aes(x = tdist, y = Elevation)) +
  # scale_color_viridis(discrete =TRUE, option = "D")+
  scale_x_continuous(name = "meters from baseline",limits = c(25,500))+
  scale_y_continuous(name = "elevation (m)", limits = c(0,7))+
  ggtitle(paste("profile", pltr ,sep = " "))
p
## ## ## ---------------------------------------------------- ## ## ##

## ## ## -- COMPARE DIFFS for error adjustments - the shrubs--## ## ##

tt1 = subset(tt, subset = Trans_ID == pltr & (variable == "elev05"))
tt2 = subset(tt, subset = Trans_ID == pltr & (variable == "elev09"))
tt3 = subset(tt, subset = Trans_ID == pltr & (variable == "elev14"))

tt12 = data.frame(value = (tt2$value - tt1$value), variable = "tt12", tdist = tt1$tdist)
tt23 = data.frame(value = (tt3$value - tt2$value), variable = "tt23", tdist = tt2$tdist)

diff = rbind(tt12, tt23)

x11()
p = ggplot()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", size = 1)+
  geom_line(data = diff, aes(x = tdist, y = value, color = variable))+
  geom_point(data = fff[fff$Type == "FDC",], aes(x = tdist, y = 0)) +
  geom_line(data = tt12, aes(x = tdist,y = tt12$value + tt23$value))+
  scale_x_continuous(name = "meters from baseline",limits = c(25,1000))+
  scale_y_continuous(name = "elevation change (m)", limits = c(-2,2))+
  ggtitle(paste("profile", pltr, "years:" , tt1$variable[1],tt2$variable[1],tt3$variable[1],sep = " "))

p


## ## ## ------------------- profile 
plotData = melt(sampleProf, id = c("OID","Trans_ID","POINT_X","POINT_Y", "tdist"))

x11()

p = ggplot() +
  geom_point(data = plotData[plotData$tdist == 0,], aes(x = POINT_X, y = POINT_Y)) +
  geom_point(data = sampleFeat[sampleFeat$Type == "FDC",], aes(x = POINT_X, y = POINT_Y, color = Year)) + 
  scale_color_viridis(discrete =TRUE, option = "D")+
  coord_fixed()
p


## ## ## ------------------ Total Fluxes
TOTsumchange = data.frame('value' = rowSums(TFTab500[,1:10]), 'variable' = 'tot500') ### summed change for transect shoreline to 500 m
sumchange = data.frame('value' = rowSums(BDFTab500[,1:10]), 'variable' = 'bd500') ### summed change for each transect and each timestep behind fd to 500 m
SLsumchange = data.frame('value' = (TOTsumchange$value - sumchange$value), 'variable' = 'slbd')

x11()
plot(BDFluxTab$Transect, sumchange$value, ylab = "total flux (2005 - 2019)", 
     main = "total flux (fd inland)", xlab = "transect", xlim = c(0,1250),ylim = c(-200,400), type = 'l', lwd = 1)
abline(h = 0, col = 'gray')
abline(v = c(136,280,402,537,763,858,1018,1147,1280), col = 'gray')
#lines(BDFluxTab$Transect, BDFluxTab$d05_19b, col = 'darkorange')
#points(BDFluxTab$Transect, BDFluxTab$d05_19b, col = 'darkorange')
#lines(BDFluxTab$Transect, sumchange)
lines(TotFluxTab500$Transect, SLsumchange$value, col = 'darkorange', lwd = 1)
legend("topleft", legend = c("dune to 500 m", "shoreline to dune"), col = c("black", "darkorange"), lwd = 1)


## ## ## -------------------HISTOGRAMS ---------------------- ## ## ##
PlotData = rbind(sumchange[TotFluxTab500$Transect<1250,], SLsumchange[TotFluxTab500$Transect<1250,])
grp.mean = data.frame(ave = c(mean(sumchange$value[TotFluxTab500$Transect<1250]), mean(SLsumchange$value[TotFluxTab500$Transect<1250])), variable = c("bd500", "slbd"))

x11()
p = ggplot(PlotData, aes(x=value, color=variable)) +
  geom_histogram(fill="white", alpha=0.5, position="identity")+
  geom_vline(data=grp.mean, aes(xintercept=ave, color=variable), linetype="dashed")+
  theme(legend.position="top")
p

## ## ## ------------------------------------------------------------------- ## ## ## Volume Distance Thing - all the points! ## ## ##
pal = wes_palette(( "Darjeeling1" ), n = 5)
pal = pal[c(4,1,5,2)] 
ECCplt = ECC[ECC$SubZone == "Rodanthe" | ECC$SubZone == "Old Sandbag",]#
ECCplt = ECC[ECC$SubZone != "Inlet Zone",]
ECCplt$Zone2 = paste(ECCplt$Zone, "2") 

x11()
ggplot(ECCplt %>% arrange(-Trans_ID))+
  geom_point(aes(x = sldist, y = cumch, col = Zone),  alpha = 0.05,  size = 1)+
  scale_color_manual(values=pal)+
  geom_smooth(aes(x = sldist, y = cumch, col = Zone2), level = 0.95, method = 'gam', formula = y~s(x, bs = "cs", k = 15) )+
  #geom_smooth(aes(x = sldist, y = cumch, col = SubZone), method = 'lm')+
  geom_hline(yintercept = 0, linetype = 'dashed')+ 
  ggtitle("cumulative elevation change at distance from shoreline")+  
  scale_x_continuous(limits = c(0,1200), breaks = pretty(ECC$sldist, n = 20))

pal = pal[c(4,5)]

x11()   ### NOT JUST TWO ZONES, all but oregon inllet
ggplot(ECC[ECC$mzb != "inoi",] %>% arrange(-Trans_ID))+
  geom_point(aes(x = sldist, y = cumch, col = Zone),  alpha = 0.5,  size = 1)+
  scale_color_manual(values=pal)+
  #geom_smooth(aes(x = sldist, y = cumch, col = Zone), level = 0.95, method = 'gam', formula = y~s(x, bs = "cs", k = 15) )+
  #geom_smooth(aes(x = sldist, y = cumch, col = SubZone), method = 'lm')+
  geom_hline(yintercept = 0, linetype = 'dashed')+ 
  ggtitle("cumulative elevation change at distance from shoreline")+  
  scale_x_continuous(limits = c(0,600), breaks = pretty(ECC$sldist, n = 20))+
  scale_y_continuous(limits = c(0,7))



## ## ## ---------------------------------------------------------------## ## ## Volume-distance box plots

x11()
p = ggplot(ECC[ECC$cumch > 0 & ECC$mzbl != "inoi",], aes(x = Zone, y = sldist, color = Zone)) +   ### m above MHW based on Duck station
  #geom_jitter(width = 100, height = 0.5, size = 0.9) +
  #geom_boxplot(aes(x = reorder(mzbl, MOW, median), y = MOW),outlier.shape = NA) +
  geom_violin(trim = FALSE)+
  geom_boxplot(outlier.shape = NA, width = 0.1) +
  #geom_point(data = ow[ow$Zone != "Inlet Zone",], aes(x = Zone, y = DeepFlux), color = 'darkgreen', size = 2)+
  #scale_y_continuous( name = "Shoreline Change Rate (m/yr)") +
  #labs(y = "Shoreline Change Rate (m/yr)", x = "regions")+
  #scale_x_date(date_breaks = "1 year", date_labels = "%Y")
  ggtitle("deposition distance")

p + scale_color_brewer(palette = "Dark2")


## ## ## Breaking distance and cumulative change into bins for better plotting  (Excluding oregon inlet (inoi))

ECC_mutate = ECC[ECC$mzbl != "inoi",] %>% mutate(DistBins = cut(sldist, breaks = seq(from = 0, to = 300, by = 10)))   
ECC_means = ECC_mutate %>% group_by(DistBins,mzbl,Zone) %>% 
  summarise(mean_cumch = mean(cumch, na.rm = TRUE), sd_cumch = sd(cumch, na.rm = TRUE))  ### all cumulative change values 
ECC_Zmeans = ECC_mutate %>% group_by(DistBins, Zone) %>% 
  summarise(mean_cumch = mean(cumch, na.rm = TRUE), sd_cumch = sd(cumch, na.rm = TRUE))  ### all cumulative change values, big ZONES
ECC_Pmeans = ECC_mutate %>% group_by(DistBins,mzbl,Zone) %>% 
  summarise(mean_cumch = mean(cumch[cumch > 0], na.rm = TRUE), sd_cumch = sd(cumch[cumch > 0], na.rm = TRUE))  ### positive cumulative change values only.
ECC_PZmeans = ECC_mutate %>% group_by(DistBins,Zone) %>% 
  summarise(mean_cumch = mean(cumch[cumch > 0], na.rm = TRUE), sd_cumch = sd(cumch[cumch > 0], na.rm = TRUE))  ### positive cumulative change values only, big zones.


#setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/FinalWriting/Figure_drafts")
#pdf("WashoverDistanceFromShorelineBar_notjustpositive_noINNI_cut300.pdf")
x11()
p = ggplot(ECC_Zmeans, aes(x=DistBins, y = (mean_cumch), fill = Zone)) +
  geom_bar(stat = "identity", position=position_dodge())+
  geom_errorbar(aes(ymin=mean_cumch - sd_cumch, ymax=mean_cumch + sd_cumch, colour = Zone), width=.2,
                position=position_dodge(.9))
p
#ggsave("WashoverDistanceFromShorelineBar_notjustpositive_noINNI_cut300.pdf", width = 5, height = 5, units = "in", dpi = 300)
#dev.off()





## ## ## ---------------------------------------------------------------## ## ## Volume-distance summed & contributing portions plots 


pal = wes_palette("FantasticFox1", 100, type = "continuous") 
pal = pal[c(1,95,10,85,45,20,55)]

# zone.sums = ECC[ECC$mzbl != "inoi",] %>% 
#   mutate(DistBins = cut(sldist, breaks = seq(from = 0, to = 300, by = 10))) %>%
#   group_by(DistBins) %>% 
#   summarise(sum_cumch = sum(cumch[cumch > 0], na.rm = TRUE))

ECC[ECC$mzbl != "inoi",] %>%            ### CUMULATIVE VOLUME DEP IN SLIVERS ###
  mutate(DistBins = cut(sldist, breaks = seq(from = 0, to = 300, by = 10))) %>%
  group_by(DistBins,mzbl.name, mzbl.order) %>% 
  ### SUM is area change in each dist bin from 1m cells (m^2 * 20m = m^3) %>%
  summarise(sum_cumch = sum(cumch[cumch > 0], na.rm = TRUE)*20)%>%  
  subset(!is.na(DistBins))%>%
  #arrange(DistBins, mzbl.order) %>%
  #ggplot(aes(fill = mzbl, x = DistBins, y = sum_cumch)) + geom_col(position = 'stack')+
  ggplot(aes(x = DistBins, y = sum_cumch, fill = mzbl.name, group = mzbl.name)) + geom_area(position = 'stack')+
  ggtitle("Cumulative elevation change with distance from shoreline")+
  ylab("cumulative elevation change m^3")+
  xlab("distance from shoreline 10 m bins")+
  scale_fill_manual(values = pal)


  
#pal = wes_palette("FantasticFox1", 100, type = "continuous") 
#pal = pal[c(1,95,10,85,45,20,55)]
#pal = pal[c(30,40,50,60,90,70,100)]

#pal = wes_palette("Zissou1", 100, type = "continuous") 
#pal = pal[c(10,20,30,40,90,50,100)]
#pal = wes_palette("Zissou1Continuous", 7, type = "continuous") 


ECC[ECC$mzbl != "inoi",] %>%             ### NORMALIZED BY ALONGSHORE LENGTH ###
  mutate(DistBins = cut(sldist, breaks = seq(from = 0, to = 300, by = 10))) %>%
  group_by(DistBins, mzbl.order,mzbl.name) %>% 
  ### SUM is area change in each dist bin from 1m cells (=m^2) -> divide by length alongshore length to calculate average alongshore 
  summarise(sum_cumch = (sum(cumch[cumch > 0], na.rm = TRUE)*20)/length(unique(Trans_ID)), alongshore.length = length(unique(Trans_ID))) %>%  
  ### sum in 10 m bin x 20 m = cumulative volume change in 10 m cross shore bins times the 20 m transect spacing in zone. divide by transects in zone to get average 
  subset(!is.na(DistBins)) %>%
  arrange(DistBins,mzbl.order) %>%
  
  ggplot(aes(x = DistBins, y = sum_cumch, group = mzbl.order, color = mzbl.order)) +
  geom_line(size = 2,alpha=0.8) +
  ggtitle("Cumulative elevation change with distance from shoreline, alongshore normalized") +
  ylab("cumulative elevation change m^3/m") +
  xlab("distance from shoreline 10 m bins")+ 
  scale_color_manual(values = pal)+
  scale_color_manual(labels = unique(ECC$mzbl.name), values = pal)



ECC[ECC$mzbl != "inoi",] %>%             ### NORMALIZED BY ALONGSHORE LENGTH ###
  mutate(DistBins = cut(sldist, breaks = seq(from = 0, to = 300, by = 10))) %>%
  group_by(DistBins, mzbl.order,mzbl.name) %>% 
  ### SUM is area change in each dist bin from 1m cells (=m^2) -> divide by length alongshore length to calculate average alongshore 
  summarise(sum_cumch = (sum(cumch[cumch > 0], na.rm = TRUE))/(length(unique(Trans_ID))*10), alongshore.length = length(unique(Trans_ID))) %>%  
  ### to get average cumulative elevation change in these bins, don't multiply by 20 (no volume, just elev change), and divide by 10 to get ave value in 10 m bin. 
  ### only positive change!
  subset(!is.na(DistBins)) %>%
  arrange(DistBins,mzbl.order) %>%
  
  ggplot(aes(x = DistBins, y = sum_cumch, group = mzbl.order, color = mzbl.order)) +
  geom_line(size = 2,alpha=0.8) +
  ggtitle("Cumulative elevation change with distance from shoreline, alongshore normalized") +
  ylab("cumulative elevation change m^3/m") +
  xlab("distance from shoreline 10 m bins")+ 
  scale_color_manual(values = pal)+
  scale_color_manual(labels = unique(ECC$mzbl.name), values = pal)



## ## ## ------------- USEFUL LINES OF CODE -------------- ## ## ##
# install.packages("")
# ord = BDFluxPY[order(BDFluxPY$d11_12),]
# head(ord, 10)      ## 10 smallest fluxes from date range
# tail(ord, 10)      ## 10 largest fluxes from date range

# abline(h = 0)

# dev.set(4)         ## activing a specific graphics device

# df = data.frame(matrix(ncol =  , nrow = )

# hist(d, breaks = 20)

# plot(x,y,asp=1)

# locator(5)  ## locate coords for 5 points on the plot

#legend(1, 95, legend=c("Line 1", "Line 2"),
#       col=c("red", "blue"), lty=1:2, cex=0.8)
#legend("topleft", legend = c("Lidar 2018a","Lidar 2018b", "Lidar 2019a","Lidar 2019b", "Lidar 2005"), col = pal[5:9], lwd = 1)

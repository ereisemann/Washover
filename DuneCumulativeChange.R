################################################################################
### Eve Eisemann -- 02/2024
### calculating cumulative dune change and comparing with shoreline change 
###
################################################################################

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
library(tidyverse)
library(wesanderson)
library(reshape)
library(viridis)
library(data.table)
library(zoo)
library(ggplot2)
library(RColorBrewer)


setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/R_analysis/VolumeDefecit") ### pc

SLCmelt = read.table("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/R_analysis/FluxCalcs/SLCmelt.csv", header = TRUE, sep = ",") ## distance between points from each timestep ## negative is retreat, positive is growth 
## wrote new SLC melt table that includes the mzbl and other zones. 
mzoneboundaries = c(0,219, 403, 492, 819, 916, 1009, 1300, 1398) 
mzonenames = c("DD","DND","NND","ND","NND","ND","NND","inlet")
mzonenames2 =  c("a_DD","b_DND","c_NND","d_ND","e_NND","f_ND","g_NND","h_inlet")
i = (1:length(mzonenames))
mzone = c(rep(mzonenames[i],(mzoneboundaries[i+1]-mzoneboundaries[i])))
mzone2 = c(rep(mzonenames2[i],(mzoneboundaries[i+1]-mzoneboundaries[i])))

mzoneLboundaries = c(0,250, 403, 492, 819, 916, 1300, 1398)    ## lumping f and g together hummocky/small/no veg dunes in NWR
mzoneLnames =  c("aa","bb","c_NND","d_ND","e_NND","f_g","h_inlet")
i = (1:length(mzoneLnames))
mzoneLUMP = c(rep(mzoneLnames[i], (mzoneLboundaries[i+1]-mzoneLboundaries[i])))

### OPTION B LUMPED ZONES
mzbl.boundaries = c(0,250, 403, 492, 819, 851, 1009, 1272, 1398)    ## m zone boundaries lumped: more morpho behavior designations. lumping f and g together hummocky/small/no veg dunes in NWR
mzbl.label = c("tnat", "tman","insc","nnat", "inni", "nmas","nman","inoi" )
i = (1:length(mzbl.label))
mzbl = c(rep(mzbl.label[i], (mzbl.boundaries[i+1]-mzbl.boundaries[i])))
### aa :     tnat: 'natural' behavior in towns- shoreline retreat, dunes shrink, shoreline prograde, dunes grow
### bb :     tman: 'managed' behavior in towns- no relationship between shoreline change and dune growth
### c_NND:   insc: inlet zone -  s-curves ephemeral inlet
### d_ND:    nnat: 'natural' behavior in NWR (stable zone, dune ridges)
### e_NND:   inni: inlet zone - new inlet
### f_g:     nman: 'managed' behavior in NWR (Ponds North) (highly managed dunes, road close to shoreline)
### f_g:     nmas: 'managed' behavior south of visitor's center ponds (Ponds South)
### h_inlet: inoi: inlet zone - Oregon Inlet
names = c("Waves & Salvo", "Rodanthe", "S-Curves Inlet", "NWR C", "New Inlet", "NWR B", "NWR A", "Oregon Inlet")
mzbl.name = c(rep(names[i], (mzbl.boundaries[i+1]-mzbl.boundaries[i])))

mzbl.ord = c(1:8)
i = (1:length(mzbl.ord))
mzbl.order = c(rep(mzbl.ord[i], (mzbl.boundaries[i+1]-mzbl.boundaries[i])))

# SLCmelt$mzone = mzone    ## adding all the above zones to slcmelt
# SLCmelt$mzone2 = mzone2
# SLCmelt$mzoneLUMP = mzoneLUMP
# SLCmelt$mzbl = mzbl
SLCmelt$mzbl.order = mzbl.order


Fe = read.table("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/R_analysis/FluxCalcs/Features.txt", header = TRUE, sep = ',')            ## dune features 
Feclip = Fe[Fe$Type == "HD",c(2,7,11)]    ## slimming down to just HD elevation, transect ID, and year
FeCast = dcast(data = Feclip, formula = Trans_ID~Year, fun.aggregate = sum, value.var = 'Elevation')

HDave = aggregate(Fe$Elevation[Fe$Type == "HD"] ~ Fe$Trans_ID[Fe$Type == "HD"], data = Fe, mean)   ### mean highest dune values  - STILL ELEVATION 


QOW_MOW = read.table("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/R_analysis/FluxCalcs/QOW_MOW.txt", header = TRUE, sep = ",")              ## EQUILIBRIUM OVERWASH FLUX (QOW) based on Carruthers et al., and Measured Overwash Flux from lidar differencing
QOW_MOW$mzoneLUMP = mzoneLUMP
QOW_MOW$mzbl = mzbl
QOW_MOW$mzbl.order = mzbl.order
QOW_MOW$mzbl.name = mzbl.name




################################################################################
###                                                                          ###
###        calculating total change over the time period in HD               ###
###                                                                          ###
################################################################################

yrs = c("d2005_","d2009_","d2011_","d2012_","d2014_","d2016_","d2017_","d2018a","d2018b","d2019a","d2019b")

colnames(FeCast) = c("Trans_ID", yrs )

fe_cumul_ch = FeCast$d2019b - FeCast$d2005_   ## total change over the time period - from 2005 to 2019

### a few transects are skipped at the end in the features - 1384 is the last one that matches up nicely
### only using 1:1384 for the plotting below. 



################################################################################
###                                                                          ###
###        adding linear dune elevation change rate                          ###
###                                                                          ###
################################################################################

sdates = as.Date(c('24-8-2005','10-08-2009','28-8-2011','5-11-2012','6-01-2014','21-11-2016','24-6-2017',
                   '24-8-2018','2-10-2018','18-6-2019','26-9-2019'), '%d-%m-%Y')
i = (1:length(sdates))
dd = difftime(sdates[i+1],sdates[i],units='days')
diffdates = as.numeric(c(dd[1:max(i)-1], difftime(sdates[max(i)], sdates[1], units = 'days'))/365)

# for(i in c(1:10)){   ## adding decimal years from first survey to date
#   SLCmelt$date[SLCmelt$variable == Labs[i]] = cumsum(c(diffdates))[i]
#   SLCmelt$Rate[SLCmelt$variable == Labs[i]] = SLCmelt$value[SLCmelt$variable == Labs[i]]/diffdates[i]
# }
CN = c("d05_09","d09_11","d11_12","d12_14","d14_16","d16_17","d17_18a","d18a_18b","d18b_19a","d19a_19b") 

FeFd = Fe[Fe$Type == "FD",]
FeChange = data.frame(matrix(ncol = length(CN)+1, nrow = max(FeFd$Trans_ID)))      ### Feature Change (FeChange) calculated in the same way as ShorelineChange
colnames(FeChange) = c("Trans_ID", CN[1:length(CN)])
FeChange$Trans_ID = 1: max(Fe$Trans_ID)

for (i in 1:(length(CN))){
  YRS = unique(FeFd$Year)
  yr1 = FeFd %>% filter(Year == YRS[i]) 
  yr2 = FeFd %>% filter(Year == YRS[i+1]) 
  
  for(ii in 1:(length(yr1$Trans_ID))){
    FeChange[ii,i+1] = yr2$Elevation[ii] - yr1$Elevation[ii]    ## negative is shorter, positive is higher
  }
}


#define function to extract overall p-value of model
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

FeCmelt = melt(FeChange, id = c('Trans_ID'))  ## melting feature change data frame 

for(i in c(1:10)){   ## adding decimal years from first survey to date
  FeCmelt$date[FeCmelt$variable == CN[i]] = cumsum(c(diffdates))[i]
 # SLCmelt$Rate[SLCmelt$variable == Labs[i]] = SLCmelt$value[SLCmelt$variable == Labs[i]]/diffdates[i]
}

for(i in 1:length(FeCmelt$Trans_ID)){   ## and linear regression dune change
  FeCmelt$CumuDuneCh[FeCmelt$Trans_ID == i] = cumsum(FeCmelt$value[FeCmelt$Trans_ID == i])
  lengthna = sum(is.na(FeCmelt$value[FeCmelt$Trans_ID == i]))  ## if there are fewer than three dune change values na  (value = shoreline positions in SLCmelt)
  if((length(FeCmelt$value[FeCmelt$Trans_ID == i])-lengthna) < 3){
    FeCmelt$DlinearRate[FeCmelt$Trans_ID == i] = NA
    FeCmelt$DrSquared[FeCmelt$Trans_ID == i] = NA
  }else{  ## if more than three slc values, create linear regression and save slope and r sq
    mod  = lm(FeCmelt$CumuDuneCh[FeCmelt$Trans_ID == i] ~ 0 + FeCmelt$date[FeCmelt$Trans_ID == i])    ### Add 0 b/c values are all position based on initial 0  
    coeff = coef(mod)
    FeCmelt$DlinearRate[FeCmelt$Trans_ID == i] = coeff[1]   ## meters per year
    FeCmelt$DrSquared[FeCmelt$Trans_ID == i] = summary(mod)$adj.r.squared
    FeCmelt$Pval[FeCmelt$Trans_ID == i] = overall_p(mod)

  }
}

### SAVING - foredune height change ###
setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/R_analysis/FluxCalcs") ### pc
write.table(FeCmelt, file = "FDchangerate.txt", sep = ",", col.names = TRUE)                   ## all profile points and elevs, includes zones and sub zones now (9/23)


FeCmelt = read.table("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/R_analysis/FluxCalcs/FDchangerate.txt", header = TRUE, sep = ",") ## distance between points from each timestep ## negative is retreat, positive is growth 

#########################################################################################################################
###                                                                                                                   ###
###        setting up dataframe including SLC, Dune heights, dune change rate, and cumulative change, aves, MOW.      ###
###                                                                                                                   ###
#########################################################################################################################

slc = SLCmelt[SLCmelt$date > 14,] ### shoreline change rate                   ### only cumulative retreat value for the last time step - total 
fec = FeCmelt[FeCmelt$date > 14,] ### dune ht change rate
df = data.frame(Transect = slc$Trans_ID[slc$Trans_ID < 1385], fe_cumul_ch = fec$CumuDuneCh[fec$Trans_ID <1385], slc = slc$cumulativeRetreat[slc$Trans_ID < 1385], Zone = slc$Zone[slc$Trans_ID < 1385], 
                slcR = slc$linearRate[slc$Trans_ID < 1385], slcR2 =  slc$rSquared[slc$Trans_ID < 1385] , mzone2 = slc$mzone2[slc$Trans_ID < 1385], mzoneLUMP = slc$mzoneLUMP[slc$Trans_ID < 1385], mzbl = slc$mzbl[slc$Trans_ID<1385],
                fecR = fec$DlinearRate[slc$Trans_ID < 1385], fecR2 = fec$DrSquared[slc$Trans_ID < 1385], mzbl.order = slc$mzbl.order[slc$Trans_ID<1385])
df$slc[df$slc < -200] = NaN                         ### erroneous shoreline change points caused by edge effects on the southern edge of the study area
df$HD_mean = HDave$`Fe$Elevation[Fe$Type == "HD"]`[HDave$`Fe$Trans_ID[Fe$Type == "HD"]`<1385]   ### adding mean of HD elevations 
df$MOW = QOW_MOW$MOW[QOW_MOW$Trans_ID < 1385]
df$mzbl.name = QOW_MOW$mzbl.name[QOW_MOW$Trans_ID < 1385]

################################################################################
###                                                                          ###
###        PLOTTING shoreline change and dune hts.                           ###
###                                                                          ###
################################################################################
x11() 
ggplot(df[df$mzone2 == "g_NND" | df$mzone2 == "d_ND"|df$mzone2 == "f_ND",], aes(x = slc, y = fe_cumul_ch, col = mzoneLUMP))+     ###--------------- Plotting dune elev change vs. shoreline change:NWR select zones
  geom_point()+
  #geom_vline(xintercept = mzoneboundaries)+
  scale_x_continuous(breaks = pretty(df$slc, n = 10), name = "Net Shoreline Change (m)", limits = c(-80,20))+
  scale_y_continuous(name = "Dune Elevation Change (m)", limits = c(-7,7)) + 
  theme(text=element_text(size=15), axis.text.x=element_text(size=15))+
  theme(text=element_text(size=15), axis.text.y=element_text(size=15))+
  geom_smooth(method = lm, se=TRUE, data = df[df$mzoneLUMP == "f_g",])+
  geom_smooth(method = lm, se=TRUE, data = df[df$mzone2 == "d_ND",])+
  geom_smooth(method = lm, formula = y~poly(x,3)  , se=TRUE, data = df[df$mzone2 == "d_ND",])
  #scale_color_manual(values=pal[c(3,4)])+
  #geom_text(data=mzonelabs, aes( x=trans_id-30, y=y-10, label=names), angle = 90, color = 'black', size = 5)+  ### LABELING THE ZONES! ###
  #geom_hline(yintercept = 0, linetype = 'dotted', size = 1.1)

### Testing significance of correlation
lmdat = df[df$mzoneLUMP == "f_g",]
linear_mod = lm(fe_cumul_ch ~ slc, data = lmdat)
summary(linear_mod)

lmdat = df[df$mzone2 == "d_ND",]
linear_mod = lm(fe_cumul_ch ~ slc, data = lmdat)
summary(linear_mod)

### -------------------------------------------------------------------------###
x11() 
ggplot(df[df$mzone2 == "g_NND" | df$mzone2 == "d_ND"|df$mzone2 == "f_ND",], aes(x = slcR, y = fe_cumul_ch, col = mzoneLUMP))+     ###------------ Plotting dune elev change vs. shoreline LINEAR change: NWR select zones
  geom_point()+
  #geom_vline(xintercept = mzoneboundaries)+
  scale_x_continuous(breaks = pretty(df$slcR, n = 10), name = "Linear Shoreline Change Rate (m/yr)", limits = c(-10,10))+
  scale_y_continuous(name = "Dune Elevation Change (m)", limits = c(-7,7)) + 
  theme(text=element_text(size=15),axis.text.x=element_text(size=15))+
  theme(text=element_text(size=15),axis.text.y=element_text(size=15))+
  geom_smooth(method = lm, se=TRUE, data = df[df$mzoneLUMP == "f_g",])+
  geom_smooth(method = lm, se=TRUE, data = df[df$mzone2 == "d_ND",])+
  geom_smooth(method = lm, formula = y~poly(x,2)  , se=TRUE, data = df[df$mzone2 == "d_ND",])
#scale_color_manual(values=pal[c(3,4)])+
#geom_text(data=mzonelabs, aes( x=trans_id-30, y=y-10, label=names), angle = 90, color = 'black', size = 5)+  ### LABELING THE ZONES! ###
#geom_hline(yintercept = 0, linetype = 'dotted', size = 1.1)

### Testing significance of correlation
lmdat = df[df$mzoneLUMP == "f_g",]
linear_mod = lm(fe_cumul_ch ~ slcR, data = lmdat)
summary(linear_mod)


lmdat = df[df$mzone2 == "d_ND",]
linear_mod = lm(fe_cumul_ch ~ slcR, data = lmdat)
summary(linear_mod)


### -------------------------------------------------------------------------###
x11() 
ggplot(df[(df$mzoneLUMP == "aa"|df$mzoneLUMP == "bb"),], aes(x = slc, y = fe_cumul_ch, col = mzoneLUMP))+     ###----------------------------- Plotting dune elev change vs. shoreline change: TOWNS
  geom_point()+
  #geom_vline(xintercept = mzoneboundaries)+
  scale_x_continuous(breaks = pretty(df$slc, n = 10), name = "Net Shoreline Change (m)", limits = c(-80,60))+
  scale_y_continuous(name = "Dune Elevation Change (m)", limits = c(-7,7)) + 
  theme(text=element_text(size=15), axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  geom_smooth(method = lm, se=TRUE)
#scale_color_manual(values=pal[c(3,4)])+
#geom_text(data=mzonelabs, aes( x=trans_id-30, y=y-10, label=names), angle = 90, color = 'black', size = 5)+  ### LABELING THE ZONES! ###
#geom_hline(yintercept = 0, linetype = 'dotted', size = 1.1)

### LUMP ZONE AA - waves and salvo - 'natural behavior'
lmdat = df[df$mzoneLUMP == "aa",]
linear_mod = lm(fe_cumul_ch ~ slc, data = lmdat)
summary(linear_mod)

### LUMP ZONE BB - rodanthe - no correlation between dune height and sl retreat
lmdat = df[df$mzoneLUMP == "bb",]
linear_mod = lm(fe_cumul_ch ~ slc, data = lmdat)
summary(linear_mod)

### -------------------------------------------------------------------------###
p = wes_palette(( "Darjeeling1" ), n = 5)
pal = c(p[2], p[5], p[4])
# pal = wes_palette("GrandBudapest1", n = 3)

setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/FinalWriting/Figure_drafts")
pdf("fecrVSslcr_NWR.pdf")
#x11() 
ggplot(df[(df$mzbl == "nnat"|df$mzbl == "nman"|df$mzbl == "nmas"),]%>%arrange(mzbl), aes(x = slcR, y = fecR))+     ###----------------------------- Plotting dune elev change rate vs. LINEAR shoreline change NWR
  geom_point(aes(x = slcR, y = fecR, col = mzbl, alpha = slcR2), size = 2)+
  scale_color_manual(values= pal)+
  scale_x_continuous(breaks = pretty(df$slcR, n = 10), name = "Linear Shoreline Change Rate (m/yr)", limits = c(-6,5))+
  scale_y_continuous(name = "Linear Dune Height Change (m/yr)", limits = c(-0.5,0.5)) + 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  geom_smooth(method = lm, se=TRUE, data = df[df$mzbl == "nman"|df$mzbl == "nmas",], col = pal[1])+
  geom_smooth(method = lm, se = TRUE, data = df[df$mzbl == "nnat",], col = pal[3])+
  theme(aspect.ratio=1)

ggsave("fecrVSslcr_NWR_II.pdf", width = 5, height = 5, units = "in", dpi = 300)
dev.off()

### examining linear fits 
lmdat = df[df$mzbl == "nnat",]
linear_mod = lm(fecR ~ slcR, data = lmdat)
summary(linear_mod)

### 
lmdat = df[df$mzbl == "nman" | df$mzbl == "nmas",]
linear_mod = lm(fecR ~ slcR, data = lmdat)
summary(linear_mod)

### -------------------------------------------------------------------------###

p = wes_palette(( "Darjeeling1" ), n = 5)
pal = c(p[2], p[4])

setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/FinalWriting/Figure_drafts")
pdf("fecrVSslcr_TOWNS.pdf")
#x11() 
ggplot(df[(df$mzbl == "tnat"|df$mzbl == "tman"),]%>%arrange(mzbl), aes(x = slcR, y = fecR))+     ###----------------------------- Plotting dune elev change rate vs. LINEAR shoreline change TOWNS
  geom_point(aes(x = slcR, y = fecR, col = mzbl, alpha = slcR2), size = 2)+
  scale_color_manual(values= pal)+
  scale_x_continuous(breaks = pretty(df$slcR, n = 10), name = "Linear Shoreline Change Rate (m/yr)", limits = c(-6,5))+
  scale_y_continuous(name = "Linear Dune Height Change (m/yr)", limits = c(-0.5,0.5)) + 
  theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  geom_smooth(method = lm, se=TRUE, data = df[df$mzbl == "tman",], col = pal[1])+
  geom_smooth(method = lm, se = TRUE, data = df[df$mzbl == "tnat",], col = pal[2])+
  theme(aspect.ratio=1)

ggsave("fecrVSslcr_TOWNS_II.pdf", width = 5, height = 5, units = "in", dpi = 300)
dev.off()

### examining linear fits 
lmdat = df[df$mzbl == "tnat",]
linear_mod = lm(fecR ~ slcR, data = lmdat)
summary(linear_mod)

### 
lmdat = df[df$mzbl == "tman",]
linear_mod = lm(fecR ~ slcR, data = lmdat)
summary(linear_mod)


### -------------------------------------------------------------------------###
x11()
p =ggplot(df[df$mzbl != "inoi",], aes(x = slcR, y = fecR))+     ###----------------------------- Plotting dune height change rate vs. LINEAR shoreline change - all zones TILED
  geom_vline(xintercept  = 0)+
  geom_hline(yintercept  = 0)+
  geom_point(aes(x = slcR, y = fecR, col = mzbl.name), size = 1, show.legend = FALSE)+
  facet_wrap(~reorder(mzbl.name, mzbl.order), ncol = 7)+
 # scale_color_manual(values= pal)+
  scale_x_continuous(breaks = pretty(df$slcR, n = 10), name = "Linear Shoreline Change Rate (m/yr)", limits = c(-6,5))+
  scale_y_continuous(name = "Linear Dune Height Change (m/yr)", limits = c(-0.5,0.5)) + 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  geom_smooth(method = lm, col = 'black')+
  theme(aspect.ratio=1)

p + scale_color_brewer(palette = "Dark2")

### -------------------------------------------------------------------------###


pp = wes_palette(( "Darjeeling1" ), n = 5)
pal = c(pp[2], pp[4])

x11()
p =ggplot(df[df$mzbl != "inoi",], aes(x = slcR, y = fecR))+     ###----------------------------- Plotting dune elev change rate vs. LINEAR shoreline change - all zones
  geom_vline(xintercept  = 0)+
  geom_hline(yintercept  = 0)+
  geom_point(aes(x = slcR, y = fecR, shape = Zone, col = Zone), size = 1.5)+
  #facet_wrap(~mzbl.order, ncol = 7)+
  # scale_color_manual(values= pal)+
  scale_x_continuous(breaks = pretty(df$slcR, n = 10), name = "Linear Shoreline Change Rate (m/yr)", limits = c(-6,5))+
  scale_y_continuous(name = "Linear Dune Height Change (m/yr)", limits = c(-0.5,0.5)) + 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  geom_smooth(method = lm, se=TRUE, data = df[df$Zone == "Developed",], col = "darkgreen")+
  geom_smooth(method = lm, se =TRUE, data = df[df$Zone == "NWR",], col = "darkorange3")+
  theme(aspect.ratio=1)

p + scale_color_manual(values = pal)

################################################################################
###                                                                          ###
###        PLOTTING net dune elevation change and net washover flux          ###
###                                                                          ###
################################################################################


### -------------------------------------------------------------------------###
x11()
p =ggplot(df[df$mzbl != "inoi",], aes(x = HD_mean, y = MOW))+     ###----------------------------- Plotting MOW vs. ave dune height - all zones TILED FACET
  geom_vline(xintercept  = 0)+
  geom_hline(yintercept  = 0)+
  geom_point(aes(x = HD_mean, y = MOW, col = mzbl), size = 1,  show.legend = FALSE)+
  facet_wrap(~reorder(mzbl.name, mzbl.order), ncol = 7)+
  # scale_color_manual(values= pal)+
  scale_x_continuous(breaks = pretty(df$HD_mean, n = 10), name = "Dune Elevation Average (m)", limits = c(3,12)) +
  scale_y_continuous(name = "Overwash Flux (m^3/m/yr)", limits = c(-5,25)) + 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  geom_smooth(method = lm, col = 'black')+
  theme(aspect.ratio=1)


p + scale_color_brewer(palette = "Dark2")

### -------------------------------------------------------------------------###

### -------------------------------------------------------------------------###
x11()
p =ggplot(df[df$mzbl != "inoi",], aes(x = HD_mean, y = MOW))+     ###----------------------------- Plotting MOW vs. ave dune height - all zones 
  geom_vline(xintercept  = 0)+
  geom_hline(yintercept  = 0)+
  geom_point(aes(x = HD_mean, y = MOW, col = Zone), size = 1)+
  #facet_wrap(~reorder(mzbl.name, mzbl.order), ncol = 7)+
  # scale_color_manual(values= pal)+
  scale_x_continuous(breaks = pretty(df$HD_mean, n = 10), name = "Dune Elevation Average (m)", limits = c(3,12)) +
  scale_y_continuous(name = "Overwash Flux (m^3/m/yr)", limits = c(-5,25)) + 
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  geom_smooth(method = lm, se=TRUE, data = df[df$Zone == "Developed",], col = "darkgreen")+
  geom_smooth(method = lm, se =TRUE, data = df[df$Zone == "NWR",], col = "darkorange3")+
  theme(aspect.ratio=1)


p + scale_color_manual(values = pal)

### -------------------------------------------------------------------------###


setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/FinalWriting/Figure_drafts")
pdf("mowVSfecr_TOWNS.pdf")

p = wes_palette(( "Darjeeling1" ), n = 5)
pal = c(p[2], p[4])
#x11() 
ggplot(df[((df$mzbl == "tnat" | df$mzbl == "tman")),], aes(x = HD_mean, y = MOW, col = mzbl)) +  ### --------------- Plotting MOW vs. mean dune elev TOWNS
  geom_point() +
  scale_color_manual(values= pal)+
  scale_x_continuous(breaks = pretty(df$HD_mean, n = 10), name = "Dune Elevation Average (m)", limits = c(3,12)) +
  scale_y_continuous(name = "Overwash Flux (m^3/m/yr)", limits = c(-5,25)) + 
  theme(text=element_text(size=15), axis.text.x=element_text(size=15)) +
  theme(text=element_text(size=15), axis.text.y=element_text(size=15)) +
  #geom_smooth(aes(x = HD_mean, y = MOW), method = lm, formula = y~(log(x)), se=TRUE, inherit.aes = FALSE, col = 'darkgrey')+
  geom_smooth(method = "nls", aes(x = HD_mean, y = MOW),inherit.aes = FALSE,
              formula = y ~ a * exp(b * x), 
              se =  FALSE, # this is important 
              method.args = list(start = list(a = .10, b = -0.5)), 
              col = "black")
  #geom_line(data = lmdat, aes(x = HD_mean, y = fit), color = 'blue')    ### FROM THE FIT ASSESSMENT BELOW - turns out identical 

ggsave("mowVSdunemean_TOWNS.pdf", width = 5, height = 5, units = "in", dpi = 300)
dev.off()

### MOW vs mean dune elevation (Towns combined, linear method = 'nls')
lmdat = df[((df$mzbl == "tnat" | df$mzbl == "tman")),]
mod = nls(formula = MOW ~ a * exp(b * HD_mean), data = lmdat, start = list(a = .10, b = -0.5))
#  lm(fe_cumul_ch ~ slc, data = lmdat)
summary(mod)
lmdat$fit = predict(mod, newdata = lmdat)

### -------------------------------------------------------------------------###
setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/FinalWriting/Figure_drafts")
pdf("mowVSfecr_TOWNS_plus_nman.pdf")

p = wes_palette(( "Darjeeling1" ), n = 5)
pal = c(p[2], p[4], p[5])
#x11() 
ggplot(df[(df$mzbl == "nnat"|df$mzbl == "nmas"|df$mzbl == "nman"),], aes(x = HD_mean, y = MOW, col = mzbl )) + ### ---------------------------------- Plotting MOW vs. ave dune height (NWR naturalish areas)
  geom_point() +
  scale_color_manual(values= pal)+
  scale_x_continuous(breaks = pretty(df$HD_mean, n = 10), name = "Dune Elevation Average (m)", limits = c(3, 12))+
  scale_y_continuous(name = "Overwash Flux (m^3/m/yr)", limits = c(-5,25)) + 
  theme(text=element_text(size=15), axis.text.x=element_text(size=15)) +
  theme(text=element_text(size=15), axis.text.y=element_text(size=15)) +
  #geom_smooth(aes(x = HD_mean, y = MOW), method = lm, formula = y~(log(x)), se=TRUE, inherit.aes = FALSE, col = 'black')+
  geom_smooth(method = "nls",  data = df[(df$mzbl == "nnat"|df$mzbl == "nmas"),],aes(x = HD_mean, y = MOW),inherit.aes = FALSE,
              formula = y ~ a * exp(b * x), 
              se =  FALSE, # this is important 
              method.args = list(start = list(a = .1, b = -0.2)), 
              col = "black")
ggsave("mowVSdunemean_NWR.pdf", width = 5, height = 5, units = "in", dpi = 300)
dev.off()


### MOW vs mean dune elevation (NWR zones B and C , linear method = 'nls')
lmdat = df[(df$mzbl == "nnat"|df$mzbl == "nmas"),]
mod = nls(formula = MOW ~ a * exp(b * HD_mean), data = lmdat, start = list(a = .10, b = -0.2))
#  lm(fe_cumul_ch ~ slc, data = lmdat)
summary(mod)
lmdat$fit = predict(mod, newdata = lmdat)

### -------------------------------------------------------------------------###

setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/FinalWriting/Figure_drafts")
pdf("mowVSfecr_TOWNS_NMAN.pdf")
#x11() 
ggplot(df[(df$mzbl == "nman"),], aes(x = HD_mean, y = MOW, col = mzbl )) + ### ------------------------------------------------ Plotting MOW vs. ave dune height 
  geom_point() +
  scale_x_continuous(breaks = pretty(df$HD_mean, n = 10), name = "Dune Elevation Average (m)", limits = c(3, 12))+
  scale_y_continuous(name = "Overwash Flux (m^3/m/yr)", limits = c(-5,25)) + 
  theme(text=element_text(size=15), axis.text.x=element_text(size=15)) +
  theme(text=element_text(size=15), axis.text.y=element_text(size=15)) +
  geom_smooth(aes(x = HD_mean, y = MOW), method = lm, formula = y~x, se=TRUE, inherit.aes = FALSE, col = 'black')
  # geom_smooth(method = "nls",
  #             formula = y ~ a * exp(b * x), 
  #             se =  FALSE, # this is important 
  #             method.args = list(start = list(a = .1, b = -0.2)), 
  #             col = "black")
ggsave("mowVSdunemean_NWR.pdf", width = 5, height = 5, units = "in", dpi = 300)
dev.off()


### -------------------------------------------------------------------------###
x11() 
ggplot(df[(df$mzbl == "nman"|df$mzbl == "nmas"),], aes(x = HD_mean, y = MOW, col = mzbl)) + ### --------------- Plotting MOW vs. dune height linear change (NWR managed)
  geom_point() +
  scale_x_continuous(breaks = pretty(df$HD_mean, n = 10), name = "Dune Elevation Average (m)") +
  scale_y_continuous(name = "Overwash Flux (m^3/m/yr)") + 
  theme(text=element_text(size=15), axis.text.x=element_text(size=15)) +
  theme(text=element_text(size=15), axis.text.y=element_text(size=15)) +
  geom_smooth(aes(x = HD_mean, y = MOW), method = lm, formula = y~(log(x)), se=TRUE, inherit.aes = FALSE, col = 'black')


### -------------------------------------------------------------------------###
x11() 
ggplot(df[(df$mzoneLUMP == "aa" | df$mzoneLUMP == "bb"),], aes(x = HD_mean, y = MOW, col = mzoneLUMP))+     ### --------------- Plotting MOW vs. AVE DUNE HEIGHT TOWNS
  geom_point()+
  scale_x_continuous(breaks = pretty(df$HD_mean, n = 10), name = "Dune Elevation Average (m)")+
  scale_y_continuous(name = "Overwash Flux (m^3/m/yr)") + 
  theme(text=element_text(size=15), axis.text.x=element_text(size=15))+
  theme(text=element_text(size=15), axis.text.y=element_text(size=15))+
  # geom_smooth(method = lm, se=TRUE, data = df[df$mzoneLUMP == "f_g",])+
  # geom_smooth(method = lm, se=TRUE, data = df[df$mzone2 == "d_ND",])+
 # geom_smooth(aes(x = HD_mean, y = MOW), method = lm, formula = y~poly(x,2), se=TRUE, inherit.aes = FALSE, col = 'black')
  geom_smooth(aes(x = HD_mean, y = MOW), method = lm, formula = y~(log(x)), se=TRUE, inherit.aes = FALSE, col = 'black')


### -------------------------------------------------------------------------###
x11() 
ggplot(df[(df$mzoneLUMP == "d_ND"),], aes(x = HD_mean, y = MOW, col = mzoneLUMP))+     ### --------------- Plotting MOW vs. AVE DUNE HEIGHT D_ND 
  geom_point()+
  scale_x_continuous(breaks = pretty(df$HD_mean, n = 10), name = "Dune Elevation Average (m)")+
  scale_y_continuous(name = "Overwash Flux (m^3/m/yr)") + 
  theme(text=element_text(size=15), axis.text.x=element_text(size=15))+
  theme(text=element_text(size=15), axis.text.y=element_text(size=15))+
  # geom_smooth(method = lm, se=TRUE, data = df[df$mzoneLUMP == "f_g",])+
  # geom_smooth(method = lm, se=TRUE, data = df[df$mzone2 == "d_ND",])+
  # geom_smooth(aes(x = HD_mean, y = MOW), method = 'gam', formula = y~s(log(x)), se=TRUE, inherit.aes = FALSE, col = 'black')   ### Look up what the gam method 's' function does 
  geom_smooth(aes(x = HD_mean, y = MOW), method = lm, formula = y~(log(x)), se=TRUE, inherit.aes = FALSE, col = 'black')

################################################################################
###                                                                          ###
###        PLOTTING box plots for results                                    ###
###                                                                          ###
################################################################################

### ----------------- ADDING BOX/violin PLOTS OF OW for ZONES --------------------###
setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/FinalWriting/Figure_drafts")
pdf("mow_jitterplots_geographic.pdf")
#x11()

p = ggplot(QOW_MOW[QOW_MOW$mzbl != "inoi",], aes(x = mzbl, y = MOW, color = mzbl)) +
  #geom_violin(aes(x = reorder(mzbl, mzbl.order), y = MOW),trim = FALSE)+
  geom_boxplot(aes(x = reorder(mzbl, mzbl.order), y = MOW),outlier.shape = NA, width = 0.2) +
  geom_jitter(alpha = 0.3, size = 0.9)+
  scale_y_continuous(limits = c(-10,30), name = 'sediment flux') +
  labs(y = "washover flux", x = "regions")+
  ggtitle("Washover Flux 2005 - 2019")

p + scale_color_brewer(palette = "Dark2")

ggsave("mow_jitterplots_geographic.pdf", width = 5, height = 5, units = "in", dpi = 300)
dev.off()

aggMOWmean = QOW_MOW %>% group_by(mzbl) %>% summarise(mean_MOW = mean(MOW, na.rm = TRUE))  ##mean slc values for zones
#aggMOW_meanval = aggSLCmelt %>% as.data.frame()
aggMOWsd = QOW_MOW %>% group_by(mzbl) %>% summarise(sd_MOW = sd(MOW, na.rm = TRUE))  ##mean slc values for zones


### -------------------------Welch's T Test  --------------###
for(mm in mzbl.label){
  for (m in mzbl.label[mzbl.label != mm]){
    print(c(mm,m))
    print(t.test(MOW ~ mzbl, data = (QOW_MOW[QOW_MOW$mzbl == mm | QOW_MOW$mzbl == m,])))
  }
} 

### -----------------WILCOX RANK SUM (MANN WHITNEY) non parametric  --------------###
for(mm in mzbl.label){
  for (m in mzbl.label[mzbl.label != mm]){
    wtst = wilcox.test(MOW ~ mzbl, data = (QOW_MOW[QOW_MOW$mzbl == mm | QOW_MOW$mzbl == m,]))   ### MEASURED FLUX
    if(wtst$p.value >= 0.05){     ## ONLY PRINTING ONES NOT SIG DIFF. 
      print(c(mm,m))
      print(wtst)
    }
  }
} 

### -----------------WILCOX RANK SUM (MANN WHITNEY) non parametric - TOWN VS NWR --------------###
wtst = wilcox.test(MOW ~ Zone, data = (QOW_MOW[QOW_MOW$mzbl != 'inoi',]))   ### MEASURED FLUX
print(wtst)

# data = (QOW_MOW[QOW_MOW$mzbl != 'inoi',])
# d = (mean(data$MOW[data$Zone=='Developed'], na.rm = TRUE)-mean(data$MOW[data$Zone=='NWR'], na.rm = TRUE))/
#   sqrt((var(data$MOW[data$Zone=='Developed'], na.rm = TRUE)+var(data$MOW[data$Zone=='NWR'], na.rm = TRUE))/2)
# d
# ## cohens d - not useful for non parametric tests!! 


### ----------------- ADDING BOX/violin PLOTS OF Dune Heights ------------------###
setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/FinalWriting/Figure_drafts")
pdf("duneHts_jitterplots_geographic.pdf")
#x11()

p = ggplot(df[df$mzbl != "inoi",], aes(x = mzbl, y = (HD_mean - 0.360), color = mzbl)) +   ### m above MHW based on Duck station
  #geom_violin(aes(x = reorder(mzbl, mzbl.order), y = (HD_mean - 0.360)),trim = FALSE)+
  geom_boxplot(aes(x = reorder(mzbl, mzbl.order), y = (HD_mean - 0.360)),outlier.shape = NA, width = 0.2) +
  geom_jitter(alpha = 0.3, size = 0.9)+
  scale_y_continuous(limits = c(0,15), name = "mean foredune elevation (m MHW)") +
  labs(y = "mean foredune elevation (m MHW)", x = "regions")+
  ggtitle("Mean Dune Height Per Transect 2005 - 2019")

p + scale_color_brewer(palette = "Dark2")

ggsave("duneHts_jitterplots_geographic.pdf", width = 5, height = 5, units = "in", dpi = 300)
dev.off()

### -----------------WILCOX RANK SUM (MANN WHITNEY) non parametric - TOWN VS NWR --------------###
wtst = wilcox.test(HD_mean ~ Zone, data = (df[df$mzbl != 'inoi',]))   ### dune heights
print(wtst)

# data = (df[df$mzbl != 'inoi',])
# d = (mean(data$HD_mean[data$Zone=='Developed'], na.rm = TRUE)-mean(data$HD_mean[data$Zone=='NWR'], na.rm = TRUE))/
#   sqrt((var(data$HD_mean[data$Zone=='Developed'], na.rm = TRUE)+var(data$HD_mean[data$Zone=='NWR'], na.rm = TRUE))/2)
# d
# ## Cohens D 

### ----------------- ADDING BOX/violin PLOTS of shoreline change  ------------------###
setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/FinalWriting/Figure_drafts")
pdf("ShoreLineCh_jitterplots_geographic.pdf")
#x11()

p = ggplot(df[df$mzbl != "inoi",], aes(x = mzbl, y = slcR, color = mzbl)) +   ### m above MHW based on Duck station
  #geom_violin(aes(x = reorder(mzbl, mzbl.order), y = slcR),trim = FALSE)+
  geom_boxplot(aes(x = reorder(mzbl, mzbl.order), y = slcR),outlier.shape = NA, width = 0.2) +
  geom_jitter(alpha = 0.3, size = 0.9)+
  scale_y_continuous( name = "Shoreline Change Rate (m/yr)") +
  labs(y = "Shoreline Change Rate (m/yr)", x = "regions")+
  ggtitle("Shoreline Change Rates 2005 - 2019")

p + scale_color_brewer(palette = "Dark2")

ggsave("ShoreLineCh_jitterplots_geographic.pdf", width = 5, height = 5, units = "in", dpi = 300)
dev.off()


### -----------------WILCOX RANK SUM (MANN WHITNEY) non parametric  --------------###
for(mm in mzbl.label){
  for (m in mzbl.label[mzbl.label != mm]){
    wtst = wilcox.test(slcR ~ mzbl, data = (df[df$mzbl == mm | df$mzbl == m,]))   ### shoreline change rate
    if(wtst$p.value >= 0.05){     ## ONLY PRINTING ONES NOT SIG DIFF. 
      print(c(mm,m))
      print(wtst)
    }
  }
} 

### -----------------WILCOX RANK SUM (MANN WHITNEY) non parametric - TOWN VS NWR --------------###
wtst = wilcox.test(slcR ~ Zone, data = (df[df$mzbl != 'inoi',]))   ### shoreline change rate
print(wtst)
# 
# data = (df[df$mzbl != 'inoi',])
# d = (mean(data$slcR[data$Zone=='Developed'], na.rm = TRUE)-mean(data$slcR[data$Zone=='NWR'], na.rm = TRUE))/
#   sqrt((var(data$slcR[data$Zone=='Developed'], na.rm = TRUE)+var(data$slcR[data$Zone=='NWR'], na.rm = TRUE))/2)
# d
# ## Cohens D 
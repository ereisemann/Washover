################################################################################
### Eve Eisemann -- 2024
### Function to calculate a mean profile across specified transect numbers
### range of values (5-95%) for each distance from the shoreline also included
################################################################################
################################################################################
###                                                                          ###
###                          the function                                    ###
###                                                                          ###
################################################################################
###
###   INPUTS   
###   elev  ::: (elevation vector) contains elevations in column for one time step
###   tdist :::  (transect distance) is distance along the transect with zero where desired, matching up with edf
###          function was built for 1-m, or ~ whole number increments. no interpolation is done here 
###   tid   ::: (transect id) is the transect id, matching up with edf
### 
###   RETURNS 
###   data frame with rounded tdist (to the whole number), 
###
################################################################################


AverageProfile = function(elev, tdist, tid){
  df = data.frame(elev,tdist, tdistR = round(tdist, digits = 0), tid)
  df =  df[order(df$tdistR),]
  avg_df = df %>% group_by(tdistR) %>% summarise(avg = mean(elev, na.rm = TRUE))
  sd_df = df %>% group_by(tdistR) %>% summarise(sd = sd(elev, na.rm = TRUE))
  
  return(data.frame(tdistR = avg_df$tdistR, avg = avg_df$avg, sd = sd_df$sd))
}
################################################################################


################################################################################
###                                                                          ###
###                   applying the function                                  ###
###                                                                          ###
################################################################################

library(tidyverse)
library(ggplot2)
library(viridis)
library(RColorBrewer)



setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/R_analysis/FluxCalcs") ### pc
Pr = read.table("Profiles.txt", header = TRUE, sep = ',')            ## raw profile elevation data for each lidar dataset (bias corrected)

###                                                                          ###
###                    defining morpho zones                                 ###
###                                                                          ###
#mzoneboundaries = c(0,219, 403, 492, 819, 916, 1009, 1300, 1398) 
#mzonenames = c("DD","DND","NND","ND","NND","ND","NND","inlet")
#mzonenames2 =  c("a_DD","b_DND","c_NND","d_ND","e_NND","f_ND","g_NND","h_inlet")
### a_DD: salvo/waves
### b_DND: Rodanthe
### c_NND: s-curves ephemeral inlet zone
### d_ND: 'stable zone' in NWR with dune ridges
### e_NND: new inlet
### f_ND: NWR highly managed dunes, some veg on them
### g_NND: NWR highly managed dunes, little veg on them (road closest to shoreline here)

### OPTION A Lumped zones
#mzbl = c(0,250, 403, 492, 819, 916, 1300, 1398)    ## m zone boundaries lumped: more morpho behavior designations. lumping f and g together hummocky/small/no veg dunes in NWR
#mzbl = c(0,250, 403, 492, 819, 851, 1272, 1398)    ## m zone boundaries lumped: more morpho behavior designations. lumping f and g together hummocky/small/no veg dunes in NWR
#mzoneL.names =  c("aa","bb","c_NND","d_ND","e_NND","f_g","h_inlet")
#mzbl.label = c("tnat", "tman","insc","nnat", "inni", "nman","inoi" )
### aa :     tnat: 'natural' behavior in towns- shoreline retreat, dunes shrink, shoreline prograde, dunes grow
### bb :     tman: 'managed' behavior in towns- no relationship between shoreline change and dune growth
### c_NND:   insc: inlet zone -  s-curves ephemeral inlet
### d_ND:    nnat: 'natural' behavior in NWR (stable zone, dune ridges)
### e_NND:   inni: inlet zone - new inlet
### f_g:     nman: 'managed' behavior in NWR (highly managed dunes, road close to shoreline)
### h_inlet: inoi: inlet zone - Oregon Inlet

### OPTION B LUMPED ZONES
mzbl = c(0,250, 403, 492, 819, 851, 1009, 1272, 1398)    ## m zone boundaries lumped: more morpho behavior designations. lumping f and g together hummocky/small/no veg dunes in NWR
mzbl.label = c("tnat", "tman","insc","nnat", "inni", "nmas","nman","inoi" )
names = c("Waves & Salvo", "Rodanthe", "S-Curves Inlet", "NWR C", "New Inlet", "NWR B", "NWR A", "Oregon Inlet")
ord = c(8:1)  ### 1 = oregon inlet (north to south)
#i = (1:length(ord))
#mzbl.order = c(rep(mzbl.ord[i], (mzbl.boundaries[i+1]-mzbl.boundaries[i])))

### aa :     tnat: 'natural' behavior in towns- shoreline retreat, dunes shrink, shoreline prograde, dunes grow
### bb :     tman: 'managed' behavior in towns- no relationship between shoreline change and dune growth
### c_NND:   insc: inlet zone -  s-curves ephemeral inlet
### d_ND:    nnat: 'natural' behavior in NWR (stable zone, dune ridges)
### e_NND:   inni: inlet zone - new inlet
### f_g:     nman: 'managed' behavior in NWR (highly managed dunes, road close to shoreline)
### f_g:     nmas: 'managed' behavior south ponds
### h_inlet: inoi: inlet zone - Oregon Inlet

yrs = c("d2005_","d2009_","d2011_","d2012_","d2014_","d2016_","d2017_","d2018a","d2018b","d2019a","d2019b")


### Creating Average Profiles 
adf = data.frame()
for(j in 1:length(yrs)){
  for(i in 1:length(mzbl-1)){
    Pr.mzone = Pr[Pr$Trans_ID >= mzbl[i] & Pr$Trans_ID < mzbl[i+1],]
    ya = AverageProfile(Pr.mzone[,4+j], Pr.mzone$tdist, Pr.mzone$Trans_ID)
    ya$mzone = rep(mzbl.label[i],length(ya$tdistR))
    ya$mzone.name = rep(names[i],length(ya$tdistR))
    ya$mzone.ord = rep(ord[i], length(ya$tdistR))
    ya$year = rep(yrs[j], length(ya$tdistR))
    adf = rbind(adf,ya)
  }
}
adf = adf %>% drop_na(mzone)     ## rows of NA for some reason. drop rows with NA as mzone 

################################################################################
###                                                                          ###
###                   plotting some things                                   ###
###                                                                          ###
################################################################################

 ##### PROFILES - average for each zone ####
#adf$mzone = factor(adf$mzone, levels = rev(mzbl.label))


x11()
plt = ggplot(adf[adf$mzone != "inoi",], aes(x = tdistR, y = avg, color = year))+
  geom_line(size = 1.3)+
  facet_wrap(~mzone.ord, ncol = 1, labeller = labeller(mzone.ord = c("2"= "NWR A", "3"="NWR B",
                                                                     "4"="New Inlet", "5"="NWR C", "6"="S-curves", 
                                                                     "7"="Rodanthe", "8"="Waves & Salvo")))
  
plt + 
  scale_color_viridis(discrete = TRUE)+
  labs(y = "Elevation (NAVD 88)", x = "distance along transect (m)" )+
  guides(color = guide_legend(title = "Year"))+
  theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  scale_x_continuous(limits = c(0,500), breaks = pretty(c(0:500), n = 10), name = "distance along transect (m)")



 #### DENSITY SETBACK plus PROFILES ###
DevInt = read.table("DevelopmentIntersects.csv", header = TRUE, sep = ',')   
DevInt$mzone.ord = 9- DevInt$mzbl.order

setwd("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/FinalWriting/Figure_drafts")
pdf("AvgProfiles_DevSetback.pdf")
#x11()

coeff = 0.05

plt = ggplot(data = DevInt[DevInt$mzbl.name != "Oregon Inlet",], aes(x = tdist))+
  geom_histogram(binwidth = 6, alpha = 0.6)+
  geom_line(data = adf[adf$mzone != "inoi",], aes(x = tdistR, y = avg/coeff, col = year), size = 1.3)+
  #geom_histogram(aes(y = after_stat(count/sum(count))), binwidth = 6)+
  #geom_histogram(aes(y=..density..), position='identity',binwidth=6)+
  facet_wrap(~mzone.ord, ncol = 1, labeller = labeller(mzone.ord = c("2"= "NWR A", "3"="NWR B", 
                                                                     "4"="New Inlet","5"="NWR C", "6"="S-curves", 
                                                                     "7"="Rodanthe", "8"="Waves & Salvo")))+
  scale_y_continuous(
    # Features of the first axis
    name = "Development Intersect Count",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="Elevation (meters NAVD 88)")
  ) +
  scale_x_continuous(limits = c(0,500))

ts = 15
plt + scale_color_viridis(discrete = TRUE)
#+
#  theme(axis.text.y=element_text(size=ts),axis.text.x=element_text(size=ts),
#        axis.title.x = element_text(size= ts), axis.title.y = element_text(size=ts))


ggsave("AvgProfiles_DevSetback.pdf", width = 6, height = 10, units = "in", dpi = 300)
dev.off()

### Density &  Setback Vs. Washover Flux ###

QOW_MOW = read.table("C:/Users/eveve/OneDrive - University of North Carolina at Chapel Hill/PhDprojects/DUNEX/R_analysis/FluxCalcs/QOW_MOW.txt", header = TRUE, sep = ",")              ## EQUILIBRIUM OVERWASH FLUX (QOW) based on Carruthers et al., and Measured Overwash Flux from lidar differencing
mzbl.boundaries = c(0,250, 403, 492, 819, 851, 1009, 1272, 1398)    ## m zone boundaries lumped: more morpho behavior designations. lumping f and g together hummocky/small/no veg dunes in NWR
names = c("Waves & Salvo", "Rodanthe", "S-Curves Inlet", "NWR C", "New Inlet", "NWR B", "NWR A", "Oregon Inlet")
i = (1:length(names))
mzbl.name = c(rep(names[i], (mzbl.boundaries[i+1]-mzbl.boundaries[i])))
mzbl.ord = c(1:8)
mzbl.order = c(rep(mzbl.ord[i], (mzbl.boundaries[i+1]-mzbl.boundaries[i])))
QOW_MOW$mzbl.order = mzbl.order
QOW_MOW$mzbl.name = mzbl.name

DevInt = read.table("DevelopmentIntersects.csv", header = TRUE, sep = ',')   
DevInt_Count_Setback = 
  DevInt %>%
  group_by(Trans_ID)%>%
  summarize(DevIntCount = mean(DevIntCount, na.rm = TRUE), DevInt_SB = min(tdist, na.rm = TRUE))
  
#DevInt_SetBack = 
#  DevInt %>%
#  group_by(Trans_ID)%>%
#  summarize(DevInt_SB = min(tdist, na.rm = TRUE))

DevInt_Count_Setback$MOW = rep(NA, length(DevInt_Count_Setback$Trans_ID))
DevInt_Count_Setback$mzbl.name = rep(NA, length(DevInt_Count_Setback$Trans_ID))
DevInt_Count_Setback$Zone = rep(NA, length(DevInt_Count_Setback$Trans_ID))   #added for stat analysis 8/14/24

for(i in DevInt_Count_Setback$Trans_ID){
  DevInt_Count_Setback$MOW[DevInt_Count_Setback$Trans_ID == i] = QOW_MOW$MOW[QOW_MOW$Trans_ID == i]
  DevInt_Count_Setback$mzbl.name[DevInt_Count_Setback$Trans_ID == i] = QOW_MOW$mzbl.name[QOW_MOW$Trans_ID == i]
  DevInt_Count_Setback$Zone[DevInt_Count_Setback$Trans_ID == i] =  QOW_MOW$Zone[QOW_MOW$Trans_ID == i]  ## added for stat analysis 8/14/24
}


### PLOTTING - development density vs measured over wash
DevInt_Count_Setback %>%
  filter(mzbl.name != 'Oregon Inlet') %>%
  mutate(DICbins = cut(DevIntCount, breaks = seq(from = 0, to = 21, by = 3))) %>%
  
  ggplot() +
  geom_boxplot(aes( x = DICbins, y = MOW))+
  geom_jitter(aes(x = DICbins, y = MOW, col = mzbl.name), alpha = 0.3)+
  #geom_point(alpha = 0.6)+
  scale_color_brewer(palette = "Dark2")


### PLOTTING - development density vs setback 
DevInt_Count_Setback %>%
  filter(mzbl.name != 'Oregon Inlet') %>%
  
  ggplot(aes(x = DevInt_SB, y = MOW, col = mzbl.name)) +
  geom_point(alpha = 0.6)+
  scale_color_brewer(palette = "Dark2")

### PLOTTING - development density vs setback in bins
DevInt_Count_Setback %>%
  filter(mzbl.name != 'Oregon Inlet') %>%
  mutate(SBbins = cut(DevInt_SB, breaks = seq(from = 0, to = 600, by = 50))) %>%
  
  
  ggplot() +
  geom_boxplot(aes( x = SBbins, y = MOW), na.rm=TRUE, outlier.shape = NA)+
  geom_jitter(aes(x = SBbins, y = MOW, col = mzbl.name), alpha = 0.6)
  #geom_point(alpha = 0.6)+
  #scale_color_brewer(palette = "Dark2")



### -----------------WILCOX RANK SUM (MANN WHITNEY) non parametric - TOWN VS NWR --------------###
wtst = wilcox.test(DevIntCount ~ Zone, data = (DevInt_Count_Setback[DevInt_Count_Setback$mzbl.name != 'Oregon Inlet',]))   ### development density 
print(wtst)

medians <- DevInt_Count_Setback %>%
  group_by(Zone) %>%
  summarize(median_value = median(DevIntCount))

print(medians)

means <- DevInt_Count_Setback %>%
  group_by(Zone) %>%
  summarize(mean_value = mean(DevIntCount))

print(means)


x11()
ggplot(DevInt_Count_Setback, aes(x = DevIntCount, fill = Zone)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.7) +
  geom_vline(data = means, aes(xintercept = mean_value, color = Zone),
             linetype = "dashed", size = 1) +
  labs(title = "Development Density",
       x = "intersect count",
       y = "transect count") +
  theme_minimal()


### -----------------WILCOX RANK SUM (MANN WHITNEY) non parametric - TOWN VS NWR --------------###
wtst = wilcox.test(DevInt_SB ~ Zone, data = (DevInt_Count_Setback[DevInt_Count_Setback$mzbl.name != 'Oregon Inlet',]))   ### development setback
print(wtst)

medians <- DevInt_Count_Setback %>%
  group_by(Zone) %>%
  summarize(median_value = median(DevInt_SB))

print(medians)

means <- DevInt_Count_Setback %>%
  group_by(Zone) %>%
  summarize(mean_value = mean(DevInt_SB))

print(means)


x11()
ggplot(DevInt_Count_Setback, aes(x = DevInt_SB, fill = Zone)) +
  geom_histogram(binwidth = 10, position = "identity", alpha = 0.7) +
  geom_vline(data = means, aes(xintercept = mean_value, color = Zone),
             linetype = "dashed", size = 1) +
  labs(title = "Development Setback",
       x = "setback distance (m)",
       y = "transect count") +
  theme_minimal()




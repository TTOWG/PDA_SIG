# To the Only Wise God - TTOWG


#Agbabu Bitumen Deposit Individual Porosity Data Analysis

# Importing required libraries
library(sp)
library(gstat)
library(dbscan)
library(som)
library(tidyverse)
library(ggExtra)
library(Rmisc)
library(gridExtra)
library(grid)

#Importing data from raw spreadsheet; data source: Adegoke et al (1980)
RawData = read.csv(file.choose(),header=T)
RawData$BoreHole = as.factor(RawData$BoreHole) #Ensuring R sees the Borehole numbers as ID and not as numerals

#Setting parameters
BitumenSG = 0.9679 #Falebita et al (2014)
GrainDensity = 2.58 #Ola (1991)
PoroCutoff = 0.5

#computations
WaterContentWTPCT = RawData$WET.TAR.WEIGHT..-RawData$DRY.TAR.WEIGHT..

TarWtgram = 0.01*RawData$DRY.TAR.WEIGHT..*20
WaterWTgram = 0.01*WaterContentWTPCT*20
GrainWTgram = 0.01*(100-RawData$WET.TAR.WEIGHT..)*20

TarVolCC = TarWtgram/BitumenSG
WaterVolCC = WaterWTgram/1
GrainVolCC = GrainWTgram/GrainDensity

Porosity = (TarVolCC+WaterVolCC)/(TarVolCC+WaterVolCC+GrainVolCC)

#Constructing the data table to include raw and computed data
ComputedData = cbind(RawData,WaterContentWTPCT,TarWtgram,WaterWTgram,GrainWTgram,TarVolCC,WaterVolCC,GrainVolCC,Porosity)
remove(WaterContentWTPCT,TarWtgram,WaterWTgram, GrainWTgram, TarVolCC, WaterVolCC,GrainVolCC, Porosity,pos=1)

#Excluding spurious/suspicious data: data rows with porosity above the cut off 
TrimmedData = ComputedData[ComputedData$Porosity<PoroCutoff&!ComputedData$Horizon=="", ]
ActualDepths = read.csv(file.choose(),header=T)
TrimmedData$DEPTH.m. =ActualDepths
Coordinates = round(read.csv(file.choose(),header=T),digits = 2)

SmartTrimmedData = data.frame(TrimmedData$BoreHole,TrimmedData$DEPTH.m.,round(TrimmedData$Porosity, digits = 4),Coordinates$Easting,Coordinates$Northing)



#On CUT2 data (Samples between Q1-IQR and Q3+IQR )
SpatialTrimmedDataCUT2 = SmartTrimmedData[SmartTrimmedData$Porosity<0.3994&SmartTrimmedData$Porosity>0.0826, ]

names(SpatialTrimmedDataCUT2) = c("BoreHoleID", "Depth", "Porosity", "x", "y")
coordinates(SpatialTrimmedDataCUT2) = c("x","y","Depth")


# Obtaining variogram cloud
Porosity_VariogramCloudCUT2_90_22.5 = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, alpha = 90, tol.hor = 22.5, cloud = TRUE)
ggplot(data = Porosity_VariogramCloudCUT2_90_22.5)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram")+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.06, by = 0.01))
ggsave("Porosity Variogram Cloud_CUT2_90_22.5.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

CloudMarginal = ggplot(data = Porosity_VariogramCloudCUT2_90_22.5)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram")+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.06, by = 0.01))
ggMarginal(CloudMarginal, type = "histogram", margins = "x")
ggsave("Porosity Variogram Cloud_CUT2_90_22.5_marginal.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


#intersamplepointdistances_histogram_90 = hist(Porosity_VariogramCloudCUT2_90_22.5$dist,freq = T, breaks = seq(from = 0, to = 8000, by = 500),col = "grey", ylim =  range(0,5000), xlab="Inter-sample-point Distances",ylab = "Frequency", las = 3)

write.table(Porosity_VariogramCloudCUT2_90_22.5,"Porosity_VariogramCloudCUT2_90_22.5.csv", sep = ",", row.names = F)

# Obtaining empirical variogram - conventional approach (bin-based)
Porosity_VariogramCUT2_90_22.5 = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, alpha = 90, tol.hor = 22.5)
ggplot(data = Porosity_VariogramCUT2_90_22.5)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))+
  geom_line(color = "red",aes(x = dist, y = gamma))
ggsave("Porosity Emperical Variogram_CUT2_90_22.5.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


# Fitting to Model - conventional approach (bin-based)
binbasedVariogModel_90 = fit.variogram(Porosity_VariogramCUT2_90_22.5, vgm("Sph","Exp","Gau"),fit.sills = TRUE, fit.ranges = TRUE)
binbasedFitRsquared_90 = function(expvar,modvar){ 
  SSErr = attr(modvar,"SSErr") 
  weight = expvar$np/(expvar$dist^2) #This is the weight factor for the default fit.method  
  SST = sum((weight*(expvar$gamma-mean(expvar$gamma)))^2)
  R_sq<-1-SSErr/SST
  return(R_sq)
}
binbasedFitRsqd_90 = round(binbasedFitRsquared_90(Porosity_VariogramCUT2_90_22.5,binbasedVariogModel_90),4)
binbasedRsqdposter_90 = Porosity_VariogramCUT2_90_22.5 %>% 
  summarise(dist = 1000,gamma = 0.0002,binbasedRsqdposterstring_90 = paste("Goodness of Fit, R-squared = ",binbasedFitRsqd_90,sep=""))
binbasedVariogModelFunc_90 = function(x){
  ifelse(x <= binbasedVariogModel_90$range[2], (binbasedVariogModel_90$psill[1]+(binbasedVariogModel_90$psill[2]*((1.5*(x/binbasedVariogModel_90$range[2]))-(0.5*((x/binbasedVariogModel_90$range[2])^3))))), sum(binbasedVariogModel_90$psill))
}

ggplot(data = Porosity_VariogramCUT2_90_22.5, aes(x = dist, y = gamma))+
  geom_point(aes(shape = "mu"), color = "blue")+
  labs(x = "Lag Distance, h (m)", y = expression(Variogram~gamma))+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, max(Porosity_VariogramCUT2_90_22.5$gamma,0.005)))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, max(max(Porosity_VariogramCUT2_90_22.5$gamma)+0.001,0.005), by = 0.001))+
  stat_function(aes(x = dist, color = 'jj'), fun = binbasedVariogModelFunc_90, xlim = c(0, 8000))+
  geom_text(aes(label = binbasedRsqdposterstring_90 ), data = binbasedRsqdposter_90, vjust = "bottom", hjust = "left")+
  scale_colour_manual(name = '', values =c('jj'='red'), labels = expression(model:0.0034+0.0011*Sph[4766]*(h)))+
  scale_shape_manual(name = '', values = c('mu' = 16), labels = c('empirical'))+
  theme(legend.position=c(0.6, 0.3))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
ggsave("Porosity Emperical&Fitted Variogram_CUT2_90_22.5_binbased.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis//Model Fits", dpi = 96)


# Obtaining empirical variogram - dbscan-aided approach

# Setting parameters
kNNdistplot(Porosity_VariogramCloudCUT2_90_22.5[,c(2,3)],k = 500)
abline(h=100, col = "red", lty=2)
# This yields eps = 100 and minPts = 500
porodbscan_90 = dbscan(Porosity_VariogramCloudCUT2_90_22.5[,c(2,3)], eps = 100, minPts = 500)

dbscanVariogCloud_90 = data.frame(Porosity_VariogramCloudCUT2_90_22.5$dist, Porosity_VariogramCloudCUT2_90_22.5$gamma, porodbscan_90$cluster, Direction = 90)    
dbscanVariogCloud_90$porodbscan_90.cluster = as.factor(dbscanVariogCloud_90$porodbscan_90.cluster) #Ensuring R sees the cluster numbers as classification ID and not as numerals
ggplot(data = dbscanVariogCloud_90[!dbscanVariogCloud_90$porodbscan_90.cluster=="0",])+
  geom_point(mapping = aes(x = Porosity_VariogramCloudCUT2_90_22.5.dist, y = Porosity_VariogramCloudCUT2_90_22.5.gamma, color = porodbscan_90.cluster, shape = porodbscan_90.cluster))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Cluster" )+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))+
  scale_color_manual(name = 'Cluster', values = c("#000000","#00FF00","#330033","#000099","#000033","#FF0000"))+
  scale_shape_manual(name = 'Cluster', values = c(16,8,15,17,18,4))
ggsave("Porosity Variogram Cloud_CUT2_90_22.5_colored_dbsacan.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

NumberofClusters_90 = length(levels(dbscanVariogCloud_90$porodbscan_90.cluster))
dbscanExpVariog_90 = matrix(0, NumberofClusters_90-1, 3, dimnames = list(c("1","2","3","4","5","6"),c("NumberofPairs", "AverageDistance", "ExpVariog")))
for(i in 1:NumberofClusters_90-1){
  dbscanExpVariog_90[i,1]=length(dbscanVariogCloud_90$Porosity_VariogramCloudCUT2_90_22.5.dist[dbscanVariogCloud_90$porodbscan_90.cluster==i]); 
  dbscanExpVariog_90[i,2]=mean(dbscanVariogCloud_90$Porosity_VariogramCloudCUT2_90_22.5.dist[dbscanVariogCloud_90$porodbscan_90.cluster==i]);
  dbscanExpVariog_90[i,3]=mean(dbscanVariogCloud_90$Porosity_VariogramCloudCUT2_90_22.5.gamma[dbscanVariogCloud_90$porodbscan_90.cluster==i])
}
dbscanExpVariog_90 = data.frame(dbscanExpVariog_90)
ggplot(data = dbscanExpVariog_90)+
  geom_point(mapping = aes(x = AverageDistance, y = ExpVariog), color = "blue")+
  geom_line(mapping = aes(x = AverageDistance, y = ExpVariog), color = "red")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.005))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.005, by = 0.001))
ggsave("Porosity Emperical Variogram_CUT2_90_22.5_dbscan.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

#Fitting to Model - dbscan-aided approach
dbscanExpVariog_90_for_fitting = dbscanExpVariog_90
names(dbscanExpVariog_90_for_fitting) = c("np","dist","gamma")
class(dbscanExpVariog_90_for_fitting) = c("gstatVariogram", "data.frame")
#dbscanVariogModel_90 = fit.variogram(dbscanExpVariog_90_for_fitting, vgm("Sph",nugget = 0.0024),fit.sills = FALSE, fit.ranges = TRUE)
dbscanVariogModel_90 = fit.variogram(dbscanExpVariog_90_for_fitting, vgm("Sph",psill = 0.0021, nugget = 0.0024, range = 3500),fit.sills = FALSE, fit.ranges = FALSE)
#dbscanFitRsquared_90 = function(expvar,modvar){ 
#SSErr = attr(modvar,"SSErr") 
#weight = expvar$np/(expvar$dist^2) #This is the weight factor for the default fit.method  
#SST = sum((weight*(expvar$gamma-mean(expvar$gamma)))^2)
#R_sq<-1-SSErr/SST
#return(R_sq)
#}
#dbscanFitRsqd_90 = round(dbscanFitRsquared_90(dbscanExpVariog_90_for_fitting,dbscanVariogModel_90),4)
#dbscanRsqdposter_90 = dbscanExpVariog_90_for_fitting %>% 
#  summarise(dist = 1000,gamma = 0.0002,dbscanRsqdposterstring_90 = paste("Goodness of Fit, R-squared = ",dbscanFitRsqd_90,sep=""))
dbscanVariogModelFunc_90 = function(x){
  ifelse(x <= dbscanVariogModel_90$range[2], (dbscanVariogModel_90$psill[1]+(dbscanVariogModel_90$psill[2]*((1.5*(x/dbscanVariogModel_90$range[2]))-(0.5*((x/dbscanVariogModel_90$range[2])^3))))), sum(dbscanVariogModel_90$psill))
}

ggplot(data = dbscanExpVariog_90_for_fitting, aes(x = dist, y = gamma))+
  geom_point(aes(shape = "mu"), color = "blue")+
  labs(x = "Lag Distance, h (m)", y = expression(Variogram~gamma))+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, max(dbscanExpVariog_90_for_fitting$gamma,0.005)))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, max(max(dbscanExpVariog_90_for_fitting$gamma)+0.001,0.005), by = 0.001))+
  stat_function(aes(x = dist, color = 'jj'), fun = dbscanVariogModelFunc_90, xlim = c(0, 8000))+
  #geom_text(aes(label = dbscanRsqdposterstring_90 ), data = dbscanRsqdposter_90, vjust = "bottom", hjust = "left")+
  scale_colour_manual(name = '', values =c('jj'='red'), labels = expression(model:0.0024+0.0021*Sph[3500]*(h)))+
  scale_shape_manual(name = '', values = c('mu' = 16), labels = c('empirical'))+
  theme(legend.position=c(0.6, 0.3))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
ggsave("Porosity Emperical&Fitted Variogram_CUT2_90_22.5_dbscan.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis//Model Fits", dpi = 96)



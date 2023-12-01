# To the Only Wise God - TTOWG


#Agbabu Bitumen Deposit Individual Porosity Data Analysis

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
#For export to SGeMS
SmartTrimmedData = data.frame(TrimmedData$BoreHole,TrimmedData$DEPTH.m.,round(TrimmedData$Porosity, digits = 4),Coordinates$Easting,Coordinates$Northing)
write.table(SmartTrimmedData,"SmartTrimmedData.txt", row.names = F)
SmartTrimmedData1 = data.frame(TrimmedData$BoreHole,TrimmedData$Horizon, TrimmedData$DEPTH.m.,round(TrimmedData$Porosity, digits = 4),Coordinates$Easting,Coordinates$Northing)
SmartTrimmedData_X = SmartTrimmedData1[SmartTrimmedData1$TrimmedData.Horizon=="X",]
write.table(SmartTrimmedData_X,"SmartTrimmedData_X.txt", row.names = F)
SmartTrimmedData_Y = SmartTrimmedData1[SmartTrimmedData1$TrimmedData.Horizon=="Y",]
write.table(SmartTrimmedData_Y,"SmartTrimmedData_Y.txt", row.names = F)

write.csv(SmartTrimmedData_X, "Horizon_X_Poro.csv", row.names = F)
write.csv(SmartTrimmedData_Y, "Horizon_Y_Poro.csv", row.names = F)

#Engineering and natural indexing of sample points - for PSRF use
delta_X = 400
delta_Y = 400
delta_Z = 1
nx = 40
ny = 13
nz = 100
SampleLocation = data.frame(matrix(0,408,7))
names(SampleLocation) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth")

for(i in 1:408){
  SampleLocation$x_index[i] = ceiling((SmartTrimmedData$Easting[i] - 700000)/delta_X)
  SampleLocation$y_index[i] = ceiling((SmartTrimmedData$Northing[i] - 732500)/delta_Y)
  SampleLocation$z_index[i] = ceiling((SmartTrimmedData$Depth[i] - 0)/delta_Z)
  SampleLocation$GridBlockID[i] = ((SampleLocation$z_index[i]-1)*nx*ny)+((SampleLocation$y_index[i]-1)*nx)+SampleLocation$x_index[i]
  SampleLocation$TrueEasting[i] = SmartTrimmedData$Easting[i]
  SampleLocation$TrueNorthing[i] = SmartTrimmedData$Northing[i]
  SampleLocation$TrueDepth[i] = SmartTrimmedData$Depth[i]
}
write.csv(SampleLocation, file = "C://Users//TTOWG//Documents//R trials//trials//Sample Location.csv", row.names = F)
#Trimming to mimick CUT2 data sample location so as to get the desired cloud configurations
TrimmedSampleLocation = SampleLocation[-c(as.numeric(rownames(Complement_SpatialTrimmedDataCUT2))),]
write.csv(TrimmedSampleLocation, file = "C://Users//TTOWG//Documents//R trials//trials//Trimmed Sample Location.csv", row.names = F)


#Sample points and data for uncertainty analysis
SamplePointandData = data.frame(matrix(0,408,8))
names(SamplePointandData) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth", "Porosity")

for(i in 1:408){
  SamplePointandData$x_index[i] = ceiling((SmartTrimmedData$Easting[i] - 700000)/delta_X)
  SamplePointandData$y_index[i] = ceiling((SmartTrimmedData$Northing[i] - 732500)/delta_Y)
  SamplePointandData$z_index[i] = ceiling((SmartTrimmedData$Depth[i] - 0)/delta_Z)
  SamplePointandData$GridBlockID[i] = ((SampleLocation$z_index[i]-1)*nx*ny)+((SampleLocation$y_index[i]-1)*nx)+SampleLocation$x_index[i]
  SamplePointandData$TrueEasting[i] = SmartTrimmedData$Easting[i]
  SamplePointandData$TrueNorthing[i] = SmartTrimmedData$Northing[i]
  SamplePointandData$TrueDepth[i] = SmartTrimmedData$Depth[i]
  SamplePointandData$Porosity[i] = SmartTrimmedData$Porosity[i]
}
write.csv(SamplePointandData, file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis//Sample Points and Data.csv", row.names = F)
#Trimming to mimick CUT2 data sample location so as to get the desired cloud configurations
TrimmedSamplePointandData = SamplePointandData[-c(as.numeric(rownames(Complement_SpatialTrimmedDataCUT2))),]
write.csv(TrimmedSamplePointandData, file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis//Trimmed Sample Points and Data.csv", row.names = F)


#For ArcGIS Pro use
names(SmartTrimmedData) = c("BoreHoleID", "Depth", "Porosity", "Easting", "Northing")
write.table(SmartTrimmedData,"SmartTrimmedData.csv", sep = ",", row.names = F)


library(sp)
library(gstat)
library(dbscan)
library(som)
library(tidyverse)
library(ggExtra)
library(Rmisc)
library(gridExtra)
library(grid)


#Exploratory Spatial Data Analysis (ESDA) - to detect presence of spatial correlation.

#On whole data (no CUTOFF other than Porocutoff)
SpatialTrimmedDataWhole = SmartTrimmedData
names(SpatialTrimmedDataWhole) = c("BoreHoleID", "Depth", "Porosity", "x", "y")
coordinates(SpatialTrimmedDataWhole) = c("x","y","Depth")
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Porosity Lagged Scatter Plot_whole.jpg")
hscat(Porosity~1,SpatialTrimmedDataWhole,(0:6)*1000)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Porosity Variogram Map_whole.jpg")
Porosity_Variogram_mapWhole = variogram(Porosity~1,SpatialTrimmedDataWhole, map = TRUE,cutoff=10000,width=1000)
plot(Porosity_Variogram_mapWhole)
dev.off()

Porosity_VariogramCloudWhole = variogram(Porosity~1,SpatialTrimmedDataWhole, cloud = TRUE)
ggplot(mapping = aes(x = dist, y = gamma) )+
  geom_point(data = Porosity_VariogramCloudWhole[c(-291,-38206,-38211,-39009),], color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram")+
  scale_x_continuous(breaks = seq(0, 6000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.12, by = 0.02))+
  geom_point(data = Porosity_VariogramCloudWhole[c(291,38206,38211,39009),], color = "red", shape = 12, size = 3)
ggsave("Porosity Variogram Cloud_whole2.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

PoroPointPairs=plot(Porosity_VariogramCloudWhole, identify = TRUE, digitize = TRUE, xlim = c(0,6000), ylim = c(0,0.14))
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Porosity Variogram CloudPoinPairs_whole.jpg")
plot(PoroPointPairs, SpatialTrimmedDataWhole)
dev.off()

Porosity_VariogramWhole = variogram(Porosity~1,SpatialTrimmedDataWhole)
ggplot(data = Porosity_VariogramWhole)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 6000), ylim = c(0, 0.012))+
  scale_x_continuous(breaks = seq(0, 6000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.012, by = 0.002))
ggsave("Porosity Emperical Variogram_whole.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

#On CUT2 data (Samples between Q1-IQR and Q3+IQR )
SpatialTrimmedDataCUT2 = SmartTrimmedData[SmartTrimmedData$Porosity<0.3994&SmartTrimmedData$Porosity>0.0826, ]
SpatialTrimmedDataCUT2_copy = SmartTrimmedData1[SmartTrimmedData1$Porosity<0.3994&SmartTrimmedData1$Porosity>0.0826, ]
SpatialTrimmedDataCUT2_X = SpatialTrimmedDataCUT2_copy[SpatialTrimmedDataCUT2_copy$Horizon=="X",]
SpatialTrimmedDataCUT2_Y = SpatialTrimmedDataCUT2_copy[SpatialTrimmedDataCUT2_copy$Horizon=="Y",]

write.csv(SpatialTrimmedDataCUT2_X, "Horizon_X_Poro_CUT2.csv", row.names = F)
write.csv(SpatialTrimmedDataCUT2_Y, "Horizon_Y_Poro_CUT2.csv", row.names = F)

names(SpatialTrimmedDataCUT2) = c("BoreHoleID", "Depth", "Porosity", "x", "y")
coordinates(SpatialTrimmedDataCUT2) = c("x","y","Depth")

jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Porosity Lagged Scatter Plot_CUT2.jpg")
hscat(Porosity~1,SpatialTrimmedDataCUT2,(0:6)*1000)
dev.off()

jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Porosity Variogram Map_CUT2.jpg")
Porosity_Variogram_mapCUT2 = variogram(Porosity~1,SpatialTrimmedDataCUT2, map = TRUE,cutoff=5000,width=1000)
plot(Porosity_Variogram_mapCUT2)
dev.off()

Porosity_VariogramCloudCUT2 = variogram(Porosity~1,SpatialTrimmedDataCUT2, cloud = TRUE)
ggplot(data = Porosity_VariogramCloudCUT2)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram")+
  scale_x_continuous(breaks = seq(0, 6000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.06, by = 0.01))
ggsave("Porosity Variogram Cloud_CUT2.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


PoroPointPairsCUT2=plot(Porosity_VariogramCloudCUT2, identify = TRUE, digitize = TRUE, xlim = c(0,6000), ylim = c(0,0.06))
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Porosity Variogram CloudPoinPairs_CUT2.jpg")
plot(PoroPointPairsCUT2, SpatialTrimmedDataCUT2)
dev.off()

Porosity_VariogramCloudCUT2_100 = variogram(Porosity~1,SpatialTrimmedDataCUT2, alpha = 100, tol.hor=22.5, cloud = TRUE)
ggplot(data = Porosity_VariogramCloudCUT2_100)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram")+
  scale_x_continuous(breaks = seq(0, 6000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.06, by = 0.01))
ggsave("Porosity Variogram Cloud_CUT2_100.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


PoroPointPairsCUT2_100=plot(Porosity_VariogramCloudCUT2_100, identify = TRUE, digitize = TRUE, xlim = c(0,6000), ylim = c(0,0.06))
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Porosity Variogram Cloud_CUT2_100.jpg")
plot(PoroPointPairsCUT2, SpatialTrimmedDataCUT2)
dev.off()

jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Porosity Emperical Variogram_CUT2.jpg")

Porosity_VariogramCUT2_binbased = variogram(Porosity~1,SpatialTrimmedDataCUT2)
ggplot(data = Porosity_VariogramCUT2_binbased)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))
ggsave("Porosity Emperical Variogram_CUT2_bindased.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

Porosity_VariogramCUT2_lagbased = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, boundaries = c(500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500))
ggplot(data = Porosity_VariogramCUT2_lagbased)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))
ggsave("Porosity Emperical Variogram_CUT2_lagbased.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

#All lag spacing-tolerance specifications: done in GSLIB and imported to R
Porosity_VariogramCUT2_lagbased_1000by300 = read.csv("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//GV_Agbabu Field-wide Individual Porosity_1000_300_omnidir_CUT2.csv",header=F)
names(Porosity_VariogramCUT2_lagbased_1000by300) = c("np", "dist", "gamma", "specifications")
ggplot(data = Porosity_VariogramCUT2_lagbased_1000by300)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))
ggsave("Porosity Emperical Variogram_CUT2_lagbased_1000by300.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

Porosity_VariogramCUT2_lagbased_1000by500 = read.csv("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//GV_Agbabu Field-wide Individual Porosity_1000_500_omnidir_CUT2.csv",header=F)
names(Porosity_VariogramCUT2_lagbased_1000by500) = c("np", "dist", "gamma", "specifications")
ggplot(data = Porosity_VariogramCUT2_lagbased_1000by500)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))
ggsave("Porosity Emperical Variogram_CUT2_lagbased_1000by500.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


Porosity_VariogramCUT2_lagbased_1000by700 = read.csv("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//GV_Agbabu Field-wide Individual Porosity_1000_700_omnidir_CUT2.csv",header=F)
names(Porosity_VariogramCUT2_lagbased_1000by700) = c("np", "dist", "gamma", "specifications")
ggplot(data = Porosity_VariogramCUT2_lagbased_1000by700)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))
ggsave("Porosity Emperical Variogram_CUT2_lagbased_1000by700.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

Porosity_VariogramCUT2_lagbased_500by300 = read.csv("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//GV_Agbabu Field-wide Individual Porosity_500_300_omnidir_CUT2.csv",header=F)
names(Porosity_VariogramCUT2_lagbased_500by300) = c("np", "dist", "gamma", "specifications")
ggplot(data = Porosity_VariogramCUT2_lagbased_500by300)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))
ggsave("Porosity Emperical Variogram_CUT2_lagbased_500by300.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

Porosity_VariogramCUT2_lagbased_500by500 = read.csv("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//GV_Agbabu Field-wide Individual Porosity_500_500_omnidir_CUT2.csv",header=F)
names(Porosity_VariogramCUT2_lagbased_500by500) = c("np", "dist", "gamma","specifications")
ggplot(data = Porosity_VariogramCUT2_lagbased_500by500)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))
ggsave("Porosity Emperical Variogram_CUT2_lagbased_500by500.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

Porosity_VariogramCUT2_lagbased_allspecs = rbind(Porosity_VariogramCUT2_lagbased_500by300,Porosity_VariogramCUT2_lagbased_500by500,Porosity_VariogramCUT2_lagbased_1000by300,Porosity_VariogramCUT2_lagbased_1000by500,Porosity_VariogramCUT2_lagbased_1000by700)
ggplot(data = Porosity_VariogramCUT2_lagbased_allspecs, aes(x = dist, y = gamma))+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 2000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))+
  geom_point(color = "blue") + 
  facet_wrap(~ specifications, nrow = 3)
ggsave("Porosity Emperical Variogram_CUT2_lagbased_allspecs.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


#Directionals
#45 degree
Porosity_VariogramCloudCUT2_45_22.5 = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, alpha = 45, tol.hor = 22.5, cloud = TRUE)
ggplot(data = Porosity_VariogramCloudCUT2_45_22.5)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram")+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.06, by = 0.01))
ggsave("Porosity Variogram Cloud_CUT2_45_22.5.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

write.table(Porosity_VariogramCloudCUT2_45_22.5,"Porosity_VariogramCloudCUT2_45_22.5.csv", sep = ",", row.names = F)

Porosity_VariogramCUT2_45_22.5 = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, alpha = 45, tol.hor = 22.5)
ggplot(data = Porosity_VariogramCUT2_45_22.5)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.008))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.008, by = 0.001))
ggsave("Porosity Emperical Variogram_CUT2_45_22.5.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

       
# Now considering the natural cluster trends present in the clouds in declaring distance class interval
Porosity_VariogramCUT2_45_22.5_cluster = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, alpha = 45, tol.hor = 22.5, boundaries = c(100, 2450, 2805, 4010, 4850, 5500, 6000, 7100, 8000))
ggplot(data = Porosity_VariogramCUT2_45_22.5_cluster)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))
ggsave("Porosity Emperical Variogram_CUT2_45_22.5_cluster.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

# Now using DBSCAN for the cluster analysis
kNNdistplot(Porosity_VariogramCloudCUT2_45_22.5[,c(2,3)],k = 500)
abline(h=230, col = "red", lty=2)
# This yields eps = 230 and minPts = 500
porodbscan_45 = dbscan(Porosity_VariogramCloudCUT2_45_22.5[,c(2,3)], eps = 230, minPts = 500)

dbscanVariogCloud_45 = data.frame(Porosity_VariogramCloudCUT2_45_22.5$dist, Porosity_VariogramCloudCUT2_45_22.5$gamma, porodbscan_45$cluster, Direction = 45)    
dbscanVariogCloud_45$porodbscan_45.cluster = as.factor(dbscanVariogCloud_45$porodbscan_45.cluster) #Ensuring R sees the cluster numbers as classification ID and not as numerals
ggplot(data = dbscanVariogCloud_45[!dbscanVariogCloud_45$porodbscan_45.cluster=="0",])+
  geom_point(mapping = aes(x = Porosity_VariogramCloudCUT2_45_22.5.dist, y = Porosity_VariogramCloudCUT2_45_22.5.gamma, color = porodbscan_45.cluster))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Bin" )+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))
  ggsave("Porosity Variogram Cloud_CUT2_45_22.5_colored_dbsacan.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

NumberofClusters_45 = length(levels(dbscanVariogCloud_45$porodbscan_45.cluster))
dbscanExpVariog_45 = matrix(0, NumberofClusters_45-1, 3, dimnames = list(c("1","2","3","4","5","6"),c("NumberofPairs", "AverageDistance", "ExpVariog")))
for(i in 1:NumberofClusters_45-1){
  dbscanExpVariog_45[i,1]=length(dbscanVariogCloud_45$Porosity_VariogramCloudCUT2_45_22.5.dist[dbscanVariogCloud_45$porodbscan_45.cluster==i]); 
  dbscanExpVariog_45[i,2]=mean(dbscanVariogCloud_45$Porosity_VariogramCloudCUT2_45_22.5.dist[dbscanVariogCloud_45$porodbscan_45.cluster==i]);
  dbscanExpVariog_45[i,3]=mean(dbscanVariogCloud_45$Porosity_VariogramCloudCUT2_45_22.5.gamma[dbscanVariogCloud_45$porodbscan_45.cluster==i])
}
dbscanExpVariog_45 = data.frame(dbscanExpVariog_45)
ggplot(data = dbscanExpVariog_45)+
  geom_point(mapping = aes(x = AverageDistance, y = ExpVariog), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))
  ggsave("Porosity Emperical Variogram_CUT2_45_22.5_dbscan.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

  dbscanExpVariog_45_for_fitting = dbscanExpVariog_45
  names(dbscanExpVariog_45_for_fitting) = c("np","dist","gamma")
  class(dbscanExpVariog_45_for_fitting) = c("gstatVariogram", "data.frame")
  dbscanVariogModel_45 = fit.variogram(dbscanExpVariog_45_for_fitting, vgm("Sph",psill = 0.0021, nugget = 0.0024, range = 3000),fit.sills = FALSE, fit.ranges = FALSE)
  
  dbscanVariogModelFunc_45 = function(x){
    ifelse(x <= dbscanVariogModel_45$range[2], (dbscanVariogModel_45$psill[1]+(dbscanVariogModel_45$psill[2]*((1.5*(x/dbscanVariogModel_45$range[2]))-(0.5*((x/dbscanVariogModel_45$range[2])^3))))), sum(dbscanVariogModel_45$psill))
  }
  ggplot(data = dbscanExpVariog_45_for_fitting, aes(x = dist, y = gamma))+
    geom_point(aes(shape = "mu"), color = "blue")+
    labs(x = "Lag Distance, h (m)", y = expression(Variogram~gamma))+
    coord_cartesian(xlim = c(0, 8000), ylim = c(0, max(dbscanExpVariog_45_for_fitting$gamma,0.005)))+
    scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
    scale_y_continuous(breaks = seq(0, max(max(dbscanExpVariog_45_for_fitting$gamma)+0.001,0.005), by = 0.001))+
    stat_function(aes(x = dist, color = 'jj'), fun = dbscanVariogModelFunc_45, xlim = c(0, 8000))+
    scale_colour_manual(name = '', values =c('jj'='red'), labels = expression(model:0.0024+0.0021*Sph[3000]*(h)))+
    scale_shape_manual(name = '', values = c('mu' = 16), labels = c('empirical'))+
    theme(legend.position=c(0.6, 0.3))+
    theme(legend.title = element_blank())+
    theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
  ggsave("Porosity Emperical&Fitted Variogram_CUT2_45_22.5_dbscan.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis//Model Fits", dpi = 96)
  

#For comparison; the colored cloud plot without cluster analysis:
#bin-based (gstat)
CloudCUT2_45_22.5_Pointpaircategory = matrix(0, length(Porosity_VariogramCloudCUT2_45_22.5$dist), 1)
for(i in 1:length(Porosity_VariogramCloudCUT2_45_22.5$dist)){
  if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<533.33){
    CloudCUT2_45_22.5_Pointpaircategory[i,1] = 1}
  else{
    if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<1066.67){
      CloudCUT2_45_22.5_Pointpaircategory[i,1] = 2}
    else{
      if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<1600){
        CloudCUT2_45_22.5_Pointpaircategory[i,1] = 3}
      else{
        if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<2133.33){
          CloudCUT2_45_22.5_Pointpaircategory[i,1] = 4}
        else{
          if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<2666.67){
            CloudCUT2_45_22.5_Pointpaircategory[i,1] = 5}
          else{
            if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<3200){
              CloudCUT2_45_22.5_Pointpaircategory[i,1] = 6}
            else{
              if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<3733.33){
                CloudCUT2_45_22.5_Pointpaircategory[i,1] = 7}
              else{
                if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<4266.67){
                  CloudCUT2_45_22.5_Pointpaircategory[i,1] = 8}
                else{
                  if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<4800){
                    CloudCUT2_45_22.5_Pointpaircategory[i,1] = 9}
                  else{
                    if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<5333.33){
                      CloudCUT2_45_22.5_Pointpaircategory[i,1] = 10}
                    else{
                      if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<5866.67){
                        CloudCUT2_45_22.5_Pointpaircategory[i,1] = 11}
                      else{
                        if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<6400){
                          CloudCUT2_45_22.5_Pointpaircategory[i,1] = 12}
                        else{
                          if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<6933.33){
                            CloudCUT2_45_22.5_Pointpaircategory[i,1] = 13}
                          else{
                            if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<7466.67){
                              CloudCUT2_45_22.5_Pointpaircategory[i,1] = 14}
                            else{
                              CloudCUT2_45_22.5_Pointpaircategory[i,1] = 15
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
CloudCUT2_45_22.5_Pointpaircategory = data.frame(CloudCUT2_45_22.5_Pointpaircategory)
Porosity_VariogramCloudCUT2_45_22.5_classified = data.frame(Porosity_VariogramCloudCUT2_45_22.5,CloudCUT2_45_22.5_Pointpaircategory)
names(Porosity_VariogramCloudCUT2_45_22.5_classified)[8] = "Bin"
Porosity_VariogramCloudCUT2_45_22.5_classified$Bin = as.factor(Porosity_VariogramCloudCUT2_45_22.5_classified$Bin)
ggplot(data = Porosity_VariogramCloudCUT2_45_22.5_classified)+
  geom_point(mapping = aes(x = dist, y = gamma, color = Bin))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", color = "Bin" )+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))
  ggsave("Porosity Variogram Cloud_CUT2_45_22.5_colored_gstat.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

#lag/tolerance-based (gslib)
CloudCUT2_45_22.5_Pointpaircategory_gslib = matrix(0, length(Porosity_VariogramCloudCUT2_45_22.5$dist), 1)
for(i in 1:length(Porosity_VariogramCloudCUT2_45_22.5$dist)){
  if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<500){
    CloudCUT2_45_22.5_Pointpaircategory_gslib[i,1] = 1}
  else{
    if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<1500){
      CloudCUT2_45_22.5_Pointpaircategory_gslib[i,1] = 2}
    else{
      if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<2500){
        CloudCUT2_45_22.5_Pointpaircategory_gslib[i,1] = 3}
      else{
        if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<3500){
          CloudCUT2_45_22.5_Pointpaircategory_gslib[i,1] = 4}
        else{
          if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<4500){
            CloudCUT2_45_22.5_Pointpaircategory_gslib[i,1] = 5}
          else{
            if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<5500){
              CloudCUT2_45_22.5_Pointpaircategory_gslib[i,1] = 6}
            else{
              if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<6500){
                CloudCUT2_45_22.5_Pointpaircategory_gslib[i,1] = 7}
              else{
                if(Porosity_VariogramCloudCUT2_45_22.5$dist[i]<7500){
                  CloudCUT2_45_22.5_Pointpaircategory_gslib[i,1] = 8}
                else{
                  CloudCUT2_45_22.5_Pointpaircategory_gslib[i,1] = 9
                }
              }
            }
          }
        }
      }
    }
  }
}
CloudCUT2_45_22.5_Pointpaircategory_gslib = data.frame(CloudCUT2_45_22.5_Pointpaircategory_gslib)
Porosity_VariogramCloudCUT2_45_22.5_classified_gslib = data.frame(Porosity_VariogramCloudCUT2_45_22.5,CloudCUT2_45_22.5_Pointpaircategory_gslib)
names(Porosity_VariogramCloudCUT2_45_22.5_classified_gslib)[8] = "Lag_Interval"
Porosity_VariogramCloudCUT2_45_22.5_classified_gslib$Lag_Interval = as.factor(Porosity_VariogramCloudCUT2_45_22.5_classified_gslib$Lag_Interval)
ggplot(data = Porosity_VariogramCloudCUT2_45_22.5_classified_gslib)+
  geom_point(mapping = aes(x = dist, y = gamma, color = Lag_Interval))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", color = "Lag Interval" )+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))
ggsave("Porosity Variogram Cloud_CUT2_45_22.5_colored_gslib.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

#A replica of gslib (lag/tolerance-based) experimental variogram
Porosity_VariogramCUT2_45_22.5_gslib = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, alpha = 45, tol.hor = 22.5, boundaries = c(500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500))
ggplot(data = Porosity_VariogramCUT2_45_22.5_gslib)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))
  ggsave("Porosity Emperical Variogram_CUT2_45_22.5_gslib.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

      #Fitting to Model - gslib
  
  lagbasedVariogModel_45 = fit.variogram(Porosity_VariogramCUT2_45_22.5_gslib, vgm("Sph","Exp","Gau"),fit.sills = TRUE, fit.ranges = TRUE)
  lagbasedFitRsquared_45 = function(expvar,modvar){ 
    SSErr = attr(modvar,"SSErr") 
    weight = expvar$np/(expvar$dist^2) #This is the weight factor for the default fit.method  
    SST = sum((weight*(expvar$gamma-mean(expvar$gamma)))^2)
    R_sq<-1-SSErr/SST
    return(R_sq)
  }
  #lagbasedFitRsqd_45 = round(lagbasedFitRsquared_45(Porosity_VariogramCUT2_45_22.5_gslib,lagbasedVariogModel_45),4)
  #dbscanRsqdposter_90 = dbscanExpVariog_90_for_fitting %>% 
  #  summarise(dist = 1000,gamma = 0.0002,dbscanRsqdposterstring_90 = paste("Goodness of Fit, R-squared = ",dbscanFitRsqd_90,sep=""))
  lagbasedVariogModelFunc_45 = function(x){
    ifelse(x <= lagbasedVariogModel_45$range[2], (lagbasedVariogModel_45$psill[1]+(lagbasedVariogModel_45$psill[2]*((1.5*(x/lagbasedVariogModel_45$range[2]))-(0.5*((x/lagbasedVariogModel_45$range[2])^3))))), sum(lagbasedVariogModel_45$psill))
  }
  
  ggplot(data = Porosity_VariogramCUT2_45_22.5_gslib, aes(x = dist, y = gamma))+
    geom_point(aes(shape = "mu"), color = "blue")+
    labs(x = "Lag Distance, h (m)", y = expression(Variogram~gamma))+
    coord_cartesian(xlim = c(0, 8000), ylim = c(0, max(Porosity_VariogramCUT2_45_22.5_gslib$gamma,0.005)))+
    scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
    scale_y_continuous(breaks = seq(0, max(max(Porosity_VariogramCUT2_45_22.5_gslib$gamma)+0.001,0.005), by = 0.001))+
    stat_function(aes(x = dist, color = 'jj'), fun = lagbasedVariogModelFunc_45, xlim = c(0, 8000))+
    #geom_text(aes(label = dbscanRsqdposterstring_90 ), data = dbscanRsqdposter_90, vjust = "bottom", hjust = "left")+
    scale_colour_manual(name = '', values =c('jj'='red'), labels = expression(model:0.0034+0.0008*Sph[2507]*(h)))+
    scale_shape_manual(name = '', values = c('mu' = 16), labels = c('empirical'))+
    theme(legend.position=c(0.6, 0.3))+
    theme(legend.title = element_blank())+
    theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
  ggsave("Porosity Emperical&Fitted Variogram_CUT2_45_22.5_gslib.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis//Model Fits", dpi = 96)
  
  

#90 degrees
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

Porosity_VariogramCUT2_90_22.5 = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, alpha = 90, tol.hor = 22.5)
ggplot(data = Porosity_VariogramCUT2_90_22.5)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))+
  geom_line(color = "red",aes(x = dist, y = gamma))
  ggsave("Porosity Emperical Variogram_CUT2_90_22.5.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


                    #Fitting to Model - bin-based
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
  
  
  # Now considering the natural cluster trends present in the clouds in declaring distance class interval
Porosity_VariogramCUT2_90_22.5_cluster = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, alpha = 90, tol.hor = 22.5, boundaries = c(300, 2000, 3600, 5600, 6800, 8000))
ggplot(data = Porosity_VariogramCUT2_90_22.5_cluster)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.005))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.005, by = 0.001))
  ggsave("Porosity Emperical Variogram_CUT2_90_22.5_cluster.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

# Now using DBSCAN for the cluster analysis
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

          #Fitting to Model - dbscan
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


# Now using SOM for the cluster analysis
porosom_90 = som(Porosity_VariogramCloudCUT2_90_22.5[,c(2,3)], xdim=6, ydim=1, topol="hexa", neigh="gaussian")
somneuronposition_90 = data.frame(porosom_90$code)
somclustereddata_90 = porosom_90$visual
somVariogCloud_90 = data.frame(Porosity_VariogramCloudCUT2_90_22.5$dist, Porosity_VariogramCloudCUT2_90_22.5$gamma, somclustereddata$x)
somVariogCloud_90$somclustereddata.x = as.factor(somVariogCloud_90$somclustereddata.x)

ggplot(data = somVariogCloud_90)+
  geom_point(mapping = aes(x = Porosity_VariogramCloudCUT2_90_22.5.dist, y = Porosity_VariogramCloudCUT2_90_22.5.gamma, color = somclustereddata.x, shape = somclustereddata.x))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "SOM cluster")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))+
  scale_color_manual(name = 'Cluster', values = c("#000000","#00FF00","#330033","#000099","#000033","#FF0000"))+
  scale_shape_manual(name = 'Cluster', values = c(16,8,15,17,18,4))
  ggsave("Porosity Variogram Cloud_CUT2_90_22.5_colored_som.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

  ggplot(data = somVariogCloud_90)+
    geom_point(mapping = aes(x = Porosity_VariogramCloudCUT2_90_22.5.dist, y = Porosity_VariogramCloudCUT2_90_22.5.gamma, color = somclustereddata.x, shape = somclustereddata.x))+
    labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "SOM cluster" )+
    coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
    scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
    scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))+
    scale_color_manual(name = 'Cluster', values = c("#000000","#00FF00","#330033","#000099","#000033","#FF0000"))+
    scale_shape_manual(name = 'Cluster', values = c(16,8,15,17,18,4))+
    geom_vline(xintercept =  somneuronposition_90$X1, color = "red")
    ggsave("Porosity Variogram Cloud_CUT2_90_22.5_colored_som_neuropos.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)
  
NumberofsomClusters_90 = length(levels(somVariogCloud_90$somclustereddata.x))
somExpVariog = matrix(0, NumberofsomClusters_90, 3, dimnames = list(c("1","2","3","4","5","6"),c("NumberofPairs", "AverageDistance", "ExpVariog")))
for(i in 1:NumberofsomClusters_90){
  somExpVariog[i,1]=length(somVariogCloud_90$Porosity_VariogramCloudCUT2_90_22.5.dist[somVariogCloud_90$somclustereddata.x ==i-1]); 
  somExpVariog[i,2]=mean(somVariogCloud_90$Porosity_VariogramCloudCUT2_90_22.5.dist[somVariogCloud_90$somclustereddata.x ==i-1]);
  somExpVariog[i,3]=mean(somVariogCloud_90$Porosity_VariogramCloudCUT2_90_22.5.gamma[somVariogCloud_90$somclustereddata.x ==i-1])
}
somExpVariog = data.frame(somExpVariog)

ggplot(data = somExpVariog)+
  geom_point(mapping = aes(x = AverageDistance, y = ExpVariog), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.005))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.005, by = 0.001))
ggsave("Porosity Emperical Variogram_CUT2_90_22.5_som.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

#For comparison; the colored cloud plot without cluster analysis:
       #bin-based (gstat)
CloudCUT2_90_22.5_Pointpaircategory = matrix(0, length(Porosity_VariogramCloudCUT2_90_22.5$dist), 1)
for(i in 1:length(Porosity_VariogramCloudCUT2_90_22.5$dist)){
  if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<533.33){
    CloudCUT2_90_22.5_Pointpaircategory[i,1] = 1}
  else{
    if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<1066.67){
      CloudCUT2_90_22.5_Pointpaircategory[i,1] = 2}
    else{
      if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<1600){
        CloudCUT2_90_22.5_Pointpaircategory[i,1] = 3}
      else{
        if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<2133.33){
          CloudCUT2_90_22.5_Pointpaircategory[i,1] = 4}
        else{
          if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<2666.67){
            CloudCUT2_90_22.5_Pointpaircategory[i,1] = 5}
          else{
            if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<3200){
              CloudCUT2_90_22.5_Pointpaircategory[i,1] = 6}
            else{
              if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<3733.33){
                CloudCUT2_90_22.5_Pointpaircategory[i,1] = 7}
              else{
                if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<4266.67){
                  CloudCUT2_90_22.5_Pointpaircategory[i,1] = 8}
                else{
                  if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<4800){
                    CloudCUT2_90_22.5_Pointpaircategory[i,1] = 9}
                  else{
                    if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<5333.33){
                      CloudCUT2_90_22.5_Pointpaircategory[i,1] = 10}
                    else{
                      if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<5866.67){
                        CloudCUT2_90_22.5_Pointpaircategory[i,1] = 11}
                      else{
                        if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<6400){
                          CloudCUT2_90_22.5_Pointpaircategory[i,1] = 12}
                        else{
                          if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<6933.33){
                            CloudCUT2_90_22.5_Pointpaircategory[i,1] = 13}
                          else{
                            if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<7466.67){
                              CloudCUT2_90_22.5_Pointpaircategory[i,1] = 14}
                            else{
                              CloudCUT2_90_22.5_Pointpaircategory[i,1] = 15
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
CloudCUT2_90_22.5_Pointpaircategory = data.frame(CloudCUT2_90_22.5_Pointpaircategory)
Porosity_VariogramCloudCUT2_90_22.5_classified = data.frame(Porosity_VariogramCloudCUT2_90_22.5,CloudCUT2_90_22.5_Pointpaircategory)
names(Porosity_VariogramCloudCUT2_90_22.5_classified)[8] = "Bin"
Porosity_VariogramCloudCUT2_90_22.5_classified$Bin = as.factor(Porosity_VariogramCloudCUT2_90_22.5_classified$Bin)
ggplot(data = Porosity_VariogramCloudCUT2_90_22.5_classified)+
    geom_point(mapping = aes(x = dist, y = gamma, color = Bin, shape = Bin))+
    labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Bin" )+
    coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
    scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
    scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))+
    scale_color_manual(name = 'Bin', values = c("#000000","#00FF00","#330033","#000099","#000033","#FF0000", "#99FF00", "#003300", "#330000", "#FFCC00", "#330066", "#FF0066", "#0000FF", "#CC0000"))+
    scale_shape_manual(name = 'Bin', values = c(16,8,15,17,18,4,0,1,2,5,3,6,13,9))
    ggsave("Porosity Variogram Cloud_CUT2_90_22.5_colored_gstat.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

ggplot(data = Porosity_VariogramCloudCUT2_90_22.5_classified)+
    geom_point(mapping = aes(x = dist, y = gamma, color = Bin, shape = Bin))+
    labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Bin" )+
    coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.06))+
    scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
    scale_y_continuous(breaks = seq(0, 0.06, by = 0.01))+
    scale_color_manual(name = 'Bin', values = c("#000000","#00FF00","#330033","#000099","#000033","#FF0000", "#99FF00", "#003300", "#330000", "#FFCC00", "#330066", "#FF0066", "#0000FF", "#CC0000"))+
    scale_shape_manual(name = 'Bin', values = c(16,8,15,17,18,4,0,1,2,5,3,6,13,9))+
    annotate("rect", xmin=2300, xmax=2600, ymin=0, ymax=0.05, alpha=.3,fill="blue")+
    annotate("text", x=2500, y=0.052, label="Sparsity Effect", family="serif",fontface="italic", colour="blue", size=4)+
    annotate("rect", xmin=5700, xmax=7000, ymin=0, ymax=0.05, alpha=.3,fill="blue")+
    annotate("text", x=6500, y=0.052, label="Split Effect", family="serif",fontface="italic", colour="blue", size=4)
ggsave("Porosity Variogram Cloud_CUT2_90_22.5_colored_gstat_annotated.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)
    
         #lag/tolerance-based (gslib)
CloudCUT2_90_22.5_Pointpaircategory_gslib = matrix(0, length(Porosity_VariogramCloudCUT2_90_22.5$dist), 1)
for(i in 1:length(Porosity_VariogramCloudCUT2_90_22.5$dist)){
  if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<500){
    CloudCUT2_90_22.5_Pointpaircategory_gslib[i,1] = 1}
  else{
    if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<1500){
      CloudCUT2_90_22.5_Pointpaircategory_gslib[i,1] = 2}
    else{
      if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<2500){
        CloudCUT2_90_22.5_Pointpaircategory_gslib[i,1] = 3}
      else{
        if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<3500){
          CloudCUT2_90_22.5_Pointpaircategory_gslib[i,1] = 4}
        else{
          if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<4500){
            CloudCUT2_90_22.5_Pointpaircategory_gslib[i,1] = 5}
          else{
            if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<5500){
              CloudCUT2_90_22.5_Pointpaircategory_gslib[i,1] = 6}
            else{
              if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<6500){
                CloudCUT2_90_22.5_Pointpaircategory_gslib[i,1] = 7}
              else{
                if(Porosity_VariogramCloudCUT2_90_22.5$dist[i]<7500){
                  CloudCUT2_90_22.5_Pointpaircategory_gslib[i,1] = 8}
                else{
                  CloudCUT2_90_22.5_Pointpaircategory_gslib[i,1] = 9
                }
              }
            }
          }
        }
      }
    }
  }
}
CloudCUT2_90_22.5_Pointpaircategory_gslib = data.frame(CloudCUT2_90_22.5_Pointpaircategory_gslib)
Porosity_VariogramCloudCUT2_90_22.5_classified_gslib = data.frame(Porosity_VariogramCloudCUT2_90_22.5,CloudCUT2_90_22.5_Pointpaircategory_gslib)
names(Porosity_VariogramCloudCUT2_90_22.5_classified_gslib)[8] = "Lag_Interval"
Porosity_VariogramCloudCUT2_90_22.5_classified_gslib$Lag_Interval = as.factor(Porosity_VariogramCloudCUT2_90_22.5_classified_gslib$Lag_Interval)
ggplot(data = Porosity_VariogramCloudCUT2_90_22.5_classified_gslib)+
  geom_point(mapping = aes(x = dist, y = gamma, color = Lag_Interval, shape = Lag_Interval))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Lag Interval" )+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))+
  scale_color_manual(name = 'Lag Interval', values = c("#000000","#FF0066","#330033","#000099","#000033","#FF0000", "#00FF00", "#003300", "#330000"))+
  scale_shape_manual(name = 'Lag Interval', values = c(16,8,15,17,3,4,0,1,2))
  ggsave("Porosity Variogram Cloud_CUT2_90_22.5_colored_gslib.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

ggplot(data = Porosity_VariogramCloudCUT2_90_22.5_classified_gslib)+
    geom_point(mapping = aes(x = dist, y = gamma, color = Lag_Interval, shape = Lag_Interval))+
    labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Lag Interval" )+
    coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.06))+
    scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
    scale_y_continuous(breaks = seq(0, 0.06, by = 0.01))+
    scale_color_manual(name = 'Lag Interval', values = c("#000000","#FF0066","#330033","#000099","#000033","#FF0000", "#00FF00", "#003300", "#330000"))+
    scale_shape_manual(name = 'Lag Interval', values = c(16,8,15,17,3,4,0,1,2))+
    annotate("rect", xmin=6500, xmax=7400, ymin=0, ymax=0.05, alpha=.3,fill="blue")+
    annotate("text", x=7000, y=0.055, label="Split-Straddle \nEffect", family="serif",fontface="italic", colour="blue", size=4)
ggsave("Porosity Variogram Cloud_CUT2_90_22.5_colored_gslib_annotated.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


#A replica of gslib (lag/tolerance-based) experimental variogram
Porosity_VariogramCUT2_90_22.5_gslib = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, alpha = 90, tol.hor = 22.5, boundaries = c(500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500))
ggplot(data = Porosity_VariogramCUT2_90_22.5_gslib)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.005))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.005, by = 0.001))
  ggsave("Porosity Emperical Variogram_CUT2_90_22.5_gslib.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

              #Fitting to Model - bin-based - visual (manual) fitting
  lagbasedVariogModel_90 = fit.variogram(Porosity_VariogramCUT2_90_22.5_gslib, vgm("Sph","Exp","Gau"),fit.sills = TRUE, fit.ranges = TRUE)
  lagbasedVariogModel_90a = fit.variogram(Porosity_VariogramCUT2_90_22.5_gslib, vgm(0.0012, "Sph", 2000, 0.0033),fit.sills = FALSE, fit.ranges = FALSE)
  
  lagbasedFitRsquared_90 = function(expvar,modvar){ 
    SSErr = attr(modvar,"SSErr") 
    weight = expvar$np/(expvar$dist^2) #This is the weight factor for the default fit.method  
    SST = sum((weight*(expvar$gamma-mean(expvar$gamma)))^2)
    R_sq = 1-SSErr/SST
    return(R_sq)
  }
  lagbasedFitRsqd_90 = round(lagbasedFitRsquared_90(Porosity_VariogramCUT2_90_22.5_gslib,lagbasedVariogModel_90),4)
  lagbasedRsqdposter_90 = Porosity_VariogramCUT2_90_22.5_gslib %>% 
    summarise(dist = 1000,gamma = 0.0002,lagbasedRsqdposterstring_90 = paste("Goodness of Fit, R-squared = ",lagbasedFitRsqd_90,sep=""))
  lagbasedVariogModelFunc_90 = function(x){
    ifelse(x <= lagbasedVariogModel_90a$range[2], (lagbasedVariogModel_90a$psill[1]+(lagbasedVariogModel_90a$psill[2]*((1.5*(x/lagbasedVariogModel_90a$range[2]))-(0.5*((x/lagbasedVariogModel_90a$range[2])^3))))), sum(lagbasedVariogModel_90a$psill))
  }
  
  ggplot(data = Porosity_VariogramCUT2_90_22.5_gslib, aes(x = dist, y = gamma))+
    geom_point(aes(shape = "mu"), color = "blue")+
    labs(x = "Lag Distance, h (m)", y = expression(Variogram~gamma))+
    coord_cartesian(xlim = c(0, 8000), ylim = c(0, max(Porosity_VariogramCUT2_90_22.5_gslib$gamma,0.005)))+
    scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
    scale_y_continuous(breaks = seq(0, max(max(Porosity_VariogramCUT2_90_22.5_gslib$gamma)+0.001,0.005), by = 0.001))+
    stat_function(aes(x = dist, color = 'jj'), fun = lagbasedVariogModelFunc_90, xlim = c(0, 8000))+
    geom_text(aes(label = lagbasedRsqdposterstring_90 ), data = lagbasedRsqdposter_90, vjust = "bottom", hjust = "left")+
    scale_colour_manual(name = '', values =c('jj'='red'), labels = expression(model:0.0034+0.0013*Sph[5848]*(h)))+
    scale_shape_manual(name = '', values = c('mu' = 16), labels = c('empirical'))+
    theme(legend.position=c(0.6, 0.3))+
    theme(legend.title = element_blank())+
    theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
  ggsave("Porosity Emperical&Fitted Variogram_CUT2_90_22.5_lagbased.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis//Model Fits", dpi = 96)
  
  
  
#135 degree
Porosity_VariogramCloudCUT2_135_22.5 = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, alpha = 135, tol.hor = 22.5, cloud = TRUE)
ggplot(data = Porosity_VariogramCloudCUT2_135_22.5)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram")+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))
  ggsave("Porosity Variogram Cloud_CUT2_135_22.5.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

write.table(Porosity_VariogramCloudCUT2_135_22.5,"Porosity_VariogramCloudCUT2_135_22.5.csv", sep = ",", row.names = F)

Porosity_VariogramCUT2_135_22.5 = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, alpha = 135, tol.hor = 22.5)
ggplot(data = Porosity_VariogramCUT2_135_22.5)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.008))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.008, by = 0.001))
  ggsave("Porosity Emperical Variogram_CUT2_135_22.5.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


# Now considering the natural cluster trends present in the clouds in declaring distance class interval
Porosity_VariogramCUT2_135_22.5_cluster = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, alpha = 135, tol.hor = 22.5, boundaries = c(100, 2200, 3000, 4100, 5500, 7000))
ggplot(data = Porosity_VariogramCUT2_135_22.5_cluster)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.008))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.008, by = 0.001))
  ggsave("Porosity Emperical Variogram_CUT2_135_22.5_cluster.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

# Now using DBSCAN for the cluster analysis
kNNdistplot(Porosity_VariogramCloudCUT2_135_22.5[,c(2,3)],k = 500)
abline(h=220, col = "red", lty=2)
# This yields eps = 180 and minPts = 500
porodbscan_135 = dbscan(Porosity_VariogramCloudCUT2_135_22.5[,c(2,3)], eps = 220, minPts = 500)
dbscanVariogCloud_135 = data.frame(Porosity_VariogramCloudCUT2_135_22.5$dist, Porosity_VariogramCloudCUT2_135_22.5$gamma, porodbscan_135$cluster, Direction = 135)    
dbscanVariogCloud_135$porodbscan_135.cluster = as.factor(dbscanVariogCloud_135$porodbscan_135.cluster) #Ensuring R sees the cluster numbers as classification ID and not as numerals
ggplot(data = dbscanVariogCloud_135[!dbscanVariogCloud_135$porodbscan_135.cluster=="0",])+
  geom_point(mapping = aes(x = Porosity_VariogramCloudCUT2_135_22.5.dist, y = Porosity_VariogramCloudCUT2_135_22.5.gamma, color = porodbscan_135.cluster))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Bin" )+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))
  ggsave("Porosity Variogram Cloud_CUT2_135_22.5_colored_dbsacan.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

NumberofClusters_135 = length(levels(dbscanVariogCloud_135$porodbscan_135.cluster))
dbscanExpVariog_135 = matrix(0, NumberofClusters_135-1, 3, dimnames = list(c("1","2","3","4","5","6"),c("NumberofPairs", "AverageDistance", "ExpVariog")))
for(i in 1:NumberofClusters_135-1){
  dbscanExpVariog_135[i,1]=length(dbscanVariogCloud_135$Porosity_VariogramCloudCUT2_135_22.5.dist[dbscanVariogCloud_135$porodbscan_135.cluster==i]); 
  dbscanExpVariog_135[i,2]=mean(dbscanVariogCloud_135$Porosity_VariogramCloudCUT2_135_22.5.dist[dbscanVariogCloud_135$porodbscan_135.cluster==i]);
  dbscanExpVariog_135[i,3]=mean(dbscanVariogCloud_135$Porosity_VariogramCloudCUT2_135_22.5.gamma[dbscanVariogCloud_135$porodbscan_135.cluster==i])
}
dbscanExpVariog_135 = data.frame(dbscanExpVariog_135)
ggplot(data = dbscanExpVariog_135)+
  geom_point(mapping = aes(x = AverageDistance, y = ExpVariog), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.008))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.008, by = 0.001))
  ggsave("Porosity Emperical Variogram_CUT2_135_22.5_dbscan.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

#Fitting to Model - dbscan  
  dbscanExpVariog_135_for_fitting = dbscanExpVariog_135[1:6,]
  names(dbscanExpVariog_135_for_fitting) = c("np","dist","gamma")
  class(dbscanExpVariog_135_for_fitting) = c("gstatVariogram", "data.frame")
  dbscanVariogModel_135 = fit.variogram(dbscanExpVariog_135_for_fitting, vgm("Sph",psill = 0.0021, nugget = 0.0024, range = 3000),fit.sills = FALSE, fit.ranges = FALSE)
  
  dbscanVariogModelFunc_135 = function(x){
    ifelse(x <= dbscanVariogModel_135$range[2], (dbscanVariogModel_135$psill[1]+(dbscanVariogModel_135$psill[2]*((1.5*(x/dbscanVariogModel_135$range[2]))-(0.5*((x/dbscanVariogModel_135$range[2])^3))))), sum(dbscanVariogModel_135$psill))
  }
  ggplot(data = dbscanExpVariog_135_for_fitting, aes(x = dist, y = gamma))+
    geom_point(aes(shape = "mu"), color = "blue")+
    labs(x = "Lag Distance, h (m)", y = expression(Variogram~gamma))+
    coord_cartesian(xlim = c(0, 8000), ylim = c(0, max(dbscanExpVariog_135_for_fitting$gamma,0.005)))+
    scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
    scale_y_continuous(breaks = seq(0, max(max(dbscanExpVariog_135_for_fitting$gamma)+0.001,0.005), by = 0.001))+
    stat_function(aes(x = dist, color = 'jj'), fun = dbscanVariogModelFunc_135, xlim = c(0, 8000))+
    scale_colour_manual(name = '', values =c('jj'='red'), labels = expression(model:0.0024+0.0021*Sph[3000]*(h)))+
    scale_shape_manual(name = '', values = c('mu' = 16), labels = c('empirical'))+
    theme(legend.position=c(0.6, 0.3))+
    theme(legend.title = element_blank())+
    theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
  ggsave("Porosity Emperical&Fitted Variogram_CUT2_135_22.5_dbscan.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis//Model Fits", dpi = 96)
  
  
#For comparison; the colored cloud plot without cluster analysis:
       #bin-based (gstat)
CloudCUT2_135_22.5_Pointpaircategory = matrix(0, length(Porosity_VariogramCloudCUT2_135_22.5$dist), 1)
for(i in 1:length(Porosity_VariogramCloudCUT2_135_22.5$dist)){
  if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<533.33){
    CloudCUT2_135_22.5_Pointpaircategory[i,1] = 1}
  else{
    if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<1066.67){
      CloudCUT2_135_22.5_Pointpaircategory[i,1] = 2}
    else{
      if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<1600){
        CloudCUT2_135_22.5_Pointpaircategory[i,1] = 3}
      else{
        if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<2133.33){
          CloudCUT2_135_22.5_Pointpaircategory[i,1] = 4}
        else{
          if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<2666.67){
            CloudCUT2_135_22.5_Pointpaircategory[i,1] = 5}
          else{
            if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<3200){
              CloudCUT2_135_22.5_Pointpaircategory[i,1] = 6}
            else{
              if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<3733.33){
                CloudCUT2_135_22.5_Pointpaircategory[i,1] = 7}
              else{
                if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<4266.67){
                  CloudCUT2_135_22.5_Pointpaircategory[i,1] = 8}
                else{
                  if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<4800){
                    CloudCUT2_135_22.5_Pointpaircategory[i,1] = 9}
                  else{
                    if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<5333.33){
                      CloudCUT2_135_22.5_Pointpaircategory[i,1] = 10}
                    else{
                      if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<5866.67){
                        CloudCUT2_135_22.5_Pointpaircategory[i,1] = 11}
                      else{
                        if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<6400){
                          CloudCUT2_135_22.5_Pointpaircategory[i,1] = 12}
                        else{
                          if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<6933.33){
                            CloudCUT2_135_22.5_Pointpaircategory[i,1] = 13}
                          else{
                            if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<7466.67){
                              CloudCUT2_135_22.5_Pointpaircategory[i,1] = 14}
                            else{
                              CloudCUT2_135_22.5_Pointpaircategory[i,1] = 15
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
CloudCUT2_135_22.5_Pointpaircategory = data.frame(CloudCUT2_135_22.5_Pointpaircategory)
Porosity_VariogramCloudCUT2_135_22.5_classified = data.frame(Porosity_VariogramCloudCUT2_135_22.5,CloudCUT2_135_22.5_Pointpaircategory)
names(Porosity_VariogramCloudCUT2_135_22.5_classified)[8] = "Bin"
Porosity_VariogramCloudCUT2_135_22.5_classified$Bin = as.factor(Porosity_VariogramCloudCUT2_135_22.5_classified$Bin)
ggplot(data = Porosity_VariogramCloudCUT2_135_22.5_classified)+
  geom_point(mapping = aes(x = dist, y = gamma, color = Bin))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Bin" )+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))
  ggsave("Porosity Variogram Cloud_CUT2_135_22.5_colored_gstat.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)
               
         #lag/tolerance-based (gslib)
CloudCUT2_135_22.5_Pointpaircategory_gslib = matrix(0, length(Porosity_VariogramCloudCUT2_135_22.5$dist), 1)
for(i in 1:length(Porosity_VariogramCloudCUT2_135_22.5$dist)){
  if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<500){
    CloudCUT2_135_22.5_Pointpaircategory_gslib[i,1] = 1}
  else{
    if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<1500){
      CloudCUT2_135_22.5_Pointpaircategory_gslib[i,1] = 2}
    else{
      if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<2500){
        CloudCUT2_135_22.5_Pointpaircategory_gslib[i,1] = 3}
      else{
        if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<3500){
          CloudCUT2_135_22.5_Pointpaircategory_gslib[i,1] = 4}
        else{
          if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<4500){
            CloudCUT2_135_22.5_Pointpaircategory_gslib[i,1] = 5}
          else{
            if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<5500){
              CloudCUT2_135_22.5_Pointpaircategory_gslib[i,1] = 6}
            else{
              if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<6500){
                CloudCUT2_135_22.5_Pointpaircategory_gslib[i,1] = 7}
              else{
                if(Porosity_VariogramCloudCUT2_135_22.5$dist[i]<7500){
                  CloudCUT2_135_22.5_Pointpaircategory_gslib[i,1] = 8}
                else{
                  CloudCUT2_135_22.5_Pointpaircategory_gslib[i,1] = 9
                }
              }
            }
          }
        }
      }
    }
  }
}
CloudCUT2_135_22.5_Pointpaircategory_gslib = data.frame(CloudCUT2_135_22.5_Pointpaircategory_gslib)
Porosity_VariogramCloudCUT2_135_22.5_classified_gslib = data.frame(Porosity_VariogramCloudCUT2_135_22.5,CloudCUT2_135_22.5_Pointpaircategory_gslib)
names(Porosity_VariogramCloudCUT2_135_22.5_classified_gslib)[8] = "Bin"
Porosity_VariogramCloudCUT2_135_22.5_classified_gslib$Bin = as.factor(Porosity_VariogramCloudCUT2_135_22.5_classified_gslib$Bin)
ggplot(data = Porosity_VariogramCloudCUT2_135_22.5_classified_gslib)+
  geom_point(mapping = aes(x = dist, y = gamma, color = Bin))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Bin" )+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))
  ggsave("Porosity Variogram Cloud_CUT2_135_22.5_colored_gslib.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

#A replica of gslib (lag/tolerance-based) experimental variogram
Porosity_VariogramCUT2_135_22.5_gslib = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 8000, alpha = 135, tol.hor = 22.5, boundaries = c(500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500))
ggplot(data = Porosity_VariogramCUT2_135_22.5_gslib)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.008))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.008, by = 0.001))
  ggsave("Porosity Emperical Variogram_CUT2_135_22.5_gslib.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

  
#180 degree. Note the cutoff has to be 2500 because the domain is of small dimension in this direction
Porosity_VariogramCloudCUT2_180_22.5 = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 2500, alpha = 180, tol.hor = 22.5, cloud = TRUE)
ggplot(data = Porosity_VariogramCloudCUT2_180_22.5)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram")+
  scale_x_continuous(breaks = seq(0, 2500, by = 500))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))
  ggsave("Porosity Variogram Cloud_CUT2_180_22.5.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

write.table(Porosity_VariogramCloudCUT2_180_22.5,"Porosity_VariogramCloudCUT2_180_22.5.csv", sep = ",", row.names = F)

Porosity_VariogramCUT2_180_22.5 = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 2500, alpha = 180, tol.hor = 22.5)
ggplot(data = Porosity_VariogramCUT2_180_22.5)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 2500), ylim = c(0, 0.005))+
  scale_x_continuous(breaks = seq(0, 2500, by = 500))+
  scale_y_continuous(breaks = seq(0, 0.005, by = 0.001))
  ggsave("Porosity Emperical Variogram_CUT2_180_22.5.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


# Now considering the natural cluster trends present in the clouds in declaring distance class interval
Porosity_VariogramCUT2_180_22.5_cluster = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 2500, alpha = 180, tol.hor = 22.5, boundaries = c(100, 600, 1200, 1700, 2000, 2500))
ggplot(data = Porosity_VariogramCUT2_180_22.5_cluster)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 2500), ylim = c(0, 0.005))+
  scale_x_continuous(breaks = seq(0, 2500, by = 500))+
  scale_y_continuous(breaks = seq(0, 0.005, by = 0.001))
  ggsave("Porosity Emperical Variogram_CUT2_180_22.5_cluster.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


# Now using DBSCAN for the cluster analysis
kNNdistplot(Porosity_VariogramCloudCUT2_180_22.5[,c(2,3)],k = 500)
abline(h=70, col = "red", lty=2)
# This yields eps = 70 and minPts = 400
porodbscan_180 = dbscan(Porosity_VariogramCloudCUT2_180_22.5[,c(2,3)], eps = 70, minPts = 400)
dbscanVariogCloud_180 = data.frame(Porosity_VariogramCloudCUT2_180_22.5$dist, Porosity_VariogramCloudCUT2_180_22.5$gamma, porodbscan_180$cluster, Direction = 180)    
dbscanVariogCloud_180$porodbscan_180.cluster = as.factor(dbscanVariogCloud_180$porodbscan_180.cluster) #Ensuring R sees the cluster numbers as classification ID and not as numerals
ggplot(data = dbscanVariogCloud_180[!dbscanVariogCloud_180$porodbscan_180.cluster=="0",])+
  geom_point(mapping = aes(x = Porosity_VariogramCloudCUT2_180_22.5.dist, y = Porosity_VariogramCloudCUT2_180_22.5.gamma, color = porodbscan_180.cluster))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Bin" )+
  coord_cartesian(xlim = c(0, 2500), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 2500, by = 500))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))
  ggsave("Porosity Variogram Cloud_CUT2_180_22.5_colored_dbsacan.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


NumberofClusters_180 = length(levels(dbscanVariogCloud_180$porodbscan_180.cluster))
dbscanExpVariog_180 = matrix(0, NumberofClusters_180-1, 3, dimnames = list(c("1","2","3","4","5"),c("NumberofPairs", "AverageDistance", "ExpVariog")))
for(i in 1:NumberofClusters_180-1){
  dbscanExpVariog_180[i,1]=length(dbscanVariogCloud_180$Porosity_VariogramCloudCUT2_180_22.5.dist[dbscanVariogCloud_180$porodbscan_180.cluster==i]); 
  dbscanExpVariog_180[i,2]=mean(dbscanVariogCloud_180$Porosity_VariogramCloudCUT2_180_22.5.dist[dbscanVariogCloud_180$porodbscan_180.cluster==i]);
  dbscanExpVariog_180[i,3]=mean(dbscanVariogCloud_180$Porosity_VariogramCloudCUT2_180_22.5.gamma[dbscanVariogCloud_180$porodbscan_180.cluster==i])
}
dbscanExpVariog_180 = data.frame(dbscanExpVariog_180)
ggplot(data = dbscanExpVariog_180)+
  geom_point(mapping = aes(x = AverageDistance, y = ExpVariog), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 2500), ylim = c(0, 0.005))+
  scale_x_continuous(breaks = seq(0, 2500, by = 500))+
  scale_y_continuous(breaks = seq(0, 0.005, by = 0.001))
  ggsave("Porosity Emperical Variogram_CUT2_180_22.5_dbscan.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

  
# Fitting to Model - dbscan
  dbscanExpVariog_180_for_fitting = dbscanExpVariog_180
  names(dbscanExpVariog_180_for_fitting) = c("np","dist","gamma")
  class(dbscanExpVariog_180_for_fitting) = c("gstatVariogram", "data.frame")
  dbscanVariogModel_180 = fit.variogram(dbscanExpVariog_180_for_fitting, vgm("Sph",psill = 0.0021, nugget = 0.0024, range = 1500),fit.sills = FALSE, fit.ranges = FALSE)

  
  dbscanVariogModelFunc_180 = function(x){
    ifelse(x <= dbscanVariogModel_180$range[2], (dbscanVariogModel_180$psill[1]+(dbscanVariogModel_180$psill[2]*((1.5*(x/dbscanVariogModel_180$range[2]))-(0.5*((x/dbscanVariogModel_180$range[2])^3))))), sum(dbscanVariogModel_180$psill))
  }
  ggplot(data = dbscanExpVariog_180_for_fitting, aes(x = dist, y = gamma))+
    geom_point(aes(shape = "mu"), color = "blue")+
    labs(x = "Lag Distance, h (m)", y = expression(Variogram~gamma))+
    coord_cartesian(xlim = c(0, 8000), ylim = c(0, max(dbscanExpVariog_180_for_fitting$gamma,0.005)))+
    scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
    scale_y_continuous(breaks = seq(0, max(max(dbscanExpVariog_180_for_fitting$gamma)+0.001,0.005), by = 0.001))+
    stat_function(aes(x = dist, color = 'jj'), fun = dbscanVariogModelFunc_180, xlim = c(0, 3000))+
    scale_colour_manual(name = '', values =c('jj'='red'), labels = expression(model:0.0024+0.0021*Sph[1500]*(h)))+
    scale_shape_manual(name = '', values = c('mu' = 16), labels = c('empirical'))+
    theme(legend.position=c(0.6, 0.3))+
    theme(legend.title = element_blank())+
    theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
  ggsave("Porosity Emperical&Fitted Variogram_CUT2_180_22.5_dbscan.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis//Model Fits", dpi = 96)
  
  
  

#For comparison; the colored cloud plot without cluster analysis:
      #bin-based (gstat)
CloudCUT2_180_22.5_Pointpaircategory = matrix(0, length(Porosity_VariogramCloudCUT2_180_22.5$dist), 1)
for(i in 1:length(Porosity_VariogramCloudCUT2_180_22.5$dist)){
  if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<166.66){
    CloudCUT2_180_22.5_Pointpaircategory[i,1] = 1}
  else{
    if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<333.33){
      CloudCUT2_180_22.5_Pointpaircategory[i,1] = 2}
    else{
      if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<500){
        CloudCUT2_180_22.5_Pointpaircategory[i,1] = 3}
      else{
        if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<666.67){
          CloudCUT2_180_22.5_Pointpaircategory[i,1] = 4}
        else{
          if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<833.33){
            CloudCUT2_180_22.5_Pointpaircategory[i,1] = 5}
          else{
            if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<1000){
              CloudCUT2_180_22.5_Pointpaircategory[i,1] = 6}
            else{
              if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<1166.67){
                CloudCUT2_180_22.5_Pointpaircategory[i,1] = 7}
              else{
                if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<1333.33){
                  CloudCUT2_180_22.5_Pointpaircategory[i,1] = 8}
                else{
                  if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<1500){
                    CloudCUT2_180_22.5_Pointpaircategory[i,1] = 9}
                  else{
                    if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<1666.67){
                      CloudCUT2_180_22.5_Pointpaircategory[i,1] = 10}
                    else{
                      if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<1833.33){
                        CloudCUT2_180_22.5_Pointpaircategory[i,1] = 11}
                      else{
                        if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<2000){
                          CloudCUT2_180_22.5_Pointpaircategory[i,1] = 12}
                        else{
                          if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<2166.67){
                            CloudCUT2_180_22.5_Pointpaircategory[i,1] = 13}
                          else{
                            if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<2333.33){
                              CloudCUT2_180_22.5_Pointpaircategory[i,1] = 14}
                            else{
                              CloudCUT2_180_22.5_Pointpaircategory[i,1] = 15
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
CloudCUT2_180_22.5_Pointpaircategory = data.frame(CloudCUT2_180_22.5_Pointpaircategory)
Porosity_VariogramCloudCUT2_180_22.5_classified = data.frame(Porosity_VariogramCloudCUT2_180_22.5,CloudCUT2_180_22.5_Pointpaircategory)
names(Porosity_VariogramCloudCUT2_180_22.5_classified)[8] = "Bin"
Porosity_VariogramCloudCUT2_180_22.5_classified$Bin = as.factor(Porosity_VariogramCloudCUT2_180_22.5_classified$Bin)
ggplot(data = Porosity_VariogramCloudCUT2_180_22.5_classified)+
  geom_point(mapping = aes(x = dist, y = gamma, color = Bin))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Bin" )+
  coord_cartesian(xlim = c(0, 2500), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 2500, by = 500))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))
  ggsave("Porosity Variogram Cloud_CUT2_180_22.5_colored_gstat.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

#lag/tolerance-based (gslib)
CloudCUT2_180_22.5_Pointpaircategory_gslib = matrix(0, length(Porosity_VariogramCloudCUT2_180_22.5$dist), 1)
for(i in 1:length(Porosity_VariogramCloudCUT2_180_22.5$dist)){
  if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<250){
    CloudCUT2_180_22.5_Pointpaircategory_gslib[i,1] = 1}
  else{
    if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<750){
      CloudCUT2_180_22.5_Pointpaircategory_gslib[i,1] = 2}
    else{
      if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<1250){
        CloudCUT2_180_22.5_Pointpaircategory_gslib[i,1] = 3}
      else{
        if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<1750){
          CloudCUT2_180_22.5_Pointpaircategory_gslib[i,1] = 4}
        else{
          if(Porosity_VariogramCloudCUT2_180_22.5$dist[i]<2250){
            CloudCUT2_180_22.5_Pointpaircategory_gslib[i,1] = 5}
          else{
            CloudCUT2_180_22.5_Pointpaircategory_gslib[i,1] = 6
          }
        }
      }
    }
  }
}
CloudCUT2_180_22.5_Pointpaircategory_gslib = data.frame(CloudCUT2_180_22.5_Pointpaircategory_gslib)
Porosity_VariogramCloudCUT2_180_22.5_classified_gslib = data.frame(Porosity_VariogramCloudCUT2_180_22.5,CloudCUT2_180_22.5_Pointpaircategory_gslib)
names(Porosity_VariogramCloudCUT2_180_22.5_classified_gslib)[8] = "Bin"
Porosity_VariogramCloudCUT2_180_22.5_classified_gslib$Bin = as.factor(Porosity_VariogramCloudCUT2_180_22.5_classified_gslib$Bin)
ggplot(data = Porosity_VariogramCloudCUT2_180_22.5_classified_gslib)+
  geom_point(mapping = aes(x = dist, y = gamma, color = Bin))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Bin" )+
  coord_cartesian(xlim = c(0, 2500), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 2500, by = 500))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))
  ggsave("Porosity Variogram Cloud_CUT2_180_22.5_colored__gslib.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)

#A replica of gslib (lag/tolerance-based) experimental variogram
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Porosity Variogram_CUT2_135_22.5_gslib.jpg")


Porosity_VariogramCUT2_180_22.5_gslib = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 2500, alpha = 180, tol.hor = 22.5, boundaries = c(250, 750, 1250, 1750, 2250, 2750))
ggplot(data = Porosity_VariogramCUT2_180_22.5_gslib)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 2500), ylim = c(0, 0.005))+
  scale_x_continuous(breaks = seq(0, 2500, by = 500))+
  scale_y_continuous(breaks = seq(0, 0.005, by = 0.001))
  ggsave("Porosity Emperical Variogram_CUT2_180_22.5_gslib.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)



#Agbabu Bitumen Deposit Well-mean Porosity Data Analysis

#Importing data
WellMeanPoroData = read.csv(file.choose(), header = F)
names(WellMeanPoroData) = c("x", "y", "Porosity", "BoreHoleID")
WellMeanPoroData$BoreHoleID = as.factor(WellMeanPoroData$BoreHoleID) #Ensuring R sees the Borehole numbers as ID and not as numerals

library(sp)
library(gstat)

coordinates(WellMeanPoroData) = c("x","y")

jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Well-Mean Analysis//Well-mean Porosity Lagged Scatter Plot_whole.jpg")
hscat(Porosity~1,WellMeanPoroData,(0:4)*1000)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Well-Mean Analysis//Well-mean Porosity Variogram Map_whole.jpg")
WellmeanPorosity_Variogram_mapWhole = variogram(Porosity~1,WellMeanPoroData, map = TRUE,cutoff=10000,width=1000)
plot(WellmeanPorosity_Variogram_mapWhole)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Well-Mean Analysis//Well-mean Porosity Variogram Cloud_whole.jpg")
WellmeanPorosity_VariogramCloudWhole = variogram(Porosity~1,WellMeanPoroData[WellMeanPoroData$Porosity<0.36, ], cloud = TRUE)
WellmeanPoroPointPairs=plot(WellmeanPorosity_VariogramCloudWhole, identify = TRUE, digitize = TRUE, xlim = c(0,6000), ylim = c(0,0.008))
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Well-Mean Analysis//Well-mean Porosity Variogram CloudPoinPairs_whole.jpg")
plot(WellmeanPoroPointPairs, WellMeanPoroData, xlab ="Easting", ylab ="Northing")
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Well-Mean Analysis//Well-mean Porosity Emperical Variogram_whole.jpg")
WellmeanPorosity_VariogramWhole = variogram(Porosity~1,WellMeanPoroData[WellMeanPoroData$Porosity<0.36, ])
plot(WellmeanPorosity_VariogramWhole)
dev.off()












#Agbabu Bitumen Deposit Individual Porosity Data Analysis - Horizon X

#No CUT
SpatialTrimmedDataX = SmartTrimmedData_X
names(SpatialTrimmedDataX) = c("BoreHoleID", "Horizon", "Depths", "Porosity", "x", "y")

library(sp)
library(gstat)

coordinates(SpatialTrimmedDataX) = c("x","y","Depths")

jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon X Analysis//Porosity Lagged Scatter Plot_X.jpg")
hscat(Porosity~1,SpatialTrimmedDataX,(0:6)*1000)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon X Analysis//Porosity Variogram Map_X.jpg")
Porosity_Variogram_mapX = variogram(Porosity~1,SpatialTrimmedDataX, map = TRUE,cutoff=10000,width=1000)
plot(Porosity_Variogram_mapX)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon X Analysis//Porosity Variogram Cloud_X.jpg")
Porosity_VariogramCloudX = variogram(Porosity~1,SpatialTrimmedDataX, cloud = TRUE)
PoroPointPairsX=plot(Porosity_VariogramCloudX, identify = TRUE, digitize = TRUE, xlim = c(0,6000), ylim = c(0,0.14))
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon X Analysis//Porosity Variogram CloudPoinPairs_X.jpg")
plot(PoroPointPairsX, SpatialTrimmedDataX)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon X Analysis//Porosity Emperical Variogram_X.jpg")
Porosity_VariogramX = variogram(Porosity~1,SpatialTrimmedDataX)
plot(Porosity_VariogramX)
dev.off()



#CUT2 Data points between Q1-IQR (0.1239) and Q3+IQR (0.3780)
SmartTrimmedDataXCUT2 = SmartTrimmedData_X
names(SmartTrimmedDataXCUT2) = c("BoreHoleID", "Horizon", "Depths", "Porosity", "x", "y")
SpatialTrimmedDataXCUT2 = SmartTrimmedDataXCUT2[SmartTrimmedDataXCUT2$Porosity<0.3780&SmartTrimmedDataXCUT2$Porosity>0.1239, ]


library(sp)
library(gstat)

coordinates(SpatialTrimmedDataXCUT2) = c("x","y","Depths")

jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon X Analysis//Porosity Lagged Scatter Plot_XCUT2.jpg")
hscat(Porosity~1,SpatialTrimmedDataXCUT2,(0:6)*1000)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon X Analysis//Porosity Variogram Map_XCUT2.jpg")
Porosity_Variogram_mapXCUT2 = variogram(Porosity~1,SpatialTrimmedDataXCUT2, map = TRUE,cutoff=10000,width=1000)
plot(Porosity_Variogram_mapXCUT2)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon X Analysis//Porosity Variogram Cloud_XCUT2.jpg")
Porosity_VariogramCloudXCUT2 = variogram(Porosity~1,SpatialTrimmedDataXCUT2, cloud = TRUE)
PoroPointPairsXCUT2=plot(Porosity_VariogramCloudXCUT2, identify = TRUE, digitize = TRUE, xlim = c(0,6000), ylim = c(0,0.04))
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon X Analysis//Porosity Variogram CloudPoinPairs_XCUT2.jpg")
plot(PoroPointPairsXCUT2, SpatialTrimmedDataXCUT2)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon X Analysis//Porosity Emperical Variogram_XCUT2.jpg")
Porosity_VariogramXCUT2 = variogram(Porosity~1,SpatialTrimmedDataXCUT2)
plot(Porosity_VariogramXCUT2)
dev.off()












#Agbabu Bitumen Deposit Individual Porosity Data Analysis - Horizon Y

#No CUT
SpatialTrimmedDataY = SmartTrimmedData_Y
names(SpatialTrimmedDataY) = c("BoreHoleID", "Horizon", "Depths", "Porosity", "x", "y")

library(sp)
library(gstat)

coordinates(SpatialTrimmedDataY) = c("x","y","Depths")

jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon Y Analysis//Porosity Lagged Scatter Plot_Y.jpg")
hscat(Porosity~1,SpatialTrimmedDataY,(0:6)*1000)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon Y Analysis//Porosity Variogram Map_Y.jpg")
Porosity_Variogram_mapY = variogram(Porosity~1,SpatialTrimmedDataY, map = TRUE,cutoff=10000,width=1000)
plot(Porosity_Variogram_mapY)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon Y Analysis//Porosity Variogram Cloud_Y.jpg")
Porosity_VariogramCloudY = variogram(Porosity~1,SpatialTrimmedDataY, cloud = TRUE)
PoroPointPairsY=plot(Porosity_VariogramCloudY, identify = TRUE, digitize = TRUE, xlim = c(0,6000), ylim = c(0,0.14))
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon Y Analysis//Porosity Variogram CloudPoinPairs_Y.jpg")
plot(PoroPointPairsY, SpatialTrimmedDataY)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon Y Analysis//Porosity Emperical Variogram_Y.jpg")
Porosity_VariogramY = variogram(Porosity~1,SpatialTrimmedDataY)
plot(Porosity_VariogramY)
dev.off()



#CUT2 Data points between Q1 (0.1399) and Q3 (0.3040)
SmartTrimmedDataYCUT2 = SmartTrimmedData_Y
names(SmartTrimmedDataYCUT2) = c("BoreHoleID", "Horizon", "Depths", "Porosity", "x", "y")
SpatialTrimmedDataYCUT2 = SmartTrimmedDataYCUT2[SmartTrimmedDataYCUT2$Porosity<0.3040&SmartTrimmedDataYCUT2$Porosity>0.1399, ]


library(sp)
library(gstat)

coordinates(SpatialTrimmedDataYCUT2) = c("x","y","Depths")

jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon Y Analysis//Porosity Lagged Scatter Plot_YCUT2.jpg")
hscat(Porosity~1,SpatialTrimmedDataYCUT2,(0:6)*1000)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon Y Analysis//Porosity Variogram Map_YCUT2.jpg")
Porosity_Variogram_mapYCUT2 = variogram(Porosity~1,SpatialTrimmedDataYCUT2, map = TRUE,cutoff=10000,width=1000)
plot(Porosity_Variogram_mapYCUT2)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon Y Analysis//Porosity Variogram Cloud_YCUT2.jpg")
Porosity_VariogramCloudYCUT2 = variogram(Porosity~1,SpatialTrimmedDataYCUT2, cloud = TRUE)
PoroPointPairsYCUT2=plot(Porosity_VariogramCloudYCUT2, identify = TRUE, digitize = TRUE, xlim = c(0,6000), ylim = c(0,0.015))
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon Y Analysis//Porosity Variogram CloudPoinPairs_YCUT2.jpg")
plot(PoroPointPairsYCUT2, SpatialTrimmedDataYCUT2)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Horizon Y Analysis//Porosity Emperical Variogram_YCUT2.jpg")
Porosity_VariogramYCUT2 = variogram(Porosity~1,SpatialTrimmedDataYCUT2)
plot(Porosity_VariogramYCUT2)
dev.off()






#Agbabu Bitumen Deposit Top-Thickness Spatial Data Analysis

#Horizon X

Top_ThicknessDataX = read.csv(file.choose(),header=T)
Top_ThicknessDataX$Borehole = as.factor(Top_ThicknessDataX$Borehole)

library(sp)
library(gstat)

SpatialTop_ThicknessDataX = Top_ThicknessDataX
names(SpatialTop_ThicknessDataX) = c("BoreHoleID", "TopX", "ThicknessX", "x", "y")
coordinates(SpatialTop_ThicknessDataX) = c("x","y")

#Depth to Top Analysis
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis//TopX Lagged Scatter Plot.jpg")
hscat(TopX~1,SpatialTop_ThicknessDataX,(0:6)*1000)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis//TopX Variogram Map.jpg")
TopXVariogram_map = variogram(TopX~1,SpatialTop_ThicknessDataX, map = TRUE,cutoff=10000,width=1000)
plot(TopXVariogram_map)
dev.off()

TopXVariogram_Cloud = variogram(TopX~1,SpatialTop_ThicknessDataX, cutoff = 8000, cloud = TRUE)
#TopXVariogram_Cloud = data.frame(TopXVariogram_Cloud)
ggplot(mapping = aes(x = dist, y = gamma) )+
  geom_point(data = TopXVariogram_Cloud[c(-245,-249),], color = "blue")+
  labs(x = "Lag Distance, m", y = "Depth-to-top Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 2000))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 2000, by = 250))+
  geom_point(data = TopXVariogram_Cloud[c(245,249),], color = "red", shape = 12, size = 3)
ggsave("TopX Variogram Cloud.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#Directional Variogram Cloud

#90 deg
TopXVariogram_Cloud_90 = variogram(TopX~1,SpatialTop_ThicknessDataX, alpha = 90, tol.hor = 45, cutoff = 8000, cloud = TRUE)
ggplot(data =TopXVariogram_Cloud_90, mapping = aes(x = dist, y = gamma) )+
  geom_point(color = "blue")+
  labs(x = "Lag Distance, m", y = "Depth-to-top Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 2000))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 2000, by = 250))
ggsave("TopX Variogram Cloud_90azim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#0 deg
TopXVariogram_Cloud_0 = variogram(TopX~1,SpatialTop_ThicknessDataX, alpha = 0, tol.hor = 45, cutoff = 3000, cloud = TRUE)
ggplot(data =TopXVariogram_Cloud_0, mapping = aes(x = dist, y = gamma) )+
  geom_point(color = "blue")+
  labs(x = "Lag Distance, m", y = "Depth-to-top Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 2000))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 2000, by = 250))
ggsave("TopX Variogram Cloud_0azim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#Combined Directional Cloud using facet plot
TopXvariogram_cloud_allazim = rbind(TopXVariogram_Cloud_90,TopXVariogram_Cloud_0)
ggplot(data = TopXvariogram_cloud_allazim, aes(x = dist, y = gamma))+
  labs(x = "Lag Distance, m", y = "Depth-to-top Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 2000))+
  geom_point(color = "blue") + 
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 2000, by = 500))+
  facet_wrap(~ dir.hor, nrow = 2)
ggsave("TopX Variogram Cloud_allazim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#Directional Empirical Variogram

#90 deg
TopXEmpVariogram_90 = variogram(TopX~1,SpatialTop_ThicknessDataX, alpha = 90, tol.hor = 45, cutoff = 8000, boundaries = c(2000, 4000, 5500, 7000,8000))
TopXVariogModel_90 = fit.variogram(TopXEmpVariogram_90, vgm("Sph", psill = 410),fit.sills = TRUE, fit.ranges = TRUE)
TopXVariogModelFunc_90 = function(x){
  ifelse(x <= TopXVariogModel_90$range, (TopXVariogModel_90$psill*((1.5*(x/TopXVariogModel_90$range))-(0.5*((x/TopXVariogModel_90$range)^3)))), TopXVariogModel_90$psill)
}

#0 deg
TopXEmpVariogram_0 = variogram(TopX~1,SpatialTop_ThicknessDataX, alpha = 0, tol.hor = 45, cutoff = 2500, boundaries = c(750, 1500, 2500))
TopXVariogModel_0 = fit.variogram(TopXEmpVariogram_0, vgm("Sph", psill = 422.36),fit.sills = FALSE, fit.ranges = TRUE)
TopXVariogModelFunc_0 = function(x){
  ifelse(x <= TopXVariogModel_0$range, (TopXVariogModel_0$psill*((1.5*(x/TopXVariogModel_0$range))-(0.5*((x/TopXVariogModel_0$range)^3)))), TopXVariogModel_0$psill)
}

#Combined Directional Empirical and model Variogram plots using facet


TopXEmpVariogram_allazim = rbind(TopXEmpVariogram_0,TopXEmpVariogram_90)

TopXdatarange_0 = seq(0,4000,1)
TopXVariogModelVal_0 = data.frame(TopXdatarange_0,TopXVariogModelFunc_0(TopXdatarange_0),c(rep(0,length(TopXdatarange_0))))
names(TopXVariogModelVal_0) = c("Moddist", "Modgamma", "dir.hor")

TopXdatarange_90 = seq(0,8000,1)
TopXVariogModelVal_90 = data.frame(TopXdatarange_90,TopXVariogModelFunc_90(TopXdatarange_90),c(rep(90,length(TopXdatarange_90))))
names(TopXVariogModelVal_90) = c("Moddist", "Modgamma", "dir.hor")

TopXVariogModelVal_allazim = rbind(TopXVariogModelVal_0,TopXVariogModelVal_90)

ModelDisplay_TopX = data.frame(
label = c("model: 422.36*Sph[2377.69]*(h)", "model: 422.36*Sph[5268.64]*(h)"),
dir.hor   = c(0, 90), parse=T
)


ggplot(data = TopXEmpVariogram_allazim, aes(x = dist, y = gamma))+
  labs(x = "Lag Distance, m", y = "Depth-to-top Variogram, m-sq")+
  geom_point(aes(shape = "mu"),color = "blue") +
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 600))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 600, by = 100))+
  geom_line(data = TopXVariogModelVal_allazim,  aes(x=Moddist, y=Modgamma, color = 'jj'), size = 0.75)+
  geom_label(data    = eval(ModelDisplay_TopX), mapping = aes(x = 4000, y = -Inf, label = label),hjust   = -0.1,vjust   = -1, parse = T, color = "blue", fill = "light grey", label.size = 0)+
  scale_colour_manual(name = '', values =c('jj'='red'), labels = c('model'))+
  scale_shape_manual(name = '', values = c('mu' = 16), labels = c('empirical'))+
  #theme(legend.position=c(0.5, 0.5))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))+
  facet_wrap(~ dir.hor, nrow = 2)
ggsave("TopX Emperical&Fitted Variogram_allallazim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)





#Thickness Analysis
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis//ThicknessX Lagged Scatter Plot.jpg")
hscat(ThicknessX~1,SpatialTop_ThicknessDataX,(0:6)*1000)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis//ThicknessX Variogram Map.jpg")
ThicknessXVariogram_map = variogram(ThicknessX~1,SpatialTop_ThicknessDataX, map = TRUE,cutoff=4000,width=800)
plot(ThicknessXVariogram_map)
dev.off()

ThicknessXVariogram_Cloud = variogram(ThicknessX~1,SpatialTop_ThicknessDataX, cutoff = 8000, cloud = TRUE)
#ThicknessXVariogram_Cloud = data.frame(ThicknessXVariogram_Cloud)
ggplot(data = ThicknessXVariogram_Cloud,mapping = aes(x = dist, y = gamma) )+
  geom_point(color = "blue")+
  labs(x = "Lag Distance, m", y = "Thickness Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 150))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 150, by = 25))
ggsave("ThicknessX Variogram Cloud.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#Directional Variogram Cloud

#45 deg
ThicknessXVariogram_Cloud_45 = variogram(ThicknessX~1,SpatialTop_ThicknessDataX, alpha = 45, tol.hor = 45, cutoff = 8000, cloud = TRUE)
ggplot(data =ThicknessXVariogram_Cloud_45, mapping = aes(x = dist, y = gamma) )+
  geom_point(color = "blue")+
  labs(x = "Lag Distance, m", y = "Thickness Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 150))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 2000, by = 25))
ggsave("ThicknessX Variogram Cloud_45azim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#0 deg
ThicknessXVariogram_Cloud_135 = variogram(ThicknessX~1,SpatialTop_ThicknessDataX, alpha = 135, tol.hor = 45, cutoff = 8000, cloud = TRUE)
ggplot(data =ThicknessXVariogram_Cloud_135, mapping = aes(x = dist, y = gamma) )+
  geom_point(color = "blue")+
  labs(x = "Lag Distance, m", y = "Thickness Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 150))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 150, by = 25))
ggsave("ThicknessX Variogram Cloud_135azim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#Combined Directional Cloud using facet plot
ThicknessXvariogram_cloud_allazim = rbind(ThicknessXVariogram_Cloud_45,ThicknessXVariogram_Cloud_135)
ggplot(data = ThicknessXvariogram_cloud_allazim, aes(x = dist, y = gamma))+
  labs(x = "Lag Distance, m", y = "Thickness Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 150))+
  geom_point(color = "blue") + 
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 150, by = 25))+
  facet_wrap(~ dir.hor, nrow = 2)
ggsave("ThicknessX Variogram Cloud_allazim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#Directional Empirical Variogram

#45 deg
ThicknessXEmpVariogram_45 = variogram(ThicknessX~1,SpatialTop_ThicknessDataX, alpha = 45, tol.hor = 45, cutoff = 8000, boundaries = c(2000, 4000, 6000,8000))
ThicknessXVariogModel_45 = fit.variogram(ThicknessXEmpVariogram_45, vgm("Sph", psill = 20),fit.sills = TRUE, fit.ranges = TRUE)
ThicknessXVariogModelFunc_45 = function(x){
  ifelse(x <= ThicknessXVariogModel_45$range, (ThicknessXVariogModel_45$psill*((1.5*(x/ThicknessXVariogModel_45$range))-(0.5*((x/ThicknessXVariogModel_45$range)^3)))), ThicknessXVariogModel_45$psill)
}

#135 deg
ThicknessXEmpVariogram_135 = variogram(ThicknessX~1,SpatialTop_ThicknessDataX, alpha = 135, tol.hor = 45, cutoff = 8000, boundaries = c(2000, 4000, 6000,8000))
ThicknessXVariogModel_135 = fit.variogram(ThicknessXEmpVariogram_135, vgm("Sph", psill = 23.13),fit.sills = FALSE, fit.ranges = TRUE)
ThicknessXVariogModelFunc_135 = function(x){
  ifelse(x <= ThicknessXVariogModel_135$range, (ThicknessXVariogModel_135$psill*((1.5*(x/ThicknessXVariogModel_135$range))-(0.5*((x/ThicknessXVariogModel_135$range)^3)))), ThicknessXVariogModel_135$psill)
}

#Combined Directional Empirical and model Variogram plots using facet


ThicknessXEmpVariogram_allazim = rbind(ThicknessXEmpVariogram_45,ThicknessXEmpVariogram_135)

ThicknessXdatarange_45 = seq(0,8000,1)
ThicknessXVariogModelVal_45 = data.frame(ThicknessXdatarange_45,ThicknessXVariogModelFunc_45(ThicknessXdatarange_45),c(rep(45,length(ThicknessXdatarange_45))))
names(ThicknessXVariogModelVal_45) = c("Moddist", "Modgamma", "dir.hor")

ThicknessXdatarange_135 = seq(0,8000,1)
ThicknessXVariogModelVal_135 = data.frame(ThicknessXdatarange_135,ThicknessXVariogModelFunc_135(ThicknessXdatarange_135),c(rep(135,length(ThicknessXdatarange_135))))
names(ThicknessXVariogModelVal_135) = c("Moddist", "Modgamma", "dir.hor")

ThicknessXVariogModelVal_allazim = rbind(ThicknessXVariogModelVal_45,ThicknessXVariogModelVal_135)

ModelDisplay_ThicknessX = data.frame(
  label = c("model: 23.13*Sph[4883.9]*(h)", "model: 23.13*Sph[1658.4]*(h)"),
  dir.hor   = c(45, 135), parse=T
)


ggplot(data = ThicknessXEmpVariogram_allazim, aes(x = dist, y = gamma))+
  labs(x = "Lag Distance, m", y = "Thickness Variogram, m-sq")+
  geom_point(aes(shape = "mu"),color = "blue") +
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 35))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 35, by = 5))+
  geom_line(data = ThicknessXVariogModelVal_allazim,  aes(x=Moddist, y=Modgamma, color = 'jj'), size = 0.75)+
  geom_label(data    = eval(ModelDisplay_ThicknessX), mapping = aes(x = 4000, y = -Inf, label = label),hjust   = -0.1,vjust   = -1, parse = T, color = "blue", fill = "light grey", label.size = 0)+
  scale_colour_manual(name = '', values =c('jj'='red'), labels = c('model'))+
  scale_shape_manual(name = '', values = c('mu' = 16), labels = c('empirical'))+
  #theme(legend.position=c(0.5, 0.5))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))+
  facet_wrap(~ dir.hor, nrow = 2)
ggsave("ThicknessX Emperical&Fitted Variogram_allallazim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)






#Horizon Y

Top_ThicknessDataY = read.csv(file.choose(),header=T)
Top_ThicknessDataY$Borehole = as.factor(Top_ThicknessDataY$Borehole)

library(sp)
library(gstat)

SpatialTop_ThicknessDataY = Top_ThicknessDataY
names(SpatialTop_ThicknessDataY) = c("BoreHoleID", "TopY", "ThicknessY", "x", "y")
coordinates(SpatialTop_ThicknessDataY) = c("x","y")

#Depth to Top Analysis
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis//TopY Lagged Scatter Plot.jpg")
hscat(TopY~1,SpatialTop_ThicknessDataY,(0:6)*1000)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis//TopY Variogram Map.jpg")
TopYVariogram_map = variogram(TopY~1,SpatialTop_ThicknessDataY, map = TRUE,cutoff=10000,width=1000)
plot(TopYVariogram_map)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis//TopY Variogram Cloud.jpg")
TopYVariogram_Cloud = variogram(TopY~1,SpatialTop_ThicknessDataY,cutoff = 8000, cloud = TRUE)
#TopYVariogram_Cloud = data.frame(TopYVariogram_Cloud)
ggplot(mapping = aes(x = dist, y = gamma) )+
  geom_point(data = TopYVariogram_Cloud[c(-186,-191),], color = "blue")+
  labs(x = "Lag Distance, m", y = "Depth-to-top Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 2250))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 2250, by = 250))+
  geom_point(data = TopYVariogram_Cloud[c(186,191),], color = "red", shape = 12, size = 3)
ggsave("TopY Variogram Cloud.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#Directional Variogram Cloud

#90 deg
TopYVariogram_Cloud_90 = variogram(TopY~1,SpatialTop_ThicknessDataY, alpha = 90, tol.hor = 45, cutoff = 8000, cloud = TRUE)
ggplot(data =TopYVariogram_Cloud_90, mapping = aes(x = dist, y = gamma) )+
  geom_point(color = "blue")+
  labs(x = "Lag Distance, m", y = "Depth-to-top Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 2000))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 2000, by = 250))
ggsave("TopY Variogram Cloud_90azim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#0 deg
TopYVariogram_Cloud_0 = variogram(TopY~1,SpatialTop_ThicknessDataY, alpha = 0, tol.hor = 45, cutoff = 2500, cloud = TRUE)
ggplot(data =TopYVariogram_Cloud_0, mapping = aes(x = dist, y = gamma) )+
  geom_point(color = "blue")+
  labs(x = "Lag Distance, m", y = "Depth-to-top Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 2000))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 2000, by = 250))
ggsave("TopY Variogram Cloud_0azim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#Combined Directional Cloud using facet plot
TopYvariogram_cloud_allazim = rbind(TopYVariogram_Cloud_90,TopYVariogram_Cloud_0)
ggplot(data = TopYvariogram_cloud_allazim, aes(x = dist, y = gamma))+
  labs(x = "Lag Distance, m", y = "Depth-to-top Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 2000))+
  geom_point(color = "blue") + 
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 2000, by = 500))+
  facet_wrap(~ dir.hor, nrow = 2)
ggsave("TopY Variogram Cloud_allazim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)



#Directional Empirical Variogram

#90 deg
TopYEmpVariogram_90 = variogram(TopY~1,SpatialTop_ThicknessDataY, alpha = 90, tol.hor = 45, cutoff = 8000, boundaries = c(2000, 4000, 5500, 7000,8000))
TopYVariogModel_90 = fit.variogram(TopYEmpVariogram_90, vgm("Sph", psill = 410),fit.sills = FALSE, fit.ranges = TRUE)
TopYVariogModelFunc_90 = function(x){
  ifelse(x <= TopYVariogModel_90$range, (TopYVariogModel_90$psill*((1.5*(x/TopYVariogModel_90$range))-(0.5*((x/TopYVariogModel_90$range)^3)))), TopYVariogModel_90$psill)
}

#0 deg
TopYEmpVariogram_0 = variogram(TopY~1,SpatialTop_ThicknessDataY, alpha = 0, tol.hor = 45, cutoff = 2500, boundaries = c(750, 1500, 2500))
TopYVariogModel_0 = fit.variogram(TopYEmpVariogram_0, vgm("Sph", psill = 410),fit.sills = FALSE, fit.ranges = TRUE)
TopYVariogModelFunc_0 = function(x){
  ifelse(x <= TopYVariogModel_0$range, (TopYVariogModel_0$psill*((1.5*(x/TopYVariogModel_0$range))-(0.5*((x/TopYVariogModel_0$range)^3)))), TopYVariogModel_0$psill)
}

#Combined Directional Empirical and model Variogram plots using facet

TopYEmpVariogram_allazim = rbind(TopYEmpVariogram_0,TopYEmpVariogram_90)

TopYdatarange_0 = seq(0,4000,1)
TopYVariogModelVal_0 = data.frame(TopYdatarange_0,TopYVariogModelFunc_0(TopYdatarange_0),c(rep(0,length(TopYdatarange_0))))
names(TopYVariogModelVal_0) = c("Moddist", "Modgamma", "dir.hor")

TopYdatarange_90 = seq(0,8000,1)
TopYVariogModelVal_90 = data.frame(TopYdatarange_90,TopYVariogModelFunc_90(TopYdatarange_90),c(rep(90,length(TopYdatarange_90))))
names(TopYVariogModelVal_90) = c("Moddist", "Modgamma", "dir.hor")

TopYVariogModelVal_allazim = rbind(TopYVariogModelVal_0,TopYVariogModelVal_90)

ModelDisplay_TopY = data.frame(
  label = c("model: 410*Sph[2283.65]*(h)", "model: 410*Sph[4619.43]*(h)"),
  dir.hor   = c(0, 90), parse=T
)


ggplot(data = TopYEmpVariogram_allazim, aes(x = dist, y = gamma))+
  labs(x = "Lag Distance, m", y = "Depth-to-top Variogram, m-sq")+
  geom_point(aes(shape = "mu"),color = "blue") +
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 600))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 600, by = 100))+
  geom_line(data = TopYVariogModelVal_allazim,  aes(x=Moddist, y=Modgamma, color = 'jj'), size = 0.75)+
  geom_label(data    = eval(ModelDisplay_TopY), mapping = aes(x = 4000, y = -Inf, label = label),hjust   = -0.1,vjust   = -1, parse = T, color = "blue", fill = "light grey", label.size = 0)+
  scale_colour_manual(name = '', values =c('jj'='red'), labels = c('model'))+
  scale_shape_manual(name = '', values = c('mu' = 16), labels = c('empirical'))+
  #theme(legend.position=c(0.5, 0.5))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))+
  facet_wrap(~ dir.hor, nrow = 2)
ggsave("TopY Emperical&Fitted Variogram_allallazim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)



#Thickness Analysis
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis//ThicknessY Lagged Scatter Plot.jpg")
hscat(ThicknessY~1,SpatialTop_ThicknessDataY,(0:6)*1000)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis//ThicknessY Variogram Map.jpg")
ThicknessYVariogram_map = variogram(ThicknessY~1,SpatialTop_ThicknessDataY, map = TRUE,cutoff=10000,width=1000)
plot(ThicknessYVariogram_map)
dev.off()
jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis//ThicknessY Variogram Cloud.jpg")
ThicknessYVariogram_Cloud = variogram(ThicknessY~1,SpatialTop_ThicknessDataY,cutoff = 8000, cloud = TRUE)
#ThicknessYVariogram_Cloud = data.frame(ThicknessYVariogram_Cloud)
ggplot(data = ThicknessYVariogram_Cloud, mapping = aes(x = dist, y = gamma) )+
  geom_point(color = "blue")+
  labs(x = "Lag Distance, m", y = "Thickness Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 400))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 400, by = 50))
ggsave("ThicknessY Variogram Cloud.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#Directional Variogram Cloud

#45 deg
ThicknessYVariogram_Cloud_45 = variogram(ThicknessY~1,SpatialTop_ThicknessDataY, alpha = 45, tol.hor = 45, cutoff = 8000, cloud = TRUE)
ggplot(data =ThicknessYVariogram_Cloud_45, mapping = aes(x = dist, y = gamma) )+
  geom_point(color = "blue")+
  labs(x = "Lag Distance, m", y = "Thickness Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 400))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 400, by = 50))
ggsave("ThicknessY Variogram Cloud_45azim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#135 deg
ThicknessYVariogram_Cloud_135 = variogram(ThicknessY~1,SpatialTop_ThicknessDataY, alpha = 135, tol.hor = 45, cutoff = 8000, cloud = TRUE)
ggplot(data =ThicknessYVariogram_Cloud_135, mapping = aes(x = dist, y = gamma) )+
  geom_point(color = "blue")+
  labs(x = "Lag Distance, m", y = "Thickness Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 400))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 400, by = 50))
ggsave("ThicknessY Variogram Cloud_135azim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)

#Combined Directional Cloud using facet plot
ThicknessYvariogram_cloud_allazim = rbind(ThicknessYVariogram_Cloud_45,ThicknessYVariogram_Cloud_135)
ggplot(data = ThicknessYvariogram_cloud_allazim, aes(x = dist, y = gamma))+
  labs(x = "Lag Distance, m", y = "Thickness Pair Variogram, m-sq")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 400))+
  geom_point(color = "blue") + 
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 400, by = 50))+
  facet_wrap(~ dir.hor, nrow = 2)
ggsave("ThicknessY Variogram Cloud_allazim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)



#Directional Empirical Variogram

#45 deg
ThicknessYEmpVariogram_45 = variogram(ThicknessY~1,SpatialTop_ThicknessDataY, alpha = 45, tol.hor = 45, cutoff = 8000, boundaries = c(2000, 3500, 5700, 8000))
ThicknessYVariogModel_45 = fit.variogram(ThicknessYEmpVariogram_45, vgm("Sph", psill = 62),fit.sills = TRUE, fit.ranges = TRUE)
ThicknessYVariogModelFunc_45 = function(x){
  ifelse(x <= ThicknessYVariogModel_45$range, (ThicknessYVariogModel_45$psill*((1.5*(x/ThicknessYVariogModel_45$range))-(0.5*((x/ThicknessYVariogModel_45$range)^3)))), ThicknessYVariogModel_45$psill)
}

#135 deg
ThicknessYEmpVariogram_135 = variogram(ThicknessY~1,SpatialTop_ThicknessDataY, alpha = 135, tol.hor = 45, cutoff = 8000, boundaries = c(2000, 4000, 8000))
ThicknessYVariogModel_135 = fit.variogram(ThicknessYEmpVariogram_135, vgm("Sph", psill = 63.42),fit.sills = FALSE, fit.ranges = TRUE)
ThicknessYVariogModelFunc_135 = function(x){
  ifelse(x <= ThicknessYVariogModel_135$range, (ThicknessYVariogModel_135$psill*((1.5*(x/ThicknessYVariogModel_135$range))-(0.5*((x/ThicknessYVariogModel_135$range)^3)))), ThicknessYVariogModel_135$psill)
}

#Combined Directional Empirical and model Variogram plots using facet

ThicknessYEmpVariogram_allazim = rbind(ThicknessYEmpVariogram_45,ThicknessYEmpVariogram_135)

ThicknessYdatarange_45 = seq(0,8000,1)
ThicknessYVariogModelVal_45 = data.frame(ThicknessYdatarange_45,ThicknessYVariogModelFunc_45(ThicknessYdatarange_45),c(rep(45,length(ThicknessYdatarange_45))))
names(ThicknessYVariogModelVal_45) = c("Moddist", "Modgamma", "dir.hor")

ThicknessYdatarange_135 = seq(0,8000,1)
ThicknessYVariogModelVal_135 = data.frame(ThicknessYdatarange_135,ThicknessYVariogModelFunc_135(ThicknessYdatarange_135),c(rep(135,length(ThicknessYdatarange_135))))
names(ThicknessYVariogModelVal_135) = c("Moddist", "Modgamma", "dir.hor")

ThicknessYVariogModelVal_allazim = rbind(ThicknessYVariogModelVal_45,ThicknessYVariogModelVal_135)

ModelDisplay_ThicknessY = data.frame(
  label = c("model: 63.42*Sph[2806.10]*(h)", "model: 63.42*Sph[1867.00]*(h)"),
  dir.hor   = c(45, 135), parse=T
)


ggplot(data = ThicknessYEmpVariogram_allazim, aes(x = dist, y = gamma))+
  labs(x = "Lag Distance, m", y = "Thickness Variogram, m-sq")+
  geom_point(aes(shape = "mu"),color = "blue") +
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 80))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 80, by = 10))+
  geom_line(data = ThicknessYVariogModelVal_allazim,  aes(x=Moddist, y=Modgamma, color = 'jj'), size = 0.75)+
  geom_label(data    = eval(ModelDisplay_ThicknessY), mapping = aes(x = 4000, y = -Inf, label = label),hjust   = -0.1,vjust   = -1, parse = T, color = "blue", fill = "light grey", label.size = 0)+
  scale_colour_manual(name = '', values =c('jj'='red'), labels = c('model'))+
  scale_shape_manual(name = '', values = c('mu' = 16), labels = c('empirical'))+
  #theme(legend.position=c(0.5, 0.5))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))+
  facet_wrap(~ dir.hor, nrow = 2)
ggsave("ThicknessY Emperical&Fitted Variogram_allallazim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Top-Thickness Analysis", dpi = 96)



#Vertical Porosity Variogram

Porosity_VariogramCloudCUT2_Vertical = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 40,  alpha = 90, tol.hor = 0, beta = 90, tol.ver = 0, cloud = TRUE)
plot(Porosity_VariogramCloudCUT2_Vertical)
ggplot(data = Porosity_VariogramCloudCUT2_Vertical)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram")+
  scale_x_continuous(breaks = seq(0, 40, by = 5))+
  scale_y_continuous(breaks = seq(0, 0.04, by = 0.01))
ggsave("Porosity Vertical Variogram Cloud.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


Porosity_VariogramCUT2_Vertical = variogram(Porosity~1,SpatialTrimmedDataCUT2, cutoff = 40, boundaries = seq(from = 2, to = 40, by = 4),  alpha = 0, tol.hor = 0, beta = 90, tol.ver = 0)
vertmodel = fit.variogram(Porosity_VariogramCUT2_Vertical, vgm(model = "Sph", range = 15, nugget = 0.0024, psill = 0.0016), fit.sills = FALSE, fit.ranges = FALSE)

VertVariogModelFunc = function(x){
  ifelse(x <= vertmodel$range[2], (vertmodel$psill[1]+(vertmodel$psill[2]*((1.5*(x/vertmodel$range[2]))-(0.5*((x/vertmodel$range[2])^3))))), sum(vertmodel$psill))
}
ggplot(data = Porosity_VariogramCUT2_Vertical, aes(x = dist, y = gamma))+
  geom_point(aes(shape = "mu"), color = "blue")+
  labs(x = "Lag Distance, h (m)", y = "Porosity Variogram")+
  coord_cartesian(xlim = c(0, 40), ylim = c(0, max(Porosity_VariogramCUT2_Vertical$gamma,0.006)))+
  scale_x_continuous(breaks = seq(0, 40, by = 5))+
  scale_y_continuous(breaks = seq(0, max(max(Porosity_VariogramCUT2_Vertical$gamma)+0.001,0.006), by = 0.001))+
  stat_function(aes(x = dist, color = 'jj'), fun = VertVariogModelFunc, xlim = c(0, 40), size = 1)+
  scale_colour_manual(name = '', values =c('jj'='red'), labels = expression(model:0.0024+0.0016*Sph[15]*(h)))+
  scale_shape_manual(name = '', values = c('mu' = 16), labels = c('empirical'))+
  theme(legend.position=c(0.6, 0.3))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
ggsave("Porosity Emperical&Fitted Vertical Variogram.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis//Model Fits", dpi = 96)


#Integrated 3-D Porosity Variogram Modeling - TTOWG!!!!!!!!!!!!!
zonalanisomodel = vgm(psill = 0.0005, "Sph", range = 1000000000, anis = c(0,90,0,3.5e-6,1.5e-6))
Integrated3Dmodel = vgm(nugget = 0.0024, psill = 0.0016, range = 3500, model = "Sph", anis=c(90,0,0,0.428571428,4.28571428e-3), add.to = zonalanisomodel)






#Generating plots with ggplot2, for presentation
ggplotPorosity_VariogramCloudWhole = ggplot(data = Porosity_VariogramCloudWhole)+
                                      geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
                                      labs(x = "Lag Distance, m", y = "Porosity Pair Variogram")+
                                      scale_x_continuous(breaks = seq(0, 6000, by = 1000))+
                                      scale_y_continuous(breaks = seq(0, 0.12, by = 0.02))
                                      

ggplotPorosity_VariogramWhole = ggplot(data = Porosity_VariogramWhole)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "red")+
  labs(x = "Lag Distance, m", y = "Porosity Experimental Variogram")+
  coord_cartesian(xlim = c(0, 6000), ylim = c(0, 0.012))+
  scale_x_continuous(breaks = seq(0, 6000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.012, by = 0.002))+
  geom_smooth(mapping = aes(x = dist, y = gamma)) # a preliminary variogram modeling attempt.
  
     #combining cloud with experimental variogram
ggplotPorosity_VariogramCloudWhole+geom_point(data = Porosity_VariogramWhole,mapping = aes(x = dist, y = gamma), color = "red")


allazimvariogram_binbased = rbind(Porosity_VariogramCUT2_45_22.5,Porosity_VariogramCUT2_90_22.5,Porosity_VariogramCUT2_135_22.5,Porosity_VariogramCUT2_180_22.5)
  ggplot(data = allazimvariogram, aes(x = dist, y = gamma))+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  geom_point(color = "blue") +
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.008))+
  scale_x_continuous(breaks = seq(0, 8000, by = 2000))+
  scale_y_continuous(breaks = seq(0, 0.008, by = 0.001))+
  facet_wrap(~ dir.hor, nrow = 2)
ggsave("Porosity Emperical Variogram_CUT2_binbased_allallazim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


   #All azimuths Empirical and Fitted Variogram Model Facet Plot - gslib
allazimvariogram_lagbased = rbind(Porosity_VariogramCUT2_45_22.5_gslib,Porosity_VariogramCUT2_90_22.5_gslib,Porosity_VariogramCUT2_135_22.5_gslib,Porosity_VariogramCUT2_180_22.5_gslib)
lagbasedVariogModelFunc_45_visual = function(x){
  ifelse(x <= 2000, (0.003+(0.0012*((1.5*(x/2000))-(0.5*((x/2000)^3))))), 0.0042)
}
lagbasedVariogModelVal_45_visual = data.frame(Porosity_VariogramCUT2_45_22.5_gslib$dist,lagbasedVariogModelFunc_45_visual(Porosity_VariogramCUT2_45_22.5_gslib$dist))
names(lagbasedVariogModelVal_45_visual) = c("Moddist", "Modgamma")
lagbasedVariogModelFunc_90_visual = function(x){
  ifelse(x <= 1700, (0.003+(0.0015*((1.5*(x/1700))-(0.5*((x/1700)^3))))), 0.0045)
}
lagbasedVariogModelVal_90_visual = data.frame(Porosity_VariogramCUT2_90_22.5_gslib$dist,lagbasedVariogModelFunc_90_visual(Porosity_VariogramCUT2_90_22.5_gslib$dist))
names(lagbasedVariogModelVal_90_visual) = c("Moddist", "Modgamma")
lagbasedVariogModelFunc_135_visual = function(x){
  ifelse(x <= 2000, (0.003+(0.001*((1.5*(x/2000))-(0.5*((x/2000)^3))))), 0.004)
}
lagbasedVariogModelVal_135_visual = data.frame(Porosity_VariogramCUT2_135_22.5_gslib$dist,lagbasedVariogModelFunc_135_visual(Porosity_VariogramCUT2_135_22.5_gslib$dist))
names(lagbasedVariogModelVal_135_visual) = c("Moddist", "Modgamma")
lagbasedVariogModelFunc_180_visual = function(x){
  ifelse(x <= 1000, (0.0034+(0.0007*((1.5*(x/1000))-(0.5*((x/1000)^3))))), 0.0041)
}
lagbasedVariogModelVal_180_visual = data.frame(Porosity_VariogramCUT2_180_22.5_gslib$dist,lagbasedVariogModelFunc_180_visual(Porosity_VariogramCUT2_180_22.5_gslib$dist))
names(lagbasedVariogModelVal_180_visual) = c("Moddist", "Modgamma")

allazimvariogramModelVal_lagbased = rbind(lagbasedVariogModelVal_45_visual,lagbasedVariogModelVal_90_visual,lagbasedVariogModelVal_135_visual,lagbasedVariogModelVal_180_visual)
allazimvariogram_EmpandMod_lagbased = data.frame(allazimvariogram_lagbased,allazimvariogramModelVal_lagbased)

ModelDisplay = data.frame(
  label = c("model:0.003+0.0012*Sph[2000]*(h)", "model:0.003+0.0015*Sph[1700]*(h)", "model:0.003+0.001*Sph[2000]*(h)","model:0.0034+0.0007*Sph[1000]*(h)"),
  dir.hor   = c(45, 90, 135, 180), parse=T
)

ggplot(data = allazimvariogram_EmpandMod_lagbased, aes(x = dist, y = gamma))+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  geom_point(color = "blue") +
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.008))+
  scale_x_continuous(breaks = seq(0, 8000, by = 2000))+
  scale_y_continuous(breaks = seq(0, 0.008, by = 0.001))+
  geom_smooth(aes(x=dist, y=Modgamma), se = FALSE, colour = "red")+
  geom_label(data    = eval(ModelDisplay), mapping = aes(x = -Inf, y = -Inf, label = label),hjust   = -0.1,vjust   = -1, parse = T, color = "blue", fill = "light grey", label.size = 0)+
  facet_wrap(~ dir.hor, nrow = 2)
ggsave("Porosity Emperical&Fitted Variogram_CUT2_lagbased_allallazim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)


allazimvariogramcloud = rbind(Porosity_VariogramCloudCUT2_45_22.5,Porosity_VariogramCloudCUT2_90_22.5,Porosity_VariogramCloudCUT2_135_22.5,Porosity_VariogramCloudCUT2_180_22.5)
ggplotPorosity_VariogramCloudCUT2_allazim = ggplot(data = allazimvariogramcloud, aes(x = dist, y = gamma))+
                                            labs(x = "Lag Distance, m", y = "Porosity Pair Variogram")+
                                            geom_point(color = "blue") + 
                                            facet_wrap(~ dir.hor, nrow = 2)
                                       

Porosity_VariogramCloudCUT2_90_22.5_classified = data.frame(Porosity_VariogramCloudCUT2_90_22.5,CloudCUT2_90_22.5_Pointpaircategory)
Porosity_VariogramCloudCUT2_90_22.5_classified$CloudCUT2_90_22.5_Pointpaircategory = as.factor(Porosity_VariogramCloudCUT2_90_22.5_classified$CloudCUT2_90_22.5_Pointpaircategory)
ggplotPorosity_VariogramCloudCUT2_90_colored = ggplot(data = Porosity_VariogramCloudCUT2_90_22.5_classified)+
                                               geom_point(mapping = aes(x = dist, y = gamma, color = CloudCUT2_90_22.5_Pointpaircategory))+
                                               labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Bin" )+
                                               coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
                                               scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
                                               scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))




ggplotPorosity_VariogramCloudCUT2_90_dbscancolored = ggplot(data = dbscanVariogCloud[!dbscanVariogCloud$porodbscan.cluster=="0",])+
  geom_point(mapping = aes(x = Porosity_VariogramCloudCUT2_90_22.5.dist, y = Porosity_VariogramCloudCUT2_90_22.5.gamma, color = porodbscan.cluster))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Bin" )+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))



names(dbscanVariogCloud_45) = c("dist", "gamma", "clusterID", "direction")
names(dbscanVariogCloud_90) = c("dist", "gamma", "clusterID", "direction")
names(dbscanVariogCloud_135) = c("dist", "gamma", "clusterID", "direction")
names(dbscanVariogCloud_180) = c("dist", "gamma", "clusterID", "direction")

allazimdbscanVariogCloud = rbind(dbscanVariogCloud_45, dbscanVariogCloud_90, dbscanVariogCloud_135, dbscanVariogCloud_180)
allazimdbscanVariogCloud$clusterID = as.factor(allazimdbscanVariogCloud$clusterID) #Ensuring R sees the cluster numbers as classification ID and not as numerals
ggplot(data = allazimdbscanVariogCloud[!allazimdbscanVariogCloud$clusterID=="0",])+
  geom_point(mapping = aes(x = dist, y = gamma, color = clusterID, shape = clusterID))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Cluster" )+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
  scale_x_continuous(breaks = seq(0, 8000, by = 2000))+
  scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))+
  scale_color_manual(name = 'Cluster', values = c("#000000","#00FF00","#330033","#000099","#000033","#FF0000"))+
  scale_shape_manual(name = 'Cluster', values = c(16,8,15,17,18,4))+
  facet_wrap(~direction, nrow = 2)
ggsave("Porosity VariogramCloud_CUT2_dbscan_colored_allallazim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)




dbscanExpVariog_90deg = data.frame(dbscanExpVariog,c(rep(90,6)))
names(dbscanExpVariog_90deg)[4] = "Direction"
dbscanExpVariog_45deg = data.frame(dbscanExpVariog_45,c(rep(45,6)))
names(dbscanExpVariog_45deg)[4] = "Direction"
dbscanExpVariog_135deg = data.frame(dbscanExpVariog_135,c(rep(135,6)))
names(dbscanExpVariog_135deg)[4] = "Direction"
dbscanExpVariog_180deg = data.frame(dbscanExpVariog_180,c(rep(180,5)))
names(dbscanExpVariog_180deg)[4] = "Direction"
allazimdbscanvariogram = rbind(dbscanExpVariog_45deg,dbscanExpVariog_90deg,dbscanExpVariog_135deg,dbscanExpVariog_180deg)


datarange_45 = seq(0,(max(dbscanExpVariog_45deg$AverageDistance)+1),1)
dbscanVariogModelVal_45 = data.frame(datarange_45,dbscanVariogModelFunc_45(datarange_45),c(rep(45,length(datarange_45))))
names(dbscanVariogModelVal_45) = c("Moddist", "Modgamma", "Direction")

datarange_90 = seq(0,(max(dbscanExpVariog_90deg$AverageDistance)+1),1)
dbscanVariogModelVal_90 = data.frame(datarange_90,dbscanVariogModelFunc_90(datarange_90),c(rep(90,length(datarange_90))))
names(dbscanVariogModelVal_90) = c("Moddist", "Modgamma", "Direction")

datarange_135 = seq(0,(max(dbscanExpVariog_135deg$AverageDistance)+1),1)
dbscanVariogModelVal_135 = data.frame(datarange_135,dbscanVariogModelFunc_135(datarange_135),c(rep(135,length(datarange_135))))
names(dbscanVariogModelVal_135) = c("Moddist", "Modgamma", "Direction")

datarange_180 = seq(0,(max(dbscanExpVariog_180deg$AverageDistance)+1),1)
dbscanVariogModelVal_180 = data.frame(datarange_180,dbscanVariogModelFunc_180(datarange_180),c(rep(180,length(datarange_180))))
names(dbscanVariogModelVal_180) = c("Moddist", "Modgamma", "Direction")

allazimvariogramModelVal_dbscan = rbind(dbscanVariogModelVal_45,dbscanVariogModelVal_90,dbscanVariogModelVal_135,dbscanVariogModelVal_180)
#allazimvariogram_EmpandMod_dbscan = data.frame(allazimdbscanvariogram,allazimvariogramModelVal_dbscan)

#ModelDisplay_dbscan = data.frame(
#  label = c("model:0.0034+0.0008*Sph[2946]*(h)", "model:0.0034+0.0011*Sph[4215]*(h)", "model:0.0034+0.0008*Sph[2946]*(h)","model:0.0033+0.0007*Sph[743]*(h)"),
#  Direction   = c(45, 90, 135, 180), parse=T
#)


ggplot(data = allazimdbscanvariogram, aes(x = AverageDistance, y = ExpVariog))+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  geom_point(aes(shape = "mu"),color = "blue") +
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.007))+
  scale_x_continuous(breaks = seq(0, 8000, by = 2000))+
  scale_y_continuous(breaks = seq(0, 0.007, by = 0.001))+
  geom_line(data = allazimvariogramModelVal_dbscan,  aes(x=Moddist, y=Modgamma, color = 'jj'), size = 0.75)+
  #geom_label(data    = eval(ModelDisplay_dbscan), mapping = aes(x = -Inf, y = -Inf, label = label),hjust   = -0.1,vjust   = -1, parse = T, color = "blue", fill = "light grey", label.size = 0)+
  scale_colour_manual(name = '', values =c('jj'='red'), labels = c('model'))+
  scale_shape_manual(name = '', values = c('mu' = 16), labels = c('empirical'))+
  #theme(legend.position=c(0.5, 0.5))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))+
  facet_wrap(~ Direction, nrow = 2)
ggsave("Updated Porosity Emperical&Fitted Variogram_CUT2_dbscan_allallazim.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)





wa = ggplot(data = Porosity_VariogramCloudCUT2_90_22.5_classified)+
  geom_point(mapping = aes(x = dist, y = gamma, color = Bin, shape = Bin))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Bin" )+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.06))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.06, by = 0.01))+
  scale_color_manual(name = 'Bin', values = c("#000000","#00FF00","#330033","#000099","#000033","#FF0000", "#99FF00", "#003300", "#330000", "#FFCC00", "#330066", "#FF0066", "#0000FF", "#CC0000"))+
  scale_shape_manual(name = 'Bin', values = c(16,8,15,17,18,4,0,1,2,5,3,6,13,9))+
  annotate("rect", xmin=2300, xmax=2600, ymin=0, ymax=0.05, alpha=.3,fill="blue")+
  annotate("text", x=2500, y=0.052, label="Sparsity Effect", family="serif",fontface="italic", colour="blue", size=4)+
  annotate("rect", xmin=5700, xmax=7000, ymin=0, ymax=0.05, alpha=.3,fill="blue")+
  annotate("text", x=6500, y=0.052, label="Split Effect", family="serif",fontface="italic", colour="blue", size=4)

wo = ggplot(data = Porosity_VariogramCloudCUT2_90_22.5_classified_gslib)+
  geom_point(mapping = aes(x = dist, y = gamma, color = Lag_Interval, shape = Lag_Interval))+
  labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Lag Interval" )+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.06))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.06, by = 0.01))+
  scale_color_manual(name = 'Lag Interval', values = c("#000000","#FF0066","#330033","#000099","#000033","#FF0000", "#00FF00", "#003300", "#330000"))+
  scale_shape_manual(name = 'Lag Interval', values = c(16,8,15,17,3,4,0,1,2))+
  annotate("rect", xmin=6500, xmax=7400, ymin=0, ymax=0.05, alpha=.3,fill="blue")+
  annotate("text", x=7000, y=0.055, label="Split-Straddle \nEffect", family="serif",fontface="italic", colour="blue", size=4)


MultipleVarCloudPlot = arrangeGrob(wa, wo, nrow = 2)
grid.draw(MultipleVarCloudPlot, recording = TRUE)
ggsave("MultipleVarCloudPlot.jpg", MultipleVarCloudPlot, path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 48)

ggplot(data = Porosity_VariogramCUT2_90_22.5)+
  geom_point(mapping = aes(x = dist, y = gamma), color = "blue")+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))+
  geom_line(color = "red",aes(x = dist, y = gamma))+
  annotate("text", x=2494, y=0.0060, label="ep1", family="serif",fontface="italic", colour="blue", size=4)+
  annotate("text", x=3816, y=0.0027, label="ep2", family="serif",fontface="italic", colour="blue", size=4)+
  annotate("text", x=5509, y=0.0057, label="ep3", family="serif",fontface="italic", colour="blue", size=4)+
  annotate("text", x=7334, y=0.0061, label="ep4", family="serif",fontface="italic", colour="blue", size=4)+
  annotate("text", x=1650, y=0.0035, label="sp1", family="serif",fontface="italic", colour="blue", size=4)+
  annotate("text", x=1450, y=0.0048, label="sp2", family="serif",fontface="italic", colour="blue", size=4)
ggsave("Porosity Emperical Variogram_CUT2_90_22.5 annotated.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Plots//Field-wide Porosity Analysis", dpi = 96)






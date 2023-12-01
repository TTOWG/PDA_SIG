# To the Only Wise God - TTOWG


library(gstat)
library(sp)
library(dbscan)
library(tidyverse)
library(hydroGOF)
library(MLmetrics)
library(ggExtra)
library(gridExtra)

#Importing Sample Data from Mosobalaje et al (2019)
SamplePointsandData = read.csv("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis//Trimmed Sample Points and Data.csv",header=T)
names(SamplePointsandData)[5:7] = c("BlockEasting_big", "BlockNorthing_big", "BlockDepth_big") 


#Discretization
delta_X_big = 400
delta_Y_big = 400
delta_Z_big = 1
nx = 40
ny = 13
nz = 100
VarCutoff = 8000
grd_big=expand.grid(1:nx,1:ny,1:nz)
names(grd_big) = c("x","y","z")
BlockEasting_big = (delta_X_big*(grd_big$x-0.5))+700000 #This is the x-coordinate of the block center.
BlockNorthing_big = (delta_Y_big*(grd_big$y-0.5))+732500 #This is the y-coordinate of the block center.
BlockDepth_big = (delta_Z_big*(grd_big$z-0.5))+0 #This is the z-coordinate of the block center.
Scaledgrd_big = data.frame(BlockEasting_big, BlockNorthing_big, BlockDepth_big) 

#Simulations
zonalanisomodel = vgm(psill = 0.0005, "Sph", range = 1000000000, anis = c(0,90,0,3.5e-6,1.5e-6))
Integrated3Dmodel = vgm(nugget = 0.0024, psill = 0.0016, range = 3500, model = "Sph", anis=c(90,0,0,0.428571428,4.28571428e-3), add.to = zonalanisomodel)
#Integrated3Dmodel = vgm(nugget = 0.0024, psill = 0.0021, range = 3500, model = "Sph", anis=c(90,0,0,0.428571428,4.28571428e-3))
SimPar_big = gstat(formula = Porosity~1, locations = ~BlockEasting_big+BlockNorthing_big+BlockDepth_big, data = SamplePointsandData[,5:8], model = Integrated3Dmodel, nmin = 20, nmax = 50 )
number_of_anisorealizations = 78
Sim_big = predict(SimPar_big, newdata = Scaledgrd_big, nsim = number_of_anisorealizations,debug=-1)
write.table(Sim_big,"Sim_big_aniso.dat", sep = " ", row.names = F)

#Variogram reproduction using exhaustive data
GrandVarSim_90azim = data.frame(matrix(0,(VarCutoff/delta_X_big),3*number_of_anisorealizations))
GrandMeanAbsolutePercentDev_90azim = data.frame(matrix(0,number_of_anisorealizations,2))
names(GrandMeanAbsolutePercentDev_90azim) = c("SIMID","SIMMAPE")
for(i in 1:number_of_anisorealizations){
  SimID = paste("Sim",i,sep = " ")
  SimSubSet = data.frame(Sim_big[,1:3],Sim_big[,i+3])
  names(SimSubSet)[4] = "sim"
  #First, histogram - to ascertain gaussianity
  #Histogram
  ggplot(SimSubSet, aes(x= sim ))+
    geom_histogram(bins = 30,color = "blue", fill = "lightgreen")+
    labs(x = "Porosity", y = "Frequency")
  ggsave(paste(SimID,"Histogram.jpg",sep = " "), path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis//Exhaustive Data Plots", dpi = 96)
  coordinates(SimSubSet) = ~BlockEasting_big+BlockNorthing_big+BlockDepth_big
  #90-deg azimuth
  VarSim_90azim = variogram(sim~1,SimSubSet, alpha = 90, tol.hor = 0, beta=0, tol.ver = 0, cutoff = VarCutoff, boundaries = seq((0.5*delta_X_big), (VarCutoff+(0.5*delta_X_big)), delta_X_big))
  destinationcol_begin_90azim = (3*i)-2
  detsinationcol_end_90azim = 3*i
  names(GrandVarSim_90azim)[destinationcol_begin_90azim:detsinationcol_end_90azim] = c(paste("np",i,sep = ""),paste("dist",i,sep = ""),paste("gamma",i,sep = ""))
  GrandVarSim_90azim[destinationcol_begin_90azim:detsinationcol_end_90azim] = data.frame(VarSim_90azim$np,VarSim_90azim$dist,VarSim_90azim$gamma)
  VarModValues_90 = variogramLine(Integrated3Dmodel, VarCutoff, (VarCutoff/delta_X_big), dir=c(1,0,0), dist_vector = seq(from = 1, to = 8000, by = 1))
  MeanAbsolutePercentDev_90azim = round(MAPE(VarSim_90azim$gamma,VarModValues_90$gamma),2)
  GrandMeanAbsolutePercentDev_90azim[i,1] = SimID
  GrandMeanAbsolutePercentDev_90azim[i,2] = MeanAbsolutePercentDev_90azim
  poster_90azim = VarSim_90azim %>% 
    summarise(dist = 1000,gamma = 0.0005,posterstring_90 = paste("Mean Absolute Percentage Deviation = ",MeanAbsolutePercentDev_90azim,"%",sep=""))
  ggplot(data = VarSim_90azim, aes(x = dist, y = gamma))+
    geom_point(aes(shape = 'mu'), color = "blue")+
    labs(x = "Lag Distance, m", y = "Porosity Variogram")+
    coord_cartesian(xlim = c(0, VarCutoff), ylim = c(0, max(VarSim_90azim$gamma,0.006)))+
    scale_x_continuous(breaks = seq(0, VarCutoff, by = 1000))+
    scale_y_continuous(breaks = seq(0, max(VarSim_90azim$gamma,0.006), by = 0.001))+
    geom_line(data = VarModValues_90, aes(color = 'jj'), size = 0.75)+
    geom_label(aes(label = posterstring_90 ), data = poster_90azim, vjust = "bottom", hjust = "left")+
    scale_colour_manual(name = '', values =c("jj"='red'), labels = c("io" = 'empirical',"jj" = 'model'))+
    scale_shape_manual(name = '', values = c('mu' = 16), labels = c("mu" = 'empirical'))+
    theme(legend.position=c(0.8, 0.4))+
    theme(legend.title = element_blank())+
    theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
  ggsave(paste(SimID,"90deg Variogram Plot.jpg",sep = " "), path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis//Exhaustive Data Plots", dpi = 96)

}

ChoiceRealization = Sim_big[4]
  
  
  #Creating suite of sample points shiftings, and sampling the choice realization for the each shiftings

# Number of configuration (i.e. shifts) possible

SampleLocation = SamplePointsandData[1:7]
names(SampleLocation) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth")


Eastwardshifts = nx-max(SampleLocation$x_index)
Northwardshifts = ny-max(SampleLocation$y_index)
Downwardshifts = nz-max(SampleLocation$z_index)
Westwardshifts = min(SampleLocation$x_index)-0
Southwradshifts = min(SampleLocation$y_index)-0
Upwardshifts = min(SampleLocation$z_index)-0
TotalPossibleshifts =(Eastwardshifts+Westwardshifts)*(Northwardshifts+Southwradshifts)*(Downwardshifts+Upwardshifts)  

# a simple case to test and to illustrate
#nx = 5
#ny = 3
#nz = 3
#SampleLocation = read.csv("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis//Test Sample Location.csv",header=T)
#ChoiceRealization = data.frame(Sim_big[1:45,4])
#Eastwardshifts = nx-max(SampleLocation$x_index)
#Northwardshifts = ny-max(SampleLocation$y_index)
#Downwardshifts = nz-max(SampleLocation$z_index)
#Westwardshifts = min(SampleLocation$x_index)-0
#Southwradshifts = min(SampleLocation$y_index)-0
#Upwardshifts = min(SampleLocation$z_index)-0
#TotalPossibleshifts =(Eastwardshifts+Westwardshifts)*(Northwardshifts+Southwradshifts)*(Downwardshifts+Upwardshifts)  


  # Now the shifting and sampling; and variogram computation

# A placeholder for the variogram values across all configurations
GrandVarConfig_90azim = data.frame(matrix(0,9,3*TotalPossibleshifts))

# Original Sample Location plus its Westward shifts
SampleConfig_Position_ends_WW = ((Upwardshifts-1)*(Northwardshifts+Southwradshifts)*(Eastwardshifts+Westwardshifts))+((Southwradshifts-1)*(Eastwardshifts+Westwardshifts))+(Westwardshifts-1)+1 
SampleConfig_Position_begins_WW = SampleConfig_Position_ends_WW-Westwardshifts+1
Planarshifts = data.frame(matrix(0,362,(7*(Northwardshifts+Southwradshifts)*(Westwardshifts+Eastwardshifts))))
Horizontalshifts = data.frame(matrix(0,362,(7*(Westwardshifts+Eastwardshifts))))
for(r in 1:Westwardshifts){
  SampleLocPosition = SampleConfig_Position_ends_WW-(r-1)
  ShiftedSampleLocation = data.frame(SampleLocation$x_index-(r-1),SampleLocation$y_index,SampleLocation$z_index,SampleLocation$GridBlockID-(r-1),SampleLocation$TrueEasting-((r-1)*delta_X_big),SampleLocation$TrueNorthing,SampleLocation$TrueDepth)
  names(ShiftedSampleLocation) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth")
  Planarshifts[,((7*(Southwradshifts-1)*(Westwardshifts+Eastwardshifts))+(7*(Westwardshifts-r))+1):((7*(Southwradshifts-1)*(Westwardshifts+Eastwardshifts))+(7*(Westwardshifts-(r-1))))] = ShiftedSampleLocation
  names(Planarshifts)[((7*(Southwradshifts-1)*(Westwardshifts+Eastwardshifts))+(7*(Westwardshifts-r))+1):((7*(Southwradshifts-1)*(Westwardshifts+Eastwardshifts))+(7*(Westwardshifts-(r-1))))] = c(paste("x_index",SampleLocPosition, sep = ""),paste("y_index",SampleLocPosition, sep = ""),paste("z_index",SampleLocPosition, sep = ""),paste("GridBlockID",SampleLocPosition, sep = ""),paste("TrueEasting",SampleLocPosition, sep = ""),paste("TrueNorthing",SampleLocPosition, sep = ""),paste("TrueDepth",SampleLocPosition, sep = ""))
  Horizontalshifts[,(7*(Westwardshifts-r)+1):(7*(Westwardshifts-(r-1)))] = ShiftedSampleLocation 
  SampleData = data.frame(matrix(0,362,1)) 
    for (s in 1:362){
    SampleData[s,] = ChoiceRealization[as.numeric(rownames(ChoiceRealization))==ShiftedSampleLocation$GridBlockID[s],]
    }
  names(SampleData) = c("SampleData")
  assign(paste("SampleConfig",SampleLocPosition,sep = ""), data.frame(ShiftedSampleLocation,SampleData))
  ShiftedSampleLocanddata = data.frame(ShiftedSampleLocation,SampleData)
  coordinates(ShiftedSampleLocanddata) = ~TrueEasting+TrueNorthing+TrueDepth
  VarConfig_90azim = variogram(SampleData~1,ShiftedSampleLocanddata, alpha = 90, tol.hor = 22.5, cutoff = VarCutoff, boundaries = c(500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500))
  Configdestinationcol_begin_90azim = (3*SampleLocPosition)-2
  Configdestinationcol_end_90azim = 3*SampleLocPosition
  GrandVarConfig_90azim[Configdestinationcol_begin_90azim:Configdestinationcol_end_90azim] = data.frame(VarConfig_90azim$np,VarConfig_90azim$dist,VarConfig_90azim$gamma)
  names(GrandVarConfig_90azim)[Configdestinationcol_begin_90azim:Configdestinationcol_end_90azim] = c(paste("np",SampleLocPosition,sep = ""),paste("dist",SampleLocPosition,sep = ""),paste("gamma",SampleLocPosition,sep = ""))
}

# Eastward shifts of original
if(Eastwardshifts>0){
SampleConfig_Position_begins_EW = ((Upwardshifts-1)*(Northwardshifts+Southwradshifts)*(Eastwardshifts+Westwardshifts))+((Southwradshifts-1)*(Eastwardshifts+Westwardshifts))+(Westwardshifts-1)+1     +1
SampleConfig_Position_ends_EW = SampleConfig_Position_begins_EW+Eastwardshifts-1
for(t in 1:Eastwardshifts){
  SampleLocPosition = SampleConfig_Position_begins_EW+(t-1)
  ShiftedSampleLocation = data.frame(SampleLocation$x_index+t,SampleLocation$y_index,SampleLocation$z_index,SampleLocation$GridBlockID+t,SampleLocation$TrueEasting+(t*delta_X_big),SampleLocation$TrueNorthing,SampleLocation$TrueDepth)
  names(ShiftedSampleLocation) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth")
  Planarshifts[,((7*(Southwradshifts-1)*(Westwardshifts+Eastwardshifts))+(7*(Westwardshifts+(t-1)))+1):((7*(Southwradshifts-1)*(Westwardshifts+Eastwardshifts))+(7*(Westwardshifts+t)))] = ShiftedSampleLocation
  names(Planarshifts)[((7*(Southwradshifts-1)*(Westwardshifts+Eastwardshifts))+(7*(Westwardshifts+(t-1)))+1):((7*(Southwradshifts-1)*(Westwardshifts+Eastwardshifts))+(7*(Westwardshifts+t)))] = c(paste("x_index",SampleLocPosition, sep = ""),paste("y_index",SampleLocPosition, sep = ""),paste("z_index",SampleLocPosition, sep = ""),paste("GridBlockID",SampleLocPosition, sep = ""),paste("TrueEasting",SampleLocPosition, sep = ""),paste("TrueNorthing",SampleLocPosition, sep = ""),paste("TrueDepth",SampleLocPosition, sep = ""))
  Horizontalshifts[,(7*(Westwardshifts+t-1)+1):(7*(Westwardshifts+t))] = ShiftedSampleLocation 
  SampleData = data.frame(matrix(0,362,1)) 
  for (u in 1:362){
    SampleData[u,] = ChoiceRealization[as.numeric(rownames(ChoiceRealization))==ShiftedSampleLocation$GridBlockID[u],]
  }
  names(SampleData) = c("SampleData")
  assign(paste("SampleConfig",SampleLocPosition,sep = ""), data.frame(ShiftedSampleLocation,SampleData))
  ShiftedSampleLocanddata = data.frame(ShiftedSampleLocation,SampleData)
  coordinates(ShiftedSampleLocanddata) = ~TrueEasting+TrueNorthing+TrueDepth
  VarConfig_90azim = variogram(SampleData~1,ShiftedSampleLocanddata, alpha = 90, tol.hor = 22.5, cutoff = VarCutoff, boundaries = c(500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500))
  Configdestinationcol_begin_90azim = (3*SampleLocPosition)-2
  Configdestinationcol_end_90azim = 3*SampleLocPosition
  GrandVarConfig_90azim[Configdestinationcol_begin_90azim:Configdestinationcol_end_90azim] = data.frame(VarConfig_90azim$np,VarConfig_90azim$dist,VarConfig_90azim$gamma)
  names(GrandVarConfig_90azim)[Configdestinationcol_begin_90azim:Configdestinationcol_end_90azim] = c(paste("np",SampleLocPosition,sep = ""),paste("dist",SampleLocPosition,sep = ""),paste("gamma",SampleLocPosition,sep = ""))
  
 }
}


# Southward shifts of all horizontal shifts
if(Southwradshifts>1){
SampleConfig_Position_ends_SW = ((Upwardshifts-1)*(Northwardshifts+Southwradshifts)*(Eastwardshifts+Westwardshifts))+((Southwradshifts-1)*(Eastwardshifts+Westwardshifts)) 
SampleConfig_Position_begins_SW = SampleConfig_Position_ends_SW-(((Southwradshifts-1)*(Eastwardshifts+Westwardshifts))-1)
for(v in 1:(Southwradshifts-1)){
  for(w in 1:(Westwardshifts+Eastwardshifts)){
  SampleLocPosition = SampleConfig_Position_ends_SW-(((v-1)*(Westwardshifts+Eastwardshifts))+(w-1))
  HorizontalExtract = Horizontalshifts[,(7*(Westwardshifts+Eastwardshifts-w)+1):(7*(Westwardshifts+Eastwardshifts-(w-1)))] 
  names(HorizontalExtract) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth")
  ShiftedSampleLocation = data.frame(HorizontalExtract$x_index,HorizontalExtract$y_index-v,HorizontalExtract$z_index,HorizontalExtract$GridBlockID-nx,HorizontalExtract$TrueEasting,HorizontalExtract$TrueNorthing-(v*delta_Y_big),HorizontalExtract$TrueDepth)
  names(ShiftedSampleLocation) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth")
  Planarshifts[,((7*(((Southwradshifts-1-(v-1))*(Westwardshifts+Eastwardshifts))-w)+1)):(7*(((Southwradshifts-1-(v-1))*(Westwardshifts+Eastwardshifts))-(w-1)))] = ShiftedSampleLocation
  names(Planarshifts)[((7*(((Southwradshifts-1-(v-1))*(Westwardshifts+Eastwardshifts))-w)+1)):(7*(((Southwradshifts-1-(v-1))*(Westwardshifts+Eastwardshifts))-(w-1)))] = c(paste("x_index",SampleLocPosition, sep = ""),paste("y_index",SampleLocPosition, sep = ""),paste("z_index",SampleLocPosition, sep = ""),paste("GridBlockID",SampleLocPosition, sep = ""),paste("TrueEasting",SampleLocPosition, sep = ""),paste("TrueNorthing",SampleLocPosition, sep = ""),paste("TrueDepth",SampleLocPosition, sep = ""))
  SampleData = data.frame(matrix(0,362,1)) 
  for (a in 1:362){
    SampleData[a,] = ChoiceRealization[as.numeric(rownames(ChoiceRealization))==ShiftedSampleLocation$GridBlockID[a],]
   }
  names(SampleData) = c("SampleData")
  assign(paste("SampleConfig",SampleLocPosition,sep = ""), data.frame(ShiftedSampleLocation,SampleData))
  ShiftedSampleLocanddata = data.frame(ShiftedSampleLocation,SampleData)
  coordinates(ShiftedSampleLocanddata) = ~TrueEasting+TrueNorthing+TrueDepth
  VarConfig_90azim = variogram(SampleData~1,ShiftedSampleLocanddata, alpha = 90, tol.hor = 22.5, cutoff = VarCutoff, boundaries = c(500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500))
  Configdestinationcol_begin_90azim = (3*SampleLocPosition)-2
  Configdestinationcol_end_90azim = 3*SampleLocPosition
  GrandVarConfig_90azim[Configdestinationcol_begin_90azim:Configdestinationcol_end_90azim] = data.frame(VarConfig_90azim$np,VarConfig_90azim$dist,VarConfig_90azim$gamma)
  names(GrandVarConfig_90azim)[Configdestinationcol_begin_90azim:Configdestinationcol_end_90azim] = c(paste("np",SampleLocPosition,sep = ""),paste("dist",SampleLocPosition,sep = ""),paste("gamma",SampleLocPosition,sep = ""))
  }
 }
}


# Northward shifts of all horizontal shifts
if(Northwardshifts>0){
  SampleConfig_Position_begins_NW = ((Upwardshifts-1)*(Northwardshifts+Southwradshifts)*(Eastwardshifts+Westwardshifts))+((Southwradshifts-1)*(Eastwardshifts+Westwardshifts))+(Westwardshifts-1)+1   +1   +Eastwardshifts-1   +1 
  SampleConfig_Position_ends_NW = SampleConfig_Position_begins_NW+(Northwardshifts*(Westwardshifts+Eastwardshifts))-1
  for(b in 1:Northwardshifts){
    for(c in 1:(Westwardshifts+Eastwardshifts)){
    SampleLocPosition = SampleConfig_Position_begins_NW+(((b-1)*(Westwardshifts+Eastwardshifts))+(c-1))
    HorizontalExtract = Horizontalshifts[,(7*(c-1)+1):(7*c)] 
    names(HorizontalExtract) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth")
    ShiftedSampleLocation = data.frame(HorizontalExtract$x_index,HorizontalExtract$y_index+b,HorizontalExtract$z_index,HorizontalExtract$GridBlockID+nx,HorizontalExtract$TrueEasting,HorizontalExtract$TrueNorthing+(b*delta_Y_big),HorizontalExtract$TrueDepth)
    names(ShiftedSampleLocation) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth")
    Planarshifts[,((7*((Southwradshifts-1)   +1   +(b-1))*(Westwardshifts+Eastwardshifts))+    (7*(c-1))   +1):((7*((Southwradshifts-1)   +1   +(b-1))*(Westwardshifts+Eastwardshifts))+    (7*c))] = ShiftedSampleLocation
    names(Planarshifts)[((7*((Southwradshifts-1)   +1   +(b-1))*(Westwardshifts+Eastwardshifts))+    (7*(c-1))   +1):((7*((Southwradshifts-1)   +1   +(b-1))*(Westwardshifts+Eastwardshifts))+    (7*c))] = c(paste("x_index",SampleLocPosition, sep = ""),paste("y_index",SampleLocPosition, sep = ""),paste("z_index",SampleLocPosition, sep = ""),paste("GridBlockID",SampleLocPosition, sep = ""),paste("TrueEasting",SampleLocPosition, sep = ""),paste("TrueNorthing",SampleLocPosition, sep = ""),paste("TrueDepth",SampleLocPosition, sep = ""))
    SampleData = data.frame(matrix(0,362,1)) 
    for (d in 1:362){
      SampleData[d,] = ChoiceRealization[as.numeric(rownames(ChoiceRealization))==ShiftedSampleLocation$GridBlockID[d],]
    }
    names(SampleData) = c("SampleData")
    assign(paste("SampleConfig",SampleLocPosition,sep = ""), data.frame(ShiftedSampleLocation,SampleData))
    ShiftedSampleLocanddata = data.frame(ShiftedSampleLocation,SampleData)
    coordinates(ShiftedSampleLocanddata) = ~TrueEasting+TrueNorthing+TrueDepth
    VarConfig_90azim = variogram(SampleData~1,ShiftedSampleLocanddata, alpha = 90, tol.hor = 22.5, cutoff = VarCutoff, boundaries = c(500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500))
    Configdestinationcol_begin_90azim = (3*SampleLocPosition)-2
    Configdestinationcol_end_90azim = 3*SampleLocPosition
    GrandVarConfig_90azim[Configdestinationcol_begin_90azim:Configdestinationcol_end_90azim] = data.frame(VarConfig_90azim$np,VarConfig_90azim$dist,VarConfig_90azim$gamma)
    names(GrandVarConfig_90azim)[Configdestinationcol_begin_90azim:Configdestinationcol_end_90azim] = c(paste("np",SampleLocPosition,sep = ""),paste("dist",SampleLocPosition,sep = ""),paste("gamma",SampleLocPosition,sep = ""))
    }
  }
}

# Upward shifts of all Planar shifts
if(Upwardshifts>1){
  VerticalUpShifts = data.frame(matrix(0,362,(7*(Upwardshifts-1)*(Northwardshifts+Southwradshifts)*(Westwardshifts+Eastwardshifts))))
  SampleConfig_Position_ends_UW = ((Upwardshifts-1)*(Northwardshifts+Southwradshifts)*(Eastwardshifts+Westwardshifts)) 
  SampleConfig_Position_begins_UW = SampleConfig_Position_ends_UW-(((Upwardshifts-1)*(Northwardshifts+Southwradshifts)*(Eastwardshifts+Westwardshifts))-1)
  for(d in 1:(Upwardshifts-1)){
    for(e in 1:(Southwradshifts+Northwardshifts)){
      for(f in 1:(Westwardshifts+Eastwardshifts)){
      SampleLocPosition = SampleConfig_Position_ends_UW-(((d-1)*(Southwradshifts+Northwardshifts)*(Eastwardshifts+Westwardshifts))+((e-1)*(Westwardshifts+Eastwardshifts))  +(f-1))
      PlanarExtract = Planarshifts[,(7*(((Westwardshifts+Eastwardshifts)*(Southwradshifts+Northwardshifts))-((e-1)*(Westwardshifts+Eastwardshifts))-f)+1):(7*(((Westwardshifts+Eastwardshifts)*(Southwradshifts+Northwardshifts))-((e-1)*(Westwardshifts+Eastwardshifts))-(f-1)))] 
      names(PlanarExtract) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth")
      ShiftedSampleLocation = data.frame(PlanarExtract$x_index,PlanarExtract$y_index,PlanarExtract$z_index-d,PlanarExtract$GridBlockID-(nx*ny),PlanarExtract$TrueEasting,PlanarExtract$TrueNorthing,PlanarExtract$TrueDepth-(d*delta_Z_big))
      names(ShiftedSampleLocation) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth")
      VerticalUpShifts[,(7*(((Upwardshifts-1)*(Westwardshifts+Eastwardshifts)*(Southwradshifts+Northwardshifts))-((d-1)*(Westwardshifts+Eastwardshifts)*(Southwradshifts+Northwardshifts))-((e-1)*(Westwardshifts+Eastwardshifts))-f)+1):(7*(((Upwardshifts-1)*(Westwardshifts+Eastwardshifts)*(Southwradshifts+Northwardshifts))-((d-1)*(Westwardshifts+Eastwardshifts)*(Southwradshifts+Northwardshifts))-((e-1)*(Westwardshifts+Eastwardshifts))-(f-1)))] = ShiftedSampleLocation
      names(VerticalUpShifts)[(7*(((Upwardshifts-1)*(Westwardshifts+Eastwardshifts)*(Southwradshifts+Northwardshifts))-((d-1)*(Westwardshifts+Eastwardshifts)*(Southwradshifts+Northwardshifts))-((e-1)*(Westwardshifts+Eastwardshifts))-f)+1):(7*(((Upwardshifts-1)*(Westwardshifts+Eastwardshifts)*(Southwradshifts+Northwardshifts))-((d-1)*(Westwardshifts+Eastwardshifts)*(Southwradshifts+Northwardshifts))-((e-1)*(Westwardshifts+Eastwardshifts))-(f-1)))] = c(paste("x_index",SampleLocPosition, sep = ""),paste("y_index",SampleLocPosition, sep = ""),paste("z_index",SampleLocPosition, sep = ""),paste("GridBlockID",SampleLocPosition, sep = ""),paste("TrueEasting",SampleLocPosition, sep = ""),paste("TrueNorthing",SampleLocPosition, sep = ""),paste("TrueDepth",SampleLocPosition, sep = ""))
      SampleData = data.frame(matrix(0,362,1)) 
      for (g in 1:362){
        SampleData[g,] = ChoiceRealization[as.numeric(rownames(ChoiceRealization))==ShiftedSampleLocation$GridBlockID[g],]
      }
      names(SampleData) = c("SampleData")
      assign(paste("SampleConfig",SampleLocPosition,sep = ""), data.frame(ShiftedSampleLocation,SampleData))
      ShiftedSampleLocanddata = data.frame(ShiftedSampleLocation,SampleData)
      coordinates(ShiftedSampleLocanddata) = ~TrueEasting+TrueNorthing+TrueDepth
      VarConfig_90azim = variogram(SampleData~1,ShiftedSampleLocanddata, alpha = 90, tol.hor = 22.5, cutoff = VarCutoff, boundaries = c(500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500))
      Configdestinationcol_begin_90azim = (3*SampleLocPosition)-2
      Configdestinationcol_end_90azim = 3*SampleLocPosition
      GrandVarConfig_90azim[Configdestinationcol_begin_90azim:Configdestinationcol_end_90azim] = data.frame(VarConfig_90azim$np,VarConfig_90azim$dist,VarConfig_90azim$gamma)
      names(GrandVarConfig_90azim)[Configdestinationcol_begin_90azim:Configdestinationcol_end_90azim] = c(paste("np",SampleLocPosition,sep = ""),paste("dist",SampleLocPosition,sep = ""),paste("gamma",SampleLocPosition,sep = ""))
      }
   }
  }
}  


# Downward shifts of all Planar shifts
if(Downwardshifts>0){
  VerticalDownShifts = data.frame(matrix(0,362,(7*Downwardshifts*(Northwardshifts+Southwradshifts)*(Westwardshifts+Eastwardshifts))))
  SampleConfig_Position_begins_DW = ((Upwardshifts-1)*(Northwardshifts+Southwradshifts)*(Eastwardshifts+Westwardshifts))+((Southwradshifts-1)*(Eastwardshifts+Westwardshifts))+(Westwardshifts-1)+1   +1   +Eastwardshifts-1   +1     +(Northwardshifts*(Eastwardshifts+Westwardshifts))-1    +1 
  SampleConfig_Position_ends_DW = SampleConfig_Position_begins_DW+(Downwardshifts*(Northwardshifts+Southwradshifts)*(Westwardshifts+Eastwardshifts))-1
  for(h in 1:Downwardshifts){
    for(m in 1:(Southwradshifts+Northwardshifts)){
      for(n in 1:(Westwardshifts+Eastwardshifts)){
        SampleLocPosition = SampleConfig_Position_begins_DW+(((h-1)*(Southwradshifts+Northwardshifts)*(Eastwardshifts+Westwardshifts))+((m-1)*(Westwardshifts+Eastwardshifts))  +(n-1))
        PlanarExtract = Planarshifts[,(7*(((m-1)*(Westwardshifts+Eastwardshifts))+(n-1))+1):(7*(((m-1)*(Westwardshifts+Eastwardshifts))+n))] 
        names(PlanarExtract) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth")
        ShiftedSampleLocation = data.frame(PlanarExtract$x_index,PlanarExtract$y_index,PlanarExtract$z_index+h,PlanarExtract$GridBlockID+(nx*ny),PlanarExtract$TrueEasting,PlanarExtract$TrueNorthing,PlanarExtract$TrueDepth+(h*delta_Z_big))
        names(ShiftedSampleLocation) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth")
        VerticalDownShifts[,(7*(((h-1)*(Northwardshifts+Southwradshifts)*(Eastwardshifts+Westwardshifts))+((m-1)*(Westwardshifts+Eastwardshifts))+(n-1))+1):(7*(((h-1)*(Northwardshifts+Southwradshifts)*(Eastwardshifts+Westwardshifts))+((m-1)*(Westwardshifts+Eastwardshifts))+n))] = ShiftedSampleLocation
        names(VerticalDownShifts)[(7*(((h-1)*(Northwardshifts+Southwradshifts)*(Eastwardshifts+Westwardshifts))+((m-1)*(Westwardshifts+Eastwardshifts))+(n-1))+1):(7*(((h-1)*(Northwardshifts+Southwradshifts)*(Eastwardshifts+Westwardshifts))+((m-1)*(Westwardshifts+Eastwardshifts))+n))] = c(paste("x_index",SampleLocPosition, sep = ""),paste("y_index",SampleLocPosition, sep = ""),paste("z_index",SampleLocPosition, sep = ""),paste("GridBlockID",SampleLocPosition, sep = ""),paste("TrueEasting",SampleLocPosition, sep = ""),paste("TrueNorthing",SampleLocPosition, sep = ""),paste("TrueDepth",SampleLocPosition, sep = ""))
        SampleData = data.frame(matrix(0,362,1)) 
        for (p in 1:362){
          SampleData[p,] = ChoiceRealization[as.numeric(rownames(ChoiceRealization))==ShiftedSampleLocation$GridBlockID[p],]
        }
        names(SampleData) = c("SampleData")
        assign(paste("SampleConfig",SampleLocPosition,sep = ""), data.frame(ShiftedSampleLocation,SampleData))
        ShiftedSampleLocanddata = data.frame(ShiftedSampleLocation,SampleData)
        coordinates(ShiftedSampleLocanddata) = ~TrueEasting+TrueNorthing+TrueDepth
        VarConfig_90azim = variogram(SampleData~1,ShiftedSampleLocanddata, alpha = 90, tol.hor = 22.5, cutoff = VarCutoff, boundaries = c(500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500))
        Configdestinationcol_begin_90azim = (3*SampleLocPosition)-2
        Configdestinationcol_end_90azim = 3*SampleLocPosition
        GrandVarConfig_90azim[Configdestinationcol_begin_90azim:Configdestinationcol_end_90azim] = data.frame(VarConfig_90azim$np,VarConfig_90azim$dist,VarConfig_90azim$gamma)
        names(GrandVarConfig_90azim)[Configdestinationcol_begin_90azim:Configdestinationcol_end_90azim] = c(paste("np",SampleLocPosition,sep = ""),paste("dist",SampleLocPosition,sep = ""),paste("gamma",SampleLocPosition,sep = ""))
        }
    }
  }
}  

#Ploting the variogram family of curves
VargFamilyPlot = ggplot(data = GrandVarConfig_90azim)+
  labs(x = "Lag Distance, m", y = "Porosity Empirical Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.005))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.005, by = 0.001))
for(q in 1:TotalPossibleshifts){
  VargFamilyPlot = VargFamilyPlot+geom_line(aes_string(x = GrandVarConfig_90azim[,((3*q)-1)], y = GrandVarConfig_90azim[,3*q]), color = "red")
}
print(VargFamilyPlot)
ggsave("VargFamilyPlot.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)


for(ab in 1:nrow(GrandVarConfig_90azim)){
  LagVarg = data.frame(t(GrandVarConfig_90azim[ab,seq(from=3, to = 234, by = 3)])) 
  names(LagVarg) = c("Variog")
  ggplot(LagVarg, aes(x= Variog ))+
    geom_histogram(color = "blue", fill = "lightgreen")+
    labs(x = "Porosity Empirical Variogram", y = "Frequency")
  ggsave(paste(ab,"Histogram Plot.jpg",sep = " "), path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
}



#Aggregating all sample points configurations

AllConfig = data.frame(VerticalUpShifts,Planarshifts,VerticalDownShifts)
write.csv(AllConfig, file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis//All Configurations.csv", row.names = F)



# A smarter approach (April, 2019) - not really smarter (February, 2020)

#The idea here is to utilize the AllConfig data frame.
#Some nomenclature here:
#The entire possibilities is refered to as a Configuration.
#Each possibility is refered to as a Shifting. Shifting 1, Shifting 2 etc

#The approach is simply to tear off the relevant portion (7 columns) for each Shifting,
#and sampling the chosen realization by matching GridblockID of the Shifting to rownames of chosen realization;
#finally, variogram of the Shifting is calculated and stored to the relevant columns of the overall holder - formerly called GrandVarConfig_90azim


#Here we go


#AllShiftings = data.frame(VerticalUpShifts,Planarshifts,VerticalDownShifts)
AllShiftings = read.csv("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis//Grand Shifts.csv",header=T)


#Initializing placeholders (one for lag-interval approach, one for dbscan approach)
#for the outputs of variogram call across all possible Shiftings
AllShiftingsVariogramHolder_90azim_lagapproach = data.frame(matrix(0,9,3*TotalPossibleshifts)) #9 rows because we would set 9 boundaries; 3*Numberofpossibleshifts because we would store np, dist and gamma for each shifts
AllShiftingsVariogramHolder_90azim_dbscanapproach = data.frame(matrix(0,6,3*TotalPossibleshifts)) #6 rows because there would be 6 clusters; 3*Numberofpossibleshifts because we would store np, dist and gamma for each shifts
AllShiftingsmeanporoHolder= data.frame(matrix(0,TotalPossibleshifts,2))
names(AllShiftingsmeanporoHolder) = c("Shifting", "MeanPoro")

#Now looping for each shifting
for(jj in 1:TotalPossibleshifts){
  #First, tearing off the appropriate portion (7 columns) for each shifting
  ShiftSamplepoint = AllShiftings[,((7*jj)-6):(7*jj)]
  names(ShiftSamplepoint) = c("x_index","y_index","z_index","GridBlockID","TrueEasting","TrueNorthing", "TrueDepth")
  #Drawing sample from the best realization; made possible by selecting rows of ChosenRealization (Line 104) whose row names are same as GridBlockID of ShiftsamplePoint (Line 123)
  ShiftSampleData = data.frame(matrix(0,362,1)) 
  for (kk in 1:362){
    ShiftSampleData[kk,] = ChoiceRealization[as.numeric(rownames(ChoiceRealization))==ShiftSamplepoint$GridBlockID[kk],]
  }
  names(ShiftSampleData) = c("Porosity")
  #Append ShiftSampleData to ShiftSamplePoint to make a single dataframe
  ShiftSamplePointandData = data.frame(ShiftSamplepoint,ShiftSampleData)
  
  
    #Empirical variogram computations - lag-interval approach
  coordinates(ShiftSamplePointandData) = ~TrueEasting+TrueNorthing+TrueDepth
  ShiftVariog_90azim_lagapproach = variogram(Porosity~1,ShiftSamplePointandData, alpha = 90, tol.hor = 22.5, cutoff = VarCutoff, boundaries = c(500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500))
  
  #Send the output to the overall holder
  
  AllShiftingsVariogramHolder_90azim_lagapproach[((3*jj)-2):(3*jj)] = data.frame(ShiftVariog_90azim_lagapproach$np,ShiftVariog_90azim_lagapproach$dist,ShiftVariog_90azim_lagapproach$gamma)
  names(AllShiftingsVariogramHolder_90azim_lagapproach)[((3*jj)-2):(3*jj)] = c(paste("np",jj,sep = ""),paste("dist",jj,sep = ""),paste("gamma",jj,sep = ""))

  #Empirical variogram computations - dbscan approach
  ShiftingVariogramCloud_90_22.5 = variogram(Porosity~1,ShiftSamplePointandData, alpha = 90, tol.hor = 22.5, cutoff = VarCutoff, cloud = TRUE)
  Shiftingdbscan_90 = dbscan(ShiftingVariogramCloud_90_22.5[,c(2,3)], eps = 100, minPts = 500)
  dbscanShiftingVariogramCloud_90_22.5 = data.frame(ShiftingVariogramCloud_90_22.5$dist, ShiftingVariogramCloud_90_22.5$gamma, Shiftingdbscan_90$cluster, Direction = 90)    
  dbscanShiftingVariogramCloud_90_22.5$Shiftingdbscan_90.cluster = as.factor(dbscanShiftingVariogramCloud_90_22.5$Shiftingdbscan_90.cluster) #Ensuring R sees the cluster numbers as classification ID and not as numerals
  
  #ggplot(data = dbscanShiftingVariogramCloud_90_22.5[!dbscanShiftingVariogramCloud_90_22.5$Shiftingdbscan_90.cluster=="0",])+
    #geom_point(mapping = aes(x = ShiftingVariogramCloud_90_22.5.dist, y = ShiftingVariogramCloud_90_22.5.gamma, color = Shiftingdbscan_90.cluster, shape = Shiftingdbscan_90.cluster))+
    #labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Cluster" )+
    #coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
    #scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
    #scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))+
    #scale_color_manual(name = 'Cluster', values = c("#000000","#00FF00","#330033","#000099","#000033","#FF0000"))+
    #scale_shape_manual(name = 'Cluster', values = c(16,8,15,17,18,4))

  
  NumberofClusters_90 = length(levels(dbscanShiftingVariogramCloud_90_22.5$Shiftingdbscan_90.cluster))
  ShiftVariog_90azim_dbscanapproach = data.frame(matrix(0, NumberofClusters_90-1, 3, dimnames = list(c("1","2","3","4","5","6"),c("np", "dist", "gamma"))))
  for(mm in 1:NumberofClusters_90-1){
    ShiftVariog_90azim_dbscanapproach[mm,1]=length(dbscanShiftingVariogramCloud_90_22.5$ShiftingVariogramCloud_90_22.5.dist[dbscanShiftingVariogramCloud_90_22.5$Shiftingdbscan_90.cluster==mm]); 
    ShiftVariog_90azim_dbscanapproach[mm,2]=mean(dbscanShiftingVariogramCloud_90_22.5$ShiftingVariogramCloud_90_22.5.dist[dbscanShiftingVariogramCloud_90_22.5$Shiftingdbscan_90.cluster==mm]);
    ShiftVariog_90azim_dbscanapproach[mm,3]=mean(dbscanShiftingVariogramCloud_90_22.5$ShiftingVariogramCloud_90_22.5.gamma[dbscanShiftingVariogramCloud_90_22.5$Shiftingdbscan_90.cluster==mm])
  }

  #Send the output to the overall holder
  
  AllShiftingsVariogramHolder_90azim_dbscanapproach[((3*jj)-2):(3*jj)] = data.frame(ShiftVariog_90azim_dbscanapproach$np,ShiftVariog_90azim_dbscanapproach$dist,ShiftVariog_90azim_dbscanapproach$gamma)
  names(AllShiftingsVariogramHolder_90azim_dbscanapproach)[((3*jj)-2):(3*jj)] = c(paste("np",jj,sep = ""),paste("dist",jj,sep = ""),paste("gamma",jj,sep = ""))
  
}

#Variogram Uncertainty (mean and Variance) as a funtion of lag/cluster: computations and plots
MeanandVarianceofLagVariogram = data.frame(matrix(0,nrow(AllShiftingsVariogramHolder_90azim_lagapproach),3))
names(MeanandVarianceofLagVariogram) = c("dist", "meanofgamma", "varianceofgamma")
MeanandVarianceofClusterVariogram = data.frame(matrix(0,nrow(AllShiftingsVariogramHolder_90azim_dbscanapproach),3))
names(MeanandVarianceofClusterVariogram) = c("dist", "meanofgamma", "varianceofgamma")
for(mo in 1:nrow(AllShiftingsVariogramHolder_90azim_lagapproach)){
  MeanandVarianceofLagVariogram[mo,1] = AllShiftingsVariogramHolder_90azim_lagapproach[mo,2] 
  MeanandVarianceofLagVariogram[mo,2] = mean(t(AllShiftingsVariogramHolder_90azim_lagapproach[mo,seq(from=3, to = 234, by = 3)])) 
  MeanandVarianceofLagVariogram[mo,3] = var(t(AllShiftingsVariogramHolder_90azim_lagapproach[mo,seq(from=3, to = 234, by = 3)]))
}
for(so in 1:nrow(AllShiftingsVariogramHolder_90azim_dbscanapproach)){
  MeanandVarianceofClusterVariogram[so,1] = AllShiftingsVariogramHolder_90azim_dbscanapproach[so,2] 
  MeanandVarianceofClusterVariogram[so,2] = mean(t(AllShiftingsVariogramHolder_90azim_dbscanapproach[so,seq(from=3, to = 234, by = 3)])) 
  MeanandVarianceofClusterVariogram[so,3] = var(t(AllShiftingsVariogramHolder_90azim_dbscanapproach[so,seq(from=3, to = 234, by = 3)]))
}


#Ploting the variogram family of curves: lag approach
VarModelValues_90 = variogramLine(Integrated3Dmodel, VarCutoff, (VarCutoff/delta_X_big), dir=c(1,0,0), dist_vector = seq(1,8000,1))
VargFamilyPlot_lagapproach = ggplot(data = AllShiftingsVariogramHolder_90azim_lagapproach)+
  geom_point(data = MeanandVarianceofLagVariogram, aes(x = dist, y = meanofgamma, shape = "dc", size = "bi"), color = "blue")+
  geom_line(data = VarModelValues_90,aes(x = dist, y = gamma, color = "cu"), size = 0.75)+
  geom_point(aes(x = dist1, y = gamma1, shape = "pe", size = "sm"),color = "red")+
  labs(x = "Lag Distance, m", y = "Porosity Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))+
  scale_colour_manual(name = '', values =c("cu" = "black"), labels = c("cu" = "Model"))+
  scale_shape_manual(name = '', values = c('dc' = 8,'pe' = 16), labels = c('dc' = 'Mean','pe' = 'Estimates' ))+
  scale_size_manual(name = '', values = c('bi' = 4,'sm' = 1), labels = c('bi' = 'Mean','sm' = 'Estimates' ))+
  theme(legend.position=c(0.6, 0.2))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
for(nn in 1:TotalPossibleshifts){
  VargFamilyPlot_lagapproach = VargFamilyPlot_lagapproach+geom_point(aes_string(x = AllShiftingsVariogramHolder_90azim_lagapproach[,((3*nn)-1)], y = AllShiftingsVariogramHolder_90azim_lagapproach[,3*nn]), color = "red")
}
VargFamilyPlot_lagapproach
ggsave("VargFamilyPlot_lagapproach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)

#Histogram, Normal QQ plots - lag approach
NumberofLags = nrow(AllShiftingsVariogramHolder_90azim_lagapproach)
AllLagsVariog = data.frame(matrix(0,(NumberofLags*TotalPossibleshifts),2))
for(pp in 1:NumberofLags){
  LagVariogram = data.frame(t(AllShiftingsVariogramHolder_90azim_lagapproach[pp,seq(from=3, to = 234, by = 3)])) 
  names(LagVariogram) = c("Variog_Lagapproach")
  FreqVec = hist(LagVariogram$Variog_Lagapproach,plot=F)
  MaxFreq = max(FreqVec$counts)
  ggplot(LagVariogram, aes(x= Variog_Lagapproach ))+
    geom_histogram(bins = 15,color = "blue", fill = "lightgreen")+
    geom_vline(xintercept = mean(LagVariogram$Variog_Lagapproach),color = "red", size = 0.75)+
    annotate("text", x=(mean(LagVariogram$Variog_Lagapproach)+0.00005), y=MaxFreq, label="Mean", family="serif",fontface="italic", colour="blue", size=4)+
    labs(x = "Porosity Empirical Variogram", y = "Frequency")
    ggsave(paste("Lag",pp,"Histogram Plot.jpg",sep = " "), path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
    ggplot(LagVariogram)+
    geom_qq(aes(sample = Variog_Lagapproach), distribution = stats::qnorm, geom = "point",  color = "blue")+
    geom_qq_line(aes(sample = Variog_Lagapproach), distribution = stats::qnorm, line.p = c(0.25, 0.75), geom = "path",  color = "red")+
    labs(x = "Theoretical Normal Quantiles", y = "Porosity Variogram Quantiles")
  ggsave(paste("Lag",pp,"Normal QQ Plot.jpg",sep = " "), path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
  
  #for facet plot purposes
  AllLagsVariog[(((pp-1)*TotalPossibleshifts)+1):(pp*TotalPossibleshifts),1] = LagVariogram
  AllLagsVariog[(((pp-1)*TotalPossibleshifts)+1):(pp*TotalPossibleshifts),2] = paste("Lag", pp,sep = " ")
  names(AllLagsVariog) = c("Variog_Lagapproach","LagID")
}

#Faceting histogram and qqplots
ggplot(AllLagsVariog, aes(x= Variog_Lagapproach ))+
  geom_histogram(bins = 25,color = "blue", fill = "lightgreen")+
  #geom_vline(xintercept = mean(AllLagsVariog$Variog_Lagapproach),color = "red", size = 0.75)+
  labs(x = "Porosity Empirical Variogram", y = "Frequency")+
  facet_wrap(~ LagID, nrow = 3)
ggsave("All Non-ergodic Histograms - Lag approach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
ggplot(AllLagsVariog)+
  geom_qq(aes(sample = Variog_Lagapproach), distribution = stats::qnorm, geom = "point",  color = "blue")+
  geom_qq_line(aes(sample = Variog_Lagapproach), distribution = stats::qnorm, line.p = c(0.25, 0.75), geom = "path",  color = "red")+
  labs(x = "Theoretical Normal Quantiles", y = "Porosity Variogram Quantiles")+
  facet_wrap(~ LagID, nrow = 3)
ggsave("All Non-ergodic QQPlot - Lag approach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
  


#Ploting the variogram family of curves: cluster approach
VarModelValues_90 = variogramLine(Integrated3Dmodel, VarCutoff, (VarCutoff/delta_X_big), dir=c(1,0,0), dist_vector = seq(1,8000,1))
VargFamilyPlot_clusterapproach = ggplot(data = AllShiftingsVariogramHolder_90azim_dbscanapproach)+
  geom_point(data = MeanandVarianceofClusterVariogram, aes(x = dist, y = meanofgamma, shape = "dc", size = "bi"),color = "blue")+
  geom_line(data = VarModelValues_90,aes(x = dist, y = gamma, color = "cu"), size = 0.75)+
  geom_point(aes(x = dist1, y = gamma1, shape = "pe", size = "sm"),color = "red")+
  labs(x = "Lag Distance, m", y = "Porosity Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))+
  scale_colour_manual(name = '', values =c("cu" = "black"), labels = c("cu" = "Model"))+
  scale_shape_manual(name = '', values = c('dc' = 8,'pe' = 16), labels = c('dc' = 'Mean','pe' = 'Estimates' ))+
  scale_size_manual(name = '', values = c('bi' = 4,'sm' = 1), labels = c('bi' = 'Mean','sm' = 'Estimates' ))+
  theme(legend.position=c(0.6, 0.3))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
for(ss in 1:TotalPossibleshifts){
  VargFamilyPlot_clusterapproach = VargFamilyPlot_clusterapproach+geom_point(aes_string(x = AllShiftingsVariogramHolder_90azim_dbscanapproach[,((3*ss)-1)], y = AllShiftingsVariogramHolder_90azim_dbscanapproach[,3*ss]), color = "red")
}
VargFamilyPlot_clusterapproach
ggsave("VargFamilyPlot_clusterapproach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)


 #Histogram, Normal QQ plots - cluster approach
NumberofClusters = nrow(AllShiftingsVariogramHolder_90azim_dbscanapproach)
AllClustersVariog = data.frame(matrix(0,(NumberofClusters*TotalPossibleshifts),2))
for(rr in 1:NumberofClusters){
  ClusterVariogram = data.frame(t(AllShiftingsVariogramHolder_90azim_dbscanapproach[rr,seq(from=3, to = 234, by = 3)])) 
  names(ClusterVariogram) = c("Variog_dbscanapproach")
  FreqVec = hist(ClusterVariogram$Variog_dbscanapproach,plot=F)
  MaxFreq = max(FreqVec$counts)
  ggplot(ClusterVariogram, aes(x= Variog_dbscanapproach ))+
    geom_histogram(bins = 15, color = "blue", fill = "lightgreen")+
    geom_vline(xintercept = mean(ClusterVariogram$Variog_dbscanapproach),color = "red", size = 0.75)+
    annotate("text", x=(mean(ClusterVariogram$Variog_dbscanapproach)+0.00005), y=MaxFreq, label="Mean", family="serif",fontface="italic", colour="blue", size=4)+
    labs(x = "Porosity Empirical Variogram", y = "Frequency")
  ggsave(paste("Cluster",rr,"Histogram Plot.jpg",sep = " "), path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
  ggplot(ClusterVariogram)+
    geom_qq(aes(sample = Variog_dbscanapproach), distribution = stats::qnorm, geom = "point",  color = "blue")+
    geom_qq_line(aes(sample = Variog_dbscanapproach), distribution = stats::qnorm,line.p = c(0.25, 0.75), geom = "path",  color = "red")+
    labs(x = "Theoretical Normal Quantiles", y = "Porosity Variogram Quantiles")
  ggsave(paste("Cluster",rr,"Normal QQ Plot.jpg",sep = " "), path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
  
  #for facet plot purposes
  AllClustersVariog[(((rr-1)*TotalPossibleshifts)+1):(rr*TotalPossibleshifts),1] = ClusterVariogram
  AllClustersVariog[(((rr-1)*TotalPossibleshifts)+1):(rr*TotalPossibleshifts),2] = paste("Cluster", rr,sep = " ")
  names(AllClustersVariog) = c("Variog_dbscanapproach","ClusterID")
}
#Faceting histogram and qqplots
ggplot(AllClustersVariog, aes(x= Variog_dbscanapproach ))+
  geom_histogram(bins = 25,color = "blue", fill = "lightgreen")+
  #geom_vline(xintercept = mean(AllClustersVariog$Variog_dbscanapproach),color = "red", size = 0.75)+
  labs(x = "Porosity Empirical Variogram", y = "Frequency")+
  facet_wrap(~ ClusterID, nrow = 3)
ggsave("All Non-ergodic Histograms - DBSCAN approach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
ggplot(AllClustersVariog)+
  geom_qq(aes(sample = Variog_dbscanapproach), distribution = stats::qnorm, geom = "point",  color = "blue")+
  geom_qq_line(aes(sample = Variog_dbscanapproach), distribution = stats::qnorm, line.p = c(0.25, 0.75), geom = "path",  color = "red")+
  labs(x = "Theoretical Normal Quantiles", y = "Porosity Variogram Quantiles")+
  facet_wrap(~ ClusterID, nrow = 3)
ggsave("All Non-ergodic QQPlot - DBSCAN approach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)



#Plotting Variogram Uncertainty; i.e. variance of variogram
ggplot(data = MeanandVarianceofLagVariogram, aes(x = dist, y = varianceofgamma))+
  geom_line(aes(color = "ol"),size = 0.75)+
  geom_line(data = MeanandVarianceofClusterVariogram, aes(color = "at"), size = 0.75)+
  labs(x = "Lag Distance, m", y = "Variance of Porosity Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 1.0e-06))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 1.6e-06,by = 2.0e-07))+
  scale_colour_manual(name = '', values =c('ol'='red', "at" = "blue"), labels = c("ol" = "Lag interval approach","at" = "DBSCAN approach"))+
  theme(legend.position=c(0.6, 0.8))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
ggsave("Variance of Porosity Variogram.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)









#Now investigating the impact of data values (across realizations for the 
# original sample point configuration) on the variogram uncertainty

AllRealizationData = data.frame(matrix(0,362,number_of_anisorealizations))
for (un in 1:362){
AllRealizationData[un,] = Sim_big[as.numeric(rownames(Sim_big))==SampleLocation$GridBlockID[un],4:(3+number_of_anisorealizations)]
}


#Initializing placeholders (one for lag-interval approach, one for dbscan approach)
#for the outputs of variogram call across all realizations
AllRealizationsVariogramHolder_90azim_lagapproach = data.frame(matrix(0,9,3*number_of_anisorealizations)) #9 rows because we would set 9 boundaries; 3*Numberofpossibleshifts because we would store np, dist and gamma for each shifts
AllRealizationsVariogramHolder_90azim_dbscanapproach = data.frame(matrix(0,6,3*number_of_anisorealizations)) #6 rows because there would be 6 clusters; 3*Numberofpossibleshifts because we would store np, dist and gamma for each shifts


#Now looping for each realization
for(de in 1:number_of_anisorealizations){
  #First, tearing off the appropriate portion (1 column) of the data, for each realization
  RealizationSampleData = data.frame(AllRealizationData[,de])
  names(RealizationSampleData) = c("Porosity")
  
  #Append RealizationSampleData to Location coordinates (from SampleLocation) to make a single dataframe
  RealizationSamplePointandData = data.frame(SampleLocation[,5:7],RealizationSampleData)
  
  #Empirical variogram computations - lag-interval approach
  coordinates(RealizationSamplePointandData) = ~TrueEasting+TrueNorthing+TrueDepth
  RealizationVariog_90azim_lagapproach = variogram(Porosity~1,RealizationSamplePointandData, alpha = 90, tol.hor = 22.5, cutoff = VarCutoff, boundaries = c(500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500))
  
  #Send the output to the overall holder
  
  AllRealizationsVariogramHolder_90azim_lagapproach[((3*de)-2):(3*de)] = data.frame(RealizationVariog_90azim_lagapproach$np,RealizationVariog_90azim_lagapproach$dist,RealizationVariog_90azim_lagapproach$gamma)
  names(AllRealizationsVariogramHolder_90azim_lagapproach)[((3*de)-2):(3*de)] = c(paste("np",de,sep = ""),paste("dist",de,sep = ""),paste("gamma",de,sep = ""))
  
  #Empirical variogram computations - dbscan approach
  RealizationVariogramCloud_90_22.5 = variogram(Porosity~1,RealizationSamplePointandData, alpha = 90, tol.hor = 22.5, cutoff = VarCutoff, cloud = TRUE)
  Realizationdbscan_90 = dbscan(RealizationVariogramCloud_90_22.5[,c(2,3)], eps = 100, minPts = 500)
  dbscanRealizationVariogramCloud_90_22.5 = data.frame(RealizationVariogramCloud_90_22.5$dist, RealizationVariogramCloud_90_22.5$gamma, Realizationdbscan_90$cluster, Direction = 90)    
  dbscanRealizationVariogramCloud_90_22.5$Realizationdbscan_90.cluster = as.factor(dbscanRealizationVariogramCloud_90_22.5$Realizationdbscan_90.cluster) #Ensuring R sees the cluster numbers as classification ID and not as numerals
  
  #ggplot(data = dbscanRealizationVariogramCloud_90_22.5[!dbscanRealizationVariogramCloud_90_22.5$Realizationdbscan_90.cluster=="0",])+
  #geom_point(mapping = aes(x = RealizationVariogramCloud_90_22.5.dist, y = RealizationVariogramCloud_90_22.5.gamma, color = Realizationdbscan_90.cluster, shape = Realizationdbscan_90.cluster))+
  #labs(x = "Lag Distance, m", y = "Porosity Pair Variogram", colour = "Cluster" )+
  #coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.05))+
  #scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  #scale_y_continuous(breaks = seq(0, 0.05, by = 0.01))+
  #scale_color_manual(name = 'Cluster', values = c("#000000","#00FF00","#330033","#000099","#000033","#FF0000"))+
  #scale_shape_manual(name = 'Cluster', values = c(16,8,15,17,18,4))
  
  
  NumberofClusters_90 = length(levels(dbscanRealizationVariogramCloud_90_22.5$Realizationdbscan_90.cluster))
  RealizationVariog_90azim_dbscanapproach = data.frame(matrix(0, NumberofClusters_90-1, 3, dimnames = list(c("1","2","3","4","5","6"),c("np", "dist", "gamma"))))
  for(pr in 1:NumberofClusters_90-1){
    RealizationVariog_90azim_dbscanapproach[pr,1]=length(dbscanRealizationVariogramCloud_90_22.5$RealizationVariogramCloud_90_22.5.dist[dbscanRealizationVariogramCloud_90_22.5$Realizationdbscan_90.cluster==pr]); 
    RealizationVariog_90azim_dbscanapproach[pr,2]=mean(dbscanRealizationVariogramCloud_90_22.5$RealizationVariogramCloud_90_22.5.dist[dbscanRealizationVariogramCloud_90_22.5$Realizationdbscan_90.cluster==pr]);
    RealizationVariog_90azim_dbscanapproach[pr,3]=mean(dbscanRealizationVariogramCloud_90_22.5$RealizationVariogramCloud_90_22.5.gamma[dbscanRealizationVariogramCloud_90_22.5$Realizationdbscan_90.cluster==pr])
  }
 
  #Send the output to the overall holder
  
  AllRealizationsVariogramHolder_90azim_dbscanapproach[((3*de)-2):(3*de)] = data.frame(RealizationVariog_90azim_dbscanapproach$np,RealizationVariog_90azim_dbscanapproach$dist,RealizationVariog_90azim_dbscanapproach$gamma)
  names(AllRealizationsVariogramHolder_90azim_dbscanapproach)[((3*de)-2):(3*de)] = c(paste("np",de,sep = ""),paste("dist",de,sep = ""),paste("gamma",de,sep = ""))
  
}

#Variogram Uncertainty (mean and Variance) as a funtion of lag/cluster: computations and plots
Realization_MeanandVarianceofLagVariogram = data.frame(matrix(0,nrow(AllRealizationsVariogramHolder_90azim_lagapproach),3))
names(Realization_MeanandVarianceofLagVariogram) = c("dist", "meanofgamma", "varianceofgamma")
Realization_MeanandVarianceofClusterVariogram = data.frame(matrix(0,nrow(AllRealizationsVariogramHolder_90azim_dbscanapproach),3))
names(Realization_MeanandVarianceofClusterVariogram) = c("dist", "meanofgamma", "varianceofgamma")
for(ai in 1:nrow(AllRealizationsVariogramHolder_90azim_lagapproach)){
  Realization_MeanandVarianceofLagVariogram[ai,1] = AllRealizationsVariogramHolder_90azim_lagapproach[ai,2] 
  Realization_MeanandVarianceofLagVariogram[ai,2] = mean(t(AllRealizationsVariogramHolder_90azim_lagapproach[ai,seq(from=3, to = 3*number_of_anisorealizations, by = 3)])) 
  Realization_MeanandVarianceofLagVariogram[ai,3] = var(t(AllRealizationsVariogramHolder_90azim_lagapproach[ai,seq(from=3, to = 3*number_of_anisorealizations, by = 3)]))
}
for(se in 1:nrow(AllRealizationsVariogramHolder_90azim_dbscanapproach)){
  Realization_MeanandVarianceofClusterVariogram[se,1] = AllRealizationsVariogramHolder_90azim_dbscanapproach[se,2] 
  Realization_MeanandVarianceofClusterVariogram[se,2] = mean(t(AllRealizationsVariogramHolder_90azim_dbscanapproach[se,seq(from=3, to = 3*number_of_anisorealizations, by = 3)])) 
  Realization_MeanandVarianceofClusterVariogram[se,3] = var(t(AllRealizationsVariogramHolder_90azim_dbscanapproach[se,seq(from=3, to = 3*number_of_anisorealizations, by = 3)]))
}


#Ploting the variogram family of curves: lag approach
VarModelValues_90 = variogramLine(Integrated3Dmodel, VarCutoff, (VarCutoff/delta_X_big), dir=c(1,0,0), dist_vector = seq(1,8000,1))
Realization_VargFamilyPlot_lagapproach = ggplot(data = AllRealizationsVariogramHolder_90azim_lagapproach)+
  geom_point(data = Realization_MeanandVarianceofLagVariogram, aes(x = dist, y = meanofgamma, shape = "dc", size = "bi"), color = "blue")+
  geom_line(data = VarModelValues_90,aes(x = dist, y = gamma, color = "cu"), size = 0.75)+
  geom_point(aes(x = dist1, y = gamma1, shape = "pe", size = "sm"),color = "red")+
  labs(x = "Lag Distance, m", y = "Porosity Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))+
  scale_colour_manual(name = '', values =c("cu" = "black"), labels = c("cu" = "Model"))+
  scale_shape_manual(name = '', values = c('dc' = 8,'pe' = 16), labels = c('dc' = 'Mean','pe' = 'Estimates' ))+
  scale_size_manual(name = '', values = c('bi' = 4,'sm' = 1), labels = c('bi' = 'Mean','sm' = 'Estimates' ))+
  theme(legend.position=c(0.6, 0.3))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
for(yo in 1:number_of_anisorealizations){
  Realization_VargFamilyPlot_lagapproach = Realization_VargFamilyPlot_lagapproach+geom_point(aes_string(x = AllRealizationsVariogramHolder_90azim_lagapproach[,((3*yo)-1)], y = AllRealizationsVariogramHolder_90azim_lagapproach[,3*yo]), color = "red")
}
Realization_VargFamilyPlot_lagapproach
ggsave("Realization_VargFamilyPlot_lagapproach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)

#Histogram, Normal QQ plots - lag approach
NumberofLags = nrow(AllRealizationsVariogramHolder_90azim_lagapproach)
Realization_AllLagsVariog = data.frame(matrix(0,(NumberofLags*number_of_anisorealizations),2))
for(nd in 1:NumberofLags){
  Realization_LagVariogram = data.frame(t(AllRealizationsVariogramHolder_90azim_lagapproach[nd,seq(from=3, to = 3*number_of_anisorealizations, by = 3)])) 
  names(Realization_LagVariogram) = c("Variog_Lagapproach")
  FreqVec = hist(Realization_LagVariogram$Variog_Lagapproach,plot=F)
  MaxFreq = max(FreqVec$counts)
  ggplot(Realization_LagVariogram, aes(x= Variog_Lagapproach ))+
    geom_histogram(bins = 15,color = "blue", fill = "lightgreen")+
    geom_vline(xintercept = mean(Realization_LagVariogram$Variog_Lagapproach),color = "red", size = 0.75)+
    annotate("text", x=mean(Realization_LagVariogram$Variog_Lagapproach), y=MaxFreq, label="Mean", family="serif",fontface="italic", colour="blue", size=4)+
    labs(x = "Porosity Empirical Variogram", y = "Frequency")
  ggsave(paste("Realization_Lag",nd,"Histogram Plot.jpg",sep = " "), path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
  ggplot(Realization_LagVariogram)+
    geom_qq(aes(sample = Variog_Lagapproach), distribution = stats::qnorm, geom = "point",  color = "blue")+
    geom_qq_line(aes(sample = Variog_Lagapproach), distribution = stats::qnorm, line.p = c(0.25, 0.75), geom = "path",  color = "red")+
    labs(x = "Theoretical Normal Quantiles", y = "Porosity Variogram Quantiles")
  ggsave(paste("Realization_Lag",nd,"Normal QQ Plot.jpg",sep = " "), path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
  
  #for facet plot purposes
  Realization_AllLagsVariog[(((nd-1)*number_of_anisorealizations)+1):(nd*number_of_anisorealizations),1] = Realization_LagVariogram
  Realization_AllLagsVariog[(((nd-1)*number_of_anisorealizations)+1):(nd*number_of_anisorealizations),2] = paste("Lag", nd,sep = " ")
  names(Realization_AllLagsVariog) = c("Variog_Lagapproach","LagID") 
}

#Faceting histogram and qqplots
ggplot(Realization_AllLagsVariog, aes(x= Variog_Lagapproach ))+
  geom_histogram(bins = 25,color = "blue", fill = "lightgreen")+
  #geom_vline(xintercept = mean(Realization_AllLagsVariog$Variog_Lagapproach),color = "red", size = 0.75)+
  labs(x = "Porosity Empirical Variogram", y = "Frequency")+
  facet_wrap(~LagID, nrow = 3)
ggsave("All Ergodic Histograms - Lag approach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
ggplot(Realization_AllLagsVariog)+
  geom_qq(aes(sample = Variog_Lagapproach), distribution = stats::qnorm, geom = "point",  color = "blue")+
  geom_qq_line(aes(sample = Variog_Lagapproach), distribution = stats::qnorm, line.p = c(0.25, 0.75), geom = "path",  color = "red")+
  labs(x = "Theoretical Normal Quantiles", y = "Porosity Variogram Quantiles")+
  facet_wrap(~LagID, nrow = 3)
ggsave("All Ergodic QQPlot - Lag approach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)



#Ploting the variogram family of curves: cluster approach
VarModelValues_90 = variogramLine(Integrated3Dmodel, VarCutoff, (VarCutoff/delta_X_big), dir=c(1,0,0), dist_vector = seq(1,8000,1))
Realization_VargFamilyPlot_clusterapproach = ggplot(data = AllRealizationsVariogramHolder_90azim_dbscanapproach)+
  geom_point(data = Realization_MeanandVarianceofClusterVariogram, aes(x = dist, y = meanofgamma, shape = "dc", size = "bi"), color = "blue")+
  geom_line(data = VarModelValues_90,aes(x = dist, y = gamma, color = "cu"), size = 0.75)+
  geom_point(aes(x = dist1, y = gamma1, shape = "pe", size = "sm"),color = "red")+
  labs(x = "Lag Distance, m", y = "Porosity Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))+
  scale_colour_manual(name = '', values =c("cu" = "black"), labels = c("cu" = "Model"))+
  scale_shape_manual(name = '', values = c('dc' = 8,'pe' = 16), labels = c('dc' = 'Mean','pe' = 'Estimates' ))+
  scale_size_manual(name = '', values = c('bi' = 4,'sm' = 1), labels = c('bi' = 'Mean','sm' = 'Estimates' ))+
  theme(legend.position=c(0.6, 0.3))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
for(af in 1:number_of_anisorealizations){
  Realization_VargFamilyPlot_clusterapproach = Realization_VargFamilyPlot_clusterapproach+geom_point(aes_string(x = AllShiftingsVariogramHolder_90azim_dbscanapproach[,((3*af)-1)], y = AllShiftingsVariogramHolder_90azim_dbscanapproach[,3*af]), color = "red")
}
Realization_VargFamilyPlot_clusterapproach
ggsave("Realization_VargFamilyPlot_clusterapproach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)


#Histogram, Normal QQ plots - cluster approach
NumberofClusters = nrow(AllRealizationsVariogramHolder_90azim_dbscanapproach)
Realization_AllClustersVariog = data.frame(matrix(0,(NumberofClusters*number_of_anisorealizations),2))
for(lu in 1:nrow(AllRealizationsVariogramHolder_90azim_dbscanapproach)){
  Realization_ClusterVariogram = data.frame(t(AllRealizationsVariogramHolder_90azim_dbscanapproach[lu,seq(from=3, to = 3*number_of_anisorealizations, by = 3)])) 
  names(Realization_ClusterVariogram) = c("Variog_dbscanapproach")
  FreqVec = hist(Realization_ClusterVariogram$Variog_dbscanapproach,plot=F)
  MaxFreq = max(FreqVec$counts)
  ggplot(Realization_ClusterVariogram, aes(x= Variog_dbscanapproach ))+
    geom_histogram(bins = 15, color = "blue", fill = "lightgreen")+
    geom_vline(xintercept = mean(Realization_ClusterVariogram$Variog_dbscanapproach),color = "red", size = 0.75)+
    annotate("text", x=(mean(Realization_ClusterVariogram$Variog_dbscanapproach)+0.00005), y=MaxFreq, label="Mean", family="serif",fontface="italic", colour="blue", size=4)+
    labs(x = "Porosity Empirical Variogram", y = "Frequency")
  ggsave(paste("Realization_Cluster",lu,"Histogram Plot.jpg",sep = " "), path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
  ggplot(Realization_ClusterVariogram)+
    geom_qq(aes(sample = Variog_dbscanapproach), distribution = stats::qnorm, geom = "point",  color = "blue")+
    geom_qq_line(aes(sample = Variog_dbscanapproach), distribution = stats::qnorm,line.p = c(0.25, 0.75), geom = "path",  color = "red")+
    labs(x = "Theoretical Normal Quantiles", y = "Porosity Variogram Quantiles")
  ggsave(paste("Realization_Cluster",lu,"Normal QQ Plot.jpg",sep = " "), path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
  
  #for facet plot purposes
  Realization_AllClustersVariog[(((lu-1)*number_of_anisorealizations)+1):(lu*number_of_anisorealizations),1] = Realization_ClusterVariogram
  Realization_AllClustersVariog[(((lu-1)*number_of_anisorealizations)+1):(lu*number_of_anisorealizations),2] = paste("Cluster", lu,sep = " ")
  names(Realization_AllClustersVariog) = c("Variog_dbscanapproach","ClusterID")  
}

#Faceting histogram and qqplots
ggplot(Realization_AllClustersVariog, aes(x= Variog_dbscanapproach ))+
  geom_histogram(bins = 25,color = "blue", fill = "lightgreen")+
  #geom_vline(xintercept = mean(Realization_AllClustersVariog$Variog_dbscanapproach),color = "red", size = 0.75)+
  labs(x = "Porosity Empirical Variogram", y = "Frequency")+
  facet_wrap(~ClusterID, nrow = 3)
ggsave("All Ergodic Histograms - DBSCAN approach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)
ggplot(Realization_AllClustersVariog)+
  geom_qq(aes(sample = Variog_dbscanapproach), distribution = stats::qnorm, geom = "point",  color = "blue")+
  geom_qq_line(aes(sample = Variog_dbscanapproach), distribution = stats::qnorm, line.p = c(0.25, 0.75), geom = "path",  color = "red")+
  labs(x = "Theoretical Normal Quantiles", y = "Porosity Variogram Quantiles")+
  facet_wrap(~ ClusterID, nrow = 3)
ggsave("All Ergodic QQPlot - DBSCAN approach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)




#Plotting Variogram Uncertainty; i.e. variance of variogram
ggplot(data = Realization_MeanandVarianceofLagVariogram, aes(x = dist, y = varianceofgamma))+
  geom_line(aes(color = "ol"),size = 0.75)+
  geom_line(data = Realization_MeanandVarianceofClusterVariogram, aes(color = "at"), size = 0.75)+
  labs(x = "Lag Distance, m", y = "Variance of Porosity Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 1.6e-06))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 1.6e-06,by = 2.0e-07))+
  scale_colour_manual(name = '', values =c('ol'='red', "at" = "blue"), labels = c("ol" = "Lag interval approach","at" = "DBSCAN approach"))+
  theme(legend.position=c(0.6, 0.8))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
ggsave("Realization_Variance of Porosity Variogram.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)





##Comparing location effects (Non-ergodic) and value effects (Ergodic)

#Lag approach

ggplot(data = MeanandVarianceofLagVariogram, aes(x = dist, y = varianceofgamma))+
  geom_line(aes(color = "ol"),size = 0.75)+
  geom_point(color = "red")+
  geom_line(data = Realization_MeanandVarianceofLagVariogram, aes(color = "at"), size = 0.75)+
  geom_point(data = Realization_MeanandVarianceofLagVariogram, color = "blue")+
  labs(x = "Lag Distance, m", y = "Variance of Porosity Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 2.5e-06))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 2.5e-06,by = 5.0e-07))+
  scale_colour_manual(name = '', values =c('ol'='red', "at" = "blue"), labels = c("ol" = "Non-ergodic","at" = "Ergodic"))+
  theme(legend.position=c(0.8, 0.8))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
ggsave("Variance of Porosity Variogram_Compare_lagapproach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)



#Cluster approach
ggplot(data = MeanandVarianceofClusterVariogram, aes(x = dist, y = varianceofgamma))+
  geom_line(aes(color = "ol"),size = 0.75)+
  geom_point(color = "red")+
  geom_line(data = Realization_MeanandVarianceofClusterVariogram, aes(color = "at"), size = 0.75)+
  geom_point(data = Realization_MeanandVarianceofClusterVariogram, color = "blue")+
  labs(x = "Lag Distance, m", y = "Variance of Porosity Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 6.0e-07))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 6.0e-07,by = 2.0e-07))+
  scale_colour_manual(name = '', values =c('ol'='red', "at" = "blue"), labels = c("ol" = "Non-ergodic","at" = "Ergodic"))+
  theme(legend.position=c(0.8, 0.2))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
ggsave("Variance of Porosity Variogram_Compare_clusterapproach.jpg", path = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Spatial Analysis//Agbabu Spatial Analysis using gstat//Uncertainty Analysis", dpi = 96)


# To the Only Wise God - TTOWG


#Agbabu Bitumen Deposit Porosity Data Analysis

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



#Porosity Analysis - Whole
   #whole-individuals analysis
Mean_Whole = mean(TrimmedData$Porosity,na.rm = T) #Mean
Median_whole = median(TrimmedData$Porosity,na.rm = T) #Median
StandardDeviation_whole = sd(TrimmedData$Porosity, na.rm = T) #Standard deviation
library(moments) #To enable the computation of kurtosis and skewness
Kurtosis_whole = kurtosis(TrimmedData$Porosity, na.rm = T) #Kurtosis
Skewness_whole = skewness(TrimmedData$Porosity, na.rm = T) #Skewness
Minimum_whole = min(TrimmedData$Porosity,na.rm = T) #Minimum
First_Quartile_whole = quantile(TrimmedData$Porosity,probs = 0.25, na.rm = T) #First Quartile
Third_Quartile_whole = quantile(TrimmedData$Porosity,probs = 0.75, na.rm = T) #Third Quartile
Maximum_whole = max(TrimmedData$Porosity, na.rm = T) #Maximum
Count_whole = length(TrimmedData$Porosity) #Number of data values

Data_Attributes = c("Mean", "Median", "Standard Deviation", "Kurtosis", "Skewness", "Minimum", "First Quartile", "Third Quartile", "Maximum", "Number of data values") #Summary Table row names
Value_whole_individuals = c(Mean_Whole, Median_whole, StandardDeviation_whole, Kurtosis_whole, Skewness_whole, Minimum_whole, First_Quartile_whole, Third_Quartile_whole, Maximum_whole, Count_whole) #Summary Table row values 
SummaryTable_whole_individuals = data.frame(Data_Attributes,Value_whole_individuals) #Constructing the Summary Table
names(SummaryTable_whole_individuals) = c("Data Attribute","Attribute Value for all Porosity values")
write.table(SummaryTable_whole_individuals,file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Tables//Whole Analysis//SummaryTable_whole_individuals.csv",sep = ",",row.names = F)

jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//Whole Analysis//whole-individual boxplot.jpg")
boxplot(TrimmedData$Porosity, main="Boxplot of all Porosity values",ylab= "Porosity (fraction)",las=1)
dev.off()

jpeg("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//Whole Analysis//whole-individual RF-N-K Distributions.jpg")
porohist_whole=hist(TrimmedData$Porosity,plot=F, breaks = seq(from = 0, to = 0.5, by = 0.05))
porohist_whole$counts = porohist_whole$counts/sum(porohist_whole$counts)
porodens_whole = density(TrimmedData$Porosity,na.rm = T,from = 0, to=0.5)
porodens_whole$y=porodens_whole$y*0.05
x = TrimmedData$Porosity
poronorm_whole = dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole)
poronorm_whole = poronorm_whole*0.05
plot(porohist_whole, freq = T, col = "grey", ylim =  range(0,porohist_whole$counts,porodens_whole$y,poronorm_whole), main="Relative Frequency, Normal and Kernel Distributions of all Porosity values",xlab="Porosity (fraction)",ylab = "Probability", las=1)
#text(0.4, 0.2, paste("Mean =", round(Mean_Whole, 4), "\n Std.Dev =", round(StandardDeviation_whole, 4)))
lines(porodens_whole, col="blue", type = "l", lty="solid", lwd=2)
curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole)*0.05,add = T, col="red", type = "l", lty="dashed", lwd=2)
legend("topright",legend=c("Kernel","Normal"),
       text.col=c("blue","red"),lty=c("solid","dashed"),col=c("blue","red"))
#curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole),add = T)
dev.off()



   #Borehole-averages analysis
write.csv(TrimmedData$BoreHole,"allBoreHoleNodes.csv")
allBoreHoleNodes = read.csv("allBoreHoleNodes.csv")
ListofBorehole = levels(as.factor(allBoreHoleNodes$x))
BoreHoleavgs = rep(0,length(ListofBorehole))
BoreHolestd_dev = rep(0,length(ListofBorehole))
for(i in 1:length(ListofBorehole)) {BoreHoleavgs[i]=mean(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole[i]])
BoreHolestd_dev[i]=sd(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole[i]])}
WellAveragesTable = data.frame(ListofBorehole,BoreHoleavgs)
WellStandardDeviationTable = data.frame(ListofBorehole,BoreHolestd_dev)

Mean_Whole_wellaverages = mean(WellAveragesTable$BoreHoleavgs,na.rm = T) #Mean
Median_Whole_wellaverages = median(WellAveragesTable$BoreHoleavgs,na.rm = T) #Median
StandardDeviation_Whole_wellaverages = sd(WellAveragesTable$BoreHoleavgs, na.rm = T) #Standard deviation
Kurtosis_Whole_wellaverages = kurtosis(WellAveragesTable$BoreHoleavgs, na.rm = T) #Kurtosis
Skewness_Whole_wellaverages = skewness(WellAveragesTable$BoreHoleavgs, na.rm = T) #Skewness
Minimum_Whole_wellaverages = min(WellAveragesTable$BoreHoleavgs,na.rm = T) #Minimum
First_Quartile_Whole_wellaverages = quantile(WellAveragesTable$BoreHoleavgs,probs = 0.25, na.rm = T) #First Quartile
Third_Quartile_Whole_wellaverages = quantile(WellAveragesTable$BoreHoleavgs,probs = 0.75, na.rm = T) #Third Quartile
Maximum_Whole_wellaverages = max(WellAveragesTable$BoreHoleavgs, na.rm = T) #Maximum
Count_Whole_wellaverages = length(WellAveragesTable$BoreHoleavgs) #Number of data values

Data_Attributes_Whole_wellaverages = c("Mean", "Median", "Standard Deviation", "Kurtosis", "Skewness", "Minimum", "First Quartile", "Third Quartile", "Maximum", "Number of data values") #Summary Table row names
Value_Attributes_Whole_wellaverages = c(Mean_Whole_wellaverages, Median_Whole_wellaverages, StandardDeviation_Whole_wellaverages, Kurtosis_Whole_wellaverages, Skewness_Whole_wellaverages, Minimum_Whole_wellaverages, First_Quartile_Whole_wellaverages, Third_Quartile_Whole_wellaverages, Maximum_Whole_wellaverages, Count_Whole_wellaverages) #Summary Table row values
SummaryTableWhole_wellaverages = data.frame(Data_Attributes_Whole_wellaverages,Value_Attributes_Whole_wellaverages) #Constructing the Summary Table
names(SummaryTableWhole_wellaverages) = c("Data Attribute","Attribute Value for Borehole Porosity Averages")
write.table(SummaryTableWhole_wellaverages,file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Tables//Whole Analysis//SummaryTableWhole_wellaverages.csv",sep = ",",row.names = F)

jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//Whole Analysis//whole-wellaverages boxplot.jpg")
boxplot(WellAveragesTable$BoreHoleavgs, main="Boxplot of Borehole Porosity Averages",ylab = "Mean Porosity (fraction)",las=1)
dev.off()

jpeg("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//Whole Analysis//whole-wellaverages RF-N-K Distributions.jpg")
porohist_whole_wellaverages=hist(WellAveragesTable$BoreHoleavgs,plot=F, breaks = seq(from = 0, to = 0.5, by = 0.05))
porohist_whole_wellaverages$counts = porohist_whole_wellaverages$counts/sum(porohist_whole_wellaverages$counts)
porodens_whole_wellaverages = density(WellAveragesTable$BoreHoleavgs,na.rm = T,from = 0, to=0.5)
porodens_whole_wellaverages$y=porodens_whole_wellaverages$y*0.05
x = WellAveragesTable$BoreHoleavgs
poronorm_whole_wellaverages = dnorm(x,mean = Mean_Whole_wellaverages,sd=StandardDeviation_Whole_wellaverages)
poronorm_whole_wellaverages = poronorm_whole_wellaverages*0.05
plot(porohist_whole_wellaverages, freq = T, col = "grey", ylim =  range(0,porohist_whole_wellaverages$counts,porodens_whole_wellaverages$y,poronorm_whole_wellaverages), main="Relative Frequency, Normal and Kernel Distributions of Borehole Porosity Averages",xlab="Mean Porosity (fraction)",ylab = "Probability", las=1)
#text(0.4, 0.2, paste("Mean =", round(Mean_Whole, 4), "\n Std.Dev =", round(StandardDeviation_whole, 4)))
lines(porodens_whole_wellaverages, col="blue", type = "l", lty="solid", lwd=2)
curve(dnorm(x,mean = Mean_Whole_wellaverages,sd=StandardDeviation_Whole_wellaverages)*0.05,add = T, col="red", type = "l", lty="dashed", lwd=2)
legend("topright",legend=c("Kernel","Normal"),
       text.col=c("blue","red"),lty=c("solid","dashed"),col=c("blue","red"))
#curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole),add = T)
dev.off()




#jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data//Data - Improved Estimation//For export to R//Plots//Whole Analysis//whole-wellaverages histogram.jpg")
#hist(WellAveragesTable$BoreHoleavgs, main="Histogram of Borehole Porosity Averages",xlab="Mean Porosity",las=1)
#dev.off()

#Porosity Analysis - Horizon
    #Horizon-individuals analysis
        #Horizon X-individuals analysis
Mean_Horizon_X = mean(TrimmedData$Porosity[TrimmedData$Horizon=="X"],na.rm = T) #Mean
Median_Horizon_X = median(TrimmedData$Porosity[TrimmedData$Horizon=="X"],na.rm = T) #Median
StandardDeviation_Horizon_X = sd(TrimmedData$Porosity[TrimmedData$Horizon=="X"], na.rm = T) #Standard deviation
Kurtosis_Horizon_X = kurtosis(TrimmedData$Porosity[TrimmedData$Horizon=="X"], na.rm = T) #Kurtosis
Skewness_Horizon_X = skewness(TrimmedData$Porosity[TrimmedData$Horizon=="X"], na.rm = T) #Skewness
Minimum_Horizon_X = min(TrimmedData$Porosity[TrimmedData$Horizon=="X"],na.rm = T) #Minimum
First_Quartile_Horizon_X = quantile(TrimmedData$Porosity[TrimmedData$Horizon=="X"],probs = 0.25, na.rm = T) #First Quartile
Third_Quartile_Horizon_X = quantile(TrimmedData$Porosity[TrimmedData$Horizon=="X"],probs = 0.75, na.rm = T) #Third Quartile
Maximum_Horizon_X = max(TrimmedData$Porosity[TrimmedData$Horizon=="X"], na.rm = T) #Maximum
Count_Horizon_X = length(TrimmedData$Porosity[TrimmedData$Horizon=="X"]) #Number of data values

Data_Attributes_Horizon_X = c("Mean", "Median", "Standard Deviation", "Kurtosis", "Skewness", "Minimum", "First Quartile", "Third Quartile", "Maximum", "Number of data values") #Summary Table row names
Value_Attributes_Horizon_X = c(Mean_Horizon_X, Median_Horizon_X, StandardDeviation_Horizon_X, Kurtosis_Horizon_X, Skewness_Horizon_X, Minimum_Horizon_X, First_Quartile_Horizon_X, Third_Quartile_Horizon_X, Maximum_Horizon_X, Count_Horizon_X) #Summary Table row values 
SummaryTable_Horizon_X = data.frame(Data_Attributes_Horizon_X,Value_Attributes_Horizon_X) #Constructing the Summary Table
names(SummaryTable_Horizon_X) = c("Data Attribute","Attribute Value for Horizon X Porosity")
write.table(SummaryTable_Horizon_X,file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Tables//Horizon Analysis//SummaryTable_Horizon_X.csv",sep = ",",row.names = F)


jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//Horizon Analysis//Horizon_X boxplot.jpg")
boxplot(TrimmedData$Porosity[TrimmedData$Horizon=="X"], main="Boxplot of all Horizon X Porosity values",ylab = "Porosity (fraction)", las=1)
dev.off()


jpeg("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//Horizon Analysis//Horizon_X RF-N-K Distributions.jpg")
porohist_Horizon_X=hist(TrimmedData$Porosity[TrimmedData$Horizon=="X"],plot=F, breaks = seq(from = 0, to = 0.5, by = 0.05))
porohist_Horizon_X$counts = porohist_Horizon_X$counts/sum(porohist_Horizon_X$counts)
porodens_Horizon_X = density(TrimmedData$Porosity[TrimmedData$Horizon=="X"],na.rm = T,from = 0, to=0.5)
porodens_Horizon_X$y=porodens_Horizon_X$y*0.05
x = TrimmedData$Porosity[TrimmedData$Horizon=="X"]
poronorm_Horizon_X = dnorm(x,mean = Mean_Horizon_X,sd=StandardDeviation_Horizon_X)
poronorm_Horizon_X = poronorm_Horizon_X*0.05
plot(porohist_Horizon_X, freq = T, col = "grey", ylim =  range(0,porohist_Horizon_X$counts,porodens_Horizon_X$y,poronorm_Horizon_X), main="Relative Frequency, Normal and Kernel Distributions of all Horizon X Porosity values",xlab="Porosity (fraction)",ylab = "Probability", las=1)
#text(0.4, 0.2, paste("Mean =", round(Mean_Whole, 4), "\n Std.Dev =", round(StandardDeviation_whole, 4)))
lines(porodens_Horizon_X, col="blue", type = "l", lty="solid", lwd=2)
curve(dnorm(x,mean = Mean_Horizon_X,sd=StandardDeviation_Horizon_X)*0.05,add = T, col="red", type = "l", lty="dashed", lwd=2)
legend("topright",legend=c("Kernel","Normal"),
       text.col=c("blue","red"),lty=c("solid","dashed"),col=c("blue","red"))
#curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole),add = T)
dev.off()




#jpeg("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data//Data - Improved Estimation//For export to R//Plots//Horizon Analysis//Horizon_X histogram.jpg")
#hist(TrimmedData$Porosity[TrimmedData$Horizon=="X"], main="Histogram of all Horizon X Porosity values",xlab="Porosity",las=1)
#dev.off()

        #Horizon Y-individuals analysis
Mean_Horizon_Y = mean(TrimmedData$Porosity[TrimmedData$Horizon=="Y"],na.rm = T) #Mean
Median_Horizon_Y = median(TrimmedData$Porosity[TrimmedData$Horizon=="Y"],na.rm = T) #Median
StandardDeviation_Horizon_Y = sd(TrimmedData$Porosity[TrimmedData$Horizon=="Y"], na.rm = T) #Standard deviation
Kurtosis_Horizon_Y = kurtosis(TrimmedData$Porosity[TrimmedData$Horizon=="Y"], na.rm = T) #Kurtosis
Skewness_Horizon_Y = skewness(TrimmedData$Porosity[TrimmedData$Horizon=="Y"], na.rm = T) #Skewness
Minimum_Horizon_Y = min(TrimmedData$Porosity[TrimmedData$Horizon=="Y"],na.rm = T) #Minimum
First_Quartile_Horizon_Y = quantile(TrimmedData$Porosity[TrimmedData$Horizon=="Y"],probs = 0.25, na.rm = T) #First Quartile
Third_Quartile_Horizon_Y = quantile(TrimmedData$Porosity[TrimmedData$Horizon=="Y"],probs = 0.75, na.rm = T) #Third Quartile
Maximum_Horizon_Y = max(TrimmedData$Porosity[TrimmedData$Horizon=="Y"], na.rm = T) #Maximum
Count_Horizon_Y = length(TrimmedData$Porosity[TrimmedData$Horizon=="Y"]) #Number of data values

Data_Attributes_Horizon_Y = c("Mean", "Median", "Standard Deviation", "Kurtosis", "Skewness", "Minimum", "First Quartile", "Third Quartile", "Maximum", "Number of data values") #Summary Table row names
Value_Attributes_Horizon_Y = c(Mean_Horizon_Y, Median_Horizon_Y, StandardDeviation_Horizon_Y, Kurtosis_Horizon_Y, Skewness_Horizon_Y, Minimum_Horizon_Y, First_Quartile_Horizon_Y, Third_Quartile_Horizon_Y, Maximum_Horizon_Y, Count_Horizon_Y) #Summary Table row values 
SummaryTable_Horizon_Y = data.frame(Data_Attributes_Horizon_Y,Value_Attributes_Horizon_Y) #Constructing the Summary Table
names(SummaryTable_Horizon_Y) = c("Data Attribute","Attribute Value for Horizon Y Porosity")
write.table(SummaryTable_Horizon_Y,file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Tables//Horizon Analysis//SummaryTable_Horizon_Y.csv",sep = ",",row.names = F)

jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//Horizon Analysis//Horizon_Y boxplot.jpg")
boxplot(TrimmedData$Porosity[TrimmedData$Horizon=="Y"], main="Boxplot of all Horizon Y Porosity values",ylab = "Porosity (fraction)", las=1)
dev.off()


jpeg("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//Horizon Analysis//Horizon_Y RF-N-K Distributions.jpg")
porohist_Horizon_Y=hist(TrimmedData$Porosity[TrimmedData$Horizon=="Y"],plot=F, breaks = seq(from = 0, to = 0.5, by = 0.05))
porohist_Horizon_Y$counts = porohist_Horizon_Y$counts/sum(porohist_Horizon_Y$counts)
porodens_Horizon_Y = density(TrimmedData$Porosity[TrimmedData$Horizon=="Y"],na.rm = T,from = 0, to=0.5)
porodens_Horizon_Y$y=porodens_Horizon_Y$y*0.05
x = TrimmedData$Porosity[TrimmedData$Horizon=="Y"]
poronorm_Horizon_Y = dnorm(x,mean = Mean_Horizon_Y,sd=StandardDeviation_Horizon_Y)
poronorm_Horizon_Y = poronorm_Horizon_Y*0.05
plot(porohist_Horizon_Y, freq = T, col = "grey", ylim =  range(0,porohist_Horizon_Y$counts,porodens_Horizon_Y$y,poronorm_Horizon_Y), main="Relative Frequency, Normal and Kernel Distributions of all Horizon Y Porosity values",xlab="Porosity (fraction)",ylab = "Probability", las=1)
#text(0.4, 0.2, paste("Mean =", round(Mean_Whole, 4), "\n Std.Dev =", round(StandardDeviation_whole, 4)))
lines(porodens_Horizon_Y, col="blue", type = "l", lty="solid", lwd=2)
curve(dnorm(x,mean = Mean_Horizon_Y,sd=StandardDeviation_Horizon_Y)*0.05,add = T, col="red", type = "l", lty="dashed", lwd=2)
legend("topright",legend=c("Kernel","Normal"),
       text.col=c("blue","red"),lty=c("solid","dashed"),col=c("blue","red"))
#curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole),add = T)
dev.off()




#jpeg("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data//Data - Improved Estimation//For export to R//Plots//Horizon Analysis//Horizon_Y histogram.jpg")
#hist(TrimmedData$Porosity[TrimmedData$Horizon=="Y"], main="Histogram of all Horizon Y Porosity values",xlab="Porosity",las=1)
#dev.off()

    #Horizon-well-averages analysis
        #Horizon X-well-averages analysis
write.csv(TrimmedData$BoreHole[TrimmedData$Horizon=="X"],"allHorizonXsamples.csv")
allHorizonXsamples = read.csv("allHorizonXsamples.csv")
ListofBorehole_X= levels(as.factor(allHorizonXsamples$x))
BoreHole_Xavgs = rep(0,length(ListofBorehole_X))
BoreHole_Xstddev = rep(0,length(ListofBorehole_X))
for(i in 1:length(ListofBorehole_X)) {BoreHole_Xavgs[i]=mean(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_X[i]&TrimmedData$Horizon=="X"])
BoreHole_Xstddev[i]=sd(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_X[i]&TrimmedData$Horizon=="X"])}
Well_XAveragesTable = data.frame(ListofBorehole_X,BoreHole_Xavgs)
Well_XStandarddeviationTable = data.frame(ListofBorehole_X,BoreHole_Xstddev)

Mean_BoreHole_Xavgs = mean(Well_XAveragesTable$BoreHole_Xavgs,na.rm = T) #Mean
Median_BoreHole_Xavgs = median(Well_XAveragesTable$BoreHole_Xavgs,na.rm = T) #Median
StandardDeviation_BoreHole_Xavgs = sd(Well_XAveragesTable$BoreHole_Xavgs, na.rm = T) #Standard deviation
Kurtosis_BoreHole_Xavgs = kurtosis(Well_XAveragesTable$BoreHole_Xavgs, na.rm = T) #Kurtosis
Skewness_BoreHole_Xavgs = skewness(Well_XAveragesTable$BoreHole_Xavgs, na.rm = T) #Skewness
Minimum_BoreHole_Xavgs = min(Well_XAveragesTable$BoreHole_Xavgs,na.rm = T) #Minimum
First_Quartile_BoreHole_Xavgs = quantile(Well_XAveragesTable$BoreHole_Xavgs,probs = 0.25, na.rm = T) #First Quartile
Third_Quartile_BoreHole_Xavgs = quantile(Well_XAveragesTable$BoreHole_Xavgs,probs = 0.75, na.rm = T) #Third Quartile
Maximum_BoreHole_Xavgs = max(Well_XAveragesTable$BoreHole_Xavgs, na.rm = T) #Maximum
Count_BoreHole_Xavgs = length(Well_XAveragesTable$BoreHole_Xavgs) #Number of data values

Data_Attributes_BoreHole_Xavgs = c("Mean", "Median", "Standard Deviation", "Kurtosis", "Skewness", "Minimum", "First Quartile", "Third Quartile", "Maximum", "Number of data values") #Summary Table row names
Value_Attributes_BoreHole_Xavgs = c(Mean_BoreHole_Xavgs, Median_BoreHole_Xavgs, StandardDeviation_BoreHole_Xavgs, Kurtosis_BoreHole_Xavgs, Skewness_BoreHole_Xavgs, Minimum_BoreHole_Xavgs, First_Quartile_BoreHole_Xavgs, Third_Quartile_BoreHole_Xavgs, Maximum_BoreHole_Xavgs, Count_BoreHole_Xavgs) #Summary Table row values
SummaryTable_BoreHole_Xavgs = data.frame(Data_Attributes_BoreHole_Xavgs,Value_Attributes_BoreHole_Xavgs) #Constructing the Summary Table
names(SummaryTable_BoreHole_Xavgs) = c("Data Attribute","Attribute Value for Horizon X Borehole Porosity Averages")
write.table(SummaryTable_BoreHole_Xavgs,file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Tables//Horizon Analysis//SummaryTable_BoreHole_Xavgs.csv",sep = ",",row.names = F)


jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//Horizon Analysis//HorizonX-BoreHoleaverages boxplot.jpg")
boxplot(Well_XAveragesTable$BoreHole_Xavgs, main="Boxplot of Horizon X Borehole Porosity Averages",ylab = "Mean Porosity (fraction)",las=1)
dev.off()


jpeg("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//Horizon Analysis//HorizonX-BoreHoleaverages RF-N-K Distributions.jpg")
porohist_HorizonX_wellaverages=hist(Well_XAveragesTable$BoreHole_Xavgs,plot=F, breaks = seq(from = 0, to = 0.5, by = 0.05))
porohist_HorizonX_wellaverages$counts = porohist_HorizonX_wellaverages$counts/sum(porohist_HorizonX_wellaverages$counts)
porodens_HorizonX_wellaverages = density(Well_XAveragesTable$BoreHole_Xavgs,na.rm = T,from = 0, to=0.5)
porodens_HorizonX_wellaverages$y=porodens_HorizonX_wellaverages$y*0.05
x = Well_XAveragesTable$BoreHole_Xavgs
poronorm_HorizonX_wellaverages = dnorm(x,mean = Mean_BoreHole_Xavgs,sd=StandardDeviation_BoreHole_Xavgs)
poronorm_HorizonX_wellaverages = poronorm_HorizonX_wellaverages*0.05
plot(porohist_HorizonX_wellaverages, freq = T, col = "grey", ylim =  range(0,porohist_HorizonX_wellaverages$counts,porodens_HorizonX_wellaverages$y,poronorm_HorizonX_wellaverages), main="Relative Frequency, Normal and Kernel Distributions of Horizon X Borehole Porosity Averages",xlab="Mean Porosity (fraction)",ylab = "Probability", las=1)
#text(0.4, 0.2, paste("Mean =", round(Mean_Whole, 4), "\n Std.Dev =", round(StandardDeviation_whole, 4)))
lines(porodens_HorizonX_wellaverages, col="blue", type = "l", lty="solid", lwd=2)
curve(dnorm(x,mean = Mean_BoreHole_Xavgs,sd=StandardDeviation_BoreHole_Xavgs)*0.05,add = T, col="red", type = "l", lty="dashed", lwd=2)
legend("topright",legend=c("Kernel","Normal"),
       text.col=c("blue","red"),lty=c("solid","dashed"),col=c("blue","red"))
#curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole),add = T)
dev.off()



#jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data//Data - Improved Estimation//For export to R//Plots//Horizon Analysis//HorizonX-BoreHoleaverages histogram.jpg")
#hist(Well_XAveragesTable$BoreHole_Xavgs, main="Histogram of Horizon X Borehole Porosity Averages",xlab="Mean Porosity",las=1)
#dev.off()

        #Horizon Y-well-averages analysis
write.csv(TrimmedData$BoreHole[TrimmedData$Horizon=="Y"],"allHorizonYsamples.csv")
allHorizonYsamples = read.csv("allHorizonYsamples.csv")
ListofBorehole_Y= levels(as.factor(allHorizonYsamples$x))
BoreHole_Yavgs = rep(0,length(ListofBorehole_Y))
BoreHole_Ystddev = rep(0,length(ListofBorehole_Y))
for(i in 1:length(ListofBorehole_Y)) {BoreHole_Yavgs[i]=mean(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_Y[i]&TrimmedData$Horizon=="Y"])
BoreHole_Ystddev[i]=sd(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_Y[i]&TrimmedData$Horizon=="Y"])}
Well_YAveragesTable = data.frame(ListofBorehole_Y,BoreHole_Yavgs)
Well_YStandardDeviationTable = data.frame(ListofBorehole_Y,BoreHole_Ystddev)

Mean_BoreHole_Yavgs = mean(Well_YAveragesTable$BoreHole_Yavgs,na.rm = T) #Mean
Median_BoreHole_Yavgs = median(Well_YAveragesTable$BoreHole_Yavgs,na.rm = T) #Median
StandardDeviation_BoreHole_Yavgs = sd(Well_YAveragesTable$BoreHol                                                                                                                                                         e_Yavgs, na.rm = T) #Standard deviation
Kurtosis_BoreHole_Yavgs = kurtosis(Well_YAveragesTable$BoreHole_Yavgs, na.rm = T) #Kurtosis
Skewness_BoreHole_Yavgs = skewness(Well_YAveragesTable$BoreHole_Yavgs, na.rm = T) #Skewness
Minimum_BoreHole_Yavgs = min(Well_YAveragesTable$BoreHole_Yavgs,na.rm = T) #Minimum
First_Quartile_BoreHole_Yavgs = quantile(Well_YAveragesTable$BoreHole_Yavgs,probs = 0.25, na.rm = T) #First Quartile
Third_Quartile_BoreHole_Yavgs = quantile(Well_YAveragesTable$BoreHole_Yavgs,probs = 0.75, na.rm = T) #Third Quartile
Maximum_BoreHole_Yavgs = max(Well_YAveragesTable$BoreHole_Yavgs, na.rm = T) #Maximum
Count_BoreHole_Yavgs = length(Well_YAveragesTable$BoreHole_Yavgs) #Number of data values

Data_Attributes_BoreHole_Yavgs = c("Mean", "Median", "Standard Deviation", "Kurtosis", "Skewness", "Minimum", "First Quartile", "Third Quartile", "Maximum", "Number of data values") #Summary Table row names
Value_Attributes_BoreHole_Yavgs = c(Mean_BoreHole_Yavgs, Median_BoreHole_Yavgs, StandardDeviation_BoreHole_Yavgs, Kurtosis_BoreHole_Yavgs, Skewness_BoreHole_Yavgs, Minimum_BoreHole_Yavgs, First_Quartile_BoreHole_Yavgs, Third_Quartile_BoreHole_Yavgs, Maximum_BoreHole_Yavgs, Count_BoreHole_Yavgs) #Summary Table row values
SummaryTable_BoreHole_Yavgs = data.frame(Data_Attributes_BoreHole_Yavgs,Value_Attributes_BoreHole_Yavgs) #Constructing the Summary Table
names(SummaryTable_BoreHole_Yavgs) = c("Data Attribute","Attribute Value for Horizon Y Borehole Porosity Averages")
write.table(SummaryTable_BoreHole_Yavgs,file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Tables//Horizon Analysis//SummaryTable_BoreHole_Yavgs.csv",sep = ",",row.names = F)


jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//Horizon Analysis//HorizonY-BoreHoleaverages boxplot.jpg")
boxplot(Well_YAveragesTable$BoreHole_Yavgs, main="Boxplot of Horizon Y Borehole Porosity Averages",ylab = "Mean Porosity (fraction)",las=1)
dev.off()


jpeg("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//Horizon Analysis//HorizonY-BoreHoleaverages RF-N-K Distributions.jpg")
porohist_HorizonY_wellaverages=hist(Well_YAveragesTable$BoreHole_Yavgs,plot=F, breaks = seq(from = 0, to = 0.5, by = 0.05))
porohist_HorizonY_wellaverages$counts = porohist_HorizonY_wellaverages$counts/sum(porohist_HorizonY_wellaverages$counts)
porodens_HorizonY_wellaverages = density(Well_YAveragesTable$BoreHole_Yavgs,na.rm = T,from = 0, to=0.5)
porodens_HorizonY_wellaverages$y=porodens_HorizonY_wellaverages$y*0.05
x = Well_YAveragesTable$BoreHole_Yavgs
poronorm_HorizonY_wellaverages = dnorm(x,mean = Mean_BoreHole_Yavgs,sd=StandardDeviation_BoreHole_Yavgs)
poronorm_HorizonY_wellaverages = poronorm_HorizonY_wellaverages*0.05
plot(porohist_HorizonY_wellaverages, freq = T, col = "grey", ylim =  range(0,porohist_HorizonY_wellaverages$counts,porodens_HorizonY_wellaverages$y,poronorm_HorizonY_wellaverages), main="Relative Frequency, Normal and Kernel Distributions of Horizon Y Borehole Porosity Averages",xlab="Mean Porosity (fraction)",ylab = "Probability", las=1)
#text(0.4, 0.2, paste("Mean =", round(Mean_Whole, 4), "\n Std.Dev =", round(StandardDeviation_whole, 4)))
lines(porodens_HorizonY_wellaverages, col="blue", type = "l", lty="solid", lwd=2)
curve(dnorm(x,mean = Mean_BoreHole_Yavgs,sd=StandardDeviation_BoreHole_Yavgs)*0.05,add = T, col="red", type = "l", lty="dashed", lwd=2)
legend("topright",legend=c("Kernel","Normal"),
       text.col=c("blue","red"),lty=c("solid","dashed"),col=c("blue","red"))
#curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole),add = T)
dev.off()


#jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data//Data - Improved Estimation//For export to R//Plots//Horizon Analysis//HorizonY-BoreHoleaverages histogram.jpg")
#hist(Well_YAveragesTable$BoreHole_Yavgs, main="Histogram of Horizon Y Borehole Porosity Averages",xlab="Mean Porosity",las=1)
#dev.off()

#Porosity Analysis - Per Borehole
    #Per Borehole - whole

ListofBorehole_exc_BH11 = c(ListofBorehole[1],ListofBorehole[3:8],ListofBorehole[10:33]) #Actually, excluding BH11 and BH22; they are offending the kernel and normal distributions.
for(i in 1:length(ListofBorehole_exc_BH11)) {
BoreholeID = paste("Borehole",ListofBorehole_exc_BH11[i],sep=" ")
Mean_BoreHole_i =mean(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]],na.rm = T) #Mean 
Median_BoreHole_i = median(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]],na.rm = T) #Median
StandardDeviation_BoreHole_i = sd(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]], na.rm = T) #Standard deviation
Kurtosis_BoreHole_i = kurtosis(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]], na.rm = T) #Kurtosis
Skewness_BoreHole_i = skewness(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]], na.rm = T) #Skewness
Minimum_BoreHole_i = min(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]],na.rm = T) #Minimum
First_Quartile_BoreHole_i = quantile(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]],probs = 0.25, na.rm = T) #First Quartile
Third_Quartile_BoreHole_i = quantile(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]],probs = 0.75, na.rm = T) #Third Quartile
Maximum_BoreHole_i = max(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]], na.rm = T) #Maximum
Count_BoreHole_i = length(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]]) #Number of data values

Data_Attributes = c("Mean", "Median", "Standard Deviation", "Kurtosis", "Skewness", "Minimum", "First Quartile", "Third Quartile", "Maximum", "Number of data values") #Summary Table row names
SummaryTable_Borehole_i = data.frame(Data_Attributes,c(Mean_BoreHole_i, Median_BoreHole_i, StandardDeviation_BoreHole_i, Kurtosis_BoreHole_i, Skewness_BoreHole_i, Minimum_BoreHole_i, First_Quartile_BoreHole_i, Third_Quartile_BoreHole_i, Maximum_BoreHole_i, Count_BoreHole_i)) #Constructing the Summary Table
names(SummaryTable_Borehole_i) = c("Data Attribute",paste("Attribute Value for",BoreholeID,sep = " "))
write.table(SummaryTable_Borehole_i,file.path("C:","Users","TTOWG","645","1 karia def","2. CU","Projects","Bitumen Recovery","Data Analysis","Explorative Data Analysis","For export to R","Tables","BoreHole Analysis",paste("SummaryTable_",BoreholeID,".csv",sep = "")),sep = ",",row.names = F)


jpeg(file.path("C:","Users","TTOWG","645","1 karia def","2. CU","Projects","Bitumen Recovery","Data Analysis","Explorative Data Analysis","For export to R","Plots","BoreHole Analysis",paste(BoreholeID,"boxplot.jpg",sep = " ")))
boxplot(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]], main=paste("Boxplot of", BoreholeID,"Porosity Values",sep = " "),ylab = "Porosity (fraction)",las=1)
dev.off()


jpeg(file.path("C:","Users","TTOWG","645","1 karia def","2. CU","Projects","Bitumen Recovery","Data Analysis","Explorative Data Analysis","For export to R","Plots","BoreHole Analysis",paste(BoreholeID,"RF-N-K Distributions.jpg",sep = " ")))
porohist_BoreHole_i=hist(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]],plot=F, breaks = seq(from = 0, to = 0.5, by = 0.05))
porohist_BoreHole_i$counts = porohist_BoreHole_i$counts/sum(porohist_BoreHole_i$counts)
porodens_BoreHole_i = density(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]],na.rm = T,from = 0, to=0.5)
porodens_BoreHole_i$y=porodens_BoreHole_i$y*0.05
x = TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_BH11[i]]
poronorm_BoreHole_i = dnorm(x,mean = Mean_BoreHole_i,sd=StandardDeviation_BoreHole_i)
poronorm_BoreHole_i = poronorm_BoreHole_i*0.05
plot(porohist_BoreHole_i, freq = T, col = "grey", ylim =  range(0,porohist_BoreHole_i$counts,porodens_BoreHole_i$y,poronorm_BoreHole_i), main=paste("Relative Frequency, Normal and Kernel Distributions of", BoreholeID,"Porosity Values",sep = " "),xlab="Porosity (fraction)",ylab = "Probability", las=1)
#text(0.4, 0.2, paste("Mean =", round(Mean_Whole, 4), "\n Std.Dev =", round(StandardDeviation_whole, 4)))
lines(porodens_BoreHole_i, col="blue", type = "l", lty="solid", lwd=2)
curve(dnorm(x,mean = Mean_BoreHole_i,sd=StandardDeviation_BoreHole_i)*0.05,add = T, col="red", type = "l", lty="dashed", lwd=2)
legend("topright",legend=c("Kernel","Normal"),
       text.col=c("blue","red"),lty=c("solid","dashed"),col=c("blue","red"))
#curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole),add = T)
dev.off()


#jpeg(file.path("C:","Users","TTOWG","645","1 karia def","2. CU","Projects","Bitumen Recovery","Data","Data - Improved Estimation","For export to R","Plots","BoreHole Analysis",paste(BoreholeID,"histogram.jpg",sep = " ")))
#hist(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole[i]], main=paste("Histogram of", BoreholeID,"Porosity Values",sep = " "),xlab="Porosity",las=1)
#dev.off()
}

    #Per Borehole - Horizon
         #Horizon X
ListofBorehole_exc_offend_X = c(ListofBorehole_X[1],ListofBorehole_X[3:6],ListofBorehole_X[8:30]) #Excluding BH11 and BH22; they are offending the kernel and normal distributions.
for(j in 1:length(ListofBorehole_exc_offend_X)) {
BoreholeID_X = paste("Borehole",ListofBorehole_exc_offend_X[j],sep=" ")
Mean_BoreHole_i_X =mean(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"],na.rm = T) #Mean 
Median_BoreHole_i_X = median(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"],na.rm = T) #Median
StandardDeviation_BoreHole_i_X = sd(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"], na.rm = T) #Standard deviation
Kurtosis_BoreHole_i_X = kurtosis(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"], na.rm = T) #Kurtosis
Skewness_BoreHole_i_X = skewness(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"], na.rm = T) #Skewness
Minimum_BoreHole_i_X = min(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"],na.rm = T) #Minimum
First_Quartile_BoreHole_i_X = quantile(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"],probs = 0.25, na.rm = T) #First Quartile
Third_Quartile_BoreHole_i_X = quantile(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"],probs = 0.75, na.rm = T) #Third Quartile
Maximum_BoreHole_i_X = max(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"], na.rm = T) #Maximum
Count_BoreHole_i_X = length(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"]) #Number of data values

Data_Attributes = c("Mean", "Median", "Standard Deviation", "Kurtosis", "Skewness", "Minimum", "First Quartile", "Third Quartile", "Maximum", "Number of data values") #Summary Table row names
SummaryTable_Borehole_i_X = data.frame(Data_Attributes,c(Mean_BoreHole_i_X, Median_BoreHole_i_X, StandardDeviation_BoreHole_i_X, Kurtosis_BoreHole_i_X, Skewness_BoreHole_i_X, Minimum_BoreHole_i_X, First_Quartile_BoreHole_i_X, Third_Quartile_BoreHole_i_X, Maximum_BoreHole_i_X, Count_BoreHole_i_X)) #Constructing the Summary Table
names(SummaryTable_Borehole_i_X) = c("Data Attribute",paste("Attribute Value for",BoreholeID_X,"Horizon X",sep = " "))
write.table(SummaryTable_Borehole_i_X,file.path("C:","Users","TTOWG","645","1 karia def","2. CU","Projects","Bitumen Recovery","Data Analysis","Explorative Data Analysis","For export to R","Tables","BoreHole Analysis",paste("SummaryTable_",BoreholeID_X,"_Horizon X.csv",sep = "")),sep = ",",row.names = F)

jpeg(file.path("C:","Users","TTOWG","645","1 karia def","2. CU","Projects","Bitumen Recovery","Data Analysis","Explorative Data Analysis","For export to R","Plots","BoreHole Analysis",paste(BoreholeID_X,"Horizon X boxplot.jpg",sep = " ")))
boxplot(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"], main=paste("Boxplot of", BoreholeID_X,"Horizon X Porosity Values",sep = " "),ylab = "Porosity (fraction)",las=1)
dev.off()


jpeg(file.path("C:","Users","TTOWG","645","1 karia def","2. CU","Projects","Bitumen Recovery","Data Analysis","Explorative Data Analysis","For export to R","Plots","BoreHole Analysis",paste(BoreholeID_X,"Horizon X RF-N-K Distributions.jpg",sep = " ")))
porohist_BoreHole_i_X=hist(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"],plot=F, breaks = seq(from = 0, to = 0.5, by = 0.05))
porohist_BoreHole_i_X$counts = porohist_BoreHole_i_X$counts/sum(porohist_BoreHole_i_X$counts)
porodens_BoreHole_i_X = density(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"],na.rm = T,from = 0, to=0.5)
porodens_BoreHole_i_X$y=porodens_BoreHole_i_X$y*0.05
x = TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_X[j]&TrimmedData$Horizon=="X"]
poronorm_BoreHole_i_X = dnorm(x,mean = Mean_BoreHole_i_X,sd=StandardDeviation_BoreHole_i_X)
poronorm_BoreHole_i_X = poronorm_BoreHole_i_X*0.05
plot(porohist_BoreHole_i_X, freq = T, col = "grey", ylim =  range(0,porohist_BoreHole_i_X$counts,porodens_BoreHole_i_X$y,poronorm_BoreHole_i_X), main=paste("Relative Frequency, Normal and Kernel Distributions of", BoreholeID_X,"Horizon X Porosity Values",sep = " "),xlab="Porosity (fraction)",ylab = "Probability", las=1)
#text(0.4, 0.2, paste("Mean =", round(Mean_Whole, 4), "\n Std.Dev =", round(StandardDeviation_whole, 4)))
lines(porodens_BoreHole_i_X, col="blue", type = "l", lty="solid", lwd=2)
curve(dnorm(x,mean = Mean_BoreHole_i_X,sd=StandardDeviation_BoreHole_i_X)*0.05,add = T, col="red", type = "l", lty="dashed", lwd=2)
legend("topright",legend=c("Kernel","Normal"),
       text.col=c("blue","red"),lty=c("solid","dashed"),col=c("blue","red"))
#curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole),add = T)
dev.off()




#jpeg(file.path("C:","Users","TTOWG","645","1 karia def","2. CU","Projects","Bitumen Recovery","Data","Data - Improved Estimation","For export to R","Plots","BoreHole Analysis",paste(BoreholeID_X,"Horizon X histogram.jpg",sep = " ")))
#hist(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_X[j]&TrimmedData$Horizon=="X"], main=paste("Histogram of", BoreholeID_X,"Horizon X Porosity Values",sep = " "),xlab="Porosity",las=1)
#dev.off()
}

    #Horizon Y
ListofBorehole_exc_offend_Y = c(ListofBorehole_Y[1:6],ListofBorehole_Y[8:13],ListofBorehole_Y[15:25]) #Excluding BH11 and BH22; they are offending the kernel and normal distributions.
for(k in 1:length(ListofBorehole_exc_offend_Y)) {
  BoreholeID_Y = paste("Borehole",ListofBorehole_exc_offend_Y[k],sep=" ")
  Mean_BoreHole_i_Y =mean(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"],na.rm = T) #Mean 
  Median_BoreHole_i_Y = median(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"],na.rm = T) #Median
  StandardDeviation_BoreHole_i_Y = sd(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"], na.rm = T) #Standard deviation
  Kurtosis_BoreHole_i_Y = kurtosis(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"], na.rm = T) #Kurtosis
  Skewness_BoreHole_i_Y = skewness(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"], na.rm = T) #Skewness
  Minimum_BoreHole_i_Y = min(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"],na.rm = T) #Minimum
  First_Quartile_BoreHole_i_Y = quantile(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"],probs = 0.25, na.rm = T) #First Quartile
  Third_Quartile_BoreHole_i_Y = quantile(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"],probs = 0.75, na.rm = T) #Third Quartile
  Maximum_BoreHole_i_Y = max(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"], na.rm = T) #Maximum
  Count_BoreHole_i_Y = length(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"]) #Number of data values
  
  Data_Attributes = c("Mean", "Median", "Standard Deviation", "Kurtosis", "Skewness", "Minimum", "First Quartile", "Third Quartile", "Maximum", "Number of data values") #Summary Table row names
  SummaryTable_Borehole_i_Y = data.frame(Data_Attributes,c(Mean_BoreHole_i_Y, Median_BoreHole_i_Y, StandardDeviation_BoreHole_i_Y, Kurtosis_BoreHole_i_Y, Skewness_BoreHole_i_Y, Minimum_BoreHole_i_Y, First_Quartile_BoreHole_i_Y, Third_Quartile_BoreHole_i_Y, Maximum_BoreHole_i_Y, Count_BoreHole_i_Y)) #Constructing the Summary Table
  names(SummaryTable_Borehole_i_Y) = c("Data Attribute",paste("Attribute Value for",BoreholeID_Y,"Horizon Y",sep = " "))
  write.table(SummaryTable_Borehole_i_Y,file.path("C:","Users","TTOWG","645","1 karia def","2. CU","Projects","Bitumen Recovery","Data Analysis","Explorative Data Analysis","For export to R","Tables","BoreHole Analysis",paste("SummaryTable_",BoreholeID_Y,"_Horizon Y.csv",sep = "")),sep = ",",row.names = F)
  
  jpeg(file.path("C:","Users","TTOWG","645","1 karia def","2. CU","Projects","Bitumen Recovery","Data Analysis","Explorative Data Analysis","For export to R","Plots","BoreHole Analysis",paste(BoreholeID_Y,"Horizon Y boxplot.jpg",sep = " ")))
  boxplot(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"], main=paste("Boxplot of", BoreholeID_Y,"Horizon Y Porosity Values",sep = " "),ylab = "Porosity (fraction)",las=1)
  dev.off()
  
  
  jpeg(file.path("C:","Users","TTOWG","645","1 karia def","2. CU","Projects","Bitumen Recovery","Data Analysis","Explorative Data Analysis","For export to R","Plots","BoreHole Analysis",paste(BoreholeID_Y,"Horizon Y RF-N-K Distributions.jpg",sep = " ")))
  porohist_BoreHole_i_Y=hist(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"],plot=F, breaks = seq(from = 0, to = 0.5, by = 0.05))
  porohist_BoreHole_i_Y$counts = porohist_BoreHole_i_Y$counts/sum(porohist_BoreHole_i_Y$counts)
  porodens_BoreHole_i_Y = density(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"],na.rm = T,from = 0, to=0.5)
  porodens_BoreHole_i_Y$y=porodens_BoreHole_i_Y$y*0.05
  x = TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_exc_offend_Y[k]&TrimmedData$Horizon=="Y"]
  poronorm_BoreHole_i_Y = dnorm(x,mean = Mean_BoreHole_i_Y,sd=StandardDeviation_BoreHole_i_Y)
  poronorm_BoreHole_i_Y = poronorm_BoreHole_i_Y*0.05
  plot(porohist_BoreHole_i_Y, freq = T, col = "grey", ylim =  range(0,porohist_BoreHole_i_Y$counts,porodens_BoreHole_i_Y$y,poronorm_BoreHole_i_Y), main=paste("Relative Frequency, Normal and Kernel Distributions of", BoreholeID_Y,"Horizon Y Porosity Values",sep = " "),xlab="Porosity (fraction)",ylab = "Probability", las=1)
  #text(0.4, 0.2, paste("Mean =", round(Mean_Whole, 4), "\n Std.Dev =", round(StandardDeviation_whole, 4)))
  lines(porodens_BoreHole_i_Y, col="blue", type = "l", lty="solid", lwd=2)
  curve(dnorm(x,mean = Mean_BoreHole_i_Y,sd=StandardDeviation_BoreHole_i_Y)*0.05,add = T, col="red", type = "l", lty="dashed", lwd=2)
  legend("topright",legend=c("Kernel","Normal"),
         text.col=c("blue","red"),lty=c("solid","dashed"),col=c("blue","red"))
  #curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole),add = T)
  dev.off() 
  
  
  
  #jpeg(file.path("C:","Users","TTOWG","645","1 karia def","2. CU","Projects","Bitumen Recovery","Data","Data - Improved Estimation","For export to R","Plots","BoreHole Analysis",paste(BoreholeID_Y,"Horizon Y histogram.jpg",sep = " ")))
  #hist(TrimmedData$Porosity[TrimmedData$BoreHole==ListofBorehole_Y[k]&TrimmedData$Horizon=="Y"], main=paste("Histogram of", BoreholeID_Y,"Horizon Y Porosity Values",sep = " "),xlab="Porosity",las=1)
  #dev.off()
}



#Top-Thickness Data



Top_ThicknessData = read.csv(file.choose(),header=T)
Top_ThicknessData$Borehole = as.factor(Top_ThicknessData$Borehole)

 #Top Analysis
       #Horizon X

Mean_Top_X = mean(Top_ThicknessData$Depth.of.Top..X,na.rm = T) #Mean
Median_Top_X = median(Top_ThicknessData$Depth.of.Top..X,na.rm = T) #Median
StandardDeviation_Top_X = sd(Top_ThicknessData$Depth.of.Top..X,na.rm = T) #Standard deviation
Kurtosis_Top_X = kurtosis(Top_ThicknessData$Depth.of.Top..X,na.rm = T) #Kurtosis
Skewness_Top_X = skewness(Top_ThicknessData$Depth.of.Top..X,na.rm = T) #Skewness
Minimum_Top_X = min(Top_ThicknessData$Depth.of.Top..X,na.rm = T) #Minimum
First_Quartile_Top_X = quantile(Top_ThicknessData$Depth.of.Top..X,na.rm = T,probs = 0.25) #First Quartile
Third_Quartile_Top_X = quantile(Top_ThicknessData$Depth.of.Top..X,na.rm = T,probs = 0.75) #Third Quartile
Maximum_Top_X = max(Top_ThicknessData$Depth.of.Top..X,na.rm = T) #Maximum
Count_Top_X = length(Top_ThicknessData$Depth.of.Top..X) #Number of data values

Data_Attributes_Top_X = c("Mean", "Median", "Standard Deviation", "Kurtosis", "Skewness", "Minimum", "First Quartile", "Third Quartile", "Maximum", "Number of data values") #Summary Table row names
Value_Attributes_Top_X = c(Mean_Top_X, Median_Top_X, StandardDeviation_Top_X, Kurtosis_Top_X, Skewness_Top_X, Minimum_Top_X, First_Quartile_Top_X, Third_Quartile_Top_X, Maximum_Top_X, Count_Top_X) #Summary Table row values 
SummaryTable_Top_X = data.frame(Data_Attributes_Top_X,Value_Attributes_Top_X) #Constructing the Summary Table
names(SummaryTable_Top_X) = c("Data Attribute","Attribute Value for Depth to Top of Horizon X")
write.table(SummaryTable_Top_X,file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Top-Thickness Analysis//Tables//SummaryTable_Top_X.csv",sep = ",",row.names = F)



jpeg("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Top-Thickness Analysis//Plots//Top_X RF-N-K Distributions.jpg")
porohist_Top_X=hist(Top_ThicknessData$Depth.of.Top..X,plot=F, breaks = 5)
porohist_Top_X$counts = porohist_Top_X$counts/sum(porohist_Top_X$counts)
porodens_Top_X = density(Top_ThicknessData$Depth.of.Top..X,na.rm = T,from = 0, to=250)
porodens_Top_X$y=porodens_Top_X$y*50
x = Top_ThicknessData$Depth.of.Top..X
poronorm_Top_X = dnorm(x,mean = Mean_Top_X,sd=StandardDeviation_Top_X)
poronorm_Top_X = poronorm_Top_X*50
plot(porohist_Top_X, freq = T, col = "grey", main="Relative Frequency, Normal and Kernel Distributions of Depths to Top of Horizon X", ylim = c(0,0.5), xlab="Depth to Top of Horizon X (ft)",ylab = "Probability", las=1)
#text(0.4, 0.2, paste("Mean =", round(Mean_Whole, 4), "\n Std.Dev =", round(StandardDeviation_whole, 4)))
lines(porodens_Top_X, col="blue", type = "l", lty="solid", lwd=2)
curve(dnorm(x,mean = Mean_Top_X,sd=StandardDeviation_Top_X)*50,add = T, col="red", type = "l", lty="dashed", lwd=2)
legend("topright",legend=c("Kernel","Normal"),
       text.col=c("blue","red"),lty=c("solid","dashed"),col=c("blue","red"))
#curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole),add = T)
dev.off()


#Horizon Y

Mean_Top_Y = mean(Top_ThicknessData$Depth.of.Top..Y,na.rm = T) #Mean
Median_Top_Y = median(Top_ThicknessData$Depth.of.Top..Y,na.rm = T) #Median
StandardDeviation_Top_Y = sd(Top_ThicknessData$Depth.of.Top..Y,na.rm = T) #Standard deviation
Kurtosis_Top_Y = kurtosis(Top_ThicknessData$Depth.of.Top..Y,na.rm = T) #Kurtosis
Skewness_Top_Y = skewness(Top_ThicknessData$Depth.of.Top..Y,na.rm = T) #Skewness
Minimum_Top_Y = min(Top_ThicknessData$Depth.of.Top..Y,na.rm = T) #Minimum
First_Quartile_Top_Y = quantile(Top_ThicknessData$Depth.of.Top..Y,na.rm = T,probs = 0.25) #First Quartile
Third_Quartile_Top_Y = quantile(Top_ThicknessData$Depth.of.Top..Y,na.rm = T,probs = 0.75) #Third Quartile
Maximum_Top_Y = max(Top_ThicknessData$Depth.of.Top..Y,na.rm = T) #Maximum
Count_Top_Y = length(Top_ThicknessData$Depth.of.Top..Y) #Number of data values

Data_Attributes_Top_Y = c("Mean", "Median", "Standard Deviation", "Kurtosis", "Skewness", "Minimum", "First Quartile", "Third Quartile", "Maximum", "Number of data values") #Summary Table row names
Value_Attributes_Top_Y = c(Mean_Top_Y, Median_Top_Y, StandardDeviation_Top_Y, Kurtosis_Top_Y, Skewness_Top_Y, Minimum_Top_Y, First_Quartile_Top_Y, Third_Quartile_Top_Y, Maximum_Top_Y, Count_Top_Y) #Summary Table row values 
SummaryTable_Top_Y = data.frame(Data_Attributes_Top_Y,Value_Attributes_Top_Y) #Constructing the Summary Table
names(SummaryTable_Top_Y) = c("Data Attribute","Attribute Value for Depth to Top of Horizon Y")
write.table(SummaryTable_Top_Y,file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Top-Thickness Analysis//Tables//SummaryTable_Top_Y.csv",sep = ",",row.names = F)



jpeg("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Top-Thickness Analysis//Plots//Top_Y RF-N-K Distributions.jpg")
porohist_Top_Y=hist(Top_ThicknessData$Depth.of.Top..Y,plot=F, breaks = 5)
porohist_Top_Y$counts = porohist_Top_Y$counts/sum(porohist_Top_Y$counts)
porodens_Top_Y = density(Top_ThicknessData$Depth.of.Top..Y,na.rm = T,from = 50, to=300)
porodens_Top_Y$y=porodens_Top_Y$y*50
x = Top_ThicknessData$Depth.of.Top..Y
poronorm_Top_Y = dnorm(x,mean = Mean_Top_Y,sd=StandardDeviation_Top_Y)
poronorm_Top_Y = poronorm_Top_Y*50
plot(porohist_Top_Y, freq = T, col = "grey", main="Relative Frequency, Normal and Kernel Distributions of Depths to Top of Horizon Y", xlim = c(50,300), ylim = c(0,0.5), xlab="Depth to Top of Horizon Y (ft)", ylab = "Probability", las=1)
#text(0.4, 0.2, paste("Mean =", round(Mean_Whole, 4), "\n Std.Dev =", round(StandardDeviation_whole, 4)))
lines(porodens_Top_Y, col="blue", type = "l", lty="solid", lwd=2)
curve(dnorm(x,mean = Mean_Top_Y,sd=StandardDeviation_Top_Y)*50,add = T, col="red", type = "l", lty="dashed", lwd=2)
legend("topright",legend=c("Kernel","Normal"),
       text.col=c("blue","red"),lty=c("solid","dashed"),col=c("blue","red"))
#curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole),add = T)
dev.off()


  #Thickness Analysis

       #Horizon X

Mean_Thickness_X = mean(Top_ThicknessData$Thickness....X,na.rm = T) #Mean
Median_Thickness_X = median(Top_ThicknessData$Thickness....X,na.rm = T) #Median
StandardDeviation_Thickness_X = sd(Top_ThicknessData$Thickness....X,na.rm = T) #Standard deviation
Kurtosis_Thickness_X = kurtosis(Top_ThicknessData$Thickness....X,na.rm = T) #Kurtosis
Skewness_Thickness_X = skewness(Top_ThicknessData$Thickness....X,na.rm = T) #Skewness
Minimum_Thickness_X = min(Top_ThicknessData$Thickness....X,na.rm = T) #Minimum
First_Quartile_Thickness_X = quantile(Top_ThicknessData$Thickness....X,na.rm = T,probs = 0.25) #First Quartile
Third_Quartile_Thickness_X = quantile(Top_ThicknessData$Thickness....X,na.rm = T,probs = 0.75) #Third Quartile
Maximum_Thickness_X = max(Top_ThicknessData$Thickness....X,na.rm = T) #Maximum
Count_Thickness_X = length(Top_ThicknessData$Thickness....X) #Number of data values

Data_Attributes_Thickness_X = c("Mean", "Median", "Standard Deviation", "Kurtosis", "Skewness", "Minimum", "First Quartile", "Third Quartile", "Maximum", "Number of data values") #Summary Table row names
Value_Attributes_Thickness_X = c(Mean_Thickness_X, Median_Thickness_X, StandardDeviation_Thickness_X, Kurtosis_Thickness_X, Skewness_Thickness_X, Minimum_Thickness_X, First_Quartile_Thickness_X, Third_Quartile_Thickness_X, Maximum_Thickness_X, Count_Thickness_X) #Summary Table row values 
SummaryTable_Thickness_X = data.frame(Data_Attributes_Thickness_X,Value_Attributes_Thickness_X) #Constructing the Summary Table
names(SummaryTable_Thickness_X) = c("Data Attribute","Attribute Value for Thickness of Horizon X")
write.table(SummaryTable_Thickness_X,file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Top-Thickness Analysis//Tables//SummaryTable_Thickness_X.csv",sep = ",",row.names = F)



jpeg("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Top-Thickness Analysis//Plots//Thickness_X RF-N-K Distributions.jpg")
porohist_Thickness_X=hist(Top_ThicknessData$Thickness....X,plot=F)
porohist_Thickness_X$counts = porohist_Thickness_X$counts/sum(porohist_Thickness_X$counts)
porodens_Thickness_X = density(Top_ThicknessData$Thickness....X,na.rm = T,from = 10, to=80)
porodens_Thickness_X$y=porodens_Thickness_X$y*10
x = Top_ThicknessData$Thickness....X
poronorm_Thickness_X = dnorm(x,mean = Mean_Thickness_X,sd=StandardDeviation_Thickness_X)
poronorm_Thickness_X = poronorm_Thickness_X*10
plot(porohist_Thickness_X, freq = T, col = "grey", main="Relative Frequency, Normal and Kernel Distributions of Thickness of Horizon X",ylim=c(0,0.5), xlab="Thickness of Horizon X (ft)",ylab = "Probability", las=1)
#text(0.4, 0.2, paste("Mean =", round(Mean_Whole, 4), "\n Std.Dev =", round(StandardDeviation_whole, 4)))
lines(porodens_Thickness_X, col="blue", type = "l", lty="solid", lwd=2)
curve(dnorm(x,mean = Mean_Thickness_X,sd=StandardDeviation_Thickness_X)*10,add = T, col="red", type = "l", lty="dashed", lwd=2)
legend("topright",legend=c("Kernel","Normal"),
       text.col=c("blue","red"),lty=c("solid","dashed"),col=c("blue","red"))
#curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole),add = T)
dev.off()



#Horizon Y

Mean_Thickness_Y = mean(Top_ThicknessData$Thickness....Y,na.rm = T) #Mean
Median_Thickness_Y = median(Top_ThicknessData$Thickness....Y,na.rm = T) #Median
StandardDeviation_Thickness_Y = sd(Top_ThicknessData$Thickness....Y,na.rm = T) #Standard deviation
Kurtosis_Thickness_Y = kurtosis(Top_ThicknessData$Thickness....Y,na.rm = T) #Kurtosis
Skewness_Thickness_Y = skewness(Top_ThicknessData$Thickness....Y,na.rm = T) #Skewness
Minimum_Thickness_Y = min(Top_ThicknessData$Thickness....Y,na.rm = T) #Minimum
First_Quartile_Thickness_Y = quantile(Top_ThicknessData$Thickness....Y,na.rm = T,probs = 0.25) #First Quartile
Third_Quartile_Thickness_Y = quantile(Top_ThicknessData$Thickness....Y,na.rm = T,probs = 0.75) #Third Quartile
Maximum_Thickness_Y = max(Top_ThicknessData$Thickness....Y,na.rm = T) #Maximum
Count_Thickness_Y = length(Top_ThicknessData$Thickness....Y) #Number of data values

Data_Attributes_Thickness_Y = c("Mean", "Median", "Standard Deviation", "Kurtosis", "Skewness", "Minimum", "First Quartile", "Third Quartile", "Maximum", "Number of data values") #Summary Table row names
Value_Attributes_Thickness_Y = c(Mean_Thickness_Y, Median_Thickness_Y, StandardDeviation_Thickness_Y, Kurtosis_Thickness_Y, Skewness_Thickness_Y, Minimum_Thickness_Y, First_Quartile_Thickness_Y, Third_Quartile_Thickness_Y, Maximum_Thickness_Y, Count_Thickness_Y) #Summary Table row values 
SummaryTable_Thickness_Y = data.frame(Data_Attributes_Thickness_Y,Value_Attributes_Thickness_Y) #Constructing the Summary Table
names(SummaryTable_Thickness_Y) = c("Data Attribute","Attribute Value for Thickness of Horizon Y")
write.table(SummaryTable_Thickness_Y,file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Top-Thickness Analysis//Tables//SummaryTable_Thickness_Y.csv",sep = ",",row.names = F)



jpeg("C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Top-Thickness Analysis//Plots//Thickness_Y RF-N-K Distributions.jpg")
porohist_Thickness_Y=hist(Top_ThicknessData$Thickness....Y,plot=F, breaks = 9)
porohist_Thickness_Y$counts = porohist_Thickness_Y$counts/sum(porohist_Thickness_Y$counts)
porodens_Thickness_Y = density(Top_ThicknessData$Thickness....Y,na.rm = T,from = 0, to=100)
porodens_Thickness_Y$y=porodens_Thickness_Y$y*10
x = Top_ThicknessData$Thickness....Y
poronorm_Thickness_Y = dnorm(x,mean = Mean_Thickness_Y,sd=StandardDeviation_Thickness_Y)
poronorm_Thickness_Y = poronorm_Thickness_Y*10
plot(porohist_Thickness_Y, freq = T, col = "grey", main="Relative Frequency, Normal and Kernel Distributions of Thickness of Horizon Y",ylim=c(0,0.5), xlim = c(0,100), xlab="Thickness of Horizon Y (ft)",ylab = "Probability", las=1)
#text(0.4, 0.2, paste("Mean =", round(Mean_Whole, 4), "\n Std.Dev =", round(StandardDeviation_whole, 4)))
lines(porodens_Thickness_Y, col="blue", type = "l", lty="solid", lwd=2)
curve(dnorm(x,mean = Mean_Thickness_Y,sd=StandardDeviation_Thickness_Y)*10,add = T, col="red", type = "l", lty="dashed", lwd=2)
legend("topright",legend=c("Kernel","Normal"),
       text.col=c("blue","red"),lty=c("solid","dashed"),col=c("blue","red"))
#curve(dnorm(x,mean = Mean_Whole,sd=StandardDeviation_whole),add = T)
dev.off()







jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//compare whole individual and whole averages.jpg")
boxplot(TrimmedData$Porosity,WellAveragesTable$BoreHoleavgs, ylab="Porosity (fraction)", las=1)
dev.off()

jpeg(file = "C://Users//TTOWG//645//1 karia def//2. CU//Projects//Bitumen Recovery//Data Analysis//Explorative Data Analysis//For export to R//Plots//compare X individual Y individual.jpg")
boxplot(TrimmedData$Porosity[TrimmedData$Horizon=="X"],TrimmedData$Porosity[TrimmedData$Horizon=="Y"], ylab="Porosity (fraction)", las=1)
dev.off()





# Lateral Profile Porosity Analysis
ProfileF = TrimmedData$Porosity[TrimmedData$BoreHole=="13"|TrimmedData$BoreHole=="16"|TrimmedData$BoreHole=="50"|TrimmedData$BoreHole=="51"|TrimmedData$BoreHole=="52"]
ProfileG = TrimmedData$Porosity[TrimmedData$BoreHole=="48"|TrimmedData$BoreHole=="17"|TrimmedData$BoreHole=="49"|TrimmedData$BoreHole=="30"|TrimmedData$BoreHole=="31"|TrimmedData$BoreHole=="47"|TrimmedData$BoreHole=="36"]
ProfileH = TrimmedData$Porosity[TrimmedData$BoreHole=="48"|TrimmedData$BoreHole=="10"|TrimmedData$BoreHole=="29"|TrimmedData$BoreHole=="27"]
ProfileI = TrimmedData$Porosity[TrimmedData$BoreHole=="28"|TrimmedData$BoreHole=="15"|TrimmedData$BoreHole=="26"|TrimmedData$BoreHole=="55"|TrimmedData$BoreHole=="27"|TrimmedData$BoreHole=="32"|TrimmedData$BoreHole=="33"|TrimmedData$BoreHole=="37"]
ProfileJ = TrimmedData$Porosity[TrimmedData$BoreHole=="23"|TrimmedData$BoreHole=="24"|TrimmedData$BoreHole=="25"|TrimmedData$BoreHole=="19"]
ProfileK = TrimmedData$Porosity[TrimmedData$BoreHole=="11"|TrimmedData$BoreHole=="20"|TrimmedData$BoreHole=="19"|TrimmedData$BoreHole=="34"|TrimmedData$BoreHole=="35"|TrimmedData$BoreHole=="22"]

ProfileSD = c(sd(ProfileF), sd(ProfileG), sd(ProfileH), sd(ProfileI), sd(ProfileJ), sd(ProfileK))

Profilemeans = c(mean(ProfileF), mean(ProfileG), mean(ProfileH), mean(ProfileI), mean(ProfileJ), mean(ProfileK))




WestData = TrimmedData$Porosity[TrimmedData$BoreHole=="21"|TrimmedData$BoreHole=="49"|TrimmedData$BoreHole=="23"|TrimmedData$BoreHole=="54"|TrimmedData$BoreHole=="28"|TrimmedData$BoreHole=="48"|TrimmedData$BoreHole=="11"|TrimmedData$BoreHole=="24"|TrimmedData$BoreHole=="26"|TrimmedData$BoreHole=="10"|TrimmedData$BoreHole=="20"|TrimmedData$BoreHole=="25"|TrimmedData$BoreHole=="55"|TrimmedData$BoreHole=="29"|TrimmedData$BoreHole=="13"|TrimmedData$BoreHole=="19"|TrimmedData$BoreHole=="45"|TrimmedData$BoreHole=="16"]

EastData = TrimmedData$Porosity[TrimmedData$BoreHole=="34"|TrimmedData$BoreHole=="32"|TrimmedData$BoreHole=="30"|TrimmedData$BoreHole=="50"|TrimmedData$BoreHole=="17"|TrimmedData$BoreHole=="35"|TrimmedData$BoreHole=="46"|TrimmedData$BoreHole=="33"|TrimmedData$BoreHole=="31"|TrimmedData$BoreHole=="51"|TrimmedData$BoreHole=="22"|TrimmedData$BoreHole=="47"|TrimmedData$BoreHole=="52"|TrimmedData$BoreHole=="37"|TrimmedData$BoreHole=="36"]
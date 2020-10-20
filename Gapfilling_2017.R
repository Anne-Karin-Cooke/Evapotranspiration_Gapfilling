
# path = "/run/media/acooke/Elements/data/OSU/GEK/Flux_tower/2016/lvl_combined/" 
#path= '/run/media/acooke/Elements/flux_tower_data/2017/'
path= '/run/media/acooke/Elements/flux_tower_data/2017/combined/'
# path = "/home/acooke/codes/pflotran/flux_tower_data/2018/"
setwd(path)
out.file<-c("")
file.names <- dir(path, pattern =".dat")
for(i in 1:length(file.names)){
  file <- read.table(file.names[i],skip=4, sep=",")
  df <- file
  # df$DateTime <- format(as.POSIXct(df[,1]) ,format = "%Y-%m-%d")
  # df$Time <- format(as.POSIXct(df[,1]) ,format = "%H:%M:%S")
  df$DateTime <- format(as.POSIXct(df[,1]) ,format = "%Y-%m-%d %H:%M:%S")
  
  file <- data.frame(cbind(df$DateTime,df[,78]))
  
  out.file <- rbind(out.file,file)
}

head(out.file)

out.file1 <- out.file
# file.names



out.file1 <- out.file1[-1,]
# out.file1[,2] <- as.character(out.file1[,2])
# out.file1[,2] <- as.numeric(out.file1[,2])
out.file1[,2]<- as.character(out.file1[,2])

out.file1[,2]<- as.numeric(gsub("," ,".", out.file1[,2]))

out.file1[,2]
typeof(out.file1[,2])



for(i in 2:nrow(out.file1) ) {
  if(out.file1[i,2] < 0) {
    out.file1[i,2] <- 0
  }
}
head(out.file1)

Rg <- out.file1


# results_2017 <- read.csv("/run/media/acooke/Elements/data/OSU/GEK/Flux_tower/concat_raw_2017/concat_raw_new/results/eddypro_GEK_full_output_2019-06-19T190839.csv",  skip=2, stringsAsFactors = F)

results_2017 <- read.csv("/run/media/acooke/Elements/data/OSU/GEK/Flux_tower/concat_raw_2017/results/results_new/eddypro_GEK_full_output_2019-06-25T191944",  skip=2, stringsAsFactors = F)
head(results_2017)


# results_2017 <- read.csv("/run/media/acooke/Elements/data/OSU/GEK/Flux_tower/2017_results_june/eddypro_GEK_full_output_2019-06-18T183512.csv",  skip=2, stringsAsFactors = F)
# head(results_2017)

# plot(results_2017[,40]*24, ylim=c(-5, 5))


data_readj <- data.frame(matrix(0, nrow=nrow(results_2017)+1, ncol=14))
# data_readj
colnames(data_readj) <- c("Year", "DoY", "Hour", "NEE","qcNEE" ,"LE","qcLE","H", "qcH","Tair", "Tsoil", "rH", "VPD", "Ustar")
data_readj[1,] <- c("-", "-", "-", "umolm-2s-1","-",	"Wm-2","-",	"Wm-2","-",	"degC",	"degC",	"%",	"hPa", "ms-1")
data_readj

# data_readj$Year <- 
data_readj$DoY[2:nrow(data_readj)] <- as.numeric(round(results_2017[,4]))

data_readj$NEE[2:nrow(data_readj)] <- as.numeric(results_2017[,14])
data_readj$qcNEE[2:nrow(data_readj)] <- as.numeric(results_2017[,15])
data_readj$LE[2:nrow(data_readj)] <- as.numeric(results_2017[,12])
data_readj$qcLE[2:nrow(data_readj)] <- as.numeric(results_2017[,13])
data_readj$H[2:nrow(data_readj)] <- as.numeric(results_2017[,10])
data_readj$qcH[2:nrow(data_readj)] <- as.numeric(results_2017[,11])
data_readj$rH[2:nrow(data_readj)] <- as.numeric(results_2017[,45])
data_readj$Ustar[2:nrow(data_readj)] <- as.numeric(results_2017[,60])
data_readj$VPD[2:nrow(data_readj)] <- as.numeric(results_2017[,46]/100) # transform Pa into hPa

data_readj$Tair[2:nrow(data_readj)] <- as.numeric(results_2017[,35])-273.15
# data_readj$Rg[2:nrow(data_readj)] <- "NA"
data_readj$Tsoil[2:nrow(data_readj)] <- "NA"

data_readj$Year[2:nrow(data_readj)] <- as.numeric(2017)
data_readj$Hour[2:nrow(data_readj)] <- as.numeric(substr(as.character(results_2017[,3]), 1,2))
data_readj$Hour[seq(2, length(data_readj$Hour), 2)] <- as.numeric(data_readj$Hour[seq(2, length(data_readj$Hour), 2)])+0.5  
# data_readj$Hour[1] <- "-"


data_readj$qcLE[2:nrow(data_readj)] <- as.numeric(data_readj$qcLE[2:nrow(data_readj)])
data_readj$qcH[2:nrow(data_readj)] <- as.numeric(data_readj$qcH[2:nrow(data_readj)])
data_readj$qcNEE[2:nrow(data_readj)] <- as.numeric(data_readj$qcNEE[2:nrow(data_readj)])

head(data_readj)

for(i in 2:nrow(data_readj)) {
  
  if(data_readj$qcLE[i] > 1)
  {
    data_readj <- data_readj[-i,]
  }
}

for(i in 2:nrow(data_readj)) {
  
  if(data_readj$qcH[i] > 1)
  {
    data_readj <- data_readj[-i,]
  }
}

for(i in 2:nrow(data_readj)) {
  
  if(data_readj$qcNEE[i] > 1)
  {
    data_readj <- data_readj[-i,]
  }
}



data_readj <- data_readj[,-c(5,7,9)]


# plot(data_readj$LE, ylim=c(-1000,1000))


data_readj <- data.frame(data_readj)
# data_readj <- data_readj[45:nrow(data_readj),]
data_readj[1,] <- c("-", "-", "-", "umolm-2s-1",	"Wm-2",	"Wm-2",	"degC",	"degC",	"%",	"hPa", "ms-1")
head(data_readj)
nrow(data_readj)

data_readj <- data.frame(data_readj)

write.table(data_readj, "/home/acooke/codes/pflotran/flux_tower_data/2017/2017_new.txt", row.names = F)

typeof(data_readj)


# REddyProc package
library("REddyProc")
# EddyData.F <- if (length(data_readj)) fLoadTXTIntoDataframe(data_readj) else
# or use example dataset in RData format provided with REddyProc
EddyData.F <- fLoadTXTIntoDataframe("/home/acooke/codes/pflotran/flux_tower_data/2017/2017_new.txt")
EddyData.F

typeof(EddyData.F$Hour)

# EddyData.F <- data.frame(EddyData.F)

EddyData.F$Hour <- as.character(EddyData.F$Hour)

EddyData.F$Hour <- as.numeric(gsub("," ,".", EddyData.F$Hour))

#+++ Add time stamp in POSIX time format
EddyDataWithPosix.F <- fConvertTimeToPosix(EddyData.F, 'YDH',Year = 'Year'
                                           ,Day = 'DoY',Hour = 'Hour')


seqdt<- seq(
  from=as.POSIXct("2017-01-01 0","%Y-%m-%d %H"),
  to=as.POSIXct("2017-12-31 23", "%Y-%m-%d %H"),
  by="30 min"
)

filling2017 <- data.frame(matrix(NA, nrow=length(seqdt)+1, ncol=12))
filling2017[1,] <- c("-", "-", "-","-", "umolm-2s-1",	"Wm-2",	"Wm-2",		"degC",	"degC",	"%",	"hPa", "ms-1")
colnames(filling2017) <- c("DateTime","Year", "DoY", "Hour", "NEE", "LE","H",  "Tair", "Tsoil", "rH", "VPD", "Ustar")

head(filling2017)
filling2017$Year[2:nrow(filling2017)] <- as.numeric(2017)
filling2017$DateTime[2:nrow(filling2017)] <- as.character(seqdt)
filling2017$Hour[2:nrow(filling2017)] <- round(as.numeric((format(strptime(filling2017$DateTime[2:nrow(filling2017)],"%Y-%m-%d %H"),'%H'))))
filling2017$Hour[seq(3, length(filling2017$Hour), 2)] <- as.numeric(filling2017$Hour[seq(3, length(filling2017$Hour), 2)])+0.5  
filling2017$DoY[2:nrow(filling2017)] <-round(as.numeric((format(strptime(filling2017$DateTime[2:nrow(filling2017)],"%Y-%m-%d %H"),'%j'))))
filling2017$Year[2:nrow(filling2017)] <- as.numeric(2017)
head(filling2017)

filling2017[1,]<- c("-", "-", "-", "-", "umolm-2s-1","Wm-2",	"Wm-2",		"degC",	"degC",	"%",	"hPa", "ms-1")

write.table(filling2017, "/home/acooke/codes/pflotran/flux_tower_data/2017/filling2017_new.txt", row.names = F)


filling2017.F <- fLoadTXTIntoDataframe("/home/acooke/codes/pflotran/flux_tower_data/2017/filling2017_new.txt")

filling2017.F$Hour <- as.character(filling2017.F$Hour)

filling2017.F$Hour <- as.numeric(gsub("," ,".", filling2017.F$Hour))


filling2017.F <- fConvertTimeToPosix(filling2017.F,'YDH',Year = 'Year'
                                     ,Day = 'DoY',Hour = 'Hour')
head(filling2017.F)

filling2017.F <- filling2017.F[,-2]


# ts.tmp = merge(EddyDataWithPosix.F , filling.F, by="Ustar", all=T)

# ts.tmp <-ts.tmp[,1:13]

ts.tmp = rbind(EddyDataWithPosix.F , filling2017.F)
# by="DateTime", all=F

ts.tmp.sorted <-  ts.tmp[order(ts.tmp$DateTime), ]
ts.tmp.sorted
# colnames(ts.tmp) <- c("DateTime","Year", "DoY", "Hour", "NEE", "LE","H", "Rg", "Tair", "Tsoil", "rH", "VPD", "Ustar")

head(ts.tmp)

df$DateTime <- format(as.POSIXct(df[,1]) ,format = "%Y-%m-%d %H:%M:%S") 
ts.tmp.sorted

for(i in 2:(nrow(ts.tmp.sorted)-1)) {
  if(ts.tmp.sorted$Hour[i] == ts.tmp.sorted$Hour[i+1]) {
    ts.tmp.sorted <- ts.tmp.sorted[-(i+1),]
  }
}

colnames(Rg) <- c("DateTime", "Rg")
wRg = merge(ts.tmp.sorted, Rg, by="DateTime", all=T)



for (i in 1:(nrow(wRg)-1)){
  wRg$Rg[i] <- wRg$Rg[i+1]
}


for (i in 1: nrow(wRg)){
  
  if(is.na(wRg$Hour[i])) {
    wRg <- wRg[-i,]
  }
}


wRg1 <- wRg

tail(wRg)

tail(Rg)
# typeof(wRg$Hour)

my_data <- wRg[, c(1,2,3,4,5,6,7,13,8,9,10,11,12)]

tail(my_data1)


# my_data[ 552:556,]
my_data1 <- my_data

for (i in 1:(nrow(my_data1)-1)) {
  if (my_data1$DateTime[i] == my_data1$DateTime[i+1]) {
    if(is.na(my_data1$Rg[i+1])){
    }
    else
    {
      my_data1$Rg[i] <- my_data1$Rg[i+1]
    }
    my_data1 <- my_data1[-(i+1),]
    
  }
  
}
# i
# # 
# my_data1[17485,]
# 
# [17410:17520,]
# 
# 
# 
# my_data1[ 552:556,]
# 
# my_data1[4127:4131,]
# nrow(my_data)
# is.na(my_data1$Hour)
# my_data[1320:1323,]
# my_data$Rg[1321] <- my_data$Rg[1322]
# my_data <- my_data[-1322,]
# # # 
# my_data[4033:4036,]
# my_data <- my_data[-4034,] 
# my_data <- my_data[-4035,] 

# my_data1 <- my_data1[-(17521:nrow(my_data1)),]

# my_data1[17520:17533,]
# 
# 
# 
# seqfix<- seq(
#   from=as.POSIXct("2017-12-31 00:30:00","%Y-%m-%d %H"),
#   to=as.POSIXct("2017-12-31 22:00:00", "%Y-%m-%d %H"),
#   by="30 min"
# )
# fix<- data.frame(matrix(NA, nrow=length(seqfix), ncol=13))
# 
# colnames(fix) <- c("DateTime","Year", "DoY", "Hour", "NEE", "LE","H",  "Rg", "Tair", "Tsoil", "rH", "VPD", "Ustar")

# rep2017 <- rep(2017, length(seqfix))
# repdoy <- rep(366, length(seqfix))
# 
# hourfix <- seq(0.5,22, by=.5)
# fix$DateTime <- seqfix
# fix$DateTime
# fix$Year <- rep2017
# fix$DoY <- repdoy
# fix$Hour <- hourfix
# fix
# 
# data2 <- my_data1
# 
# data3 <- rbind(data2, fix)
# 
# ncol(fix)
# ncol(data2)
# 
# tail(data3)
# head(fix)
# data3.sorted <-  data3[order(data3$DateTime), ]
# tail(data3)
# 
# head(data3.sorted)
# tail(data3.sorted)


# 17521, 17522, 17524, 17526, 17527, 17528, 17530, 17531, 17532

EddyProc.C <- sEddyProc$new('2017',my_data1, 
                            c('NEE','LE','Rg','Tair','VPD', 'Ustar'))


# EddyProc.C$sPlotFingerprintY('NEE', Year.i = 2017)

#
# uStarTh <- EddyProc.C$sEstUstarThresholdDistribution(
#   nSample = 100L, probs = c(0.05, 0.5, 0.95)) 

EddyProc.C$sEstimateUstarScenarios(
  nSample = 100L, probs = c(0.05, 0.5, 0.95))
EddyProc.C$sGetEstimatedUstarThresholdDistribution()


EddyProc.C$sGetUstarScenarios()
#filter(uStarTh, aggregationMode == "year")
# select(uStarTh, -seasonYear)
#

# uStarThAnnual <- usGetAnnualSeasonUStarMap(uStarTh)[-2]
# uStarSuffixes <- colnames(uStarThAnnual)[-1]
# print(uStarThAnnual)
# # Example_DETha98

EddyProc.C$sMDSGapFillUStarScens('LE')
# UstarThres.df = uStarThAnnual,
# UstarSuffix.V.s = uStarSuffixes,
# FillAll = TRUE
# EddyProc.C$sMDSGapFillAfterUStarDistr('LE',
#                                       UstarThres.df = uStarThAnnual,
#                                       UstarSuffix.V.s = uStarSuffixes,
#                                       FillAll = TRUE                                      

# )

# getwd()

grep("LE_.*_f$",names(EddyProc.C$sExportResults()), value = TRUE)
grep("LE_.*_fsd$",names(EddyProc.C$sExportResults()), value = TRUE)

EddyProc.C$sPlotFingerprintY('LE_U50_f', Year.i = 2017)

FilledEddyData.F <- EddyProc.C$sExportResults()
# CombinedData.F <- cbind(EddyData.F, FilledEddyData.F)
# fWriteDataframeToFile(CombinedData.F, 'DE-Tha-Results.txt', Dir.s = tempdir())

# EddyProc.C$sExportResults()
# FilledEddyData.F$LE_U05_f

data_filled <- my_data1
data_filled$LE <- FilledEddyData.F$LE_U05_f

# EddyProc.C$sMDSGapFillAfterUStarDistr('Tair',
#                                       UstarThres.df = uStarThAnnual,
#                                       UstarSuffix.V.s = uStarSuffixes,
#                                       FillAll = TRUE
# )
EddyProc.C$sMDSGapFillUStarScens('Tair')



# getwd()

grep("Tair_.*_f$",names(EddyProc.C$sExportResults()), value = TRUE)
grep("Tair_.*_fsd$",names(EddyProc.C$sExportResults()), value = TRUE)

FilledEddyData.F <- EddyProc.C$sExportResults()

data_filled$Tair <- FilledEddyData.F$Tair_U50_f

plot(data_filled$Tair)
head(data_filled)
######

ET <- data.frame(matrix(0, nrow=nrow(data_filled), ncol=3))


colnames(ET) <- c("DateTime","Tair_degC","LE_Wm.2" )

ET$DateTime <- data_filled$DateTime
ET$Tair_degC <- data_filled$Tair
ET$LE_Wm.2 <- data_filled$LE

ET_2017 <- ET

head(ET)
for (i in 1:nrow(ET)) {
  ET$lambda[i] <- 1000*(3147.5-(2.37*(as.numeric(ET$Tair_degC[i])+273.15)))
}

for (i in 1:nrow(ET)) {
  ET$ETp[i] <- as.numeric(ET$LE_Wm.2[i])/as.numeric(ET$lambda[i])
}

head(ET)

plot(ET$ETp)

# ET2


ET$Date <- as.Date(ET$DateTime)
ET <- data.frame(ET)




ET_daily <- aggregate(ETp~Date, ET, sum)
ET_daily$ETp <- ET_daily$ETp/1000 # divided by volumetric mass of water
ET_daily$ETp <- (ET_daily$ETp)*86400
plot(ET_daily$Date,ET_daily$ETp, type="b", ylab="ETa [mm]", xlab="2017", main="Flux tower after Gapfilling")
abline(h=0, col="red")


### ET BEFORE GAPFILLING


ET1 <- data.frame(matrix(0, nrow=nrow(data_filled), ncol=3))


colnames(ET1) <- c("DateTime","Tair_degC","LE_Wm.2" )

# data_readj <- fConvertTimeToPosix(data_readj,'YDH',Year = 'Year'
#                                      ,Day = 'DoY',Hour = 'Hour')
# 
# head(data_readj)
# as.POSIXct(data_readj$Year,data_readj$,"%Y-%m-%d %H"),

ET1$DateTime <- my_data1$DateTime
ET1$Tair_degC <- my_data1$Tair
ET1$LE_Wm.2 <- my_data1$LE

ET1_2018 <- ET1

head(ET1)
for (i in 1:nrow(ET1)) {
  ET1$lambda[i] <- 1000*(3147.5-(2.37*(as.numeric(ET1$Tair_degC[i])+273.15)))
}

for (i in 1:nrow(ET1)) {
  ET1$ETp[i] <- as.numeric(ET1$LE_Wm.2[i])/as.numeric(ET1$lambda[i])
}

head(ET1)

plot(ET1$ETp)

# ET2


ET1$Date <- as.Date(ET1$DateTime)
ET1 <- data.frame(ET1)


min(Rg$Rg)


ET_daily1 <- aggregate(ETp~Date, ET1, sum)
ET_daily1$ETp <- ET_daily1$ETp/1000 # divided by volumetric mass of water
ET_daily1$ETp <- (ET_daily1$ETp)*86400
plot(ET_daily1$Date,ET_daily1$ETp, type="l", ylab="ETa [mm]", xlab="2017", main="Flux tower after/before Gapfilling", ylim=c(-0.15, 0.4))
abline(h=0, col="red")
lines(ET_daily$Date,ET_daily$ETp, col="blue")
legend("topright",legend=c("after", "before"), 
       text.col=c("blue","black"),cex=c(0.7),pch=c(15),col=c("blue","black"),
       bty = "n")
# 
# 
# data_readj$lambda <- c(0)
# data_readj$ETp <- c(0)
# data_readj <- data_readj[-2,]
# head(data_readj)
# 
# for (i in 1:nrow(data_readj)) {
#   data_readj$lambda[i] <- 1000*(3147.5-(2.37*(as.numeric(data_readj$Tair[i])+273.15)))
# }
# 
# for (i in 1:nrow(data_readj)) {
#   data_readj$ETp[i] <- as.numeric(data_readj$LE[i])/as.numeric(data_readj$lambda[i])
# }
# 
# plot(data_readj$ETp)
# 
# data_readj$ETp <- data_readj$ETp/1000 # divided by volumetric mass of water
# data_readj$ETp <- (data_readj$ETp)*86400
# 
# 
# plot(data_readj$ETp)
# 
# 
# ET_daily <- aggregate(ETp~Date, ET, sum)
# ET_daily$ETp <- ET_daily$ETp/1000 # divided by volumetric mass of water
# ET_daily$ETp <- (ET_daily$ETp)*86400
# plot(ET_daily$Date,ET_daily$ETp, type="l")
# abline(h=0, col="red")




# 
# data_readj$Tair
# 
# boxplot(ET_daily$ET)
# ET_daily


# library(readr)
# pluvio_jass_pub_21_09 <- read_csv("/home/acooke/data/OSU/RainfallGEK20172017/pluvio_l2_day_jasse_27_adj.csv")
# pluvio <- pluvio_jass_pub_21_09
# pluvio <- as.data.frame(pluvio)
# pluvio
# library("Hmisc")
# names(pluvio) <- c("Date", "daily_rainfall")
# pluvio$date <- as.character(pluvio$Date)
# 
# # pluvio <- pluvio[which(pluvio$date > "2017-12-01"),]
# pluvio[,1] <- as.Date(pluvio[,1])
# ET_daily[,1] <- as.Date(ET_daily[,1])
# 
# 
# 
# final_eff_rf <- merge(ET_daily, pluvio, by="Date")
# 
# final_eff_rf$eff_rf <- final_eff_rf$daily_rainfall + final_eff_rf$ET
# plot(final_eff_rf$Date,final_eff_rf$eff_rf, type="l")
# 
# abline(h=0, col="red")
# head(final_eff_rf)
# nrow(final_eff_rf)
# 
# plot(cumsum(final_eff_rf$eff_rf))
# 
# # final_eff_rf_2017 <- final_eff_rf
# # ET_daily_2017 <- ET_daily
# 
# # write.table(ET_daily_2017, "ET_daily_2017.txt", row.names = F)
# # write.table(my_data, "data_2017.txt", row.names = F)
# # write.table(data_filled, "data_filled_2017.txt", row.names = F)
# 
# 
# ET_daily_2017



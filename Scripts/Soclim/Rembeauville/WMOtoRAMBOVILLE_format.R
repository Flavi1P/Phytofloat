###########################################################################################
# Routine developpee par C. Schmechtig Antoine Poteau et Henry Bittig
#
# March 2016 : Read Coriolis Data 
#############################################################################################
library(marmap)
require(oce)
require(RNetCDF)
library(akima)
library(zoo)




uf=commandArgs()

floats  <- uf[2]
WMO  <- uf[3]
DAC  <- uf[4]


#WMO <- c('6901474')
#DAC <- c("coriolis")
#floats = c('006b')


#WMO <- c('6901585','6901583','6901004','6902739','6902738')
#DAC <- c("coriolis","coriolis","coriolis","coriolis","coriolis")
#floats = c('049b','036b','037c','107c','104c')

#WMO <- c('6901585')
#DAC <- ("coriolis")
#floats = c('049b')



ii <- 0
listWMO <- unique(WMO)
for (WMOWMO in listWMO) {
print(WMOWMO)
AA <- NULL

ii <- ii + 1
GDAC <- DAC[which(  WMO  == WMOWMO)  [1]]
print(WMOWMO)

prof_bio_core=paste("/home/admt/GDAC/",GDAC,"/",WMOWMO,"/",WMOWMO,"_Mprof.nc",sep="")
print(WMOWMO)

if (file.exists(prof_bio_core) == TRUE) {  
prof_bio_core=open.nc(prof_bio_core,readunlim=FALSE) 
PARAM <-  var.get.nc(prof_bio_core,"PARAMETER")
PARAM_MODE <-  var.get.nc(prof_bio_core,"PARAMETER_DATA_MODE")
PRES=var.get.nc(prof_bio_core,"PRES")
sans <- NULL
for (i in seq(1,dim(PRES)[2])){
sansNA <- length(which(!is.na(PRES[,i]==TRUE)))
if (length(unique(PRES[,i])) < 5)  sans <- c(sans,i)
if (sansNA < 100)  sans <- c(sans,i)
}


print(WMOWMO)

CHLA =var.get.nc(prof_bio_core,"CHLA")
QCC_CHLA <- as.numeric(unlist(strsplit(var.get.nc(prof_bio_core,"CHLA_QC"), split = "") ))

print(QCC_CHLA)
CHLA[which(QCC_CHLA==4)]  <- NaN

BBP700 =var.get.nc(prof_bio_core,"BBP700")
QCC_BBP700 <- as.numeric(unlist(strsplit(var.get.nc(prof_bio_core,"BBP700_QC"), split = "") ))
BBP700[which(QCC_BBP700==4)]  <- NaN

CP660 =var.get.nc(prof_bio_core,"CP660")
QCC_CP660 <- as.numeric(unlist(strsplit(var.get.nc(prof_bio_core,"CP660_QC"), split = "") ))
CP660[which(QCC_CP660==4)]  <- NaN

ED490 =var.get.nc(prof_bio_core,"DOWN_IRRADIANCE490")
QCC_ED490 <- as.numeric(unlist(strsplit(var.get.nc(prof_bio_core,"DOWN_IRRADIANCE490_QC"), split = "") ))
ED490[which(QCC_ED490==4)]  <- NaN


DOWNWELLING_PAR =var.get.nc(prof_bio_core,"DOWNWELLING_PAR")
QCC_DOWNWELLING_PAR <- as.numeric(unlist(strsplit(var.get.nc(prof_bio_core,"DOWNWELLING_PAR_QC"), split = "") ))
DOWNWELLING_PAR[which(QCC_DOWNWELLING_PAR==4)]  <- NaN


LAT=var.get.nc(prof_bio_core,"LATITUDE")[-sans]
LON=var.get.nc(prof_bio_core,"LONGITUDE")[-sans]

JULD=var.get.nc(prof_bio_core,"JULD")[-sans]
PRES=var.get.nc(prof_bio_core,"PRES")[,-sans]
PSAL=try(var.get.nc(prof_bio_core,"PSAL")[,-sans])
CYCLE =var.get.nc(prof_bio_core,"CYCLE_NUMBER")[-sans]
CHL=CHLA[,-sans]

BBP700=BBP700[,-sans]
TEMP=try(var.get.nc(prof_bio_core,"TEMP")[,-sans])
OX=try(var.get.nc(prof_bio_core,"DOXY")[,-sans])
CP=CP660[,-sans]

PAR=DOWNWELLING_PAR[,-sans]

ED490=ED490[,-sans]

 


print('ffff')


for (ivv in seq(1:dim(CHL)[2])){
VVV <- which(is.na(CHL[,ivv])==FALSE)
try(CHL[VVV,ivv] <- rollmedian(CHL[VVV,ivv],5,na.pad=TRUE))
try(BBP700[VVV,ivv] <- rollmedian(BBP700[VVV,ivv],5,na.pad=TRUE))

VVV <- which(is.na(CHL[,ivv])==FALSE)
try(CHL[VVV,ivv] <- rollmedian(CHL[VVV,ivv],7,na.pad=TRUE))
try(BBP700[VVV,ivv] <- rollmedian(BBP700[VVV,ivv],7,na.pad=TRUE))

VVV <- which(is.na(ED490[,ivv])==FALSE)
try(ED490[VVV,ivv] <- rollmedian(ED490[VVV,ivv],5,na.pad=TRUE))
try(PAR[VVV,ivv] <- rollmedian(PAR[VVV,ivv],5,na.pad=TRUE))

VVV <- which(is.na(CHL[,ivv])==FALSE)
try(ED490[VVV,ivv] <- rollmedian(ED490[VVV,ivv],7,na.pad=TRUE))
try(PAR[VVV,ivv] <- rollmedian(PAR[VVV,ivv],7,na.pad=TRUE))

VVV <- which(is.na(CP[,ivv])==FALSE)
try(CP[VVV,ivv] <- rollmedian(CP[VVV,ivv],5,na.pad=TRUE))

VVV <- which(is.na(CP[,ivv])==FALSE)
try(CP[VVV,ivv] <- rollmedian(CP[VVV,ivv],7,na.pad=TRUE))
}





PRES_INT <- seq(1,1000)
CHL_INT <- matrix(ncol=dim(CHL)[2],nrow=length(PRES_INT))
PSAL_INT <- matrix(ncol=dim(TEMP)[2],nrow=length(PRES_INT))
TEMP_INT <- matrix(ncol=dim(TEMP)[2],nrow=length(PRES_INT))
OX_INT <- matrix(ncol=dim(OX)[2],nrow=length(PRES_INT))
CP_INT <- matrix(ncol=dim(CP)[2],nrow=length(PRES_INT))
PAR_INT <- matrix(ncol=dim(PAR)[2],nrow=length(PRES_INT))
ED490_INT <- matrix(ncol=dim(ED490)[2],nrow=length(PRES_INT))
BBP_INT <- matrix(ncol=dim(CHL)[2],nrow=length(PRES_INT))

for (i in seq(1,dim(CHL)[2])) {
if (length(na.omit(CHL[,i]))>100) CHL_INT[,i] <-  as.double(approx(PRES[,i],CHL[,i],PRES_INT, rule = 2)$y)
if (length(na.omit(PSAL[,i]) )>100) PSAL_INT[,i] <-  as.double(approx(PRES[,i],PSAL[,i],PRES_INT, rule = 2)$y)
if (length(na.omit(TEMP[,i]))>100) TEMP_INT[,i] <-  as.double(approx(PRES[,i],TEMP[,i],PRES_INT, rule = 2)$y)
if (length(na.omit(OX[,i]) )>100) OX_INT[,i] <-  as.double(approx(PRES[,i],OX[,i],PRES_INT, rule = 2)$y)
if (length(na.omit(CP[,i]))>100) CP_INT[,i] <-  as.double(approx(PRES[,i],CP[,i],PRES_INT, rule = 2)$y)
if (length(na.omit(PAR[,i]) )>100) PAR_INT[,i] <-  as.double(approx(PRES[,i],PAR[,i],PRES_INT, rule = 2)$y)
if (length(na.omit(ED490[,i]))>100) ED490_INT[,i] <-  as.double(approx(PRES[,i],ED490[,i],PRES_INT, rule = 2)$y)
if (length(na.omit(BBP700[,i]) )>100) BBP_INT[,i] <-  as.double(approx(PRES[,i],BBP700[,i],PRES_INT, rule = 2)$y)
}



PSAL_2 <- NULL
CHL_2 <- NULL
BBP_2 <- NULL
CP_2 <- NULL
PAR_2 <- NULL
ED_2 <- NULL
TEMP_2 <- NULL
TIME_2  <- NULL
LAT_2 <- NULL
LON_2 <- NULL
OX_2<- NULL


print('ffff')

iiii <- 0
for (iiv in (unique(JULD))){
iiii <- iiii + 1

if (length(which(JULD == iiv & is.na(colMeans(CP,na.rm=T))==FALSE)) > 0 & length(which(JULD == iiv & is.na(colMeans(PSAL_INT,na.rm=T))==FALSE)) > 0 & length(which(JULD == iiv & is.na(colMeans(CHL_INT,na.rm=T))==FALSE)) >0 ){

TEMP_2 <-    cbind(TEMP_2,TEMP_INT[,which(JULD == iiv & is.na(colMeans(TEMP_INT,na.rm=T))==FALSE)])
PSAL_2 <-    cbind(PSAL_2,PSAL_INT[,which(JULD == iiv & is.na(colMeans(PSAL_INT,na.rm=T))==FALSE)])
CHL_2 <-    cbind(CHL_2,CHL_INT[,which(JULD == iiv & is.na(colMeans(CHL_INT,na.rm=T))==FALSE)])
BBP_2 <-  cbind(BBP_2,BBP_INT[,which(JULD == iiv & is.na(colMeans(BBP_INT,na.rm=T))==FALSE)])
CP_2 <-  cbind(CP_2,CP_INT[,which(JULD == iiv & is.na(colMeans(CP_INT,na.rm=T))==FALSE)])
PAR_2 <-  cbind(PAR_2,PAR_INT[,which(JULD == iiv & is.na(colMeans(PAR_INT,na.rm=T))==FALSE)])
ED_2 <-  cbind(ED_2,ED490_INT[,which(JULD == iiv & is.na(colMeans(ED490_INT,na.rm=T))==FALSE)])
OX_2 <-  cbind(OX_2,OX_INT[,which(JULD == iiv & is.na(colMeans(OX_INT,na.rm=T))==FALSE)])

TIME_2 <- c(TIME_2,iiv)
LAT_2 <- c(LAT_2,LAT[iiii])

LON_2 <- c(LON_2,LON[iiii])


}
}

print('ffff')


PSAL_2[is.na(PSAL_2)] <- 'nan'
CHL_2[is.na(CHL_2)] <- 'nan'
TEMP_2[is.na(TEMP_2)] <- 'nan'
OX_2[is.na(OX_2)] <- 'nan'
BBP_2[is.na(BBP_2)] <- 'nan'
CP_2[is.na(CP_2)] <- 'nan'
ED_2[is.na(ED_2)] <- 'nan'
PAR_2[is.na(PAR_2)] <- 'nan'
LON_2[is.na(LON_2)] <- 'nan'
LAT_2[is.na(LAT_2)] <- 'nan'

densiteRHO <- PSAL_2

for (ivv in seq(1:dim(TEMP_2)[2])) {
densiteRHO[,ivv] <- swSigmaTheta(as.double(PSAL_2[,ivv]),as.double(TEMP_2[,ivv]),PRES_INT)
}


CP_3 <- NULL

datadir =  paste('/home/admt/PROGRAM_AP/CODE_REMBAUVILLLE/data/',floats[ii],"/", sep='')
write.table(PSAL_2, paste(datadir,"SAL.txt",sep=""),quote = FALSE, sep = ",",row.names = FALSE,col.names = FALSE)
write.table(CHL_2, paste(datadir,"CHLA.txt",sep=""),quote = FALSE, sep = ",",row.names = FALSE,col.names = FALSE)
write.table(TEMP_2, paste(datadir,"TEMP.txt",sep=""),quote = FALSE, sep = ",",row.names = FALSE,col.names = FALSE)
write.table(OX_2, paste(datadir,"OX.txt",sep=""),quote = FALSE, sep = ",",row.names = FALSE,col.names = FALSE)
write.table(BBP_2, paste(datadir,"BBP.txt",sep=""),quote = FALSE, sep = ",",row.names = FALSE,col.names = FALSE)
write.table(CP_2, paste(datadir,"CP.txt",sep=""),quote = FALSE, sep = ",",row.names = FALSE,col.names = FALSE)
write.table(ED_2, paste(datadir,"ED490.txt",sep=""),quote = FALSE, sep = ",",row.names = FALSE,col.names = FALSE)
write.table(PAR_2, paste(datadir,"PAR.txt",sep=""),quote = FALSE, sep = ",",row.names = FALSE,col.names = FALSE)
write.table(LON_2, paste(datadir,"LON.txt",sep=""),quote = FALSE, sep = ",",row.names = FALSE,col.names = FALSE)
write.table(LAT_2, paste(datadir,"LAT.txt",sep=""),quote = FALSE, sep = ",",row.names = FALSE,col.names = FALSE)
write.table(densiteRHO, paste(datadir,"SIGMA.txt",sep=""),quote = FALSE, sep = ",",row.names = FALSE,col.names = FALSE)

timedec <- ((abs(TIME_2 - trunc(TIME_2)) *24 ) )

TIME_H <- format(strptime( paste(floor(timedec), round((timedec-floor(timedec))*60), sep=":"), format="%H:%M"),"%H")
TIME_MM <- format(strptime( paste(floor(timedec), round((timedec-floor(timedec))*60), sep=":"), format="%H:%M"),"%M")
TIME_Y <- format(as.Date(trunc(TIME_2), origin="1950-01-01 00:00:00 UTC"),"%Y")
TIME_M <- format(as.Date(trunc(TIME_2), origin="1950-01-01 00:00:00 UTC"),"%m")
TIME_D <- format(as.Date(trunc(TIME_2), origin="1950-01-01 00:00:00 UTC"),"%d")

TIME_MM[is.na(TIME_MM)] <- 0
TIME_H[is.na(TIME_H)] <- 0


print(dim(CP_2))
print(dim(PSAL_2))
print(dim(BBP_2))
print(dim(CHL_2))

write.table(paste(TIME_Y,TIME_M,TIME_D,TIME_H,TIME_MM,sep=","), paste(datadir,"TIME.txt",sep=""),quote = FALSE, sep = ",",row.names = FALSE,col.names = FALSE)




}
}

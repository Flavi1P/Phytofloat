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
require(oce)
require(fields)
require(mapdata)
require(mapproj)
require(maps)
require(ncdf4)
require(R.utils)
require(gmt)
require(plotuti)
require('oce')
require("stringr")
verbose <- Verbose(threshold=0)

library(marmap)
library(ncdf4)
library(reticulate)
uf=commandArgs()
library(RcppCNPy)

floats  <- uf[2]
WMO  <- uf[3]
DAC  <- uf[4]

WMO <- c('6901585','6901583','6901004','6902739','6902738')
floats = c('049b','036b','037c','107c','104c')
#floats = c('048b')
#floats = c('049b')
#floats <- '107c'
for (i in floats) {
pdf(file= paste("Output/",i,".pdf",sep=""),onefile=TRUE)
#jpeg(file= paste("/var/www/oao/bioargo/PHP/AP/REMBAU/COMP_",i,".jpeg",sep="")  , width = 3000, height = 3000 )
par(mfrow = c(4, 4))  # 3 rows and 2 columns

datadir =  paste('Data/Soclim/data/',i,"/", sep='')

A<-read.table(paste(datadir,"SAL.txt",sep=""),sep="," )
B<-read.table(paste(datadir,"CHLA.txt",sep=""),sep="," )
C<-read.table(paste(datadir,"TEMP.txt",sep=""),sep="," )
D<-read.table(paste(datadir,"OX.txt",sep=""),sep="," )
E<-read.table(paste(datadir,"BBP.txt",sep=""),sep="," )
F<-read.table(paste(datadir,"CP.txt",sep=""),sep="," )
G<-read.table(paste(datadir,"ED490.txt",sep=""),sep="," )
H<-read.table(paste(datadir,"PAR.txt",sep=""),sep="," )
I<-read.table(paste(datadir,"LON.txt",sep=""),sep="," )
J<-read.table(paste(datadir,"LAT.txt",sep=""),sep="," )

AA<-npyLoad(paste(datadir,"SAL.npy",sep=""), dotranspose=FALSE)
AB<-npyLoad(paste(datadir,"CHLA.npy",sep=""), dotranspose=FALSE)
AC<-npyLoad(paste(datadir,"TEMP.npy",sep=""), dotranspose=FALSE)
AD<-npyLoad(paste(datadir,"OX.npy",sep=""), dotranspose=FALSE)
AE<-npyLoad(paste(datadir,"BBP.npy",sep=""), dotranspose=FALSE)
AF<-npyLoad(paste(datadir,"CP.npy",sep=""), dotranspose=FALSE)
AG<-npyLoad(paste(datadir,"ED490.npy",sep=""), dotranspose=FALSE)
AH<-npyLoad(paste(datadir,"PAR.npy",sep=""), dotranspose=FALSE)
AI<-npyLoad(paste(datadir,"LON.npy",sep=""), dotranspose=FALSE)
AJ<-npyLoad(paste(datadir,"LAT.npy",sep=""), dotranspose=FALSE)

#for (iVV in seq(1,dim(B)[2])) plot(B[,iVV],col="red")      
#for (iVV in seq(1,dim(AB)[2])) plot(AB[,iVV])


plot(B,col="red")
plot(AB/2)


dev.off()


}

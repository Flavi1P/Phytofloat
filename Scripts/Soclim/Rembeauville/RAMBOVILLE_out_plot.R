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


uf=commandArgs()

floats  <- uf[2]
WMO  <- uf[3]
DAC  <- uf[4]



#WMO <- c('6901585','6901583','6901004','6902739','6902738')
#floats = c('049b','036b','037c','107c','104c')
#floats = c('048b')
#floats = c('049b')

floats <- '049b'
i <- floats

for (i in floats) {
print(i)
bact <- read.table(paste("Data/Soclim/DATA_TEXT_OUT/bact",i,".txt",sep=""),sep = ",")
bbp <- read.table(paste("Data/Soclim/DATA_TEXT_OUT/bbp",i,".txt",sep=""),sep = ",")
chl <- read.table(paste("Data/Soclim/DATA_TEXT_OUT/chl",i,".txt",sep=""),sep = ",")
cp <- read.table(paste("Data/Soclim/DATA_TEXT_OUT/cp",i,".txt",sep=""),sep = ",")
depth <- read.table(paste("Data/Soclim/DATA_TEXT_OUT/depth",i,".txt",sep=""),sep = ",")[,1]
diat <- read.table(paste("Data/Soclim/DATA_TEXT_OUT/diat",i,".txt",sep=""),sep = ",")
nano <- read.table(paste("Data/Soclim/DATA_TEXT_OUT/nano",i,".txt",sep=""),sep = ",")
pico <- read.table(paste("Data/Soclim/DATA_TEXT_OUT/pico",i,".txt",sep=""),sep = ",")
time <- read.table(paste("Data/Soclim/DATA_TEXT_OUT/time",i,".txt",sep=""),sep = ",")[,1]


time <- as.Date(time, origin="01-01-01 00:00:00 UTC")

diat <- as.double(unlist(diat))
bact <- as.double(unlist(bact))
pico <- as.double(unlist(pico))
nano <- as.double(unlist(nano))
chl <- as.double(unlist(chl))
bbp <- as.double(unlist(bbp))
cp <- as.double(unlist(cp))
 



#pdf(file= paste("/home/admt/PROGRAM_AP/CODE_REMBAUVILLLE/DATA_TEXT_OUT/",i,".pdf",sep="")  , width = 2000, height = 2000 )
jpeg(file= paste("Output/Soclim/Plots/",i,"/",i,".jpeg",sep="")  , width = 2000, height = 2000 )


par(mfrow = c(4, 2))  # 3 rows and 2 columns

icolor <- 0
par(omi=c(0,.5,.5,0))
par(mgp=c(3,1,0))
par(mai=  c(1,1,1,1))


 time2 <- NULL
 for (iiv in time) time2 <- c(time2,rep(iiv,length(depth)))

 depth2 <- NULL
 depth2 <- rep(depth,length(time))


time <-  time2
depth <-  depth2



#for (DATAplot in c("chl")){
DATAplot1 <- 'chl'

for (DATAplot1 in c("chl","bbp","cp","bact", "pico","nano","diat")){
print(DATAplot1)

Lim_depth <- 250

DATAplotALL <-  ((diat)) + ((bact)) +  ((pico)) + ((nano))

#DATAplotALL[which(DATAplotALL == NaN) ] <- 0

if (DATAplot1 == "chl") DATAplot <-  ((chl))
if (DATAplot1 == "bbp") DATAplot <-  ((bbp))
if (DATAplot1 == "cp") DATAplot <-   ((cp))
minQC <- quantile(na.omit(DATAplot),0.001)
macQC <- quantile(na.omit(DATAplot),0.999)
DATAplot[DATAplot <= minQC ]<- minQC
DATAplot[DATAplot >= macQC ]<- macQC

Rrange <- range(DATAplot, na.rm = T)
breaks <- pretty(DATAplot, n = 50)


if (DATAplot1 == "diat") DATAplot <-  (as.double(unlist(diat))/DATAplotALL)
if (DATAplot1 == "nano") DATAplot <-  (as.double(unlist(nano))/DATAplotALL)
if (DATAplot1 == "pico") DATAplot <-  (as.double(unlist(pico))/DATAplotALL)
if (DATAplot1 == "bact") DATAplot <-  (as.double(unlist(bact))/DATAplotALL)




if (DATAplot1 == "diat" |  DATAplot1 == "nano" |  DATAplot1 == "pico"  |  DATAplot1 == "bact") {

minQC <- quantile(na.omit(DATAplot),0.001)
macQC <- quantile(na.omit(DATAplot),0.999)
DATAplot[DATAplot <= minQC ]<- minQC
DATAplot[DATAplot >= macQC ]<- macQC

DATAplot[DATAplot <= 0 ]<- 0.001
DATAplot[DATAplot >= 1 ]<- 0.999

Rrange <- range(seq(0,1,.01), na.rm = T)
breaks <- pretty(seq(0,1,.01), n = 100)

}

cs <- list(cols = tim.colors(length(breaks) - 1),breaks = breaks,name = "",unit = "",labels = seq(1,length(breaks), 5))
cols <- cs.use(DATAplot, cs)


plot(time , depth,ylim = c(Lim_depth,0), xlab = "", ylab = "", pch = 20, lwd="2", cex.lab=1.0, col = cols,font.axis=2,las=1,cex.axis = 1.5,cex.main =4  ,  yaxt="n",  xaxt="n")

title(DATAplot1, cex.main=6  ,line=3 )

abline(h=c(0))
DatePro <- time  
axis(2, col = "black" , col.axis = "black", cex.axis=2, tck=-.02,las=2, font.axis=2 )
mtext("Depth (m)", side=2, line=4.5, cex=3,las=3, col="black")
	summ.dte <- summary(DatePro)
	summ.dte2 <- summary(DatePro)
	min.dte <- min(DatePro,na.rm=T)
	max.dte <- max(DatePro,na.rm=T)

axis(1, as.Date(seq(min.dte, max.dte , length=10)), labels = FALSE)

axis.Date(1, DatePro, at = as.Date(seq(min.dte, max.dte , length=10)), format ="%d/%m/%y",col = "black", col.axis="black", cex.axis=3, tck=-.03,las=2, font.axis=2)
	
digits <- 5
cs.draw(cs, name = "", unit = "",
length =1, width = 0.02, horiz = F, pos = 1.02,
side = 1,  offset = 0, cex = 1,
border = NA, lty = NULL, lwd = par("lwd"), xpd = T, digits = digits )

}
dev.off()
}













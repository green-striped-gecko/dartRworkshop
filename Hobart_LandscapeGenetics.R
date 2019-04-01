## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(tidy = TRUE)
library(tufte)
library(knitr)
#knitr::opts_chunk$set(prompt = TRUE)


## You can find an eletronic pdf version of this script at the following github page:


## ------------------------------------------------------------------------
library(dartR)
library(PopGenReport)
library(raster)



## ---- fig.height=4-------------------------------------------------------

iso <- gl.ibd(possums.gl, projected = TRUE)



## ---- echo=T-------------------------------------------------------------

landscape.sim <- readRDS(system.file("extdata","landscape.sim.rdata", package="dartR"))



## ---- fig.height=4-------------------------------------------------------

xs <- tapply(possums.gl@other$latlong[,"lon"], pop(possums.gl), mean)
ys <- tapply(possums.gl@other$latlong[,"lat"], pop(possums.gl), mean)

plot(landscape.sim)
points(xs, ys, pch=16, col="black", cex=2)
text(xs+10, ys, popNames(possums.gl), col="black", cex=1.5)

coords <- cbind(xs, ys)



## ------------------------------------------------------------------------
eucl <- as.matrix(dist(coords))


## 

## ------------------------------------------------------------------------
cost <- gl.costdistances(landscape = landscape.sim, locs = coords, method = "leastcost", NN=8)


## 

## ------------------------------------------------------------------------
gd <-as.matrix(as.dist(gl.fst.pop(possums.gl, nboots=1)))


## ------------------------------------------------------------------------
library(PopGenReport)
wassermann(gen.mat = gd, eucl.mat = eucl, cost.mats = list(cost=cost),plot = F)

lgrMMRR(gen.mat=gd, cost.mats = list(cost=cost), eucl.mat = eucl)



## ---- fig.height=4-------------------------------------------------------

glc <- gl.genleastcost(possums.gl, fric.raster = landscape.sim, gen.distance = "D", NN = 8, pathtype = "leastcost")
wassermann(gen.mat =gd , eucl.mat = eucl, cost.mats = list(cost=cost),plot = F)


## 

## ---- fig.height=4-------------------------------------------------------

fp <- "https://raw.githubusercontent.com/green-striped-gecko/dartRworkshop/master/data"

poss <- gl.read.dart(file.path(fp,"Possums_SNP.CSV"), file.path(fp, "Possum_Covariate_DArT_180619.csv"), probar = FALSE)

gl.map.interactive(poss)

#subset to only NZ data using longitude
poss.nz <- poss[poss@other$latlong$lon>160,]



## ---- fig.height=6-------------------------------------------------------
ts <- raster(file.path(fp, "TreeScrub.tif"))
rb <- raster(file.path(fp, "RiverBridges.tif"))

par(mfrow=c(2,1))
plot(ts)
plot(rb)


## ---- fig.height=4-------------------------------------------------------

proj4string(rb)
proj4string(ts)




## ------------------------------------------------------------------------

library(rgdal)

nz.pr4 <- proj4string(rb)


xy <- project(as.matrix(poss.nz@other$latlong[,2:1]), proj = nz.pr4)

poss.nz@other$xy <- data.frame(x=xy[,1], y=xy[,2])





## ---- fig.height=4-------------------------------------------------------

#only pupulation 3
possub.nz <- poss.nz[pop(poss.nz)==3,]

set.seed(5) #leave it during the first run
#only 11 individuals from that population
possub.nz <- possub.nz[sample(1:nInd(possub.nz),11, replace = F),1:1000]

gl.map.interactive(possub.nz)


## ------------------------------------------------------------------------
ts40 <- ts
values(ts40) <- (values(ts40))*39+1


## ---- fig.height=4-------------------------------------------------------

system.time(glc <- gl.genleastcost(possub.nz, fric.raster =ts40, gen.distance = "Kosman", NN = 8, pathtype = "leastcost"))

wassermann(eucl.mat = glc$eucl.mat, cost.mat = glc$cost.mats,  gen.mat = glc$gen.mat)



## 

## 

## ---- eval=FALSE---------------------------------------------------------
## 
## #river/bridge resistance
## rb50 <- rb
## values(rb50) <- values(rb50)*200+1
## 
## #combine rivers and forest
## tcrb <- ts
## values(tcrb) <- (values(ts))*39 + values(rb)*99 + 1
## mats <- stack(ts40, rb50, tcrb)
## 
## names(mats)<- c("ts40", "rb50", "tsrb")
## 
## 
## system.time(glc <- gl.genleastcost(possub.nz, fric.raster =mats, gen.distance = "Kosman", NN = 8, pathtype = "leastcost"))
## 
## 
## lgrMMRR(gen.mat = glc$gen.mat, cost.mats = glc$cost.mats, eucl.mat = glc$eucl.mat)
## 
## 


## 

## ---- eval=FALSE---------------------------------------------------------
## 
## library(GGally)
## ss <- sapply(glc$cost.mats, function(x) as.dist(x))
## ggpairs(data.frame(ss))
## 


## ---- eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE----------------
library(secr)
library(PopGenReport)#should already be installed  
library(raster) #should already be installed 



## ---- eval=FALSE, echo=TRUE----------------------------------------------
## install.packages("secr")
## library(secr)
## library(PopGenReport)#should already be installed
## library(raster) #should already be installed
## 


## ------------------------------------------------------------------------
source("https://raw.github.com/green-striped-gecko/lgfun/master/lgfuns.r")


## ---- fig.height=4-------------------------------------------------------
r <- create.resistance(nx = 50, ny = 50, p = 0.5, A=0.5, resVal = 10)



## 

## ---- echo=TRUE, fig.height=4--------------------------------------------
locs <-create.pops(n=8, mindist = 3, landscape = r, plot = TRUE)
locs


## ------------------------------------------------------------------------
para<- list()
#Define populations (dynamics)
para$n.pops=8
para$n.ind=50

para$sex.ratio <- 0.5
#age distribution....

para$n.cov <- 3 
#number of covariates (before the loci in the data.frame, do not change this!!)


## ------------------------------------------------------------------------

#reproduction
para$n.offspring = 2

#migration
para$mig.rate <- 0.2 

#dispersal: exponential dispersal with maximal distance in map units
para$disp.max=100   #average  dispersal of an individual in meters
para$disp.rate = 0.1 #proportion of dispersing individuals

#Define genetics [create 100 SNPs]
para$n.allels <- 2
para$n.loci <- 200  
para$mut.rate <- 0.00001


## ---- fig.height=4-------------------------------------------------------
para$method <- "leastcost" #rSPDdistance, commute
para$NN <- 8  #number of neighbours for the cost distance method

# Initialize simulation of populations from scratch

 landscape<- r  #<-raster(system.file("external/rlogo.grd", package="raster"))

# Define x and y locations
 para$locs <- locs
 #give the population some names 
 rownames(para$locs) <- LETTERS[1:para$n.pops]
  
  
# Create a costdistance matrix 
 
  cost.mat <- gl.costdistances(landscape, para$locs, 
                                          para$method, para$NN) 
  #needed for the simulation
  eucl.mat <- as.matrix(dist(para$locs))  #needed for the analysis later

# Plot your landscape with the populations....
  
  plot(landscape)
  points(para$locs[,1], para$locs[,2], pch=16, cex=2, col="orange")
  text(para$locs[,1],para$locs[,2], row.names(para$locs), cex=1.5)
  
# Check the parameter list
  
  para



## ------------------------------------------------------------------------
simpops <- init.popgensim(para$n.pops, para$n.ind, para$sex.ratio, 
                            para$n.loci, para$n.allels, para$locs, para$n.cov )  


## ------------------------------------------------------------------------
names(simpops)  #the names of the subpopulations


## ------------------------------------------------------------------------
glsp <- pops2gl(simpops, locs =para$locs)
glsp #check the genlight object
summary(glsp)  
gen.mat <- as.matrix(as.dist(gl.fst.pop(glsp, nboots = 1)))
#overall pairwise fst
mean(as.dist(gen.mat))


## 

## ------------------------------------------------------------------------
simpops <- run.popgensim(simpops, steps=3, cost.mat, 
                         n.offspring=para$n.offspring,n.ind=para$n.ind,
                         para$mig.rate, para$disp.max, para$disp.rate, 
                         para$n.allels, para$mut.rate,
                         n.cov=para$n.cov, rec="none")


## 

## ---- eval=FALSE---------------------------------------------------------
## #initialise
## simpops <- init.popgensim(para$n.pops, para$n.ind, para$sex.ratio,
##                             para$n.loci, para$n.allels, para$locs, para$n.cov )
## #Calculate overall Fsts
## glsp <- pops2gl(simpops, locs =para$locs)
## gen.mat <- as.matrix(as.dist(gl.fst.pop(glsp, nboots = 1)))
## #overall pairwise fst
## tempfst <- mean(as.dist(gen.mat))
## 
## glc <- gl.genleastcost(glsp, fric.raster = landscape, gen.distance = "D", NN = 8, pathtype = "leastcost")
## temp.mantel <- wassermann(gen.mat =glc$gen.mat , eucl.mat = glc$eucl.mat,  cost.mats = glc$cost.mats,plot = F)
## 
## #find the p value in temp
## mantelp <- as.numeric(temp.mantel$mantel.tab["1","p"])
## 
## res <- data.frame(generation=0, fst=tempfst, mantelp)
## for (i in 1:10)
## {
## simpops <- run.popgensim(simpops, steps=2, cost.mat,
##                          n.offspring=para$n.offspring,n.ind=para$n.ind,
##                          para$mig.rate, para$disp.max, para$disp.rate,
##                          para$n.allels, para$mut.rate,
##                          n.cov=para$n.cov, rec="none")
## 
## glsp <- pops2gl(simpops)
## gensim.mat <- as.matrix(as.dist(gl.fst.pop(glsp, nboots = 1)))
## #overall pairwise fst
## tempfst <- mean(as.dist(gensim.mat))
## temp.mantel <- wassermann(gen.mat =gensim.mat , eucl.mat = glc$eucl.mat,  cost.mats = glc$cost.mats,plot = F)
## 
## mantelp <- as.numeric(temp.mantel$mantel.tab["1","p"])
## 
## 
## res[i,] <- c(i*2, tempfst, mantelp)
## 
## cat(paste("Generation:", i*2,"\n"))
## 
##  }
## 
## plot(res$generation, res$fst, col=(res$mantelp<0.05)+1, pch=16)
## 
## res
## 


## 

## ---- fig.height=4, message=FALSE----------------------------------------
#install via: library(devtools)
#install_github("dkahle/ggmap")
library(ggmap) 
library(dartR)
library(osmdata)
#also load some helper functions

source("https://raw.github.com/green-striped-gecko/dartRworkshop/master/code/lgfuns.r")



## define your genlight data set 
## here we use the bandicoots in NSW
glNSW <- bandicoot.gl[pop(bandicoot.gl)=="NSW",1:100]

## bounding box in lat/lon
ll <-glNSW@other$latlong
bb <-  matrix(c(min(ll$lon),min(ll$lat),max(ll$lon), max(ll$lat) ),nrow=2,ncol=2) 

#or you can use geocode functions
#bb <- getbb("South Australia")

#get a map
# you may need to play with zoom to get a suitable map
# also you could try terrain-lines to get roads if your study is close 
# to a city
map <- get_stamenmap(bb,zoom=6, maptype = "terrain-background")
plot(map)



## ---- fig.height=4-------------------------------------------------------
## converts the image to a raster (only the red values are returned)
map.raster <- ggmap_to_raster(map)[[1]]
map.raster
plot(map.raster)


## ---- fig.height=4-------------------------------------------------------
reslay <- map.raster #copy to keep the original map
reslay[values(reslay)==153] <- 1000
plot(reslay)
points(glNSW@other$latlong, pch=16)


## ------------------------------------------------------------------------
glNSW2 <- glNSW[glNSW$other$latlong$lon<152 & glNSW$other$latlong$lon > 143  & glNSW$other$latlong$lat< -30 & glNSW$other$latlong$lat > -36,]


## ---- fig.height=4-------------------------------------------------------
glc <- gl.genleastcost(glNSW2, fric.raster = reslay, gen.distance = "propShared", NN = 8, pathtype = "leastcost")


## ------------------------------------------------------------------------
wassermann(gen.mat =glc$gen.mat , eucl.mat = glc$eucl.mat,  cost.mats = glc$cost.mats,plot = F)


## 

## 

## ---- fig.height=4-------------------------------------------------------
r <- raster(xmn=142, xmx=153, ymn=-37, ymx=-29)
r[] <- 1
r
rbarrier <- r
rbarrier[,180:190]<- 10
plot(rbarrier, axes=FALSE)
points(glNSW@other$latlong, pch=16)



## ---- message=FALSE, warning=FALSE---------------------------------------
library(rgdal)
canberra <- c(149.13, -35.2809)
hobart <- c(147.3272, -42.8821)

ll <- rbind(canberra, hobart)
ll.sp <- SpatialPoints(ll, CRS("+proj=longlat +ellps=WGS84"))

xy <-spTransform(ll.sp, CRS( "+proj=utm +zone=55 +south +ellps=WGS84 +datum=WGS84") )
coordinates(xy)



## ------------------------------------------------------------------------
utms <- SpatialPoints(xy,  CRS("+proj=utm +zone=55 +south +ellps=WGS84 +datum=WGS84"))
ll.sp2 <- spTransform(utms,CRS( "+proj=longlat +ellps=WGS84"))
coordinates(ll.sp2)
coordinates(ll.sp)


## ---- fig.height=4-------------------------------------------------------
#reload the raster layer (to be used as the template)
ts <- raster("./data/TreeScrub.tif")

#elevation layer from the source mentioned above
res <- raster("./data/srtm_72_20.tif")

#use poss.nz and its projection NZGD2000 to reproject the layer
#resolution should be as ts, but too big here (takes to long)
res.proj <- projectRaster(res, crs = proj4string(ts), res=500)

# define the extent of the new resistance layer by samples
xmin <- min(poss.nz@other$xy$x-1000)
xmax <- max(poss.nz@other$xy$x+1000)
ymin <- min(poss.nz@other$xy$y-1000)
ymax <- max(poss.nz@other$xy$y+1000)

#define extent by other resistance layers
#xmin <- extent(ts)@xmin
#xmax <- extent(ts)@xmax
#ymin <- extent(ts)@ymin
#ymax <- extent(ts)@ymax

res.temp <- raster( xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=proj4string(ts))


res.proj.crop <- crop(res.proj, res.temp)
plot(res.proj.crop)
points(poss.nz@other$xy, pch=16)

res.proj.crop


## ---- fig.height=4-------------------------------------------------------


data(possums.gl)
library(raster)  #needed for that example
landscape.sim <- readRDS(system.file("extdata","landscape.sim.rdata", package="dartR"))

plot(landscape.sim)
points(possums.gl@other$xy, pch=16)


### new gl.gdm function (not tested yet)
#install.packages("gdm")
library(gdm)

gl.gdm <- function(x, glc, geo=TRUE, splines=NULL, knots=NULL)
{
  gdis <- data.frame(site=colnames(glc$gen.mat), glc$gen.mat)
  
  #recover coordinates
  fac <- pop(x)
  xys <- apply(x@other$xy, 2, function(q) tapply(q, fac, mean))
  predtab <- data.frame(site=rownames(xys), x=xys[,1], y=xys[,2] )
  costmats <- lapply(glc$cost.mats, function (x) data.frame(site=rownames(xys), x )
  )
  gtab <- formatsitepair(bioData = gdis, bioFormat = 3,XColumn = "x", YColumn = "y", pred=predtab, siteColumn = "site", distPreds = costmats)
  gdmmodell <-gdm(gtab, geo=geo, splines=splines, knots=knots) 
  return(gdmmodell)
}
  
glc <- gl.genleastcost(x=possums.gl,fric.raster=landscape.sim , gen.distance = "D", NN=8, pathtype = "leastcost",plotpath = TRUE)

gdm1 <- gl.gdm(possums.gl, glc)
summary(gdm1)
plot(gdm1)



## 

## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(tidy = TRUE)
library(tufte)
library(knitr)
#knitr::opts_chunk$set(prompt = TRUE)


## ---- warning=FALSE, message=FALSE---------------------------------------
#run withour error
library(adegenet)
library(Rcpp)
library(raster)
library(rgdal)
library(PopGenReport)
library(dartR)  
library(Sunder)
library(ecodist)
library(gdm)
library(GGally)
library(lme4)




## ---- echo=TRUE----------------------------------------------------------
source("./WEEG/Prac1_Mon/code/helper functions landscape genetics.R")


## ---- eval=TRUE----------------------------------------------------------
dir("./WEEG/Prac1_Mon/data/")


## 
## **Is there an effect of roads, eucalypt density and elevation on the population structure of koalas in the Katoomba area?**

## 
## The data set consists of a sample of 20 animals of koalas sampled around Katoomba. The samples have been genotyped and produced 30000 genetic markers (SNPs). For each sample the coordinates were recorded (using Map Grid Australia 94 (=UTM) as the coordinate system).

## 
## Using your GIS skills you were able to source three different maps of the Katoomba area. Each map is a raster data set and covers 1000 x 1000 pixel, whereas each pixel has a dimension of 16x16m. Luckily the coordinate system between your samples and the maps match [otherwise you would need to know how to reproject your data sets. At the end of this tutorial you find some short scripts that explain how to do that].

## 
## The first map is an digital elevation map of the area in meters [from 110-1160 m]. The second map represents the road network and has different values for highway, normal roads and path/tracks. The third and final map shows the eucalypt density of of koala food trees of the area for each pixel [from 0 to 1].

## 
## Your task in this tutorial is to find out which 'landscape' has an effect on the population structure of the koala population.

## 
## 

## ------------------------------------------------------------------------
snps.data <- read.csv("./WEEG/Prac1_Mon/data/snptable.csv")
kable(snps.data)




## ------------------------------------------------------------------------

snps.gl <- new("genlight", gen=snps.data[,-1], ind.names=snps.data[,1], loc.names=colnames(snps.data)[-1], ploidy=rep(2, nrow(snps.data)))



## ---- fig.height=4, warning=FALSE----------------------------------------
snps.gl  #gives an overview of the genlight object
plot(snps.gl)
gl.report.hwe(snps.gl)


## 
## **Task 1**

## 
## Here comes the first task. In the 'data' folder you find the data set 'koalas_snps.csv'. Load this data set and convert it into a genlight object called ```koalas```. The genlight object should have 20 individuals and 30000 loci (SNPs).

## 
## In case you get stuck, you can use the hint() function (which is available for every task during the tutorial). Simply type:

## 
## ```hint(1) ```

## 
## for this task and you should get some help. In case you want to have the full solution, type:

## 
## ```solution(1)```

## 
## 
## Below you can see the output if you solve it, if you type:

## 
## koalas

## 
## into the console.

## 
## 

## ---- echo=FALSE---------------------------------------------------------
snps.data <- read.csv("./WEEG/Prac1_Mon/data/koalas_snps.csv")

koalas <- new("genlight", gen=snps.data[,-1], ind.names=snps.data[,1], loc.names=colnames(snps.data)[-1], ploidy=rep(2, nrow(snps.data)))

koalas #check the genlight object



## ------------------------------------------------------------------------
#Genetic distance matrix
Gdis <- 1-gl.propShared(koalas)


## ---- fig.height=3.5-----------------------------------------------------
image (Gdis)
table.value(Gdis, csize=0.4)


## **Extra**

## Try to understand the structure of the genlight object (koalas). For example nInd(koalas) returns the number of individuals.

## 

## ------------------------------------------------------------------------
latlongs <- read.csv("./WEEG/Prac1_Mon/data/koalas_locs.csv")
head(latlongs)






## ------------------------------------------------------------------------
#check if labels match
sum(indNames(koalas)==latlongs$id)



## ------------------------------------------------------------------------
koalas@other$latlongs <- latlongs[,2:3]
head(koalas@other$latlongs)


## ----eval=TRUE-----------------------------------------------------------
pop(koalas)<- 1:20
gl.map.interactive(koalas)


## **Task 2**

## Find via google the proj4 string for the projection:

## 
## **MGA96 Zone 56** [and compare it to UTM 56 South]

## 
## at [spatialreference.org](spatialreference.org)

## 
## 
## Then use the ```project``` function to reproject your latlongs into a new object called ```xy```. Be aware that you need to convert the latlongs into a matrix using the ```as.matrix``` function.

## 
## Again you can type ```hint(2)``` or ```solution(2)``` if you get stuck.

## 
## Compare your coordinates with the output below.

## 

## ---- echo=FALSE---------------------------------------------------------
xy <- project(as.matrix(koalas@other$latlongs), proj = '+proj=utm +zone=56 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')


## ------------------------------------------------------------------------
colnames(xy)<- c("x","y")
head(xy) 


## ------------------------------------------------------------------------
koalas@other$xy <- data.frame(xy)
koalas


## ------------------------------------------------------------------------
par(mfrow=c(1,2), mai=c(0.5,0.5,0,0))
plot(koalas@other$xy, pch=16, asp=1)
text(koalas@other$xy+500, labels=1:20)

plot(koalas@other$latlongs, pch=16, asp=1)



## ------------------------------------------------------------------------
Edis <- as.matrix(dist(koalas@other$xy))
dim(Edis)


## ------------------------------------------------------------------------
table.value(Edis, col.labels = 1:20, csize = 0.3)


## 
## **Task 3**

## 
## a) Find the largest pairwise distance within Edis and check with the plots on the coordinates above.

## 
## b) Find the smallest pairwise distance (ignoring the diagonal)

## 
## c) Which individual is on average the "most" isolated individual

## 
## 

## ---- fig.height=3-------------------------------------------------------
Edis.vec <- lower(Edis)
Gdis.vec <- lower(Gdis)
length(Edis.vec)  #Why 190?

plot(Gdis.vec ~ Edis.vec)




## **Task 4**

## 
## Run a simple linear regression of Gdis (response) against Edis (predictor).

## Check the regression coefficient r and the $R^2$-value.

## 
## 

## ------------------------------------------------------------------------
ecodist::mantel(Gdis.vec ~ Edis.vec)


## ---- fig.height=3-------------------------------------------------------
roads <- raster("./WEEG/Prac1_Mon/data/roads.tif")
roads
plot(roads)
points(koalas@other$xy, pch=16, col="orange")


## ------------------------------------------------------------------------
extent(roads)
range(koalas@other$xy$x)
range(koalas@other$xy$y)




## **Extra task**

## 
## You can try to find a way to "formally" test if all locations are within the extent (there are GIS function in the raster package to do so, e.g. extract, intersect).


## ------------------------------------------------------------------------
table(values(roads)) #returns a count on all values in 


## ------------------------------------------------------------------------
crs(roads)


## **Task 5**

## 

## Load the second may "eucs.tif" and store it under the name "eucs". Plot 'eucs' and check the extent, projection and the "coding". The map is showing the cover of eucalpytus trees within a pixel, so a value of 1 represents a cell that has 100% eucalypt cover and a value of 0 means there are no euclyptus tree in that cell. The idea is that koalas like eucalypts for dispersal and that areas with a high density of eucalypts might promote disperal, whereas no tree limit the dispersal ability of individuals.

## 
## 
## 

## ---- echo=FALSE, fig.height=4-------------------------------------------
eucs <- raster("./WEEG/Prac1_Mon/data/eucs.tif")
plot(eucs)


## ------------------------------------------------------------------------
ele <- raster("./WEEG/Prac1_Mon/data/elevation.asc")
ele
plot(ele)
points(koalas@other$xy)
range(values(ele)) #no missing data!!
crs(ele)



## ------------------------------------------------------------------------
crs(ele) <- crs(roads)
ele


## ---- fig.height=2-------------------------------------------------------

rec.roads <- roads #first copy the original layer
rec.roads[values(rec.roads)==3] <- 150 # a strong effect for motorways
rec.roads[values(rec.roads)==2] <-80  #intermediate for roads
rec.roads[values(rec.roads)==1] <- 1 # a very small effect of tracks
rec.roads[values(rec.roads)==0] <- 0 # no road has no resistance


#check if recoding worked as expected
par(mai=c(0,0,0,0), mfrow=c(1,2))
plot(roads)
plot(rec.roads)

table(values(roads))
table(values(rec.roads))



## ---- eval=FALSE---------------------------------------------------------
## CD <- list() #creates and empty list
## 
## #be aware that takes about 2 minutes
## system.time(CD$roads <- gl.costdistances(landscape = rec.roads+1, locs = koalas@other$xy, method = "commute", NN = 8))


## ---- eval=TRUE, echo=FALSE----------------------------------------------
CD <- readRDS("./WEEG/Prac1_Mon/data/CD.rdata")


## **Task 6**

## a) Plot Gdis (genetic distances) against CD$roads.

## b) Run a mantel test between the roads cost distance and genetic distances. Do you think there is an effect of roads?

## 

## ------------------------------------------------------------------------
rec.eucs <- (1-eucs)* 50 
hist(values(rec.eucs))



## ---- eval=FALSE---------------------------------------------------------
## #another two minutes....
## CD$eucs <- gl.costdistances(landscape = rec.eucs+1,locs = koalas@other$xy, method = "commute", NN=8 )


## ------------------------------------------------------------------------
rec.ele <- ele

rec.ele <- (rec.ele-min(values(rec.ele)))/diff(range(values(rec.ele)))*50
hist(values(rec.ele))



## ---- eval=FALSE---------------------------------------------------------
## 
## CD$ele <- gl.costdistances(rec.ele+1, locs = koalas@other$xy, method = "commute", NN=8)


## 
## 
## !!! Be aware each cost distance calculation takes about 2-3 minutes. Therefore we prepared a list with all cost distances we came up with and save it under 'CD.rdata' as a file. So if you do not want to wait for the commands below to finish, you can simply type:

## 
## CD <- readRDS("./WEEG/Prac1_Mon/data/CD.rdata")

## 
## 

## ---- eval=FALSE---------------------------------------------------------
## 
## CD <- list()
## CD$roads <- gl.costdistances(rec.roads+1, locs = koalas@other$xy, method = "commute", NN=8)
## 
## CD$eucs <- gl.costdistances(rec.eucs+1, locs = koalas@other$xy, method = "commute", NN=8)
## CD$ele <- gl.costdistances(rec.ele+1, locs = koalas@other$xy, method = "commute", NN=8)
## 
## CD$eucs.roads <- gl.costdistances(rec.eucs+rec.roads+1, locs = koalas@other$xy, method = "commute", NN=8)
## 
## CD$eucs.ele <- gl.costdistances(rec.eucs+rec.ele+1, locs = koalas@other$xy, method = "commute", NN=8)
## 
## CD$ele.roads <- gl.costdistances(rec.ele+rec.roads+1, locs = koalas@other$xy, method = "commute", NN=8)
## 
## 
## 
## 


## ------------------------------------------------------------------------
names(CD)

## ---- eval=FALSE---------------------------------------------------------
## 
## ggpairs(data.frame(lapply(CD[c("roads","eucs","ele")], lower)))
## 
## 


## **Task 7**

## 
## Create some resistance matrices

## 
## a) Your hypothesis is that eucs and roads both have an effect. Therefore you can simply add the two resistance values together: rec.eucs+rec.roads.

## 
## b) Your next hypothesis is that all three resistance layers have an effect. So you add them together, but the maximum cost value should be 100 and minumum cost 0.

## 
## c) Create a random resistance (uniform random values between 0 and 50)

## 
## d) Create your own (univariate or multivariate) resistance layer from the three available landscape layers.

## 
## Once you finished you can calculate the cost distance matrices from the resistance layers using the gl.costdistances() function and store them in the CD object, but read the hint below first, before you spent all afternoon on cost distances.

## 

## ------------------------------------------------------------------------
all.mantels(Gdis, Edis, CD)


## 
## Which predictors are good ones?

## 
## Which would you keep?

## 
## Do we have IBD?

## 
## Discuss with your neighbour!!!!

## 
## Be aware your answers might be different as you may have created different cost matrices.

## 

## ------------------------------------------------------------------------
Alldis <- CD #copy all costdistances
Alldis$"_GDis" <- Gdis  #add our genetic distance matrix
Alldis$Edis <- Edis   #add Euclidean distance matrix

#sort to have _GDis first [underscore is used to make _GDis first]
Alldis <- Alldis[order(names(Alldis))]

names(Alldis)





## ---- eval=FALSE---------------------------------------------------------
## ggpairs(data.frame(lapply(Alldis, lower)))


## ------------------------------------------------------------------------
wassermann(Gdis, CD["eucs"], Edis, plot = FALSE)
wassermann(Gdis, CD["roads"], Edis, plot = FALSE)
wassermann(Gdis, CD[c("ele")], Edis, plot = FALSE)


## 
## Again which of the base resistance matrices turn out to be important to describe the population structure of koalas?

## 

## ------------------------------------------------------------------------
wassermann(Gdis, CD[c("eucs","roads","ele")], Edis, plot = FALSE)


## ------------------------------------------------------------------------
lgrMMRR(Gdis, CD[c("eucs","roads","ele")], Edis )



## ------------------------------------------------------------------------
lgrMMRR(Gdis, CD[c("eucs.roads","ele")], Edis )


## ------------------------------------------------------------------------
lgrMMRR(Gdis, CD[], Edis )


## ------------------------------------------------------------------------
names(Alldis)
## let us run the base resistances (this time by numbers)
#_Gdis, Edis, ele, eucs, roads
ca.out<-CAmrdm(Alldis[c(1,2,3,5, 8)], regrnperm = 999, bootn = 100, bootprop = 0.25)
ca.out



## ------------------------------------------------------------------------
ids <- To.From.ID(nInd(koalas))
head(ids)


## ------------------------------------------------------------------------
df <- data.frame(lapply(Alldis, lower))
names(df)

df$pop <- ids$pop1




## ------------------------------------------------------------------------
mlpe.full <- mlpe_rga(formula = X_GDis ~ Edis+  roads + eucs + ele +(1|pop), data=df)
summary(mlpe.full)
AIC(mlpe.full)


## ------------------------------------------------------------------------
mlpe.roads <- mlpe_rga(formula = X_GDis ~ Edis+  roads+(1|pop), data=df)
summary(mlpe.roads)
mlpe.ele <-  mlpe_rga(formula = X_GDis ~ Edis+  ele+(1|pop), data=df)
summary(mlpe.ele)


## ------------------------------------------------------------------------
AIC(mlpe.roads,  mlpe.ele, mlpe.full)


## ---- eval=FALSE---------------------------------------------------------
## 
## ## Warning takes "forever"
## #Bayesian approach of Botta et al. 2015, which is implemented in package "Sunder":
## 
## 
## #Running roads
## allel.counts <- gl2sunderarray(koalas)
## # setting parameters
## nit <- 10^2 ## just for the example, should be much larger, e.g. 50000
## run  <- c(1,1,1)
## thinning  <- 1 #  just for the example, should be larger, e.g. max(nit/10^3,1)
## ud   <- c(0,1,1,0,0)
## theta.init <- c(1,2,1,1,0.01)
## n.validation.set  <- dim(allel.counts)[1]*dim(allel.counts)[2]/10
## theta.max  <- c(10,10*max(Edis),10*max(CD$roads),1,0.01)
## plot  <- FALSE
## trace <- FALSE
## 
## # now running the method
## # this is where we use the allele counts calculated above (koalas@other$a.counts)
## # the straight line distance ('geo')
## # and the effective distance created from he resistance surface based on roads ('koalas@other$eff.roads') - make sure to adjust this input for your own effective distances!
## 
## sunder.out<-MCMCCV(allel.counts, D_G = Edis, D_E = CD$roads,
##                    nit, thinning, theta.max, theta.init,
##                    run, ud, n.validation.set, print.pct=TRUE)
## 
## # these calculations will take a bit of time...
## 
## sunder.out$mod.lik
## 


## 
## That is the end. Well done you finished.

## Feel free to have another go or have a well deserved beverage!!!

## 

## ---- fig.height=4-------------------------------------------------------
gdm <- gl.gdm(xy=koalas@other$xy, gdis = Gdis, cdis = CD$roads, geo = T)



## ----eval=FALSE----------------------------------------------------------
## #only the first 5 locations (for testing). The whole map takes about an hour to run.
## CS <- gl.runCS(landscape = rec.roads+rec.eucs+1, outpath = tempdir(), locs = koalas@other$xy[1:5,] )
## 


## ---- eval=FALSE---------------------------------------------------------
## plot(log(CS$map), col=heat.colors(255))


## ---- eval=FALSE, fig.height=4-------------------------------------------
## #takes some minutes to run!!!
## 
## glc <- gl.genleastcost(koalas[1:5,], rec.roads+1, "propShared", NN=8, pathtype = "leastcost")
## lgrMMRR(gen.mat = glc$gen.mat, cost.mats = glc$cost.mats,  eucl.mat = glc$eucl.mat)
## 


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


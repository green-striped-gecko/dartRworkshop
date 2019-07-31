###########################################################
# These are some "helper functions"
# for the landscape genetics tutorial
# by Bernd Gruber and Niko Balkenhol
# at the Workshop on Ecological and Evolutionary Genomics 
# held from 4-9th August 2019
# in Katoomba, Blue Mountains, NSW, Australia
###########################################################

################################################
# FUNCTION 1: proportion of shared alleles
################################################
# nice helper function provided by Bernd Gruber to calculate the "proportion of shared alleles" (PSA) fast!!!!

# start function++++++++++++++++++++++++++++++++++++++++++++
gl.propShared <- function(x) {
  xx <- as.matrix(x)
  cppFunction('NumericMatrix glpropSharedC(NumericMatrix x) {
              int nrow = x.nrow();
              NumericMatrix out(nrow,nrow);
              for (int i=0; i<(nrow-1); i++) {
              for (int j=(i+1); j<nrow; j++) {
              out(j,i) = 1-mean( abs( na_omit(x(i,_)- x(j,_) ))/2);
              }
              }
              return out;
}')
  res <- glpropSharedC(xx)
  res <- as.matrix(as.dist(res))
  diag(res)<- 1
  res
}
# END function++++++++++++++++++++++++++++++++++++++++++++


################################################
# FUNCTION 2: allele counts
################################################
# Another nice helper function provided by Bernd Gruber
# for getting allele counts from a genlight object.
# The result can be used in the package "Sunder"

# start function++++++++++++++++++++++++++++++++++++++++++++
gl2sunderarray <- function(x) #for individuals
{
  
  malloc <- 2  #diploid only
  A <- array(0, dim=c(nInd(x), nLoc(x), malloc))
  xx <- as.matrix(x)
  A[,,1] <- abs(2-xx)
  A[,,2] <- xx
  return (A)
}
# END function++++++++++++++++++++++++++++++++++++++++++++++++



################################################
# FUNCTION 3: Mantel and partial Mantel tests for preliminary # analyses
################################################

all.mantels<-function(gen.dist, geo.dist, eff.dist.l, n.perm=299) {

	# gen.dist = pairwise matrix of genetic distances among sampling locations
	# geo.dist = pairwise matrix straight line distances
	# eff.dist.l = list holding one or several effective 			pairwise distances
	# n.perm = number of permutations, default is 299
	
	nperm<-n.perm
	gen.dist.l<-lower(as.data.frame((gen.dist)))
	geo.dist.l<-lower(as.data.frame((geo.dist)))
	n.dist<-length(eff.dist.l)
	mantel.out <- matrix(NA, ncol=5, nrow = n.dist+1)
	colnames(mantel.out)<-c('Variable', 'Mantel_r', 'Mantel_p', 'partMantel_r', 'partMantel_p')

	for (i in 1:n.dist) 	{

		eff.l<-lower(as.data.frame(eff.dist.l[i]))
		m.test.geo<-mantel(gen.dist.l~geo.dist.l, nperm=nperm)
		m.test<-mantel(gen.dist.l~eff.l, nperm=nperm)
		pm.test<-mantel(gen.dist.l~eff.l + geo.dist.l, nperm=nperm)
		
		var.name<-names(eff.dist.l[i])
		m.r<-m.test['mantelr']
		m.p<-m.test['pval1']
		pm.r<-pm.test['mantelr']		
		pm.p<-pm.test['pval1']
		
		mantel.out[1,1]<-'Geo'
		mantel.out[1,2]<-round(as.numeric(m.test.geo['mantelr']), digits=3)
		mantel.out[1,3]<-round(as.numeric(m.test.geo['pval1']), digits=3)
		
		j<-i+1
		mantel.out[j,1]<-var.name
		mantel.out[j,2]<-round(as.numeric(m.r), digits=3)
		mantel.out[j,3]<-round(as.numeric(m.p), digits=3)
		mantel.out[j,4]<-round(as.numeric(pm.r), digits=3)
		mantel.out[j,5]<-round(as.numeric(pm.p), digits=3)

				} # end for loop
		
	print(noquote(mantel.out))

} # END function 'all.mantels'




################################################
# FUNCTION 4: commonality analysis
################################################

# This is code developed by Jerome Prunier
# Make sure to check out his website and papers at
# http://jeromeprunier.fr


#######################################################
CAmrdm = function(mydata, regrnperm, bootn, bootprop){
  #######################################################
  
  #'Prunier JG, Colyn M, Legendre X, Nimon KF, Flamand MC (2015) 
  #'Multicollinearity in spatial genetics: Separating the wheat from the chaff using commonality analyses. 
  #'Molecular Ecology 24, 263-283.
  #'
  #'OUTPUT:
  #'A matrix of pairwise zero-order correlations among predictors and VIF for each predictor
  #'A model fit and its significance level
  #'A matrix providing typical MRDM results and additional parameters derived from regression commonality analysis
  #'A matrix of commonalities
  #'A plot of beta weights with 95% confidence intervals
  #'A plot of commonalities with 95% confidence intervals
  #'
  #'INPUT:
  #'mydata : a list object containing all original squared matrices of pairwise genetic (dependent variable) and environnemental (predictors) differences. 
  #'    The dependent matrix (Y) is to be the first attribute of the list.
  #'    You can create lists directly using the list() function. For instance:
  #'        Y=read.table('Y.txt',h=F)
  #'        X1=read.table('X1.txt',h=F)
  #'        X2=read.table('X2.txt',h=F)
  #'        mydata=list(GD=y,F1=X1,F2=X2)
  #'        names(mydata)
  #'            [1] "GD" "F1" "F2"
  #'optionsel: an option for further analyses
  #'    optionsel == 1 triggers Multiple Regression on Distance Matrices (MRDM)
  #'    optionsel == 2 triggers Logistic Regression on Distance Matrices (LRDM)
  #'    optionsel == 3 offers an interactive choice between MRDM and LRDM
  #'regrnperm: number of permutations for tests of significance
  #'bootn: number of bootstrap iterations for assessing parameters' robustness to the removal of a random subset of objects
  #'bootprop: proportion of sampled objects (populations or individuals) to be randomly removed at each bootstrap iteration'
  
  
  #loading libraries
  require(fmsb)
  require(yhat)
  require(ecodist)
  
  #extraction of useful parameters
  npop=dim(mydata[[1]])[1]    
  nbvar=length(mydata)-1
  listvar=c(1:nbvar)
  
  #vectorisation and Z-transformation of matrices
  prep = function(x){ x=data.matrix(x)
  x=x[upper.tri(x, diag = FALSE)]
  x=(x-mean(x))/sqrt(var(x)) }
  mydata2=lapply(mydata, prep)     
  names(mydata2) = paste0("X", seq_along(mydata2)-1)
  names(mydata2)[1]="Y"
  mydata2=as.data.frame(mydata2)
  
  #Computation of zero-order correlations among predictors and VIF
  rp=cor(mydata2,method="pearson")
  VIF=data.frame(VIF=rep(0,nbvar))
  for (i in c(1:nbvar)){  xnam = paste("X", listvar[-(i)], sep="")
  fmla = as.formula(paste(paste("X", listvar[i], sep="")," ~ ", paste(xnam, collapse= "+")))
  VIF$VIF[(i)]=VIF(lm(fmla,data=mydata2))[1] }
  cormat=cbind(rp[c(2: dim(rp)[1]),],VIF) 
  
  #MRDM and CA
  xnamtot = paste("X", listvar, sep="")
  fmlatot = as.formula(paste(paste("Y ~ ", paste(xnamtot, collapse= "+"))))
  ca=regr(lm(fmlatot,data=mydata2))
  mrm=MRM(fmlatot,data=mydata2,nperm=regrnperm)
  mylm=lm(fmlatot,data=mydata2)
  rs=calc.yhat(mylm)$PredictorMetrics[c(1:nbvar),4]
  residuals=linear2matrix(as.data.frame(mylm$residuals),npop)
  #normalitytest=shapiro.test(mylm$residuals)
  regroutput=cbind(rp[c(2:(dim(rp)[2])),1],mrm$coef[c(2:(nbvar+1)),],ca$Commonality_Data$CCTotalbyVar)
  suppression=cbind((abs(sum((ca$Commonality_Data$CC[,1])[ca$Commonality_Data$CC[,1]<0]))),(abs(sum((ca$Commonality_Data$CC[,1])[ca$Commonality_Data$CC[,1]<0])/ca$Commonality_Data$CC[dim(ca$Commonality_Data$CC)[1],1])))
  suppression=suppression*100
  colnames(suppression)=c('%Y','%modelfit')
  modelfit=mrm$r.squared
  
  #Computation of parameters' robustness to random removal of populations through bootstrapping
  nbcomm=2^(nbvar)-1
  beta=matrix(data = 0, nrow = bootn, ncol = nbvar)
  boot=matrix(data = 0, nrow = bootn, ncol = nbcomm)
  subselnpop=npop-round(bootprop*npop)
  resBeta=data.frame(n=rep(0,nbvar),o=rep(0,nbvar),l=rep(0,nbvar),u=rep(0,nbvar))
  resBootstrap=data.frame(n=rep(0,nbcomm),o=rep(0,nbcomm),l=rep(0,nbcomm),u=rep(0,nbcomm))
  mydata3=mydata2[c(1:((subselnpop*(subselnpop-1))/2)),]
  for (i in 1:bootn){ print(i)
    rarray=sort(sample(npop,subselnpop,replace=F))
    for (j in 1:(nbvar+1)){ mydata3[,j]=prep(mydata[[j]][rarray,rarray])}
    comm=regr(lm(fmlatot,data=mydata3))
    beta[i,]=comm$LM_Output$coefficients[c(2:(nbvar+1)),1]
    boot[i,]=comm$Commonality_Data$CC[c(1:nbcomm),1] }
  for (i in 1:nbcomm){  q=quantile(boot[,i], c(.025,.975))
  resBootstrap[i,1]=i
  resBootstrap[i,3]=q[1]
  resBootstrap[i,4]=q[2] }
  for (i in 1:nbvar){ q=quantile(beta[,i], c(.025,.975))
  resBeta[i,1]=i
  resBeta[i,3]=q[1]
  resBeta[i,4]=q[2]}    
  resBeta[,2]=ca$LM_Output$coefficients[c(2:(nbvar+1)),1]
  resBootstrap[,2]=ca$Commonality_Data$CC[c(1:nbcomm),1]
  
  namecomm=names(ca$Commonality_Data$CC[c(1:nbcomm),1])
  for (i in 1:nbcomm){namecomm[i]=substr(names(ca$Commonality_Data$CC[c(1:nbcomm),1])[i],11,max(nchar(names(ca$Commonality_Data$CC[c(1:nbcomm),1]))))}
  for (i in 1:nbcomm){temp=namecomm[i]
  convar= unique(na.omit(as.numeric(unlist(strsplit(unlist(temp), "[^0-9]+")))))
  namecomm[i]=paste(names(mydata)[convar+1], collapse= ", ")}                 
  
  #plots
  
  close.screen(all.screens = T)
  split.screen(figs=c(1,2))
  screen(n=1)
  plot(resBootstrap[,2],-resBootstrap[,1],arrows(resBootstrap[,3],-resBootstrap[,1],resBootstrap[,4],-resBootstrap[,1],angle=90,code=3,length=0.05), xlim=c(min(resBootstrap[,3])-0.02, max(resBootstrap[,4])+0.02), ylab="Sets of predictors",xlab="% of explained variance",font=5,lab=c(nbcomm, 7, 1),xaxt="n",yaxt="n",cex.lab=0.7) 
  abline(v=0,lty=2)
  points(resBootstrap[,2],-resBootstrap[,1],pch=21,bg="white") 
  axis(4,cex.axis=0.6,at=-c(1:nbcomm),labels=namecomm,las=2)
  axis(1,cex.axis=0.7)
  title(main=" Plot of commonalities")
  
  close.screen(all.screens = T)
  plot(resBeta[,1],resBeta[,2],arrows(resBeta[,1],resBeta[,3],resBeta[,1],resBeta[,4],angle=90,code=3,length=0.05), ylim=c(min(resBeta[,3])-0.02, max(resBeta[,4])+0.02),xlim=c(0.5,nbvar+0.5), xlab="Predictors", ylab=names(mydata)[1],xaxt="n",yaxt="n",cex.lab=0.7) 
  abline(a = 0, b = 0,lty=2)
  points(resBeta[,1],resBeta[,2],pch=21,bg="white") 
  axis(1,cex.axis=0.7,at=c(1:nbvar),labels=names(mydata[c(2:(nbvar+1))]))
  axis(2,cex.axis=0.7)
  title(main=" Plot of beta weights")
  
  
  close.screen(all.screens = T)
  minTOT=min(mydata2[,-1])
  maxTOT=max(mydata2[,-1])
  pok=seq(minTOT,maxTOT,0.1)
  p=as.data.frame(matrix(data=0,ncol=nbvar,nrow=length(pok)))
  for (i in c(1:nbvar)){
    p[,i]=mean(mydata2[,i+1])
  }
  plot.new=T  
  par(mfrow=c(1,1))
  
  prGLMtemp=as.data.frame(matrix(data=0,ncol=nbvar,nrow=length(pok)))
  for (i in c(1:nbvar)){
    Newdata = p
    Newdata[,i]=pok
    names(Newdata)=names(mydata2)[-1]
    prGLMtemp[,i] = predict(mylm,newdata=Newdata,se=TRUE)$fit}
  
  for (i in c(1:nbvar)){
    Newdata = p
    Newdata[,i]=pok
    names(Newdata)=names(mydata2)[-1]
    prGLM = predict(mylm,newdata=Newdata,se=TRUE)
    if (i==1){plot(pok,(prGLM$fit),type="l",col=i,lty=1,lwd=1,bg = "transparent",xlab="Z-scores",ylab="Dependent variable",xlim=c(minTOT, maxTOT), ylim=c(min(prGLMtemp), max(prGLMtemp)))
    } else {lines(pok,prGLM$fit,type="l",col=i,lty=1,lwd=1,bg = "transparent")}   
    par(new=T)}
  legend("bottomright", legend=names(mydata)[-1], text.col=c(1:nbvar),cex=0.5)
  
  #outputs
  names(cormat)[c(1:(nbvar+1))]=names(mydata)
  rownames(cormat)=names(mydata[c(2:(nbvar+1))])
  regroutput2=as.data.frame(c(as.data.frame(regroutput[,1]),as.data.frame(rs),as.data.frame(regroutput[,2]),resBeta[,c(3:4)],as.data.frame(regroutput[,c(3:6)])))
  colnames(regroutput2)[c(1:5)]=c("r","rs","betas","CIinf","CIsup")
  rownames(regroutput2)=names(mydata[c(2:(nbvar+1))])
  for (j in 1:nbcomm){rownames(ca$Commonality_Data$CC)[j]=namecomm[j]}
  output=list(Matrixcorrelation=cormat,modelfit=modelfit,suppression=suppression,regroutput=regroutput2,commonalities=ca$Commonality_Data$C) #ResidualsNormality=NA)
  print(output)
}


################################################################# 
linear2matrix = function(data,dimm){
  ################################################################# 
  mat=list()
  for (k in c(1:dim(data)[2])){ 
    extrait=data[,k]
    mat2=matrix(data=0,ncol=dimm,nrow=dimm)
    for (i in c(1:(dimm-1))){
      if (length(extrait)>2) {
        temp=extrait[-c(1:(length(extrait)-(dimm-i)))]
        mat2[c(1:(dimm-i)),(dimm-i+1)]=temp
        mat2[(dimm-i+1),c(1:(dimm-i))]=temp
        extrait=extrait[c(1:(length(extrait)-(dimm-i)))]}
      else{temp=extrait
      mat2[c(1:(dimm-i)),(dimm-i+1)]=temp
      mat2[(dimm-i+1),c(1:(dimm-i))]=temp
      }
    }
    mat[[k]]=mat2
  }
  return(mat2)
}



hint <- function(hint=NULL)
{
  if (is.null(hint)) stop("You need to provide the task number you want to get a hint from (e.g. hint(1) for a hint for task 1.")
  
  if (hint==1)
  {
    cat("Look carefully how the sample data set was loaded. \nReplace the file name by the file name of the koala snp data set")
  }
  
  if (hint==2) {
    cat("Within google search type: 'proj4 MGA94 Zone 56'\n")
    cat("Then use the project function: \nproject(as.matrix(latlongs), proj=proj4string)")
  }
  
  if (hint==3)
  {
    cat("Use the functions min, max and colMeans in combination with which.\nBe aware that you need to take care for the zero entries along the diagonal to answer b) and c).\n ")
  }
  
  if (hint==4) cat("Use the lm function in the form:\n
                   lm(y ~ x)\n")
  if(hint==5) cat("Use the raster function to load the eucs.tif file. Refer to the roads example. To check the values the range(), summary() or hist functions() are a good idea.\n")
  
  if (hint==6) {
    cat("You need to convert Gdis and CD$roads into vectors via the lower() function and plot it\n")
    cat("Then use the mantel() function on the vectors for the mantel test") 
  }
  if(hint==7){
    cat("a)  Not sure what to say here. \n")
    cat("b)  Add all of them together and then standardise as elevation above between zero and one and then multiply by 100\n")
    cat("c) Find out how many random value you need and then use runif() \n")
    cat("d) There is no wrong answer...")
  }  
  
  
  
    
} 
  
  
solution <- function(solution=NULL)
{
  if (is.null(solution)) stop("You need to provide the task number you want to get the solution for (e.g. solution(1) prints the solution for task 1.")
  
  if (solution==1)
  {
    cat('snps.data <- read.csv("./WEEG/Day2/data/koala_snps.csv")\n')
    
    cat('koalas <- new("genlight", gen=snps.data[,-1], ind.names=snps.data[,1], loc.names=colnames(snps.data)[-1], ploidy=rep(2, nrow(snps.data)))\n')
    
    cat('koalas #check the genlight object\n')
    
  }
  
  if (solution==2) cat ("xy <- project(as.matrix(koalas@other$latlongs), proj = '+proj=utm +zone=56 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')\n")
  
  if (solution==3)
  {
    cat("a)\nwhich(Edis==max(Edis), arr.ind=T)\n")
    
    cat("b)\nee <- Edis #we want to keep Edis for later\ndiag(ee)<- NA\nwhich(ee==min(ee, na.rm=T), arr.ind=T)\n")
    cat("c)\ncolMeans(ee, na.rm=T)\nwhich.max(colMeans(ee, na.rm=T))")
  }
  if (solution==4) {
    cat("summary(lm(Gdis.vec ~ Edis.vec))\n")
    cat("cor(Gdis.vec , Edis.vec)")
  }
  if (solution==5) {
    cat('eucs <- raster("./WEEG/Day2/data/eucs.tif")')
    cat()
  }
  if (solution==6) {
    cat("plot(lower(Gdis)~ lower(CD$roads))\n")
    cat("mantel(lower(Gdis)~ lower(CD$roads))\n")
  }
  if(solution==7)
  {
cat("a)\n rec.eucs.roads <- rec.eucs+rec.roads\n")
        cat("b)\n rec.all <- rec.eucs+rec.roads+rec.ele\n")
    cat("rec.all <- (rec.all-min(values(rec.all)))/diff(range(values(rec.all)))*100\n")

  cat("c) \n rec.random <- roads #copy an layer\n values(rec.random)  <- runif(1000*1000)*50\n")
  cat("d) \n rec.whatever <- rec.roads/(rec.eucs+1)*7\n")
  
}  
}



library(gdm)

gl.gdm <- function(x=NULL, xy=NULL, gdis=NULL,cdis=NULL, geo=TRUE, splines=NULL, knots=NULL)
{
  if (is.null(gdis))
    gdis <- data.frame(site=colnames(x@other$Gdis), x@other$Gdis) else gdis <- data.frame(site=colnames(gdis),gdis)
    
    #recover coordinates
    if (!is.null(xy)) xys <- xy else xys <- x@other$xy
    predtab <- data.frame(site=rownames(xys), x=xys[,1], y=xys[,2] )
    if (is.null(cdis))   costmats <- list(data.frame(site=rownames(xys), x@other$Cdis)) else costmats <- list(data.frame(site=rownames(xys), list(cdis)))
    gtab <- formatsitepair(bioData = gdis, bioFormat = 3,XColumn = "x", YColumn = "y", pred=predtab, siteColumn = "site", distPreds = costmats)
    gdmmodell <-gdm(gtab, geo=geo, splines=splines, knots=knots) 
    plot(gdmmodell)
    print(summary(gdmmodell))
    return(gdmmodell)
}


##################################
### Run Circuitscape (Windows) ###
##################################

# load raster library
library(raster)


gl.runCS <- function(landscape, locs,outpath=tempdir(), plot=FALSE, CS_exe = 'C:/"Program Files"/Circuitscape/cs_run.exe' )
{
  # Cost surface
  cost <- landscape
  # Locs
  sites <- locs
  
  # Plot it
  if (plot)
  {
    plot(cost)
    points(sites,pch=19,col=2)
  }
  # Rasterize points using the cost extent
  
  sites <- rasterize(x = sites,y = cost)
  # Write rasters to your working directory
  
  writeRaster(sites,file.path(outpath, "sites_rast.asc"),overwrite=TRUE)
  writeRaster(cost,file.path(outpath,"cost_rast.asc"),overwrite=TRUE)
  
  # Make an .ini file
  CS_ini <- c("[circuitscape options]",            
              "data_type = raster",
              "scenario = pairwise",
              "write_cur_maps = 1",
              paste(c("point_file =",
                      "habitat_file =",
                      "output_file ="),
                    paste(outpath,c("sites_rast.asc",
                                    "cost_rast.asc",
                                    "CS.out"),
                          sep="/")))
  
  # Write it to your working directory
  writeLines(CS_ini,file.path(outpath,"myini.ini"))
  
  # Make the CS run cmd
  CS_run <- paste(CS_exe, file.path(outpath,"myini.ini"))
  
  # Run the command
  system(CS_run)
  
  # Import the effective resistance
  CSdis <- as.dist(read.csv(file.path(outpath,"CS_resistances.out"),sep=" ",row.names=1,header=1))
  
  CS_currentmap <- raster(file.path(outpath, "CS_cum_curmap.asc"))
  
  return(list(CSdis=CSdis, map=CS_currentmap))
}






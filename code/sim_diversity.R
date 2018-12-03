#initialise 
simpops <- init.popgensim(para$n.pops, para$n.ind, para$sex.ratio, 
                          para$n.loci, para$n.allels, para$locs, para$n.cov )  
#Calculate overall Fsts
glsp <- pops2gl(simpops, locs =para$locs)
gen.mat <- as.matrix(as.dist(gl.fst.pop(glsp, nboots = 1)))
#overall pairwise fst
tempfst <- mean(as.dist(gen.mat))

glc <- gl.genleastcost(glsp, fric.raster = landscape, gen.distance = "D", NN = 8, pathtype = "leastcost")
temp.mantel <- wassermann(gen.mat =glc$gen.mat , eucl.mat = glc$eucl.mat,  cost.mats = glc$cost.mats,plot = F)

#find the p value in temp
mantelp <- as.numeric(temp.mantel$mantel.tab["1","p"])

gg <- list(glsp)
nn<-15
res <- data.frame(generation=0, fst=tempfst, mantelp)
for (i in 1:nn)
{
  simpops <- run.popgensim(simpops, steps=2, cost.mat,
                           n.offspring=para$n.offspring,n.ind=para$n.ind,
                           para$mig.rate, para$disp.max, para$disp.rate,
                           para$n.allels, para$mut.rate,
                           n.cov=para$n.cov, rec="none")
  
  glsp <- pops2gl(simpops)
  gensim.mat <- as.matrix(as.dist(gl.fst.pop(glsp, nboots = 1)))
  #overall pairwise fst
  tempfst <- mean(as.dist(gensim.mat))
  temp.mantel <- wassermann(gen.mat =gensim.mat , eucl.mat = glc$eucl.mat,  cost.mats = glc$cost.mats,plot = F)
  
  mantelp <- as.numeric(temp.mantel$mantel.tab["1","p"])
  
  
  res[i+1,] <- c(i*2, tempfst, mantelp)
  gg[[i+1]] <- glsp
  
  cat(paste("Generation:", i*2,"\n"))
  
}

plot(res$generation, res$fst, col=(res$mantelp<0.05)+1, pch=16)

res
par(mfrow=c(4,3))
for (i in 1:(nn+1))
  
{
  
dd <- gl.diversity(gg[[i]], probar = F)  

tmat <- as.matrix(as.dist(dd$zero_D_beta))
temp.mantel <- wassermann(gen.mat =tmat , eucl.mat = glc$eucl.mat,  cost.mats = glc$cost.mats,plot = F)
mantelp <- as.numeric(temp.mantel$mantel.tab["1","p"])
res$nall[i]<- mantelp

tmat <- as.matrix(as.dist(dd$one_D_beta))
temp.mantel <- wassermann(gen.mat =tmat , eucl.mat = glc$eucl.mat,  cost.mats = glc$cost.mats,plot = F)
mantelp <- as.numeric(temp.mantel$mantel.tab["1","p"])
res$shannon[i]<- mantelp

tmat <- as.matrix(as.dist(dd$two_D_beta))
temp.mantel <- wassermann(gen.mat =tmat , eucl.mat = glc$eucl.mat,  cost.mats = glc$cost.mats,plot = F)
mantelp <- as.numeric(temp.mantel$mantel.tab["1","p"])
res$D[i]<- mantelp
}
  
matplot(res$generation, res[,c(4, 5,6,3)], type="b", xlab="#generations", ylab="diversity", pch=as.character(0:3))
legend("top", legend=c( "zero_H", "one_H","two_H", "Fst"), col=1:4, pch=as.character(0:3))
abline(h=0.05, col="orange", lwd=2)

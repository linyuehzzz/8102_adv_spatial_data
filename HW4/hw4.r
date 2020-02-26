# 1. The GeoR package
library(geoR)
library(maps)
library(mvtnorm) # package for simulation
rm(list=ls()) # clear workspace


# 2. Simulating a Gaussian spatial process
# 2.1 Obtain the limits of our spatial domain
ia.range = map("state", "iowa", plot=FALSE)$range
rx = ia.range[2] - ia.range[1] # range in x direction
ry = ia.range[4] - ia.range[3] # range in y direction

# 2.2 Define the spatial points (locations)
xg = seq(ia.range[1], ia.range[2], length.out=30)
yg = seq(ia.range[3], ia.range[4], length.out=30)
iagrid.locs = expand.grid(xg, yg)
map("state", "iowa")
points(iagrid.locs, pch=20)
ntot = dim(iagrid.locs)[1]

# 2.3 Get a matrix of pairwise distances
distmat = as.matrix(dist(iagrid.locs))
max_dist = max(distmat)

# 2.4 Use a true exponential covariance model
theta1 = 5
theta2 = 0.5
Sig = theta1 * exp(-distmat / theta2)
plot(seq(0, max_dist, length.out=20), 
     theta1 * exp(-seq(0, max_dist, length.out=20) / theta2),
     type="b", ylab="cov", xlab="distance")

# 2.5 Simulate a spatial random process from the Gaussian model
y = rmvnorm(1, matrix(5,ntot,1), Sig)
image(matrix(y, length(xg), length(yg)), x=xg, y=yg,
        col=rev(rainbow(100, start=0, end=.7)))
map("state", "iowa", lwd=3, add=TRUE)
box() # add a boundary box to the plot


# 3. Use a subset of the simulated process as working data
idxkeep = sort(sample(1:ntot, round(0.3*ntot)))
ymask = matrix(0, ntot, 1)
ymask[idxkeep, ] = 1
image(matrix(y, length(xg), length(yg)), x=xg, y=yg,
      col = rev(rainbow(100,start=0,end=.7)))
image(matrix(ymask, length(xg), length(yg)), x=xg, y=yg,
      col=c("white","transparent"), add=TRUE)
map("state", "iowa", add=TRUE, lwd=3)
box()


# 4. Estimate the empirical semi-variogram
# 4.1 Create a geodata object using the subset of simulated data
n = length(idxkeep)
ydat = matrix(y[idxkeep], n, 1)
dat.locs = iagrid.locs[idxkeep, ]
ygeodat = as.geodata(cbind(dat.locs, ydat))

# 4.2 Estimate the empirical semi-variogram of the sample data
distsample = as.matrix(dist(dat.locs)) # distance between sample pairs
maxd = max(distsample)/2
y.v = variog(ygeodat, max.dist=maxd) # estimate empirical semi-variogram
y.v$beta.ols # estimated trend of the sample data

# 4.3 Estimate the true semi-variogram
vt = theta1 - theta1 * exp(-y.v$u/theta2)

# 4.4 Plot the empirical semi-variogram and true semi-variogram
xmax = max(y.v$u, maxd) # xmax to set xlim for plot
ymax = max(y.v$v, vt) # ymax to set ylim for plot
plot(y.v, xlim = c(0,xmax), ylim = c(0,ymax), col="red") # empricial semivariogram
lines(y.v$u, vt, type="b", pch=20, col="blue", lty=2) # add true semivariogram
legend("bottomright", legend=c("Empricial", "True"), col=c("red", "blue"), lty=1:2, cex=0.8)

# 4.5 Plot the empirical semi-variogram in each of 4 principle directions
plot(variog4(ygeodat, max.dist=maxd))


# 5. Fit an exponential semi-variogram model 
y.v.wls=variofit(y.v, cov.model="exponential", wei="cressie")
y.v.wls
plot(y.v ,xlim=c(0,xmax), ylim=c(0,ymax), col="red") # empirical semivariogram
lines(y.v.wls, col="black", lty=2) # fitted semivariogram
lines(y.v$u,vt, type="b", pch=20, col="blue", lty=3) # true semivariogram
legend("bottomright", legend=c("Empricial", "Fitted", "True"), 
       col=c("red", "black", "blue"), lty=1:3, cex=0.8)


# 6. Ordinary kriging
pred.locs = iagrid.locs[-idxkeep, ] # predict un-sampled locations only will save time
y.krig.ok = krige.conv(ygeodat, loc=pred.locs, 
                       krige=krige.control(type.krige="ok",obj.m=y.v.wls))
ypred.ok = matrix(0, ntot, 1)
ypred.ok[idxkeep, ] = ydat
ypred.ok[-idxkeep, ] = y.krig.ok$predict
ypredvar.ok = matrix(0,ntot,1)
ypredvar.ok[-idxkeep, ] = y.krig.ok$krige.var

# 6.1 Plot simulated map
par(mar=c(4,3,2,0.5)) # set plot margin
op = par(mfrow=c(2,2), pty="s")
image(matrix(y,length(xg),length(yg)),
      x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),
      main="Full Dataset")
map("state","iowa",lwd=2,add=TRUE)
box()

# 6.2 Plot sampled map
image(matrix(y,length(xg),length(yg)),x=xg,y=yg,
      col=rev(rainbow(100,start=0,end=.7)),main="Sample Data")
image(matrix(ymask,length(xg),length(yg)),
      x=xg,y=yg,col=c("white","transparent"),add=TRUE)
map("state","iowa",lwd=2,add=TRUE)
ypredvar.ok[-idxkeep, ] = y.krig.ok$krige.var
box()

# 6.3 Plot prediction map
image(matrix(ypred.ok,length(xg),length(yg)),x=xg,y=yg,
      col=rev(rainbow(100,start=0,end=.7)),main="Prediction(OK)")
map("state","iowa",lwd=2,add=TRUE)
box()

# 6.4 Plot predicted standard error
image(matrix(sqrt(ypredvar.ok),length(xg),length(yg)),
      x=xg,y=yg,col=rev(rainbow(100,start=0, end=.7)),main="Prediction Standard Errors (OK)")
image(matrix(1-ymask,length(xg),length(yg)),
      x=xg,y=yg,col=c("white","transparent"), add=TRUE)
map("state","iowa",lwd=2,add=TRUE)
box()
par(mfrow = c(1,1))


# 7. Load data
library(geoR)
library(maptools)
rm(list=ls())
data(meuse) # data provided with sp package
?meuse # information about the dataset
coords = cbind(meuse$x,meuse$y) # read coordinates
coordinates(meuse) = coords
bubble(meuse, "lead") # plot points with lead column

data(meuse.grid)
Pb.df = SpatialPixelsDataFrame(points = meuse.grid[c("x", "y")],
                               data = meuse.grid)
coords.grid = cbind(meuse.grid$x,meuse.grid$y) # read coordinates
coordinates(meuse.grid) = coords.grid
plot(meuse.grid) # grid locations to predict


# 8. Estimate the empirical semi-variogram
Pbgeodat = as.geodata(cbind(meuse$x,meuse$y,meuse$lead))
distmat = as.matrix(dist(cbind(meuse$x,meuse$y)))
maxd = max(distmat)/2
maxd
y.v = variog(Pbgeodat,max.dist= maxd)
y.v$beta.ols

# 8.1 Plot empirical semi-variogram
xmax = max(y.v$u,maxd)
ymax = max(y.v$v)
plot(y.v,xlim = c(0,xmax), ylim = c(0,ymax))

# 8.2 Fit an Exponential Semi-Variogram Model
y.v.wls = variofit(y.v, ini.cov.pars=c(10,1000), cov.model="exponential", wei="cressie")
y.v.wls
plot(y.v, xlim = c(0,xmax), ylim = c(0,ymax))
lines(y.v.wls,col=2)

# 9. Kriging
# 9.1 Ordinary Kriging
y.krig.ok = krige.conv(Pbgeodat,loc=coords.grid,
                       krige=krige.control(type.krige="ok",obj.m=y.v.wls))
ypred.ok = y.krig.ok$pred
yVar.ok = y.krig.ok$krige.var
Pb.df$OK_pred = ypred.ok
Pb.df$OK_se = sqrt(yVar.ok)

# 9.2 Set the plot color scheme
maxPb = ceiling(max(ypred.ok)) # maxPb = 13
minPb = floor(min(ypred.ok)) # minPb = 0
bluepal = colorRampPalette(c("azure3", "blue"))
brks = seq(minPb,maxPb,2)
cols = bluepal(length(brks) - 1)
maxSe = ceiling(max(Pb.df$OK_se)) # maxSe = 4
minSe = floor(min(Pb.df$OK_se)) # minSe = 2
brks.se = seq(minSe,maxSe,0.3)
cols.se = bluepal(length(brks.se) - 1)

# 9.3 Plot the Kriging prediction
dev.new(width=8, height=4) # create new plot space with specific size
par(mfrow = c(1, 2), mar=c(1,1,2,1), pty = "s")
image(Pb.df, "OK_pred", col = cols)
symbols(coords, circles=meuse$lead*5, fg="black", inches=F, add=T)
legend("topleft", fill=cols, legend=leglabs(brks), bty="n", cex=0.8)
title(main = "Pb - Ordinary Kriging")
box()

# 9.4 Plot the Kriging standard errors:
image(Pb.df, "OK_se", col = cols.se)
symbols(coords, circles=meuse$lead*0+20, fg="black", inches=F, add=T)
legend("topleft", fill=cols.se, legend=leglabs(brks.se), bty="n", cex=0.8)
title(main = "Pb - Kriging Error(OK)")
box()
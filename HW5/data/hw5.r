# Part I: ESDA Using R
# 1. Load necessary R packages
library(spdep)
library(rgdal) # we used this in Lab 1
library('RSEIS')
library(RColorBrewer)
library(spatialreg)


# 2. Read data
ColData = readOGR(".", "columbus")
coords = coordinates(ColData) # coordinates of the census tract centroids
IDs = row.names(as(ColData, "data.frame"))


# 3. Neighborhood
# 3.1 Adjacency connections
col_nb = poly2nb(ColData, queen = TRUE)
summary(col_nb)
class(col_nb) # check the class of col_nb
help(poly2nb) # check help to see detailed information about poly2nb function

nb = unlist(col_nb)
class(nb) # check the class of nb
nb
nbc1 = table(nb) #table shows the number of neighbors for each polygon
nbc2 = table(nbc1) #table shows the distribution of neighbor connections
#i.e., the number of cases that have certain number of neighbors
nbc1
nbc2
barplot(nbc2,col=rainbow(20), xlab="number of neighbors",ylab="number of cases")

plot(ColData, border = "grey60")
plot(col_nb, coords, pch = 19, cex = 0.6, add = TRUE)
title("Adjacency Connection")
box()

# 3.2 K nearest neighbors
nbk1 = knn2nb(knearneigh(coords, k=1), row.names=IDs)
nbk2 = knn2nb(knearneigh(coords, k=2), row.names=IDs)

dev.off() # close previous plot if there is any
dev.new(width=8, height=4)
par(mfrow = c(1, 2), mar=c(1,1,1,1), pty = "s")
plot(ColData, border="grey60")
plot(nbk1, coords, add=TRUE, pch=19, cex=0.6)
text(bbox(ColData)[1,1] + 0.5, bbox(ColData)[2,2], labels="k=1")
box()
plot(ColData, border="grey60")
plot(nbk2, coords, add=TRUE, pch=19, cex=0.6)
text(bbox(ColData)[1,1] + 0.5, bbox(ColData)[2,2], labels="k=2")
box()
par(mfrow=c(1,1))

# 3.3 Distance-based neighbors
dsts = unlist(nbdists(col_nb, coords))
summary(dsts)
hist(dsts)

dist1=0.5
dist2=0.8
nb1 = dnearneigh(coords, d1=0, d2=dist1, row.names=IDs)
nb2 = dnearneigh(coords, d1=0, d2=dist2, row.names=IDs)

dev.off()
dev.new(width=8, height=4)
par(mfrow = c(1, 2), mar=c(1,1,1,1), pty = "s")
# Band width=0.5
plot(ColData, border="grey60")
plot(nb1, coords, add=TRUE, pch=19, cex=0.5)
text(bbox(ColData)[1,1] + 1, bbox(ColData)[2,2], labels="bw=0.5")
box()
# Band width=0.8
plot(ColData, border="grey60")
plot(nb2, coords, add=TRUE, pch=19, cex=0.8)
text(bbox(ColData)[1,1] + 1, bbox(ColData)[2,2], labels="bw=0.8")
box()
par(mfrow=c(1,1))


# 4. Global spatial autocorrelation
col.W = nb2listw(col_nb, style="W") # weight list
crime = ColData$CRIME
col_I = moran(crime, col.W, length(col_nb), Szero(col.W))
col_I$I # get Moran's I
dev.off()
moran.plot(crime, col.W, pch=19) # Moran scatterplot
moran.test(crime, col.W,zero.policy=T)
geary.test(crime, col.W,zero.policy=T)


# 5. Local spatial autocorrelation
# 5.1 Anselin's Local Moran's I
locI = localmoran(crime, col.W, zero.policy=T)
s = dim(locI)[1]
locI = cbind(locI,matrix(0,s,1))

idx1 = locI[,4]>=1.96 # index for High-high or Low-Low
idx2 = locI[,4]<=-1.96 # index for High-Low or Low-High
idxh = crime>=mean(crime) # high CRIME
idxl = crime<mean(crime) # low CRIME
locI[idx1&idxh,6] = 1 # high-high cluster
locI[idx1&idxl,6] = 2 # low-low cluster
locI[idx2&idxh,6] = 3 # high-low outlier
locI[idx2&idxl,6] = 4 # low-high outlier
ColData$locI = locI[,6]

dev.off()
cols.c = c("ghostwhite","red","blue","lightblue","pink")
clr.c = cols.c[ColData$locI+1]
leglabs = c("Not Significant","High-High","Low-Low","High-Low","Low-High")
plot(ColData,fill=T,col=clr.c)
legend("topleft", fill = cols.c,legend = leglabs,bty = "n", cex = 0.8)
title(main = "LISA Cluster Map(CRIME)")
box()

# 5.2 Getis-Ord's Gi and Gi*
col.W = nb2listw(col_nb, style="W") # same as we used for moran() and localmoran()
G = localG(crime, col.W) # Gi
Gc = matrix(0,s,1)
Gc[G>=1.96] = 1 # index for HH
Gc[G<=-1.96] = 2 # index for LL
ColData$G = Gc

col.W2 = nb2listw(include.self(col_nb), style="W")
G2 = localG(crime, col.W2) # Gi*
Gc2 = matrix(0,s,1)
Gc2[G2>=1.96] = 1 # index for HH
Gc2[G2<=-1.96] = 2 # index for LL
ColData$G2 = Gc2

dev.new(width=8, height=4)
par(mfrow = c(1, 2), mar=c(1,1,1.5,1), pty = "s")

cols.c = c("ghostwhite","red","blue")
clr.c = cols.c[ColData$G+1]
leglabs = c("Not Significant","High-High","Low-Low")
plot(ColData,fill=T,col=clr.c)
legend("topleft", fill = cols.c,legend = leglabs,bty = "n", cex = 0.8)
title(main=expression(paste("Values of the ", G[i], " statistic")))
box()

cols.c = c("ghostwhite","red","blue")
clr.c2 = cols.c[ColData$G2+1]
leglabs = c("Not Significant","High-High","Low-Low")
plot(ColData,fill=T,col=clr.c2)
legend("topleft", fill = cols.c,legend = leglabs,bty = "n", cex = 0.8)
title(main=expression(paste("Values of the ", G[i]^"*", " statistic")))
box()
par(mfrow = c(1, 1))


# Part II: Spatial Regression Models Using R
# 6. Residual dependence test
Col.lm=lm(CRIME ~ INC + HOVAL, data=ColData)
summary(Col.lm)
res=resid(Col.lm)
res=(res-min(res))/diff(range(res)) # standardize
plot(ColData,forcefill=FALSE,col=gray(1-res))
lm.morantest(Col.lm,col.W)


# 7. Spatial Proximity Matrices
col_nb = poly2nb(ColData, queen = TRUE) # recall that we have done this in Part I, 3(1)
col.listw=nb2listw(col_nb,style="W")

col.listw.sym=similar.listw(col.listw)

W=as.matrix(as_dgRMatrix_listw(col.listw))
W.sym=as.matrix(as_dgRMatrix_listw(col.listw.sym))

image(W) #gives image plot of sparse matrix
#(non-zero elements are lighter in color and zeros are white)
image(W.sym)


# 8. Fitting SAR and CAR models
crime.sar=spautolm(CRIME~INC+HOVAL,data= ColData,col.listw,family="SAR")
summary(crime.sar)

crime.car=spautolm(CRIME~INC+HOVAL,data= ColData,col.listw.sym,family="CAR")
summary(crime.car)


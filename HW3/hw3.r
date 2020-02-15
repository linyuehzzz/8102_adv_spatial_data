# 1. Load the R Packages
library(spatstat)
?spatstat
rm(list=ls())

# 2. Exploratory analysis of spatial point patterns
data(japanesepines)
jp = japanesepines
class(jp)
summary(jp)
plot(jp,axes=T,main=" 65 Japanese black pine saplings")

# 2.1 First-order effect 
## kernel estimation (sigma = 0.05)
?density.ppp
jp.Z.05 = density.ppp(jp, 0.05)
par(mar=c(0,0,1,1)) # set plot margin
tl.05 = expression(paste("Kernel Estimation of JP: ", sigma, " = 0.05"))
plot(jp.Z.05, main = tl.05)
plot(jp.Z.05, main="Kernel Estimation of JP: sigma = 0.05")
points(jp, pch="+", col="6")

## kernel estimation (sigma = 0.1)
jp.Z.1 = density.ppp(jp, 0.1)
tl.1 = expression(paste("Kernel Estimation of JP: ", sigma, " = 0.1"))
plot(jp.Z.1, main=tl.1)
plot(jp.Z.1, main="Kernel Estimation of JP: sigma = 0.1")
points(jp, pch="+", col="6")

# 2.2 Second-order effect
## nearest neighbor distances (event-to-nearest-event distances)
jp.ghat = Gest(jp)
g.max = max(jp.ghat$r)
plot(jp.ghat, cbind(rs,theo)~r, main="G Estimates", 
     xlab="r", ylab="G(r)")

## nearest neighbor distances (point-to-nearest-event distances)
jp.fhat = Fest(jp)
f.max = max(jp.fhat$r)
plot(jp.fhat, cbind(rs,theo)~r, main="F Estimates",
     xlab="r", ylab="F(r)")

## K function
?Kest
jp.khat = Kest(jp)
plot(jp.khat, cbind(border, theo)~r, main="K function for JP")

## L function
plot(jp.khat, sqrt(border/pi)-r ~ r, ylab="L(r)",
       main="L function for JP")
abline(h=0, lty=2, col="red")

# 3. Hypothesis testing of spatial point patterns
# 3.1 Simulating CSR
?runifpoint
N = 65 # number of points to generate
r1 = runifpoint(N) #Generate N uniform random points
par(mar=c(0,0,0,0)) # set plot margin
plot(r1, pch="+", main="65 points under CSR")

# 3.2 Testing CSR using nearest neighbor distances and Monte Carlo simulation
## plot Ghat
ghat.env = function(n, s, r, win=owin(c(0,1),c(0,1)))
{
  hold = matrix(0, s, length(r))
  for(i in 1:s)
  {
    hold[i,] = Gest(runifpoint(n, win=win), r=r)$rs
  }
  mn = apply(hold, 2, mean)
  Up = apply(hold, 2, max)
  Down = apply(hold, 2, min)
  return(data.frame(mn, Up, Down))
}
jp.ghat = Gest(jp)
jp.win = window(jp)
plot(jp.ghat, rs~r, main="G estimates")
jp.genv = ghat.env(n=jp$n, s=100, r=jp.ghat$r, win=jp.win)
lines(jp.ghat$r, jp.genv$Up, lty=5, col=2)
lines(jp.ghat$r, jp.genv$Down, lty=5, col=3)

## plot Fhat
fhat.env = function(n, s, r, win=owin(c(0,1),c(0,1)))
{
  hold = matrix(0, s, length(r))
  for(i in 1:s)
  {
    hold[i,] = Fest(runifpoint(n, win=win), r=r)$rs
  }
  mn = apply(hold, 2, mean)
  Up = apply(hold, 2, max)
  Down = apply(hold, 2, min)
  return(data.frame(mn, Up, Down))
}
jp.fhat = Fest(jp)
jp.win = window(jp)
jp.fenv = fhat.env(n=jp$n, s=100, r=jp.fhat$r, jp.win)
plot(jp.fhat,rs~r, main="F estimates")
lines(jp.fhat$r, jp.fenv$Up, lty=5, col=2)
lines(jp.fhat$r, jp.fenv$Down, lty=5, col=3)

## plot Khat
khat.env = function(n, s, r, win=owin(c(0,1),c(0,1)))
{
  hold = matrix(0, s, length(r))
  for(i in 1:s)
  {
    hold[i,] = Kest(runifpoint(n, win=win), r=r)$border
  }
  mn = apply(hold, 2, mean)
  Up = apply(hold, 2, max)
  Down = apply(hold, 2, min)
  return(data.frame(mn, Up, Down))
}
jp.khat = Kest(jp)
jp.win = window(jp)
jp.kenv = khat.env(n=jp$n, s=100, r=jp.khat$r, jp.win)
plot(jp.khat, border~r, main="K function for JP")
lines(jp.khat$r, jp.kenv$Up, lty=5, col=2)
lines(jp.khat$r, jp.kenv$Down, lty=5, col=3)

# 4. Conduct point pattern analysis with another dataset
# 4.1 redwood dataset
data(redwood)
rw = redwood
?redwood
summary(redwood)

## kernel estimation (sigma = 0.05)
?density.ppp
rw.Z.05 = density.ppp(rw, 0.05)
par(mar=c(0,0,1,1)) # set plot margin
tl.05 = expression(paste("Kernel Estimation of RW: ", sigma, " = 0.05"))
plot(rw.Z.05, main = tl.05)
plot(rw.Z.05, main="Kernel Estimation of RW: sigma = 0.05")
points(rw, pch="+", col="6")

## kernel estimation (sigma = 0.1)
rw.Z.1 = density.ppp(rw, 0.1)
tl.1 = expression(paste("Kernel Estimation of RW: ", sigma, " = 0.1"))
plot(rw.Z.1, main=tl.1)
plot(rw.Z.1, main="Kernel Estimation of RW: sigma = 0.1")
points(rw, pch="+", col="6")

## nearest neighbor distances (event-to-nearest-event distances)
rw.ghat = Gest(rw)
g.max = max(rw.ghat$r)
plot(rw.ghat, cbind(rs,theo)~r, main="G Estimates", 
     xlab="r", ylab="G(r)")

## nearest neighbor distances (point-to-nearest-event distances)
rw.fhat = Fest(rw)
f.max = max(rw.fhat$r)
plot(rw.fhat, cbind(rs,theo)~r, main="F Estimates",
     xlab="r", ylab="F(r)")

## K function
?Kest
rw.khat = Kest(rw)
plot(rw.khat, cbind(border, theo)~r, main="K function for RW")

## L function
plot(rw.khat, sqrt(border/pi)-r ~ r, ylab="L(r)",
     main="L function for RW")
abline(h=0, lty=2, col="red")

## test for CSR with simulations (Fhat)
rw.fhat = Fest(rw)
rw.win = window(rw)
rw.fenv = fhat.env(n=rw$n, s=100, r=rw.fhat$r, rw.win)
plot(rw.fhat, rs~r, main="F estimates")
lines(rw.fhat$r, rw.fenv$Up, lty=5, col=2)
lines(rw.fhat$r, rw.fenv$Down, lty=5, col=3)

## test for CSR with simulations (Khat)
rw.khat = Kest(rw)
rw.win = window(rw)
rw.kenv = khat.env(n=rw$n, s=100, r=rw.khat$r, rw.win)
plot(rw.khat, border~r, main="K function for RW")
lines(rw.khat$r, rw.kenv$Up, lty=5, col=2)
lines(rw.khat$r, rw.kenv$Down, lty=5, col=3)

# 4.2 Regular points
rp = rsyst(nx=10) #generate regular pattern
summary(rp)

## kernel estimation (sigma = 0.05)
?density.ppp
rp.Z.05 = density.ppp(rp, 0.05)
par(mar=c(0,0,1,1)) # set plot margin
tl.05 = expression(paste("Kernel Estimation of RP: ", sigma, " = 0.05"))
plot(rp.Z.05, main = tl.05)
plot(rp.Z.05, main="Kernel Estimation of RP: sigma = 0.05")
points(rp, pch="+", col="6")

## kernel estimation (sigma = 0.1)
rp.Z.1 = density.ppp(rp, 0.1)
tl.1 = expression(paste("Kernel Estimation of RP: ", sigma, " = 0.1"))
plot(rp.Z.1, main=tl.1)
plot(rp.Z.1, main="Kernel Estimation of RP: sigma = 0.1")
points(rp, pch="+", col="6")

## nearest neighbor distances (event-to-nearest-event distances)
rp.ghat = Gest(rp)
g.max = max(rp.ghat$r)
plot(rp.ghat, cbind(rs,theo)~r, main="G Estimates", 
     xlab="r", ylab="G(r)")

## nearest neighbor distances (point-to-nearest-event distances)
rp.fhat = Fest(rp)
f.max = max(rp.fhat$r)
plot(rp.fhat, cbind(rs,theo)~r, main="F Estimates",
     xlab="r", ylab="F(r)")

## K function
?Kest
rp.khat = Kest(rp)
plot(rp.khat, cbind(border, theo)~r, main="K function for RP")

## L function
plot(rp.khat, sqrt(border/pi)-r ~ r, ylab="L(r)",
     main="L function for RP")
abline(h=0, lty=2, col="red")

## test for CSR with simulations (Fhat)
rp.fhat = Fest(rp)
rp.win = window(rp)
rp.fenv = fhat.env(n=rp$n, s=100, r=rp.fhat$r, rp.win)
plot(rp.fhat, rs~r, main="F estimates")
lines(rp.fhat$r, rp.fenv$Up, lty=5, col=2)
lines(rp.fhat$r, rp.fenv$Down, lty=5, col=3)

## test for CSR with simulations (Khat)
rp.khat = Kest(rp)
rp.win = window(rp)
rp.kenv = khat.env(n=rp$n, s=100, r=rp.khat$r, rp.win)
plot(rp.khat, border~r, main="K function for RP")
lines(rp.khat$r, rp.kenv$Up, lty=5, col=2)
lines(rp.khat$r, rp.kenv$Down, lty=5, col=3)
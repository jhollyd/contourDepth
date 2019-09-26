#wrapCD.R has R functions for:
#  0) omprXY(n): return a df of n random Gaussian points named xydf
#              and also output the data as file xy<#>.csv.
#  1) omprSD(xydf): compute parallelized simplicial depth integer vector from df xydf. 
#  2) omprTD(xydf): compute parallelized Tukey depth integer vector from df xydf.
#  3) omprCD(xydf): compute both sd and td and output 6-column file xyst.csvincluding i, x, y, ang, sd, td.
#  4) plotD(): plot the data in file xyst.csv using base-r as .png.
#  5) plotDG(): plot the data in file xyst.csv using ggplot from package(tidyverse).
# 
#  Start R
#  Source this file:
#  source(file="wrapCD.R",echo=TRUE,print.eval=TRUE,max.deparse.length=2000)
#

omprXY <- function(n){
#n is the number of 2D random Gaussian data points desired
if(!is.loaded("cd")) {
  dyn.load("cd.so")
}
print("Number of points n is:")
print(n)
rst <- .C("omprXY",
as.integer(n),
x=double(length=as.integer(n)),
y=double(length=as.integer(n)))
#xydf <- data.frame(x,y)
xvec <- rst[2]
yvec <- rst[3]
xydf <- data.frame(xvec,yvec)
write.csv(xydf, "xy.csv")
print("Returning xydf, data frame with n points.")
return (xydf)
#
}


omprSD <- function(xydf)
{
if(!is.loaded("cd")) {
  dyn.load("cd.so")
}
d <- dim(xydf)
rst <- .C("omprSD",
as.integer(d[1]),
as.integer(d[2]),
xydf=as.double(unlist(xydf)),
sdepth=double(length=as.integer(d[1])))
#R makes a sub environment and puts argument xydf in it
#print(ls())
#unsorted vector of simplicial depths from R-C-R structure
sd <- rst[["sdepth"]]
#return(as.integer(sd))
return(sd)
#
}

omprTD <- function(xydf)
{
if(!is.loaded("cd")) {
  dyn.load("cd.so")
}
d <- dim(xydf)
rst <- .C("omprTD",
as.integer(d[1]),
as.integer(d[2]),
xydf=as.double(unlist(xydf)),
tdepth=double(length=as.integer(d[1])))
#R makes a sub environment and puts xydf in it
#print(ls())
#unsorted vector of Tukey depths fr R-C-R structure 
td <- rst[["tdepth"]]
#return(as.integer(td))
return(td)
#
}

omprCD <- function(xydf){
print("Create 6-column file xyst.csv for use in plots.")
if(!is.loaded("cd")) {
  dyn.load("cd.so")
}
d <- dim(xydf)
i <- c(1:d[1])
i <- as.integer(i)
#xydf[1] has dim 10, xydf$x has dim NULL typeof doubles 
#make inputs to atan2 numeric vectors 
x <- xydf$x
y <- xydf$y
#use func td (omprTD in C) to get Tukey depth as integer
td <- omprTD(xydf)
#use func sd (omprSD in C) to get simplicial depth as integer
sd <- omprSD(xydf)
#get angle for plot routine
ang <- unlist(atan2(y,x))
outdf <- data.frame(i,x,y,ang,sd,td)
print("Made a dataframe from vectors i,x,y,ang,sd,td, did write.csv to xyst.csv row.names=FALSE and returned df as CD6col.")
#print(outdf)
write.csv(outdf,file="xyst.csv",row.names=FALSE)
print("  Saved it to 6-column file xyst.csv.")
CD6col <- outdf
return(CD6col) 
#
}

plotD <- function(){
#base packages only
print("plotD() running...will read file xyst.csv")
all <- read.csv(file="xyst.csv",header=TRUE,sep=",")
d <- dim(all)
print("dimensions of input file xyst.csv:")
print(d)
#xyplot, no sorting, use columns for x, y, depth
x <- unlist(all[2])
y <- unlist(all[3])
tdlab <- as.character(unlist(all[6]))
png(filename="plotTD.png",width=480,height=480,pointsize=12,bg="white")
plot(x,y,main="Plot Tukey depths from file xyst.csv, label if numPt <- 50",xlab="x coordinate",ylab="y coordinate",ylim=c(-3,3),xlim=c(-3,3),cex=1,pch=1,col="gray")
if (d[1] <= 50) {
  text(x,y,tdlab)
} else {
  print("No labels output for TD because numPt > 50.")
}
dev.off()
print("Plot saved as plotTD.png.")
sdlab <- as.character(unlist(all[5]))
png(filename="plotSD.png",width=480,height=480,pointsize=12,bg="white")
plot(x,y,main="Plot simplicial depths from file xyst.csv, label if numPt <= 50",xlab="x coordinate",ylab="y coordinate",ylim=c(-3,3),xlim=c(-3,3),cex=1,pch=1,col="gray")
if (d[1] <= 50) {
  text(x,y,sdlab)
} else {
  print("No labels output for SD if numPt > 50.")
}
dev.off()
print("Plot saved as plotSD.png.")
#
}
# 

plotDG <- function(){
#install.packages("tidyverse")
library("ggplot2")
print("plotDG() running with ggplot2...")
df <- read.csv(file="xyst.csv",header=TRUE,sep=',')
d <- dim(df)
print("dimensions of input file xyst.csv:")
print(d)
n <- d[1]
#sort by Tukey depth and angle
df <- df[with(df,order(td,ang)),]
x <- unlist(df[2])
y <- unlist(df[3])
td <- df[,6,drop=FALSE]
vectd <- unlist(td)
sd <- df[,5,drop=FALSE]
vecsd <- unlist(sd)
#plot p1: point plot for Tukey depths, any n
ptd1 <- ggplot(df,aes(x=x,y=y,color=td)) + geom_point()
string <- paste("Tukey depths indicated by point color intensity for n=",as.character(n))
p1 <- ptd1 + labs(x="x coordinate",y="y coordinate") +ggtitle(paste(string))
ggsave("plot_td_points.png")
#plot p2: plot for Tukey contour depths, any n
#label when n < 150
if (n < 150){
ptd <- ggplot(df,aes(x=x,y=y,color=td))
p2td <- ptd + geom_path(aes(group=td))
p3td <- p2td + geom_text(aes(label=td))
string <- paste("Tukey depth contours for n=",as.character(n))
p2 <- p3td + labs(x="x coordinate",y="y coordinate",subtitle="last segment of each contour is not displayed") + ggtitle(paste(string))
ggsave("plot_td_contours.png")
} else {
ptd <- ggplot(df,aes(x=x,y=y,color=td))
p2td <- ptd + geom_path(aes(group=td))
string <- paste("Tukey depth contours for n=",as.character(n))
p2 <- p2td + labs(x="x coordinate",y="y coordinate",subtitle="last segment of each contour is not displayed") + ggtitle(paste(string))
ggsave("plot_td_contours.png")
}
#sort by simplicial depth and angle
df <- df[with(df,order(sd,ang)),]
x <- unlist(df[2])
y <- unlist(df[3])
td <- df[,6,drop=FALSE]
vectd <- unlist(td)
sd <- df[,5,drop=FALSE]
vecsd <- unlist(sd)
#plot p3: point plot for Simplicial depths, any n
psd1 <- ggplot(df,aes(x=x,y=y,color=sd)) + geom_point()
string <- paste("Simplicial depths indicated by point color intensity for n=",as.character(n))
p3 <- psd1 + labs(x="x coordinate",y="y coordinate")+ggtitle(paste(string))
ggsave("plot_sd_points.png")
#plot p4: contour plot for simplicial depth, any n
#label when n < 150
if (n < 150){
  psd <- ggplot(df,aes(x=x,y=y,color=sd)) 
  p2sd <- psd + geom_path(aes(group=sd))
  p3sd <- p2sd + geom_text(aes(label=sd))
  string <- paste("Simplicial depth contours (very few) for n=",as.character(n))
  p4 <- p3sd + labs(x="x coordinate",y="y coordinate",subtitle="last segment of each contour is not displayed") +ggtitle(paste(string))
  ggsave("plot_sd_contours.png")
} else {
  psd <- ggplot(df,aes(x=x,y=y,color=sd))
  p2sd <- psd + geom_path(aes(group=sd))
  string <- paste("Simplicial depth contours (very few) for n=",as.character(n))
  p4 <- p2sd + labs(x="x coordinate",y="y coordinate",subtitle="last segment of each contour is not displayed") + ggtitle(paste(string))
  ggsave("plot_sd_contours.png")
}
#
print("Data sorted by depth and angle and plots saved as: ")
print("  plot_td_points.png, plot_td_contours.png,")
print("  plot_sd_points.png, plot_sd_contours.png.")
return(p4)
}

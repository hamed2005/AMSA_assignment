library(plotrix)    ## for biplot
library(rpart)      ## for Tree-based Modeling
library(MASS)       ## for linear discriminant analysis
library(DMwR)      ## for outlier detection

## setting the working directory
setwd("/home/hamed/KUL/Multivar/Project/HomoTarget")
##

## reading the CSV files for positive and negative instances
colnames <- c("mirnaSeq","targetSeq", "totalScore", "seedScore", "WCPairs", "WobblePairs", "mismatches", "NumberBulges", "A", "C", "G", "U", "AU", "minFreeEnergy")
pos <- read.csv("dataset.csv.p", header = FALSE, col.names = colnames)
neg <- read.csv("dataset.csv.n", header = FALSE, col.names = colnames)

## appending two dataframes into one
mirnaDF <- rbind(pos,neg)
## omitting the sequence columns
mirna <- mirnaDF[,3:14]

## overview
summary(mirna)

## pairs plot
pairs(mirna, pch=20, col="#383838")


### outlier detection
outlier.scores <- lofactor(mirna, k = 5)    ##K = # of neighbors used in the calculation of the local outlier factors
plot(density(outlier.scores))
## pick top 5 as outliers
outliers <- order(outlier.scores, decreasing=T)[1:5]
## who are outliers
print(outliers)

n <- nrow(mirna)
labels <- 1:n
labels[-outliers] <- ""

## showing outliers in biplot
x<-mirna

xm<-apply(x,2,mean)
y<-sweep(x,2,xm)
ss<-(t(y)%*%as.matrix(y))
s<-ss/(nrow(x)-1)
d<-(diag(ss))^(-1/2)
e<-diag(d,nrow=ncol(x),ncol=ncol(x))
z<-as.matrix(y)%*%e
r<-t(z)%*%z
q<-svd(z)
gfd<-((q$d[1])+(q$d[2]))/sum(q$d)
gfz<-(((q$d[1])^2)+((q$d[2])^2))/sum((q$d)^2)
gfr<-(((q$d[1])^4)+((q$d[2])^4))/sum((q$d)^4)
l<-diag(q$d,nrow=ncol(x),ncol=ncol(x))
R.B<-q$u        #scores matrix
C.B<-q$v%*%l    #loadings

#possibility to stretch scores by a scale factor
scalefactor<-3.5
R.B<-q$u *scalefactor

par(mar=c(4,4,4,4),pty='s',oma=c(5,0,0,0),font=2)
plot(R.B[ ,1],R.B[ ,2],axes=F,xlim=c(-1,1),ylim=c(-1,1),xlab=' ',ylab=' ',cex=.8)
mtext('First component',side=1,line=3,cex=.8)
mtext('Second component',side=2,line=3,cex=.8)
axis(1,at=c(-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1),cex=.8)
axis(2,at=c(-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1),cex=.8)
box( )

###display outlier obs number on the biplot
text(R.B[,1]-.05,R.B[,2]+.05,as.character(labels),cex=0.5)
points(R.B[,1],R.B[,2],pch=".")


points(C.B[,1],C.B[,2],pch=".")
text(C.B[,1]-.05,C.B[,2]+.05,as.character(dimnames(x)[[2]]),cex=0.8)

## drawing the arrows
for (i in seq(1,nrow(C.B),by=1))
        arrows(0,0,C.B[i,1],C.B[i,2])

#Draw circle unit
draw.circle(0,0,1,border='black')

#removing outliers
mirna <- mirna[-103,]
mirna <- mirna[-94,]


## centering and normalizing the data (unit variance)
mirna.c <- scale(mirna, scale = F, center = T)
mirna.mean <- attr(mirna.c, "scaled:center")
mirna <- mirna.c[,]
mirna.n <- scale(mirna, scale = T, center = T)

## variance-covariance and correlation matrix extraction
mirna.cov <- cov(mirna)     ##same as cor matrix
mirna.cor <- cor(mirna)

### PCA
mirna.pca <- princomp(mirna, cor = T)    ##using correlation matrix
mirna.pca
screeplot(mirna.pca, type="lines")

mirna.pca$loadings

summary(mirna.pca)

par(mfrow = c(3,3))

plot(mirna.pca$scores[,1]~mirna.pca$scores[,2], xlim = c(-5,+5), ylim = c(-5,+5) ,cex=0.5)
plot(mirna.pca$scores[,2]~mirna.pca$scores[,3], xlim = c(-5,+5), ylim = c(-5,+5) ,cex=0.5)
plot(mirna.pca$scores[,3]~mirna.pca$scores[,4], xlim = c(-5,+5), ylim = c(-5,+5) ,cex=0.5)
plot(mirna.pca$scores[,4]~mirna.pca$scores[,5], xlim = c(-5,+5), ylim = c(-5,+5) ,cex=0.5)
plot(mirna.pca$scores[,5]~mirna.pca$scores[,6], xlim = c(-5,+5), ylim = c(-5,+5) ,cex=0.5)
plot(mirna.pca$scores[,6]~mirna.pca$scores[,7], xlim = c(-5,+5), ylim = c(-5,+5) ,cex=0.5)
plot(mirna.pca$scores[,7]~mirna.pca$scores[,8], xlim = c(-5,+5), ylim = c(-5,+5) ,cex=0.5)
plot(mirna.pca$scores[,8]~mirna.pca$scores[,9], xlim = c(-5,+5), ylim = c(-5,+5) ,cex=0.5)
plot(mirna.pca$scores[,9]~mirna.pca$scores[,10], xlim = c(-5,+5), ylim = c(-5,+5) ,cex=0.5)

par(mfrow =c(1,1))



### Biplot
x <- mirna
xm<-apply(x,2,mean)
y<-sweep(x,2,xm)
ss<-(t(y)%*%y)
s<-ss/(nrow(x)-1)
d<-(diag(ss))^(-1/2)
e<-diag(d,nrow=ncol(x),ncol=ncol(x))
z<-y%*%e
r<-t(z)%*%z
q<-svd(z)
gfd<-((q$d[1])+(q$d[2]))/sum(q$d)
gfz<-(((q$d[1])^2)+((q$d[2])^2))/sum((q$d)^2)
gfr<-(((q$d[1])^4)+((q$d[2])^4))/sum((q$d)^4)
l<-diag(q$d,nrow=ncol(x),ncol=ncol(x))
R.B<-q$u        #scores matrix
C.B<-q$v%*%l    #loadings

#possibility to stretch scores by a scale factor
scalefactor<-3.5
R.B<-q$u *scalefactor
        
par(mar=c(4,4,4,4),pty='s',oma=c(5,0,0,0),font=2)
plot(R.B[ ,1],R.B[ ,2],axes=F,xlim=c(-1,1),ylim=c(-1,1),xlab=' ',ylab=' ',cex=.8)
mtext('First component',side=1,line=3,cex=.8)
mtext('Second component',side=2,line=3,cex=.8)
axis(1,at=c(-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1),cex=.8)
axis(2,at=c(-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1),cex=.8)
box( )

points(R.B[,1],R.B[,2],pch=".")


points(C.B[,1],C.B[,2],pch=".")
text(C.B[,1]-.05,C.B[,2]+.05,as.character(dimnames(x)[[2]]),cex=0.8)

## drawing the arrows
for (i in seq(1,nrow(C.B),by=1))
        arrows(0,0,C.B[i,1],C.B[i,2])

#Draw circle unit
draw.circle(0,0,1,border='black')
        
mtext('PCA Biplot',side=1,outer=T,cex=1,line=3)
#PCA.biplot(mirna)



### Tree based modelling

par(mar=c(5,4,4,2) + 0.1)    ## to reset margin settings
par(oma=c(3,3,3,3))

##adding class variables
pos<-cbind(pos,'positive')
colnames(pos)[15]<-"class"
neg<-cbind(neg,'negative')
colnames(neg)[15]<-"class"
##merging pos and neg to build the full dataset
mirna <- rbind(pos,neg)[,3:15]

mirna$class <- as.factor(mirna$class)

##removing the outliers
mirna <- mirna[-103,]
mirna <- mirna[-94,]

##desicion tree construction using rpart

mirna.rDTree <- rpart(class ~ ., 
                      data = mirna , method="class",
                      parms=list(prior=c(0.5, 0.5)))

##parameters #######TO BE SET
rpart.control(minsplit=10,maxcompete=5,maxsurrogate=5,usesurrogate=2,xval=20,maxdepth=30)

#summary(mirna.rDTree)
par(xpd=NA)    ## to show the labels completely in output

plot(mirna.rDTree, main = "Original Decision Tree")
text(mirna.rDTree ,splits=T , use.n = T, all=F, cex = 0.7)


##results of cross validation 
plotcp(mirna.rDTree)

##plotting cost complexity in reation to number of splits
plot(mirna.rDTree$cptable[,2],mirna.rDTree$cptable[,1],
     xlab='Numberbof splits',ylab='Cost complexity parameter,cp')

##pruning
mirna.rDTree.pruned <- prune(mirna.rDTree, cp = 0.062)

plot(mirna.rDTree.pruned, uniform = F, main = "Pruned Decision Tree")
text(mirna.rDTree.pruned, use.n = T, all = T, cex = 0.7)

##plotting two trees side by side
par(mfrow = c(1,2), xpd = NA)
plot(mirna.rDTree, uniform = F, main = "Original Decision Tree")
text(mirna.rDTree ,splits=T ,all=F, cex = 0.7)
plot(mirna.rDTree.pruned, uniform = F, main = "Pruned Decision Tree")
text(mirna.rDTree.pruned, cex = 0.7)

print(mirna.rDTree.pruned)

##misclassification error and CV for unpruned tree
out <- predict(mirna.rDTree) # predict probabilities
pred.response <- colnames(out)[max.col(out, ties.method = c("random"))] # predict response
mean(mirna$class != pred.response) # % misclassification error
table(mirna$class, pred.response)

##misclassification error and CV for pruned tree
out <- predict(mirna.rDTree.pruned) # predict probabilities
pred.response <- colnames(out)[max.col(out, ties.method = c("random"))] # predict response
mean(mirna$class != pred.response) # % misclassification error
table(mirna$class, pred.response)

### Discriminant Analysis

mirna.lda <- lda(class~ + totalScore + seedScore + WCPairs + WobblePairs + mismatches + NumberBulges + A + C + G + U + AU + minFreeEnergy, data=mirna)
#####NOTE:: VARIABLES ARE COLLINEAR######
## LDA using first 9 PCs
mirna.lda <- lda(mirna$class~+mirna.pca$scores[,1:9])
plot(mirna.lda)
## cross-validation
mirna.lda <- lda(mirna$class~+mirna.pca$scores[,1:9], CV=T)
head(mirna.lda$posterior)
table(mirna$class, mirna.lda$class, dnn = c("From", "classified into"))

## lda using the variables used by tree
mirna.lda.tree <- lda(class ~ + totalScore + seedScore + WCPairs + WobblePairs + NumberBulges + C + G + AU + minFreeEnergy, data = mirna, CV=T)
head(mirna.lda.tree$posterior)
table(mirna$class, mirna.lda.tree$class, dnn = c("From", "classified into"))
### Clustering

dlimit <- min(apply(mirna.n,2,min)) -0.1
ulimit <- max(apply(mirna.n,2,max)) +0.1

par(mfrow = c(1,1))
# Set up blank plot for profiles
plot(c(0,12),c(dlimit, ulimit),type="n",xlab="Variables",ylab="Value",main="Profile Plot")

for (k in (1:200))
{
        points(1:12,mirna.n[k,],type="l")
}

## modified andrews plot
plot(c(-pi,pi),c(-7,6),type="n",xlab="Variables",ylab="Value",main="Modified Andrews Plot")
t <-seq(-pi,pi,length=500)
for (k in (1:nrow(mirna.n)))
{
        crseqm <- (1/sqrt(2))*(mirna.n[k,1]+mirna.n[k,2]*(sin(t)+cos(t)) + mirna.n[k,3]*(sin(t)-cos(t))+mirna.n[k,4]*(sin(2*t)+cos(2*t)))
        points(t,crseqm,type="l")
}

## stars plot
stars(mirna.n,draw.segments = F, cex=0.25, scale = T)

## hierarchical clustering

#computing the cluster object
mirna.hclust <- hclust(dist(mirna.c), method="average")  # AVARAGE LINK

# plclust plots the dendrogram for the cluster object.
plot(mirna.hclust,xlab="mirna Data",ylab="Average Link Distance", sub="", labels = F)
# cutree finds the cluster groupings for a cluster object for a 
# specifed (k) number of clusters.
mirna.hclust.cut <- cutree(mirna.hclust,k=3)
# You might want to see how many observations in each cluster.
table(mirna.hclust.cut)

# projection on the first 2 principal components
#clusplot(mirna.c, mirna.hclust.cut,stand=TRUE,labels=5,main="Average link")

# Compute pairwise plots for the data with color coded cluster ids.
par(mfrow=c(1,2))
plot(mirna$totalScore,mirna$seedScore, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$totalScore,mirna$WCPairs, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$totalScore,mirna$WobblePairs, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$totalScore,mirna$mismatches, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$totalScore,mirna$NumberBulges, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$totalScore,mirna$A, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$totalScore,mirna$C, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$totalScore,mirna$G, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$totalScore,mirna$U, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$totalScore,mirna$AU, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$totalScore,mirna$minFreeEnergy, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")

plot(mirna$minFreeEnergy,mirna$seedScore, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$minFreeEnergy,mirna$WCPairs, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$minFreeEnergy,mirna$WobblePairs, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$minFreeEnergy,mirna$mismatches, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$minFreeEnergy,mirna$NumberBulges, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$minFreeEnergy,mirna$A, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$minFreeEnergy,mirna$C, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$minFreeEnergy,mirna$G, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$minFreeEnergy,mirna$U, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$minFreeEnergy,mirna$AU, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")
plot(mirna$minFreeEnergy,mirna$totalScore, col = mirna.hclust.cut, cex=1.5,main="Average Link Clusters")

## KMeans
mirna.kmclust <- kmeans(mirna.c, 5, 20)

# Perform the pairwise plots with color coded cluster id included.

plot(mirna$minFreeEnergy,mirna$totalScore, col = mirna.kmclust$cluster)
points(mirna.kmclust$centers[,c(1,2)], col = 1:2, pch ="+", cex=2)

plot(mirna$minFreeEnergy,mirna$A, col = mirna.kmclust$cluster)
points(mirna.kmclust$centers[,c(1,3)], col = 1:2, pch ="+", cex=2)

plot(mirna$minFreeEnergy,mirna$AU, col = mirna.kmclust$cluster)
points(mirna.kmclust$centers[,c(1,4)], col = 1:2, pch ="+", cex=2)

plot(mirna$minFreeEnergy,mirna$NumberBulges, col = mirna.kmclust$cluster)
points(mirna.kmclust$centers[,c(2,3)], col = 1:2, pch ="+", cex=2)

plot(mirna$WobblePairs,mirna$AU, col = mirna.kmclust$cluster)
points(mirna.kmclust$centers[,c(2,4)], col = 1:2, pch ="+", cex=2)

# Run the kmeans clustering on the first two PCs
par(mfrow=c(1,1))
mirna.pca <- princomp(mirna.n)
mirna.pckclust <- kmeans(mirna.pca$scores, 2, 20)
plot(mirna.pca$scores[,1:2], col = mirna.pckclust$cluster)
points(mirna.pckclust$centers[,c(1,2)], col = 1:2, pch ="+", cex=2)
title ("K Means Cluster on first two PCs")
#clusplot(mirna.n,mirna.pckclust$cluster,stand=TRUE,labels=2,main="kmeans, 2 clusters")

### Factor Analysis

### It was mentioned in the paper that the dataset contains features from 3 classes.
### Lets try to retain these 3 groups as a proof of concept.
### Classes are : Structural, thermodynamic and positional.

##FA using PCA extraction --JUST FOR TESTING--
mirna.pca <- eigen(mirna.cor)
mirna.p <- mirna.pca$vectors[,1:5]
mirna.d <- diag(sqrt(mirna.pca$values[1:5]))

#factor loadings
mirna.B <- mirna.p%*%mirna.d
rownames(mirna.B) <- colnames(mirna[,1:12])
colnames(mirna.B) <- c("C1","C2","C3", "C4", "C5")    ##3 classes
#psi residuals correlation
mirna.psi <- mirna.cor - mirna.B%*%t(mirna.B)
#communalities
1-diag(mirna.psi)
#RMS
mirna.residuals <- (mirna.psi - diag(mirna.psi))^2
mirna.RMS.overall <- sqrt(sum(mirna.residuals)/length(mirna.residuals))    ##or /(ncol(mirna.residuals)*(ncol(mirna.residuals)-1))
mirna.RMS.overall


##FA using Maximum Likelihood
mirna.fa <- factanal(mirna[,1:12], 2, rotation = "promax", control = list(nstart=30))    ## nstart: number of starting points. o.w. ml wouldnt converge

mirna.fa$loadings

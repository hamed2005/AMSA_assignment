##
setwd("/home/hamed/KUL/Multivar/Project/HomoTarget")
##

## reading the CSV files for positive and negative instances
pos <- read.csv("dataset.csv.p", header = FALSE, col.names = c("mirnaSeq","targetSeq", "tatalScore", "seedScore", "WCPairs", "WobblePairs", "mismatches", "NumberBulges", "A", "C", "G", "U", "AU", "minFreeEnergy"))
neg <- read.csv("dataset.csv.n", header = FALSE, col.names = c("mirnaSeq","targetSeq", "tatalScore", "seedScore", "WCPairs", "WobblePairs", "mismatches", "NumberBulges", "A", "C", "G", "U", "AU", "minFreeEnergy"))

## appending two dataframes into one
mirnaDF <- rbind(pos,neg)

## overview
summary(mirnaDF)
pairs(mirnaDF)

mirna <- mirnaDF[,3:14]

## centering and normalizing the data (unit variance)
mirna.c <- scale(mirna, scale = F, center = T)
mirna.mean <- attr(mirna.c, "scaled:center")
mirna <- mirna.c[,]

## variance-covariance and correlation matrix extraction
mirna.cov <- cov(mirna)     ##same as cor matrix
mirna.cor <- cor(mirna)

## PCA
mirna.pca <- princomp(mirna, cor = T)    ##using correlation matrix
mirna.pca
screeplot(mirna.pca, type="lines")

mirna.pca$sdev
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


## Factor Analysis
### It was mentioned in the paper that the dataset contains features from 3 classes.
### Lets try to retain these 3 groups as a proof of concept.
### Classes are : Structural, thermodynamic and positional.

##FA using PCA extraction --JUST FOR TESTING--
mirna.pca <- eigen(mirna.cor)
mirna.p <- mirna.pca$vectors[,1:5]
mirna.d <- diag(sqrt(mirna.pca$values[1:5]))

#factor loadings
mirna.B <- mirna.p%*%mirna.d
rownames(mirna.B) <- colnames(mirna)
colnames(mirna.B) <- c("C1","C2","C3", "C4", "C5")    ##3 classes
#psi residuals correlation
mirna.psi <- mirna.cor - mirna.B%*%t(mirna.B)
#communalities
1-diag(mirna.psi)
#RMS
mirna.residuals <- (mirna.psi - diag(mirna.psi))^2
mirna.RMS.overall <- sqrt(sum(mirna.residuals)/length(mirna.residuals))
mirna.RMS.overall


##FA using Maximum Likelihood
mirna.fa <- factanal(mirna, 2, rotation = "promax", control = list(nstart=10))    ## nstart: number of starting points. o.w. ml wouldnt converge


## Biplot



## Discriminant Analysis



## Tree based modelling



## Clustering



## HICUUP!



## Multidimensional scaling



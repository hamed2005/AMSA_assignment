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
mirna.c <- scale(mirna, scale = T)
mirna.mean <- attr(mirna.c, "scaled:center")
mirna <- mirna.c[,]

## variance-covariance and correlation matrix extraction
mirna.cov <- cov(mirna)     ##same as cor matrix
mirna.cor <- cor(mirna)

## PCA
mirna.pca <- princomp(mirna)    ##using correlation matrix
mirna.pca
screeplot(mirna.pca, type="lines")

mirna.pca$sdev
mirna.pca$loadings


## Factor Analysis
### It was mentioned in the paper that the dataset contains features from 3 classes.
### Lets try to retain these 3 groups as a proof of concept.
### Classes are : Structural, thermodynamic and positional.

##FA using PCA extraction
mirna.pca <- eigen(mirna.cor)
mirna.p <- mirna.pca$vectors[,1:3]
mirna.d <- diag(sqrt(mirna.pca$values[1:3]))

#factor loadings
mirna.B <- mirna.p%*%mirna.d
rownames(mirna.B) <- colnames(mirna)
colnames(mirna.B) <- c("C1","C2","C3")    ##3 classes
#psi residuals correlation
mirna.psi <- mirna.cor - mirna.B%*%t(mirna.B)






## Biplot



## Discriminant Analysis



## Tree based modelling



## Clustering



## HICUUP!



## Multidimensional scaling



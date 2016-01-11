library(plotrix)    ## for biplot


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

PCA.biplot<- function(x) {
        #x is the matrix to biplot, x is numeric, thus variables on ratio or interval scales
        #x has dimnames(x)[[2]] defined (and eventually dimnames(x)[[1]] defined)
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
        
        #ability to plot rownames on the plot, dimnames(x)[[1]] necessary
        
        ###could be uncommented to display obs numbers on plot
        #text(R.B[,1]-.05,R.B[,2]+.05,as.character(dimnames(x)[[1]]),cex=0.5)
        points(R.B[,1],R.B[,2],pch=".")
        
        
        points(C.B[,1],C.B[,2],pch=".")
        text(C.B[,1]-.05,C.B[,2]+.05,as.character(dimnames(x)[[2]]),cex=0.8)
        
        ## drawing the arrows
        for (i in seq(1,nrow(C.B),by=1))
                arrows(0,0,C.B[i,1],C.B[i,2])
        
        #Draw circle unit
        draw.circle(0,0,1,border='black')
        
        ##printing some statistics to output
        mtext('PCA Biplot',side=1,outer=T,cex=1,line=3)
        results<-list('correlation matrix'=r,'column effects'=C.B,'row effects'=R.B)
        cat('The goodness of fit for the correlation matrix is',gfr,'for the centered, standardized design matrix',gfz,'and for the Mahalanobis distances is',gfd,' ') 
        results
}

PCA.biplot(mirna)



## Tree based modelling



## Discriminant Analysis


## Clustering



## HICUUP!



## Multidimensional scaling



pos <- read.csv("dataset.csv.0", header = FALSE, col.names = c("mirnaSeq","targetSeq", "tatalScore", "seedScore", "WCPairs", "WobblePairs", "mismatches", "NumberBulges", "A", "C", "G", "U", "AU", "minFreeEnergy"))
neg <- read.csv("dataset.csv.1", header = FALSE, col.names = c("mirnaSeq","targetSeq", "tatalScore", "seedScore", "WCPairs", "WobblePairs", "mismatches", "NumberBulges", "A", "C", "G", "U", "AU", "minFreeEnergy"))

mirnaDF <- rbind(pos,neg)

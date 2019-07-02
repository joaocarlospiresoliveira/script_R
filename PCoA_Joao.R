################################################################################
### PCoA - Análise de Coordenadas Principais 
################################################################################

# Load the required packages
# (vegan must be loaded after ade4 to avoid some conflicts)
library(ade4)
library(vegan)
library(ape)
library(recluster)


# Import the data from CSV files
spe<-read.table(file.choose(),row.names=1,header=T,sep=",")
dim(spe)
View(spe)
spe<-as.matrix(spe)
spe<-decostand(spe, 'normalize')
View(spe)
# Diagnosing possible problems with matrix:
# generate species and sites sums these will be NA is any are missing 
csum <- colSums(spe)
# check if any are missing 
any(is.na(csum)) 
[1] TRUE 
# yes, some missing, so which ones? 
which(is.na(csum)) 
# Check how many are missing 
summary(spe[, c("x")]) #change 'x' by the column name

# Removing unicates:
unicates<-apply(spe,2,sum) 
spe<-spe[,which(unicates>1)] ###removing unicates
dim(spe)
edit(spe)

# Check if any site has zero sum:
rsum <- rowSums(spe)
rsum
# If any, you should remove it (them) from your matrix:
#spe<-spe[cbind(-15),]
#dim(spe)


# PCoA on a percentage difference (Bray-Curtis) dissimilarity matrix of 
# species
# *********************************************************************

spe.bray <- vegdist(spe, 'euc')
spe.b.pcoa <- cmdscale(spe.bray, k=(nrow(spe)-1), eig=TRUE)
# Plot of the sites
dev.new(title="PCoA on fish species - Percentage difference")
ordiplot(scores(spe.b.pcoa, choices=c(1,2)), type="t", main="PCoA with species weighted averages")
abline(h=0, lty=3)
abline(v=0, lty=3)
# Add weighted average projection of species
spe.wa <- wascores(spe.b.pcoa$points[,1:2], spe)
text(spe.wa, rownames(spe.wa), cex=0.7, col="red")


# PCoA and projection of species vectors using function pcoa()
spe.h.pcoa <- pcoa(dist(spe, 'euclidean') )
# Biplots
dev.new(title="PCoA with species vectors", width=14, height=8)
par(mfrow=c(1,1))
# First biplot: Hellinger-transformed species data
biplot.pcoa(spe.h.pcoa, spe, dir.axis1=-1) 
abline(h=0, lty=3)
abline(v=0, lty=3)



# Comparison of PCoA results with Euclidean and non-Euclidean
# dissimilarity matrices
# ***********************************************************

# PCoA on a Hellinger distance matrix
is.euclid(dist(spe))
summary(spe.h.pcoa) 
spe.h.pcoa$values

# PCoA on a percentage difference (Bray-Curtis) dissimilarity matrix
is.euclid(spe.bray)
spe.bray.pcoa <- pcoa(spe.bray) 
spe.bray.pcoa$values		# Observe eigenvalues 18 and following

# PCoA on the square root of a percentage difference (Bray-Curtis)
# dissimilarity matrix
is.euclid(sqrt(spe.bray))
spe.braysq.pcoa <- pcoa(sqrt(spe.bray))
spe.braysq.pcoa$values	# Observe the eigenvalues

# PCoA on a percentage difference (Bray-Curtis) dissimilarity matrix with
# Lingoes correction
spe.brayl.pcoa <- pcoa(spe.bray, correction="lingoes")
spe.brayl.pcoa$values		# Observe the eigenvalues, col. 1 and 2

# PCoA on a percentage difference (Bray-Curtis) dissimilarity matrix with
# Cailliez correction
spe.brayc.pcoa <- pcoa(spe.bray, correction="cailliez")
spe.brayc.pcoa$values		# Observe the eigenvalues, col. 1 and 2



########### PCoA with Simpson distance #################
#Comparison of PCoA results with Simpson distance
simpdiss<- recluster.dist(spe1) #converte os dados originais em uma matriz com distância de Simpson
simpdiss <- as.matrix(simpdiss)
spe.simp.pcoa <- pcoa(simpdiss, correction="lingoes")
spe.simp.pcoa <- cmdscale(simpdiss, eig=TRUE)

treat=c(rep("Treatment1",7),rep("Treatment2",10),rep("Treatment3",6))
ordiplot(spe.simp.pcoa,type="n")
orditorp(spe.simp.pcoa,display="sites",groups=treat,draw="polygon", col=c("green4","cyan",
                                                 "gold1"))


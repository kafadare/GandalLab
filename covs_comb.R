setwd("~/project-gandalm/GandalLab")
output_dir <- ("~/project-gandalm/comb_data/processed")
plot_dir <- "/u/home/k/kafadare/project-gandalm/plots"
source("data_functions.R")
source("analysis_fcts.R")
# List of required packages
required_packages <- c("magrittr", "dplyr", "data.table", "ggplot2", "reshape2","PCAForQTL")
# Install or load missing packages
load_install_pkg(required_packages)

# Load data
datExpr <- fread(args$expr, data.table=F)
rownames(datExpr) <- datExpr$V1
datExpr <- datExpr[,-1]

meta <- read.table(args$meta, header=T, stringsAsFactors=FALSE)


std_datExpr <- standardize(t(datExpr))

# Generate covariates - PCA
prcompResult<-prcomp(std_datExpr,center=F,scale.=F) #This should take less than a minute.
PCs <-prcompResult$x 

##choosing # of PCs

#`runElbow()` implements an automatic elbow detection method (based on distance to the diagonal line). In our experience, the number of PCs chosen using this method **tends to match the visual elbow point reasonably well**. To run this function, simply input the previously obtained `prcompResult`. This method selects $K=12$ for our example data.
resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)
print(resultRunElbow)

#`runBE()` implements the BE algorithm, a permutation-based approach for choosing $K$ in PCA. Intuitively, the BE algorithm retains PCs that explain more variance in the data than by random chance and discards those that do not. In our experience, the number of PCs chosen via BE **tends to signal an upper bound of the reasonable number of PCs to choose**. To run this function, we need to input the data matrix (must be observation by feature) and may optionally specify `B`, the number of permutations (default is 20), and `alpha`, the significance level (default is $0.05$). We may also specify `mc.cores` to change how many cores are used for parallel computing (default is `B` or the number of available cores minus 1, whichever is smaller). For reproducibility, Linux and Mac users must change the random number generator (RNG) type (unless `mc.cores` is 1) and set the seed. On the other hand, Windows users must set `mc.cores=1` to avoid error and set the seed for reproducibility (no need to change the RNG type). The BE method selects $K=29$ for our example data.

RNGkind("L'Ecuyer-CMRG")
set.seed(1)
resultRunBE<-PCAForQTL::runBE(expr,B=20,alpha=0.05)
print(resultRunBE$numOfPCsChosen)

#visualize number of PCs chosen
K_elbow<-resultRunElbow #12.
K_BE<-resultRunBE$numOfPCsChosen #29.
K_GTEx<-60 #GTEx uses 60 PEER factors, and they are almost identical to the top 60 PCs.
PCAForQTL::makeScreePlot(prcompResult,labels=c("Elbow","BE","GTEx"),values=c(K_elbow,K_BE,K_GTEx),
                         titleText="Colon - Transverse")

ggplot2::ggsave("./Colon_Transverse.jpg",width=16,height=11,unit="cm")

##Filter out known covariates
identical(rownames(meta),rownames(datExpr)) #TRUE is good.
#Suppose we have decided to use the number of PCs chosen via BE for our example data.
PCsTop<-PCs[,1:K_BE] #368*29.

#Now we use `filterKnownCovariates()` to **filter out the known covariates that are captured well by the top PCs** (unadjusted $R^2\geq 0.9$ by default). This function returns the known covariates that should be kept. We use unadjusted $R^2$ instead of adjusted $R^2$ because we do not want to penalize for model complexity here. The cutoff value may be customized using the argument `unadjustedR2_cutoff`.
knownCovariatesFiltered<-PCAForQTL::filterKnownCovariates(knownCovariates,PCsTop,unadjustedR2_cutoff=0.9)

#Finally, we combine the remaining known covariates and the top PCs and use them as covariates in the QTL analysis. To avoid potential numerical inaccuracies in the QTL analysis, you may optionally scale the PCs to unit variance, though theoretically this would not change the QTL result from regression-based methods such as Matrix eQTL and FastQTL.
PCsTop<-scale(PCsTop) #Optional. Could be helpful for avoiding numerical inaccuracies.
covariatesToUse<-cbind(knownCovariatesFiltered,PCsTop)

# Generate covariates (age + sex + gPC + HCP)
n_subj <- ncol(datExpr)
cov <- matrix(nrow = 2 + num_hcp + num_gpc , ncol = n_subj + 1)
cov <- as.data.frame(cov)

colnames(cov) <- c("id",colnames(datExpr))
for (i in 1:num_gpc) {
  cov[i,1] <- paste0("PC",i)
}
cov[num_gpc + 1, 1] <- "sex"
for (i in (num_gpc + 2):(num_gpc + 1 + num_hcp)) {
  cov[i,1] <- paste0("HCP",i-1-num_gpc)
}
cov[nrow(cov),1] <- "age"

# gPC
geno_pc <- geno_pc[,-1]
rownames(geno_pc) <- geno_pc[,1]
geno_pc <- geno_pc[,-1]
# first 2504 are 1kg
geno_pc <- geno_pc[2505:nrow(geno_pc),]

for (i in 2:ncol(cov)){
  sample <- colnames(cov)[i]
  for (j in 1:num_gpc) {
    cov[j,sample] <- geno_pc[sample,j]
  }
}

# sex and age
row.names(meta) <- meta[,1]
for (i in 2:ncol(cov)){
  sample <- colnames(cov)[i]
  cov[nrow(cov),sample] <- meta[sample,"Age"]
  cov[1 + num_gpc, sample] <- meta[sample,"inferSex"]
}

# HCP
hcp_cov <- hcpMat[,2:(num_hcp+1)]
hcp_cov <- as.data.frame(hcp_cov)
for (i in 2:ncol(cov)) {
  sample <- colnames(cov)[i]
  for(j in (num_gpc + 2):(num_gpc + 1 + num_hcp)) {
    cov[j,sample] <- hcp_cov[sample,j-1-num_gpc]
  }
}

write.table(cov, paste0(output_dir, num_hcp, "hcp_cov.txt"), quote = F, row.names = F, col.names = T, sep = "\t")
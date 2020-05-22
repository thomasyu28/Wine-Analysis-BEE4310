## Predicting red wine quality using physicochemical characteristics  
## Thomas Yu
## Created 5/18/20
## Updated 5/22/20
# ------------------------------------------------------------------
library("car")
library("glmnet")

# read in and clean data
wwine.dat = read.csv2("winequality-white.csv")
rwine.dat = read.csv2("winequality-red.csv")
for (i in 1:11){
  wwine.dat[,i] = as.numeric(levels(wwine.dat[,i]))[wwine.dat[,i]]
  rwine.dat[,i] = as.numeric(levels(rwine.dat[,i]))[rwine.dat[,i]]
}
wwine.dat.scaled = as.data.frame(scale(wwine.dat))
rwine.dat.scaled = as.data.frame(scale(rwine.dat))

## data visualization

# histogram for quality variable
png("Figure 1.png", width = 4, height = 4, units = 'in', res = 500)
hist(wwine.dat$quality, breaks=5, main=NULL,xlab="Quality")
dev.off()

# table info for physicochemical variables
phchem.info = matrix(nrow=11, ncol=4)
colnames(phchem.info) = c("min", "max", "mean", "median")
rownames(phchem.info) = colnames(wwine.dat)[1:11]
for (i in 1:11){
  phchem.info[i,] = c(min(wwine.dat[,i]), max(wwine.dat[,i]), mean(wwine.dat[,i]), median(wwine.dat[,i]))
}

## OLS regression
wwine.mod.ols = lm(quality ~ ., data = wwine.dat.scaled) # scaled

# find VIF 
wwine.vif = matrix(vif(wwine.mod.ols), dimnames = list(colnames(wwine.dat)[1:11], "VIF"))

## ridge and LASSO regression

#' Create ridge and LASSO regression models for wine data
#' 
#' @param dat the wine data, scaled
#' 
#' @return A list of 2 things: the ridge regression model and the LASSO regression model
ridgeandlasso.wine = function(dat){
  # seperate data into covariates and response
  x = data.matrix(dat[,1:(ncol(dat)-1)])
  y = data.matrix(dat[,ncol(dat)]); colnames(y)=c("quality")
  # cross validation
  cv.lasso = cv.glmnet(x,y,alpha=1,lambda=seq(0,.5,length.out=100))
  cv.ridge = cv.glmnet(x,y,alpha=0,lambda=seq(0,.5,length.out=100))
  # fit regressions
  lambda.lasso = cv.lasso$lambda.min
  mod.lasso = glmnet(x,y,lambda=lambda.lasso,alpha=1)
  lambda.ridge = cv.ridge$lambda.min
  mod.ridge = glmnet(x,y,lambda=lambda.ridge,alpha=0)
  return(list(ridge = mod.ridge, LASSO = mod.lasso))
}

wwine.rl = ridgeandlasso.wine(wwine.dat.scaled)
wwine.mod.ridge = wwine.rl$ridge
wwine.mod.lasso = wwine.rl$LASSO

## PCA regression

#' Create PCA regression models for wine data
#' 
#' @param dat the wine data, scaled
#' @param trunc the number of PCs retained
#' 
#' @details Also graphs the scree plot
#' 
#' @return A list of 2 things: 
#'         - the PCA regression model
#'         - vector of coefficients for x's
pcar.wine = function(dat, trunc=3){
  x = data.matrix(dat[,1:(ncol(dat)-1)])
  y = data.matrix(dat[,ncol(dat)]); colnames(y)=c("quality")
  n = nrow(x)
  
  PCA = prcomp(x,center=TRUE,scale=TRUE)
  pcs = PCA$x
  eigenvals = PCA$sdev^2
  var.expl = eigenvals/sum(eigenvals)
  
  stdev = eigenvals*sqrt(2/n)
  eigenvals.upper = eigenvals+stdev
  eigenvals.lower = eigenvals-stdev
  
  num.eig = 10
  plot(1:num.eig,eigenvals[1:num.eig],xlab="EOF",ylab="Eigenvalues",
       ylim=c(min(eigenvals.lower),max(eigenvals.upper)),pch=19,type="b",col="blue")
  points(1:num.eig,eigenvals.upper[1:num.eig],col="blue",pch="-",cex=1.3)
  points(1:num.eig,eigenvals.lower[1:num.eig],col="blue",pch="-",cex=1.3)
  pcs.trunc = pcs[,1:trunc]
  
  pcadat = as.data.frame(cbind(pcs.trunc, y))
  mod.pca = lm(quality~.+0, data=pcadat)
  beta.z = as.matrix(mod.pca$coefficients)
  pca.beta.x = as.matrix(PCA$rotation[,1:trunc]) %*% beta.z
  return(list(PCA.model = mod.pca, beta.x = pca.beta.x))
}

png("Figure 2.png", width = 4, height = 4, units = 'in', res = 500)
wwine.pca = pcar.wine(wwine.dat.scaled)
dev.off()
wwine.pca = pcar.wine(wwine.dat.scaled)
wwine.mod.pca = wwine.pca$PCA.model
wwine.pca.beta.x = wwine.pca$beta.x

## comparison of coefficients
wwine.coefs = data.frame('OLS'=round(wwine.mod.ols$coef[-1],4),
                         'LASSO'=round(wwine.mod.lasso$beta[,1],4),
                         'ridge'=round(wwine.mod.ridge$beta[,1],4),
                         'PCA'=round(wwine.pca.beta.x[,1],4))

#' Compare out-of-sample predictions of OLS, ridge, LASSO, and PCA models for wine data
#' 
#' @param dat the wine data, scaled
#' @param n the number of runs
#' @param trunc the number of PCs retained for PCA regression
#' 
#' @return A matrix of RMSE values for each model
models.oos.pred = function(dat, n=100, trunc=3){
  rmse = array(NA, c(n, 4))
  colnames(rmse) = c("OLS","ridge","LASSO","PCA")
  # run n times
  for (i in 1:n) {
    # find indexes to randomly subset data
    subset1 = sample(1:nrow(dat),size=round(nrow(dat)/2),replace=F)
    subset2 = (1:nrow(dat))[-subset1]
    # find RMSE of OLS
    mod.ols = lm(quality~.+0, data=dat[subset1,])
    ols.pred = predict(mod.ols,newdata=dat[subset2,])
    ols.rmse = sqrt(mean((ols.pred-dat$quality[subset2])^2))
    # find RMSE of ridge and LASSO
    rl = ridgeandlasso.wine(dat[subset1,])
    mod.ridge = rl$ridge
    mod.lasso = rl$LASSO
    new.x = data.matrix(dat[subset2,-12])
    ridge.pred = predict(mod.ridge,newx=new.x)
    lasso.pred = predict(mod.lasso,newx=new.x)
    ridge.rmse = sqrt(mean((ridge.pred-dat$quality[subset2])^2))
    lasso.rmse = sqrt(mean((lasso.pred-dat$quality[subset2])^2))
    # find RMSE of PCA
    beta.x.pca = pcar.wine(dat[subset1,], trunc)$beta.x
    pca.pred = data.matrix(dat[subset2,])[,-length(dat)]%*%beta.x.pca
    pca.rmse = sqrt(mean((pca.pred-dat$quality[subset2])^2))
    # combine RMSEs 
    rmse[i,] <- c(ols.rmse,ridge.rmse,lasso.rmse,pca.rmse)
  }
  return(rmse)
}

# boxplot of RMSEs
wwine.rmse = models.oos.pred(wwine.dat.scaled, n=100) 
png("Figure 3.png", width = 4, height = 4, units = 'in', res = 500)
boxplot(wwine.rmse,ylab="RMSE")
dev.off()
boxplot(wwine.rmse,ylab="RMSE")

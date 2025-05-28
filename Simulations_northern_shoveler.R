### Load the libraries
library(glmnet)
library(rtrim)
library(missMDA)
library(glmmTMB)
library(lme4)
library(lori)
library(FactoMineR)
library(softImpute)
library(mice)
library(missForest)
library(softImpute)
library(parallel)
library(ggplot2)
source("~/Dropbox/waterbirds/review/code_review/cv_sft.R")


### Load the data
y <-
  read.csv2("~/Dropbox/waterbirds/souchet_reduit/souchet.csv",
            sep = ";")

sites <- y$site
cov_sites <-
  read.csv2("~/Dropbox/waterbirds/souchet_reduit/cov_sites.csv",
            sep = ";")

y <- y[,2:29]
y <- y[,12:28]
d <- dim(y)
n <- d[1]
p <- d[2]
cov_sites <- cov_sites[, 2:11]
cov_sites$ecosystem.1 <- 1*(cov_sites$ecosystem==1)
cov_sites$ecosystem.2 <- 1*(cov_sites$ecosystem==2)
cov_sites$ecosystem.3 <- 1*(cov_sites$ecosystem==3)
cov_sites <- cov_sites[, -9]
cov_years <-
  read.csv2("~/Dropbox/waterbirds/souchet_reduit/cov_years.csv",
            sep = ";")
cov_years <- cov_years[12:28, 2:8]
cov_sites_years <- read.csv2("~/Dropbox/waterbirds/souchet_reduit/cov_sites_years.csv",
                             sep = ";")
cov_sites_years <- cov_sites_years[which(cov_sites_years$year>2000), 3:5]

cov <-covmat(n,p,cov_sites, cov_years, cov_sites_years, center = F)
covsc <- cov
covsc[, c(3:8,13:22)] <- scale(cov[, c(3:8,13:22)])
covsc$algeria.dist_town <- covsc$algeria*covsc$dist_town
covsc$morocco.dist_town <- covsc$morocco*covsc$dist_town
cov_cat <- kmeans(covsc, centers=2)$cluster
row <- as.factor(rep(1:n, p))
cov_mice <- aggregate(covsc, by=list(row), FUN=mean)[,2:6] # covariates for the mice function (take means by row to add to predictor matrix)

parallel <- TRUE # change to FALSE to compute on a single core
if (parallel) n_cores <- detectCores() else n_cores <- 1


### Functions for the two missing data patterns
na_func_rand <- function(x, prob=0.1){
  x <- as.matrix(x)
  yp <- x
  d <- dim(x)
  n <- d[1]
  p <- d[2]
  idx <- which(c(as.matrix(x))<mean(c(as.matrix(x)), na.rm=T))
  n_na <- round(prob*sum(!is.na(x)))
  yp[sample(idx, n_na)] <- NA
  return(yp)
}

na_func_pattern <- function(x, prob=0.1){
  x <- as.matrix(x)
  yp <- x
  d <- dim(x)
  n <- d[1]
  p <- d[2]
  n_na <- round(prob*sum(!is.na(x)))
  n_na_sites <- round(n_na/n)
  n_na_years <- round(n_na/p)
  idx_sites <- sample(which(rowSums(!is.na(x))>n_na_sites+1),n/2)
  idx_years <- sample(which(colSums(!is.na(x))>n_na_years+1),p/2)
  for(i in idx_sites){
    idx <- which(x[i,]<mean(c(as.matrix(x)), na.rm=T))
    if(length(idx)<n_na_sites){
      idx <- which(!is.na(x[i,]))
    }
    yp[i, sample(idx, n_na_sites)] <- NA
  }
  for(j in idx_years){
    idx <- which(x[,j]<mean(c(as.matrix(x)), na.rm=T))
    if(length(idx)<n_na_years){
      idx <- which(!is.na(x[,j]))
    }
    yp[sample(idx, n_na_years), j] <- NA
  }
  return(yp)
}

### Simus -- 100*prob % of missing values
N <- 100 # number of replications
probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6) # proportion of missing values (change to obtain the different figures)
for(prob in probs){
  ylist_rand <- lapply(1:N, function(k)
    na_func_rand(as.matrix(y), prob = prob))
  ylist_pattern <- lapply(1:N, function(k)
    na_func_pattern(as.matrix(y), prob = prob))
  omega <- !is.na(y)


  #--LORI--#
  res.cv_rand <- cv.lori(ylist_rand[[1]],as.matrix(covsc), trace.it = T, parallel=parallel) # cross-validation to choose lori params
  res.cv_pattern <- cv.lori(ylist_pattern[[1]],as.matrix(covsc), trace.it = T, parallel=parallel) # cross-validation to choose lori params

  time_lori <- Sys.time()
  error_lori_rand <- sapply(1:N, function(i) {
    res_rand <- lori(ylist_rand[[i]], covsc, lambda1 = res.cv_rand$lambda1,
                     lambda2 = res.cv_rand$lambda2)
    norm(y[omega] - res_rand$imputed[omega], type="2")^2
  })
  error_lori_pattern <- sapply(1:N, function(i) {
    res_pattern <- lori(ylist_pattern[[i]], covsc, lambda1 = res.cv_pattern$lambda1,
                        lambda2 = res.cv_pattern$lambda2)
    norm(y[omega] - res_pattern$imputed[omega], type="2")^2
  })
  time_lori <- Sys.time() - time_lori

  #--MICE--#
  time_mice <- Sys.time()
  error_mice_rand <- sapply(1:N, function(k){
    mat <- cbind(ylist_rand[[k]], cov_mice)
    colnames(mat) <- 1:(p+ncol(cov_mice))
    mat <- as.data.frame(mat)
    colnames(mat) <- paste("V", 1:(p+ncol(cov_mice)), sep="")
    ini <- mice(mat, maxit=0)
    meth <- ini$meth
    pred <- ini$pred
    pred[(p+1):(p+ncol(cov_mice)), ] <- 0
    imp1 <- mice(mat)
    res_rand <- mice::complete(imp1)
    norm(y[omega] - res_rand[,1:p][omega], type="2")^2
  })
  error_mice_pattern <- sapply(1:N, function(k){
    mat <- cbind(ylist_pattern[[k]], cov_mice)
    colnames(mat) <- 1:(p+ncol(cov_mice))
    mat <- as.data.frame(mat)
    colnames(mat) <- paste("V", 1:(p+ncol(cov_mice)), sep="")
    ini <- mice(mat, maxit=0)
    meth <- ini$meth
    pred <- ini$pred
    pred[(p+1):(p+ncol(cov_mice)), ] <- 0
    imp1 <- mice(mat)
    res_pattern <- mice::complete(imp1)
    norm(y[omega] - res_pattern[,1:p][omega], type="2")^2
  })
  time_mice <- Sys.time() - time_mice

  #--MISSFOREST--#
  time_missForest <- Sys.time()
  error_missForest_rand <- sapply(1:N, function(k){
    res_rand <- missForest(cbind(ylist_rand[[k]], cov_mice))$ximp[,1:p]
    norm(y[omega] - res_rand[omega], type="2")^2
  })
  error_missForest_pattern <- sapply(1:N, function(k) {
    res_pattern <- missForest(cbind(ylist_pattern[[k]], cov_mice))$ximp[,1:p]
    norm(y[omega] - res_pattern[omega], type="2")^2
  })
  time_missForest <- Sys.time() - time_missForest

  #--TRIM--#

  time_trim <- Sys.time()
  error_trim_rand <- sapply(1:N, function(k) {
    dat_rand <-
      data.frame(
        site = rep(1:n, p),
        time = rep(1:p, each = n),
        count = unlist(c(ylist_rand[[k]])),
        cov = cov_cat
      )
    res_rand <- tryCatch({trim(count ~ site + time, data = dat_rand, model = 2, overdisp=T)},
                         error=function(cond) {
                           return(NA)})
    if(is.null(res_rand)) re <- NA
    else{
      rowidx <- which(!(1:nrow(y))%in%res_rand$site_id)
      omega_trim <- matrix(TRUE, nrow(y), ncol(y))
      omega_trim[rowidx,] <- FALSE
      tt <- matrix(0, nrow(y), ncol(y))
      tt[setdiff((1:nrow(y)), rowidx),] <- res_rand$imputed
      re <- norm(y[(omega)*(omega_trim)==1]-tt[omega*(omega_trim)==1], type="2")^2
    }
    re
  })
  error_trim_pattern <- sapply(1:N, function(k){
    dat_pattern <-
      data.frame(
        site = rep(1:n, p),
        time = rep(1:p, each = n),
        count = unlist(c(ylist_pattern[[k]])),
        cov =cov_cat
      )
    res_pattern <- tryCatch({trim(count ~ site + time, data = dat_pattern, model = 2, overdisp=T)},
                            error=function(cond) {
                              return(NULL)})
    if(is.null(res_pattern)) re <- NA
    else{
      rowidx <- which(!(1:nrow(y))%in%res_pattern$site_id)
      omega_trim <- matrix(TRUE, nrow(y), ncol(y))
      omega_trim[rowidx,] <- FALSE
      tt <- matrix(0, nrow(y), ncol(y))
      tt[setdiff((1:nrow(y)), rowidx),] <- res_pattern$imputed
      re <- (norm(y[(omega)*(omega_trim)==1]-tt[omega*(omega_trim)==1], type="2")^2)
    }
    re
  })
  time_trim <- Sys.time()-time_trim

  #--CA--#
  time_ca <- Sys.time()
  error_ca_rand <- sapply(1:N, function(k){
    res_rand <- imputeCA(ylist_rand[[k]], ncp=2)
    norm(y[omega]-res_rand[omega], type="2")^2
  })
  error_ca_pattern <- sapply(1:N, function(k){
    res_pattern <- imputeCA(ylist_pattern[[k]], ncp=2)
    norm(y[omega]-res_pattern[omega], type="2")^2
  })
  time_ca <- Sys.time()-time_ca

  #--softImpute--#
  lambda_sft_rand <- cv_sft(ylist_rand[[1]])$lambda
  lambda_sft_pattern <- cv_sft(ylist_pattern[[1]])$lambda

  time_softImpute <- Sys.time()
  error_softImpute_rand <- sapply(1:N, function(k){
    udv <- softImpute(ylist_rand[[k]], lambda=lambda_sft_rand)
    res_rand <- udv$u%*%diag(udv$d)%*%t(udv$v)
    norm(y[omega]-res_rand[omega], type="2")^2
  })
  error_softImpute_pattern <- sapply(1:N, function(k){
    udv <- softImpute(ylist_pattern[[k]], lambda=lambda_sft_rand)
    res_pattern <- udv$u%*%diag(udv$d)%*%t(udv$v)
    norm(y[omega]-res_pattern[omega], type="2")^2
  })
  time_softImpute <- Sys.time()-time_softImpute


  #--GLM--#
  dat_list_rand <- lapply(ylist_rand, function(yy) {
    dat <- data.frame(site = rep(1:n, p), time = rep(1:p, each=n),
                      count = c(yy))
    dat <- cbind.data.frame(dat, covsc)
    return(dat)
  })
  dat_train_rand <- lapply(1:N, function(k) dat_list_rand[[k]][!is.na(ylist_rand[[k]]),])
  dat_test_rand <- lapply(1:N, function(k) dat_list_rand[[k]][is.na(ylist_rand[[k]])*(!is.na(y))==1,])

  dat_list_pattern <- lapply(ylist_pattern, function(yy) {
    dat <- data.frame(site = rep(1:n, p), time = rep(1:p, each=n),
                      count = c(yy))
    dat <- cbind.data.frame(dat, covsc)
    return(dat)
  })
  dat_train_pattern <- lapply(1:N, function(k) dat_list_pattern[[k]][!is.na(ylist_pattern[[k]]),])
  dat_test_pattern <- lapply(1:N, function(k) dat_list_pattern[[k]][is.na(ylist_pattern[[k]])*(!is.na(y))==1,])
  time_glm <- Sys.time()
  error_glmm_rand <- sapply(1:N, function(k){
    d_rand <- dat_train_rand[[k]]
    res_rand <- glm(count~1+latitude+longitude+dist_towns+dam+economy, family="poisson", data = d_rand)
    impute_glmm_rand <- predict(res_rand, dat_test_rand[[k]], type="response")
    norm(y[is.na(ylist_rand[[k]])*(!is.na(y))==1]-impute_glmm_rand, type="2")^2
  })
  error_glmm_pattern <- sapply(1:N, function(k){
    d_pattern <- dat_train_pattern[[k]]
    res_pattern <- glm(count~1+latitude+longitude+dist_towns+dam+economy, family="poisson", data = d_pattern)
    impute_glmm_pattern <- predict(res_pattern, dat_test_pattern[[k]], type="response")
    norm(y[is.na(ylist_pattern[[k]])*(!is.na(y))==1]-impute_glmm_pattern, type="2")^2
  })
  time_glm <- Sys.time()-time_glm


  #--MEAN--#
  time_mean <- Sys.time()
  impute_moy <- lapply(1:N, function(k) {
    columnsMeans <- t(matrix(rep(colMeans(ylist_rand[[k]], na.rm = T), n), nrow = p))
    imp_rand <- ylist_rand[[k]]
    imp_rand[is.na(imp_rand)] <- columnsMeans[is.na(imp_rand)]
    imp_rand[is.na(imp_rand)] <- mean(imp_rand, na.rm=T)
    columnsMeans <- t(matrix(rep(colMeans(ylist_pattern[[k]], na.rm = T), n), nrow = p))
    imp_pattern <- ylist_pattern[[k]]
    imp_pattern[is.na(imp_pattern)] <- columnsMeans[is.na(imp_pattern)]
    imp_pattern[is.na(imp_pattern)] <- mean(imp_pattern, na.rm=T)
    return(list(res_rand=imp_rand, res_pattern=imp_pattern))
  })
  time_mean <- Sys.time()-time_mean
  error_moy_rand <-
    sapply(1:N, function(k)
      norm(y[omega]-impute_moy[[k]]$res_rand[omega], type="2")^2)

  error_moy_pattern <-
    sapply(1:N, function(k)
      norm(y[omega]-impute_moy[[k]]$res_pattern[omega], type="2")^2)

  errors <- data.frame(error=c(error_moy_rand, error_moy_pattern,
                               error_ca_rand, error_ca_pattern,
                               error_trim_rand, error_trim_pattern,
                               error_glmm_rand, error_glmm_pattern,
                               error_lori_rand, error_lori_pattern,
                               error_mice_rand, error_mice_pattern,
                               error_missForest_rand, error_missForest_pattern, error_softImpute_rand, error_softImpute_pattern),
                       method=rep(c("MEAN", "CA", "TRIM", "GLM", "LORI", "MICE", "MISSFOREST", "SOFTIMPUTE"),
                                  each = 2*N),
                       bootstrap=rep(rep(c("random", "pattern"), each=N), 8))


  errors$method <- factor(errors$method, levels = c("LORI", "MICE", "MISSFOREST", "SOFTIMPUTE", "TRIM", "CA", "MEAN", "GLM"))
  dat <- errors
  dat$error <- dat$error/(norm(y[!is.na(y)], type="2")^2)
  dat$error <- sqrt(dat$error)
  ## change file name to your path
  save(errors, file=paste("~/Dropbox/waterbirds/results_review/errors_souchet", 100*prob, ".Rdata", sep=""))
  dat <- errors
  dat$error <- dat$error/(norm(y[!is.na(y)], type="2")^2)
  dat$error <- sqrt(dat$error)
  dat <- dat[which(dat$method!="SOFTIMPUTE"),]
  ggplot(data = dat, aes(x=method, y=error)) + geom_violin(outlier.shape = NA)+
    scale_x_discrete(name="") +
    scale_y_continuous( name="Relative RMSE")+theme_bw()+facet_wrap(~bootstrap, scales = "free_y")+
    #coord_cartesian(ylim=c(0, 0.45))+#!the ylim values change with prob!#
    theme(axis.text=element_text(size=20, angle=45, vjust = 1, hjust = 1),
          axis.title=element_text(size=25,face="bold"),
          strip.text = element_text(size=25))+
    stat_summary(fun.y=median, geom="point", size=0.5, color="red")+
    stat_summary(fun.y=mean, geom="point", size=0.5, color="black",shape=23)
  times <- c(time_lori, time_mice, time_missForest, time_trim, time_ca, time_softImpute, time_glm, time_mean)/N
  names(times) <- c("LORI", "MICE", "MISSFOREST", "TRIM", "CA", "SOFTIMPUTE", "GLM", "MEAN")
  ## change file name to your path
  save(times, file=paste("~/Dropbox/waterbirds/results_review/times_souchet", 100*prob, ".Rdata", sep=""))
}

rm(list=ls())
# packages ---------------------
library(survival)
library(mice)
library(smcfcs)
library(mitools)
# setup -----------------------
setup = "nocorr"
filepath=paste0("simulations/",setup,"/")
nimp=10 
n.it = 100
npara=3
nsim=1000
nmethods = 8
# Matrix to store the results
res_mat = matrix(NA,nrow=nsim,ncol=2*npara*nmethods)

# interaction formulas 
formula = "Surv(t,d)~x+z1+z2"
sm_formula = "Surv(t,d)~x+z1+z2"
sm_formula_caco = "Surv(entertime,t,d)~x+z1+z2"

predictors_aprx = c("z1","z2","d","chaz")
predictors_rs = c("z1","z2")
#-------------------------------
#===============================
# Run analyses 
#===============================
j=1
set.seed(1001)
for(j in seq(1,nsim)){
  # print(j)
  #===============================
  # cox analysis using full cohort data
  #===============================
  # load data set
  cohort = readRDS(file=paste0(filepath,"cohort",j,".rds")) # loads cohort.ncc
  model=coxph(as.formula(formula),data=cohort)
  print(model)
  res_full = c(model$coefficients,sqrt(diag(model$var)))
  
  # size of full cohort
  n = dim(cohort)[1]
  
  #===============================
  # traditional case-control analysis 
  #===============================
  
  # load data set
  caco= readRDS(file=paste0(filepath,"caco",j,".rds")) 
  
  # fit the model
  model=cch(as.formula(formula), 
            data=caco, subcoh=~subco, id=~id, method="LinYing", cohort.size=n)
  print(model)
  res_caco = c(model$coefficients,sqrt(diag(model$var)))
  
  #===============================
  #MI-approx:  full-cohort approach
  #===============================
  cohort.caco= readRDS(file=paste0(filepath,"cohort.caco",j,".rds")) #this is 'cohort.ncc'
  
  # Compute Nelson-Aalen estimate of the cumulative hazard
  cohort.caco$chaz=nelsonaalen(cohort.caco,t,d)
  
  #predictor matrix which determines the imputation models for x1
  pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
  colnames(pred.mat)=names(cohort.caco)
  rownames(pred.mat)=names(cohort.caco)
  pred.mat["x",predictors_aprx]=1
  
  #method of imputation for x1 
  method.vec=rep("",dim(cohort.caco)[2])
  method.vec[which(colnames(cohort.caco)=="x")]="norm"
  
  #perform the imputation 
  imp<-mice(cohort.caco, m = nimp, method = method.vec, 
            predictorMatrix = pred.mat, 
            maxit = n.it, diagnostics = FALSE, printFlag = F)
  
  
  # asess convergence
  #densityplot(imp, ~x)
  #stripplot(imp, x~.imp, pch=20, cex=2)
  
  # plot(density(complete(imp)$x[is.na(cohort.ncc$x)]),col="darkred")
  # lines(density(complete(imp)$x)$x,density(complete(imp)$x)$y,col="red")
  # lines(density(cohort$x)$x,density(cohort$x)$y)
  # lines(density(cohort.ncc$x[!is.na(cohort.ncc$x)])$x,
  #       density(cohort.ncc$x[!is.na(cohort.ncc$x)])$y,col="blue")
  # lines(density(cohort$x[is.na(cohort.ncc$x)])$x,
  #              density(cohort$x[is.na(cohort.ncc$x)])$y,col="green")
  
  # Fit the analysis model in each imputed data set
  models<-with(imp,coxph(as.formula(formula)))
  
  # Combine estimates across the imputed data sets using Rubin's Rules
  summary_aprx = summary(pool(models))
  res_aprx = c(summary_aprx[,"estimate"],summary_aprx[,"std.error"])
  
  
  #===============================
  #MI-SMC: full-cohort approach
  #===============================
  cohort.caco= readRDS(file=paste0(filepath,"cohort.caco",j,".rds")) #this is 'cohort.ncc'
  
  #predictor matrix which determines the imputation models for x1 (the outcomes are not included - see smcfcs help)
  pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
  colnames(pred.mat)=names(cohort.caco)
  rownames(pred.mat)=names(cohort.caco)
  pred.mat["x",predictors_rs]=1
  
  #method of imputation for x1
  method.vec=rep("",dim(cohort.caco)[2])
  method.vec[which(colnames(cohort.caco)=="x")]="norm"
  
  #perform the imputation
  imp <- smcfcs(cohort.caco, smtype="coxph", smformula=sm_formula,method=method.vec,
                predictorMatrix=pred.mat,m = nimp, numit = n.it, rjlimit = 10000,noisy=F)
  
  # convergence
  plot(imp$smCoefIter[1,1,]) #imputation 1, fourth parameter against iteration number 
  
  # # imputed values
  # plot(density(imp$impDatasets[[1]]$x[is.na(cohort.ncc$x)]),col="darkred")
  # lines(density(imp$impDatasets[[1]]$x)$x,density(imp$impDatasets[[1]]$x)$y,col="red")
  # lines(density(cohort$x)$x,density(cohort$x)$y)
  # lines(density(cohort.ncc$x[!is.na(cohort.ncc$x)])$x,
  #       density(cohort.ncc$x[!is.na(cohort.ncc$x)])$y,col="blue")
  # lines(density(cohort$x[is.na(cohort.ncc$x)])$x,
  #       density(cohort$x[is.na(cohort.ncc$x)])$y,col="green")
  # 
  
  # obtain estimates across the imputed data sets and combine using Rubin's Rules
  impobj <- imputationList(imp$impDatasets)
  models <- with(impobj, coxph(as.formula(formula)))
  coef = MIcombine(models)$coefficients
  se = sqrt(diag(MIcombine(models)$variance))
  res_rej = c(coef,se)
  
  #===============================
  #MI-approx:  superset ncc
  #===============================
  cohort.caco = readRDS(file=paste0(filepath,"cohort.super.caco",j,".rds"))
  
  # Compute Nelson-Aalen estimate of the cumulative hazard for full cohort
  cohort$chaz=nelsonaalen(cohort,t,d)
  #add cumulative hazard into ncc data
  cohort.merge<-cohort[,c("id","chaz")] 
  cohort.caco<-merge(cohort.caco,cohort.merge,by.x="id")
  
  #predictor matrix which determines the imputation models for x1
  pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
  colnames(pred.mat)=names(cohort.caco)
  rownames(pred.mat)=names(cohort.caco)
  pred.mat["x",predictors_aprx]=1
  
  #method of imputation for x1 
  method.vec=rep("",dim(cohort.caco)[2])
  method.vec[which(colnames(cohort.caco)=="x")]="norm"
  
  #perform the imputation 
  imp<-mice(cohort.caco, m = nimp, method = method.vec, 
            predictorMatrix = pred.mat, 
            maxit = n.it, diagnostics = FALSE, printFlag = F)
  
  # Fit the analysis model in each imputed data set
  models <- vector("list", nimp)
  for (k in 1:nimp){
    model=cch(as.formula(formula),data=complete(imp,k),subcoh=~subco.super, id=~id, method="LinYing", cohort.size=n)
    models[[k]] = model
  }
  
  # Combine estimates across the imputed data sets using Rubin's Rules
  res_aprx_sup = c(MIcombine(models)$coef,sqrt(diag(MIcombine(models)$variance)))
  
  #===============================
  #MI-approx:  NAIVE superset caco
  #===============================
  cohort.caco = readRDS(file=paste0(filepath,"cohort.super.caco",j,".rds"))
  
  # Compute Nelson-Aalen estimate of the cumulative hazard for superset
  cohort.caco$chaz = nelsonaalen(cohort.caco,t,d)
  
  #predictor matrix which determines the imputation models for x1
  pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
  colnames(pred.mat)=names(cohort.caco)
  rownames(pred.mat)=names(cohort.caco)
  pred.mat["x",predictors_aprx]=1
  
  #method of imputation for x1 
  method.vec=rep("",dim(cohort.caco)[2])
  method.vec[which(colnames(cohort.caco)=="x")]="norm"
  
  #perform the imputation 
  imp<-mice(cohort.caco, m = nimp, method = method.vec, 
            predictorMatrix = pred.mat, 
            maxit = n.it, diagnostics = FALSE, printFlag = F)
  
  # Fit the analysis model in each imputed data set
  models <- vector("list", nimp)
  for (k in 1:nimp){
    model=cch(as.formula(formula),data=complete(imp,k),subcoh=~subco.super, id=~id, method="LinYing", cohort.size=n)
    models[[k]] = model
  }
  
  # Combine estimates across the imputed data sets using Rubin's Rules
  res_aprx_sup_nai = c(MIcombine(models)$coef,sqrt(diag(MIcombine(models)$variance)))
  
  #===============================
  #MI-SMC: superset ncc
  #===============================
  cohort.caco= readRDS(file=paste0(filepath,"cohort.super.caco",j,".rds"))
  
  #predictor matrix which determines the imputation models for x1 (the outcomes are not included - see smcfcs help)
  pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
  colnames(pred.mat)=names(cohort.caco)
  rownames(pred.mat)=names(cohort.caco)
  pred.mat["x",predictors_rs]=1
  
  #method of imputation for x1
  method.vec=rep("",dim(cohort.caco)[2])
  method.vec[which(colnames(cohort.caco)=="x")]="norm"
  
  my.sampfrac = (cohort.caco$subco.super==1)/n
  
  #perform the imputation
  imp <- smcfcs.casecohort(cohort.caco,smformula=sm_formula_caco,sampfrac=my.sampfrac,in.subco="subco.super",
                           method=method.vec,predictorMatrix=pred.mat,m=nimp,numit=100,rjlimit=10000,noisy=FALSE)
  
  # obtain estimates across the imputed data sets and combine using Rubin's Rules
  impobj <- imputationList(imp$impDatasets)
  models <- with(impobj, coxph(as.formula(paste0(sm_formula_caco,"+cluster(id)"))))
  coef = MIcombine(models)$coefficients
  se = sqrt(diag(MIcombine(models)$variance))
  res_rej_sup = c(coef,se)
  
  #===============================
  #MI-SMC: NAIVE superset ncc
  #===============================
  cohort.caco= readRDS(file=paste0(filepath,"cohort.super.caco",j,".rds")) #this is 'cohort.ncc'
  
  #predictor matrix which determines the imputation models for x1 (the outcomes are not included - see smcfcs help)
  pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
  colnames(pred.mat)=names(cohort.caco)
  rownames(pred.mat)=names(cohort.caco)
  pred.mat["x",predictors_rs]=1
  
  #method of imputation for x1
  method.vec=rep("",dim(cohort.caco)[2])
  method.vec[which(colnames(cohort.caco)=="x")]="norm"
  
  #perform the imputation
  imp <- smcfcs(cohort.caco, smtype="coxph", smformula=sm_formula,method=method.vec,
                predictorMatrix=pred.mat,m = nimp, numit = n.it, rjlimit = 10000,noisy=F)
  
  # Fit the analysis model in each imputed data set
  models <- vector("list", nimp)
  for (k in 1:nimp){
    model=cch(as.formula(formula),data=imp$impDatasets[[k]],subcoh=~subco.super, id=~id, method="LinYing", cohort.size=n)
    models[[k]] = model
  }
  
  # Combine estimates across the imputed data sets using Rubin's Rules
  res_rej_sup_nai = c(MIcombine(models)$coef,sqrt(diag(MIcombine(models)$variance)))
  
  #--------------------------------
  #================================
  # Add results from simulation j
  res_mat[j,] = c(res_full,res_caco,res_aprx,res_aprx_sup,res_aprx_sup_nai,
                  res_rej,res_rej_sup,res_rej_sup_nai)
}

#======== END FOR LOOP ==========
#--------------------------------

#===============================
# Performance measurements
#===============================
true_parameters = matrix(c(1,1,0.5),nrow=npara)

#rownames(true_parameters) = c("beta_x1","beta_x2","beta_z")
#colnames(true_parameters) = "true"

par_idx = c()
se_idx = c()
for(k in seq(0,nmethods-1)){
  par_idx= c(par_idx,seq(1,npara)+2*npara*k)
  se_idx = c(se_idx,seq(1+npara,2*npara)+2*npara*k)
}
par_idx
se_idx
parameters = res_mat[,par_idx]
parameter_se = res_mat[,se_idx]

# Bias
para_mean = apply(parameters,2,mean)
mean_mat = matrix(para_mean,nrow=nsim,ncol=npara*nmethods,byrow=T)
apply(parameters,2,mean)-rep(true_parameters,nmethods)
bias = apply(parameters,2,mean)-rep(true_parameters,nmethods)
bias_mat = matrix(bias,nrow=npara,ncol=nmethods,byrow=F)

# Model SE
model_se = apply(parameter_se,2,mean)
model_se_mat = matrix(model_se,nrow=npara,ncol=nmethods,byrow=F)


# Empirical se
emp_se = sqrt((1/(nsim-1))*apply((parameters-mean_mat)^2,2,sum))
emp_se_mat = matrix(emp_se,nrow=npara,ncol=nmethods,byrow=F)

#MSE
mse = apply((parameters-matrix(true_parameters,ncol=nmethods*npara,nrow=nsim,byrow=T))^2,2,mean)
mse_mat = matrix(mse,nrow=npara,ncol=nmethods,byrow=F)


# 95 percent coverage
lower = parameters -1.96*parameter_se
higher = parameters +1.96*parameter_se
in_int = matrix(true_parameters,ncol=nmethods*npara,nrow=nsim,byrow = T) >= lower & matrix(true_parameters,ncol=nmethods*npara,nrow=nsim,byrow = T) <= higher
cover_mat = matrix(apply(in_int,2,mean),nrow=npara,ncol=nmethods,byrow=F)

# relative efficience (compared to full cohort)
rel_eff= apply(matrix(parameter_se[,seq(1,npara)],nrow=nsim,ncol=nmethods*npara)^2/parameter_se^2,
               2,mean)
rel_eff = matrix(rel_eff,nrow=npara,ncol=nmethods,byrow=F)

### table of results
tab_res = rbind(bias_mat,model_se_mat,emp_se_mat,rel_eff,mse_mat,cover_mat)
tab_res2 = cbind(c("Bias","","","ModelSE","","","EmpSE","","","RelEff","","","MSE","","","Cov","",""),
                 matrix(rep(c(" $\\beta_x$"," $\\beta_{z_1}$"," $\\beta_{z_2}$"),6)),
                 round(tab_res,3))
tab_res2
print(tab_res2)
library(xtable)
xx = xtable(tab_res2)
print.xtable(xx,
             sanitize.text.function=function(x){x},
             include.rownames=F,include.colnames=F)

#=== SAVE raw results
saveRDS(res_mat,file=paste0("results/cch_",setup,nsim,".rds"))
#=====================================================================================================

# Analyse ncc data
# Code based on https://github.com/ruthkeogh/MI-CC by Ruth H.Keogh
# packages ---------------------
library(survival)
library(mice)
library(smcfcs)
library(mitools)
# setup -----------------------
setup=  "standard"
filepath=paste0("simulations/",setup,"/")
nimp= 10 
n.it = 100
npara=3
nsim= 1000
nmethods = 8
res_mat = matrix(NA,nrow=nsim,ncol=2*npara*nmethods)

# standard formulas
formula_full = "Surv(t,d)~x+z1+z2"
formula_ncc = "Surv(t,case)~x+z1+z2+strata(setno)"
sm_formula_full = "Surv(t,d)~x+z1+z2"
sm_formula_ncc = "Surv(t,case)~x+z1+z2+strata(setno)"
predictors_aprx = c("z1","z2","d","chaz")
predictors_rs = c("z1","z2")

#-------------------------------
#===============================
# Run analyses 
#===============================
set.seed(1001)
for(j in seq(1,nsim)){
  #===============================
  # cox analysis using full cohort data
  #===============================
  # # load data set
  cohort = readRDS(file=paste0(filepath,"cohort",j,".rds")) 
  model=coxph(as.formula(formula_full),data=cohort)
  res_full = c(model$coefficients,sqrt(diag(model$var)))
  
  #===============================
  # traditional analysis using nested case-control sample
  #===============================
  # load data set
  ncc = readRDS(file=paste0(filepath,"ncc",j,".rds")) 
  
  # fit the model
  model = coxph(as.formula(formula_ncc),data=ncc)
  res_ncc = c(model$coefficients,sqrt(diag(model$var)))
  
  #===============================
  # MI-approx:  full-cohort approach
  #===============================
  cohort.ncc= readRDS(file=paste0(filepath,"cohort.ncc",j,".rds")) 
  
  # Compute Nelson-Aalen estimate of the cumulative hazard
  cohort.ncc$chaz=nelsonaalen(cohort.ncc,t,d)
  
  # predictor matrix which determines the imputation models for x1
  pred.mat=matrix(0,nrow=dim(cohort.ncc)[2],ncol=dim(cohort.ncc)[2])
  colnames(pred.mat)=names(cohort.ncc)
  rownames(pred.mat)=names(cohort.ncc)
  pred.mat["x",predictors_aprx]=1
  
  # method of imputation for x1 
  method.vec=rep("",dim(cohort.ncc)[2])
  method.vec[which(colnames(cohort.ncc)=="x")]="norm"
  
  # perform the imputation 
  imp<-mice(cohort.ncc, m = nimp, method = method.vec, 
            predictorMatrix = pred.mat, 
            maxit = n.it, diagnostics = FALSE, printFlag = F)
  
  # Fit the analysis model in each imputed data set
  models<-with(imp,coxph(as.formula(formula_full)))
  
  # Combine estimates across the imputed data sets using Rubin's Rules
  summary_aprx = summary(pool(models))
  res_aprx = c(summary_aprx[,"estimate"],summary_aprx[,"std.error"])
  
  #===============================
  #MI-SMC: full-cohort approach
  #===============================
  cohort.ncc= readRDS(file=paste0(filepath,"cohort.ncc",j,".rds")) 
  
  # predictor matrix which determines the imputation models for x1 
  pred.mat=matrix(0,nrow=dim(cohort.ncc)[2],ncol=dim(cohort.ncc)[2])
  colnames(pred.mat)=names(cohort.ncc)
  rownames(pred.mat)=names(cohort.ncc)
  pred.mat["x",predictors_rs]=1
  
  # method of imputation for x1
  method.vec=rep("",dim(cohort.ncc)[2])
  method.vec[which(colnames(cohort.ncc)=="x")]="norm"
  
  # perform the imputation
  imp <- smcfcs(cohort.ncc, smtype="coxph", smformula=sm_formula_full,
                method=method.vec,predictorMatrix=pred.mat,m = nimp,
                numit =n.it, rjlimit = 10000,noisy=F)
  
  # obtain estimates from imputed data sets and combine using Rubin's Rules
  impobj <- imputationList(imp$impDatasets)
  models <- with(impobj, coxph(as.formula(formula_full)))
  coef = MIcombine(models)$coefficients
  se = sqrt(diag(MIcombine(models)$variance))
  res_rej = c(coef,se)
  #===============================
  #MI-approx:  superset A ncc
  #===============================
  cohort.ncc= readRDS(file=paste0(filepath,"cohort.super.ncc",j,".rds"))
  
  # Nelson-Aalen estimate of the cumulative hazard for full cohort
  cohort$chaz=nelsonaalen(cohort,t,d)
  
  #add cumulative hazard into superset ncc data
  cohort.merge<-cohort[,c("id","chaz")] 
  cohort.ncc<-merge(cohort.ncc,cohort.merge,by.x="id")
  
  # predictor matrix for the imputation models for x1 (not incl. outcome)
  pred.mat=matrix(0,nrow=dim(cohort.ncc)[2],ncol=dim(cohort.ncc)[2])
  colnames(pred.mat)=names(cohort.ncc)
  rownames(pred.mat)=names(cohort.ncc)
  pred.mat["x",predictors_aprx]=1
  
  # method of imputation for x1 
  method.vec=rep("",dim(cohort.ncc)[2])
  method.vec[which(colnames(cohort.ncc)=="x")]="norm"
  
  # perform the imputation 
  imp<-mice(cohort.ncc, m = nimp, method = method.vec, 
            predictorMatrix = pred.mat, 
            maxit = n.it, diagnostics = FALSE, printFlag = F)
  
  # Fit the analysis model in each imputed data set
  models<-with(imp,coxph(as.formula(formula_ncc)))
  
  # Combine estimates across the imputed data sets using Rubin's Rules
  summary_aprx = summary(pool(models))
  res_aprx_sup = c(summary_aprx[,"estimate"],summary_aprx[,"std.error"])
  
  #===============================
  # MI-approx: (naive) Superset B ncc
  #===============================
  cohort.ncc= readRDS(file=paste0(filepath,"cohort.super.ncc",j,".rds"))
  
  # Nelson-Aalen estimate of the cumulative hazard for the superset ncc
  cohort.ncc$chaz = nelsonaalen(cohort.ncc,t,d)
  
  # predictor matrix which determines the imputation models for x1
  pred.mat=matrix(0,nrow=dim(cohort.ncc)[2],ncol=dim(cohort.ncc)[2])
  colnames(pred.mat)=names(cohort.ncc)
  rownames(pred.mat)=names(cohort.ncc)
  pred.mat["x",predictors_aprx]=1
  
  #method of imputation for x1 
  method.vec=rep("",dim(cohort.ncc)[2])
  method.vec[which(colnames(cohort.ncc)=="x")]="norm"
  
  #perform the imputation 
  imp<-mice(cohort.ncc, m = nimp, method = method.vec, 
            predictorMatrix = pred.mat, 
            maxit = n.it, diagnostics = FALSE, printFlag = F)
  
  # Fit the analysis model in each imputed data set
  models<-with(imp,coxph(as.formula(formula_ncc)))
  
  # Combine estimates across the imputed data sets using Rubin's Rules
  summary_aprx = summary(pool(models))
  res_aprx_sup_nai = c(summary_aprx[,"estimate"],summary_aprx[,"std.error"])
  
  #===============================
  #MI-SMC: superset A ncc
  #===============================
  cohort.ncc= readRDS(file=paste0(filepath,"cohort.super.ncc",j,".rds")) 
  
  # Compute number at risk at each event time using the full cohort data
  nrisk.fit<-survfit(Surv(t,d)~1,data=cohort)
  ord.t.d1<-order(cohort$t[cohort$d==1])
  # number at risk at each unique event time
  numrisk<-summary(nrisk.fit,censored=F)$n.risk 
  
  # add numbers at risk time into the nested case-control data set
  cohort.ncc$numrisk<-NA
  cohort.ncc$numrisk[cohort.ncc$case==1][ord.t.d1]<-numrisk
  
  # assign number at to every individual in each set
  cohort.ncc$numrisk<-ave(cohort.ncc$numrisk, cohort.ncc$setno, 
                          FUN = function(x) sum(x, na.rm=T))
  
  #predictor matrix which determines the imputation models for x
  pred.mat=matrix(0,nrow=dim(cohort.ncc)[2],ncol=dim(cohort.ncc)[2])
  colnames(pred.mat)=names(cohort.ncc)
  rownames(pred.mat)=names(cohort.ncc)
  pred.mat["x",predictors_rs]=1
  
  #method of imputation for x1
  method.vec=rep("",dim(cohort.ncc)[2])
  method.vec[which(colnames(cohort.ncc)=="x")]="norm"
  
  #perform the imputation
  imp<-smcfcs.nestedcc(cohort.ncc,smformula=sm_formula_ncc,
                       set="setno",event="d",nrisk="numrisk",
                       method=method.vec, predictorMatrix=pred.mat,
                       m=nimp,numit=n.it,rjlimit=1000,noisy=F) 

  # obtain estimates and combine using Rubin's Rules
  impobj <- imputationList(imp$impDatasets)
  models <- with(impobj, coxph(as.formula(formula_ncc)))
  coef = MIcombine(models)$coefficients
  se = sqrt(diag(MIcombine(models)$variance))
  res_rej_sup = c(coef,se)
  
  #===============================
  # MI-SMC: (NAIVE) superset B ncc
  #===============================
  cohort.ncc= readRDS(file=paste0(filepath,"cohort.super.ncc",j,".rds")) 
  # predictor matrix which determines the imputation models for x1 
  pred.mat=matrix(0,nrow=dim(cohort.ncc)[2],ncol=dim(cohort.ncc)[2])
  colnames(pred.mat)=names(cohort.ncc)
  rownames(pred.mat)=names(cohort.ncc)
  pred.mat["x",predictors_rs]=1
  
  # method of imputation for x1
  method.vec=rep("",dim(cohort.ncc)[2])
  method.vec[which(colnames(cohort.ncc)=="x")]="norm"
  
  # perform the imputation
  imp <- smcfcs(cohort.ncc, smtype="coxph", smformula=sm_formula_full,
                method=method.vec,predictorMatrix=pred.mat,m =nimp, 
                numit = n.it, rjlimit = 10000,noisy=F)
  
  # obtain estimates from imputed data sets and combine using Rubin's Rules
  impobj <- imputationList(imp$impDatasets)
  models <- with(impobj, coxph(as.formula(formula_ncc)))
  coef = MIcombine(models)$coefficients
  se = sqrt(diag(MIcombine(models)$variance))
  res_rej_sup_nai = c(coef,se)
  #--------------------------------
  #================================
  # Add results from simulation j
  res_mat[j,] = c(res_full,res_ncc,res_aprx,res_aprx_sup,res_aprx_sup_nai,
                  res_rej,res_rej_sup,res_rej_sup_nai)
}
#======== END FOR LOOP ==========

#===============================
# Performance measurements
#===============================
true_parameters = matrix(c(1,1,0.5),nrow=3)
par_idx = c()
se_idx = c()
for(k in seq(0,nmethods-1)){
  par_idx= c(par_idx,seq(1,npara)+2*npara*k)
  se_idx = c(se_idx,seq(1+npara,2*npara)+2*npara*k)
}
parameters = res_mat[,par_idx]
parameter_se = res_mat[,se_idx]
# Bias
para_mean = apply(parameters,2,mean)
mean_mat = matrix(para_mean,nrow=nsim,ncol=npara*nmethods,byrow=T)
bias = apply(parameters,2,mean)-rep(true_parameters,nmethods)
bias_mat = matrix(bias,nrow=npara,ncol=nmethods,byrow=F)
# Model SE
model_se = apply(parameter_se,2,mean)
model_se_mat = matrix(model_se,nrow=npara,ncol=nmethods,byrow=F)
# Empirical se
emp_se = sqrt((1/(nsim-1))*apply((parameters-mean_mat)^2,2,sum))
emp_se_mat = matrix(emp_se,nrow=npara,ncol=nmethods,byrow=F)
# MSE
mse = apply((parameters-matrix(true_parameters,ncol=nmethods*npara,
                               nrow=nsim,byrow=T))^2,2,mean)
mse_mat = matrix(mse,nrow=npara,ncol=nmethods,byrow=F)
# 95 percent coverage
lower = parameters -1.96*parameter_se
higher = parameters +1.96*parameter_se
in_int = matrix(true_parameters,ncol=nmethods*npara,nrow=nsim,byrow = T) >= 
  lower & matrix(true_parameters,ncol=nmethods*npara,nrow=nsim,byrow = T) <= 
  higher
cover_mat = matrix(apply(in_int,2,mean),nrow=npara,ncol=nmethods,byrow=F)

# relative efficience (compared to full cohort)
rel_eff= apply(matrix(parameter_se[,seq(1,npara)],nrow=nsim,
                      ncol=nmethods*npara)^2/parameter_se^2,2,mean)
rel_eff = matrix(rel_eff,nrow=npara,ncol=nmethods,byrow=F)

### table of results
tab_res = rbind(bias_mat,model_se_mat,emp_se_mat,rel_eff,mse_mat,cover_mat)
tab_res2 = cbind(c("Bias","","","ModelSE","","","EmpSE","","",
                   "RelEff","","","MSE","","","Cov","",""),
                 matrix(rep(c(" $\\beta_x$"," $\\beta_{z_1}$","
                              $\\beta_{z_2}$"),6)),
                 round(tab_res,3))
print(tab_res2)
# Monte Carlo SE of estimates:
mc_bias = sqrt((1/(nsim*(nsim-1)))*apply((parameters-matrix(true_parameters,
                                      ncol=nmethods*npara,
                                      nrow=nsim,byrow=T))^2,2,sum))
mc_bias_mat = matrix(mc_bias,nrow=npara,ncol=nmethods,byrow=F)
mc_emp_se = (1/sqrt(2*(nsim-1)))*emp_se_mat
mc_mse = sqrt(apply(((parameters-matrix(true_parameters,ncol=nmethods*npara,
                                    nrow=nsim,byrow=T))^2 
                     -matrix(mse_mat,ncol=nmethods*npara,
                             nrow=nsim,byrow=T))^2,2,sum)/(nsim*(nsim-1)))
mc_cover = sqrt((cover_mat*(1-cover_mat))/nsim)
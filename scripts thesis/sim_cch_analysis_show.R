# Analyse case-cohort data
# Code based on https://github.com/ruthkeogh/MI-CC by Ruth H.Keogh
# packages ---------------------
library(survival)
library(mice)
library(smcfcs)
library(mitools)
# setup -----------------------
setup = "standard"
filepath=paste0("simulations/",setup,"/")
nimp=10 
n.it = 100
npara=3
nsim=1000
nmethods = 8
res_mat = matrix(NA,nrow=nsim,ncol=2*npara*nmethods)
# interaction formulas 
formula = "Surv(t,d)~x+z1+z2"
sm_formula = "Surv(t,d)~x+z1+z2"
sm_formula_caco = "Surv(entertime,t,d)~x+z1+z2"
predictors_aprx = c("z1","z2","d","chaz")
predictors_rs = c("z1","z2")
#===============================
# Run analyses 
#===============================
set.seed(1001)
for(j in seq(1,nsim)){
  #===============================
  # cox analysis using full cohort data
  #===============================
  # load data set
  cohort = readRDS(file=paste0(filepath,"cohort",j,".rds"))  
  model=coxph(as.formula(formula),data=cohort)
  res_full = c(model$coefficients,sqrt(diag(model$var)))
  # size of full cohort
  n = dim(cohort)[1]
  #===============================
  # traditional case-control analysis 
  #===============================
  # load data set
  caco= readRDS(file=paste0(filepath,"caco",j,".rds")) 
  # fit the model
  model=cch(as.formula(formula), data=caco, 
            subcoh=~subco, id=~id, method="LinYing", cohort.size=n)
  res_caco = c(model$coefficients,sqrt(diag(model$var)))
  #===============================
  # MI-approx:  full-cohort approach
  #===============================
  cohort.caco= readRDS(file=paste0(filepath,"cohort.caco",j,".rds"))
  # Compute Nelson-Aalen estimate of the cumulative hazard
  cohort.caco$chaz=nelsonaalen(cohort.caco,t,d)
  # predictor matrix which determines the imputation models for x1
  pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
  colnames(pred.mat)=names(cohort.caco)
  rownames(pred.mat)=names(cohort.caco)
  pred.mat["x",predictors_aprx]=1
  # method of imputation for x1 
  method.vec=rep("",dim(cohort.caco)[2])
  method.vec[which(colnames(cohort.caco)=="x")]="norm"
  # perform the imputation 
  imp<-mice(cohort.caco, m = nimp, method = method.vec, 
            predictorMatrix = pred.mat, 
            maxit = n.it, diagnostics = FALSE, printFlag = F)
  # Fit the analysis model in each imputed data set
  models<-with(imp,coxph(as.formula(formula)))
  # Combine estimates across the imputed data sets using Rubin's Rules
  summary_aprx = summary(pool(models))
  res_aprx = c(summary_aprx[,"estimate"],summary_aprx[,"std.error"])
  #===============================
  # MI-SMC: full-cohort approach
  #===============================
  cohort.caco= readRDS(file=paste0(filepath,"cohort.caco",j,".rds"))
  #predictor matrix which determines the imputation models for x1 
  pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
  colnames(pred.mat)=names(cohort.caco)
  rownames(pred.mat)=names(cohort.caco)
  pred.mat["x",predictors_rs]=1
  # method of imputation for x1
  method.vec=rep("",dim(cohort.caco)[2])
  method.vec[which(colnames(cohort.caco)=="x")]="norm"
  # perform the imputation
  imp <- smcfcs(cohort.caco, smtype="coxph", smformula=sm_formula,
                method=method.vec,predictorMatrix=pred.mat,m = nimp, 
                numit = n.it, rjlimit = 10000,noisy=F)
  # estimates from imputed data sets and combine using Rubin's Rules
  impobj <- imputationList(imp$impDatasets)
  models <- with(impobj, coxph(as.formula(formula)))
  coef = MIcombine(models)$coefficients
  se = sqrt(diag(MIcombine(models)$variance))
  res_rej = c(coef,se)
  #===============================
  # MI-approx: superset A cch
  #===============================
  cohort.caco=readRDS(file=paste0(filepath,"cohort.super.caco",j,".rds"))
  # Compute Nelson-Aalen estimate of the cumulative hazard for full cohort
  cohort$chaz=nelsonaalen(cohort,t,d)
  # add cumulative hazard into ncc data
  cohort.merge<-cohort[,c("id","chaz")] 
  cohort.caco<-merge(cohort.caco,cohort.merge,by.x="id")
  # predictor matrix which determines the imputation models for x1
  pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
  colnames(pred.mat)=names(cohort.caco)
  rownames(pred.mat)=names(cohort.caco)
  pred.mat["x",predictors_aprx]=1
  # method of imputation for x1 
  method.vec=rep("",dim(cohort.caco)[2])
  method.vec[which(colnames(cohort.caco)=="x")]="norm"
  #perform the imputation 
  imp<-mice(cohort.caco, m = nimp, method = method.vec, 
            predictorMatrix = pred.mat, 
            maxit = n.it, diagnostics = FALSE, printFlag = F)
  # Fit the analysis model in each imputed data set
  models <- vector("list", nimp)
  for (k in 1:nimp){
    model=cch(as.formula(formula),data=complete(imp,k),
              subcoh=~subco.super, id=~id, method="LinYing", cohort.size=n)
    models[[k]] = model
  }
  # Combine estimates across the imputed data sets using Rubin's Rules
  res_aprx_sup=c(MIcombine(models)$coef,
                 sqrt(diag(MIcombine(models)$variance)))
  #===============================
  #MI-approx: (naive) superset B cch
  #===============================
  cohort.caco=readRDS(file=paste0(filepath,"cohort.super.caco",j,".rds"))
  # Compute Nelson-Aalen estimate of the cumulative hazard for superset
  cohort.caco$chaz = nelsonaalen(cohort.caco,t,d)
  # predictor matrix which determines the imputation models for x1
  pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
  colnames(pred.mat)=names(cohort.caco)
  rownames(pred.mat)=names(cohort.caco)
  pred.mat["x",predictors_aprx]=1
  # method of imputation for x1 
  method.vec=rep("",dim(cohort.caco)[2])
  method.vec[which(colnames(cohort.caco)=="x")]="norm"
  #perform the imputation 
  imp<-mice(cohort.caco, m = nimp, method = method.vec, 
            predictorMatrix = pred.mat, 
            maxit = n.it, diagnostics = FALSE, printFlag = F)
  # Fit the analysis model in each imputed data set
  models <- vector("list", nimp)
  for (k in 1:nimp){
    model=cch(as.formula(formula),data=complete(imp,k),subcoh=~subco.super,
              id=~id, method="LinYing", cohort.size=n)
    models[[k]] = model
  }
  # Combine estimates across imputed data sets using Rubin's Rules
  res_aprx_sup_nai=c(MIcombine(models)$coef,
                     sqrt(diag(MIcombine(models)$variance)))
  #===============================
  #MI-SMC: superset A cch
  #===============================
  cohort.caco=readRDS(file=paste0(filepath,"cohort.super.caco",j,".rds"))
  # predictor matrix imputation models for x1 (not incl. outcomes)
  pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
  colnames(pred.mat)=names(cohort.caco)
  rownames(pred.mat)=names(cohort.caco)
  pred.mat["x",predictors_rs]=1
  # method of imputation for x1
  method.vec=rep("",dim(cohort.caco)[2])
  method.vec[which(colnames(cohort.caco)=="x")]="norm"
  # sampling fraction
  my.sampfrac = sum(cohort.caco$subco.super==1)/n 
  # perform the imputation
  imp <- smcfcs.casecohort(cohort.caco,smformula=sm_formula_caco,
                           sampfrac=my.sampfrac,in.subco="subco.super",
                           method=method.vec,predictorMatrix=pred.mat,
                           m=nimp,numit=100,rjlimit=10000,noisy=FALSE)
  # estimates from imputed data sets and combine using Rubin's Rules
  impobj <- imputationList(imp$impDatasets)
  models <- with(impobj,
                 coxph(as.formula(paste0(sm_formula_caco,"+cluster(id)"))))
  coef = MIcombine(models)$coefficients
  se = sqrt(diag(MIcombine(models)$variance))
  res_rej_sup = c(coef,se)
  #===============================
  #MI-SMC: superset B cch
  #===============================
  cohort.caco=readRDS(file=paste0(filepath,"cohort.super.caco",j,".rds")) 
  #predictor matrix which determines the imputation models for x1 
  pred.mat=matrix(0,nrow=dim(cohort.caco)[2],ncol=dim(cohort.caco)[2])
  colnames(pred.mat)=names(cohort.caco)
  rownames(pred.mat)=names(cohort.caco)
  pred.mat["x",predictors_rs]=1
  #method of imputation for x1
  method.vec=rep("",dim(cohort.caco)[2])
  method.vec[which(colnames(cohort.caco)=="x")]="norm"
  #perform the imputation
  imp <- smcfcs(cohort.caco, smtype="coxph", smformula=sm_formula,
                method=method.vec,predictorMatrix=pred.mat,m = nimp, 
                numit = n.it, rjlimit = 10000,noisy=F)
  # Fit the analysis model in each imputed data set
  models <- vector("list", nimp)
  for (k in 1:nimp){
    model=cch(as.formula(formula),data=imp$impDatasets[[k]],
              subcoh=~subco.super,id=~id, method="LinYing",cohort.size=n)
    models[[k]] = model
  }
  # Combine estimates across the imputed data sets using Rubin's Rules
  res_rej_sup_nai=c(MIcombine(models)$coef,
                    sqrt(diag(MIcombine(models)$variance)))
  #--------------------------------
  #================================
  # Add results from simulation j
  res_mat[j,] = c(res_full,res_caco,res_aprx,
                  res_aprx_sup,res_aprx_sup_nai,
                  res_rej,res_rej_sup,res_rej_sup_nai)
} #======== END FOR LOOP ==========
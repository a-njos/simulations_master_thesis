# ------------------------------
# GENERATE FULL-COHORT DATA
# ------------------------------
# Code based on https://github.com/ruthkeogh/MI-CC by Ruth H.Keogh
# ------------------------------
sim_setup = "standard"     #  "standard","aux","int","nocorr"
sim_idx = 1
filepath = paste0("simulations/",sim_setup[sim_idx],"/")
nsim = 1000
# ------------------------------
# Set parameter values, etc
n = 5000     # cohort size
close.time=15  # Maximum follow-up time

beta.x=1   # log hazard ratio for X
beta.z1=1   # log hazard ratio for  binary Z1
beta.z2=0.5     # log hazard ratio for continous Z2
beta_int = ifelse(sim_setup[sim_idx]=="int",0.5,0) # interaction X and Z1  

beta_x_drop =0    # log hazard ratio for X
beta_z1_drop =0     # log hazard ratio for Z1
beta_z2_drop  =0      # log hazard ratio for Z2

p_z1 = 0.5   # Bernoulli probability for Z1
mu_z2 = 0    # mean for normal Z2

a_x = 0      # constant for X
b_x=0.25   # influence of Z1 on X
c_x=0.25   # influence of Z2 on X

if(sim_setup[sim_idx]=="nocorr"){
  b_x=0 
  c_x=0 
}
eta_v = 0.8 #  Gaussian standard deviation for auxiliary 
lambda=NA  # Weibull baseline scale for event of interest
if(sim_setup[sim_idx]=="standard"){
  lambda=0.00000040
}else if(sim_setup[sim_idx]=="int"){
  lambda=0.00000025
}else if(sim_setup[sim_idx]=="nocorr"){
  lambda=0.00000055
}
kappa = 4  # Weibull baseline shape 
lambda_drop = 0.00002    # Weibull baseline scale dropout time
kappa_drop=4     # Weibull baseline shape for droput time
# ------------------------------
# starting data generating process
set.seed(123)
for(j in seq(1,nsim)){
  # Generate id numbers
  id=seq(1,n)
  # ------------------------------
  # Generate covariates 
  z1=rbinom(n,1,p_z1)
  z2=rnorm(n,mu_z2,1)
  x=rnorm(n,a_x +b_x*z1+c_x*z2,1)   
  v = x + rnorm(n,0,eta_v)
  # ------------------------------
  # Generate potential event times
  u=runif(n,0,1)
  t.event=(-log(u)*(1/lambda)*
             exp(-(beta.x*x+beta.z1*z1+
                     beta.z2*z2+x*z1*beta_int)))^(1/kappa)
  # ------------------------------
  # Generate potential drop-out time
  u=runif(n,0,1)
  t.drop=(-log(u)*(1/lambda_drop)*
            exp(-(beta_x_drop*x+beta_z1_drop*z1+
                    beta_z2_drop*z2)))^(1/kappa_drop)
  # ------------------------------
  # Generate time for event or drop out
  t=pmin(t.event,t.drop,close.time)
  cause=1*(t==t.event)+2*(t==t.drop)+3*(t==close.time) 
  # 1: event, 2: drop out, 3: administrative censoring
  d=ifelse(cause==1,1,0) 
  # ------------------------------
  # the full-cohort data with no missingness
  cohort=data.frame(id,t,d,x,z1,z2,v,cause)
  saveRDS(cohort,file=paste0(filepath,"cohort",j,".rds")) 
  
  # ------------------------------
  # GENERATE NESTED CASE-CONTROL DATA
  # ------------------------------
  # Generate nested-case-control superset sample
  n.controls.super = 3 # number of controls per case
  n.controls.super.ext = 7
  n.controls.ncc = 1
  ncc.super.ext = NULL
  ncc.super=NULL
  ncc=NULL
  no.sample=0
  for (i in which(cohort$d==1))
  {
    # Select control(s) for nested case-control
    possible.controls=which(cohort$t>=cohort$t[i])
    if (length(possible.controls)>=n.controls.super.ext){
      controls.super.ext=sample(possible.controls,n.controls.super.ext)
      controls.super=sample(controls.super.ext,n.controls.super)
      controls.ncc=sample(controls.super,n.controls.ncc)
      ncc.super.ext=rbind(ncc.super.ext,cohort[i,])
      ncc.super.ext=rbind(ncc.super.ext,cohort[controls.super.ext,])
      ncc.super=rbind(ncc.super,cohort[i,])
      ncc.super=rbind(ncc.super,cohort[controls.super,])
      ncc=rbind(ncc,cohort[i,])
      ncc=rbind(ncc,cohort[controls.ncc,])
      no.sample=no.sample+1}
  }
  ncc.super.ext$setno=rep(1:no.sample,each=n.controls.super.ext+1)
  ncc.super.ext$case=rep(c(1,rep(0,n.controls.super.ext)),no.sample)
  ncc.super$setno=rep(1:no.sample,each=n.controls.super+1)
  ncc.super$case=rep(c(1,rep(0,n.controls.super)),no.sample)
  ncc$setno=rep(1:no.sample,each=n.controls.ncc+1)
  ncc$case=rep(c(1,rep(0,n.controls.ncc)),no.sample)
  #-----------------------
  # generate indicator of being in the nested case-control sample
  cohort.ncc = cohort
  cohort.ncc$in.ncc <- cohort.ncc$id%in%ncc$id
  
  cohort.super = ncc.super
  cohort.super$in.ncc<-cohort.super$id%in%ncc$id
  
  cohort.super.ext = ncc.super.ext
  cohort.super.ext$in.ncc<-cohort.super.ext$id%in%ncc$id
  
  #-----------------------
  # make x missing in those outside the nested case-control sample
  cohort.ncc$x<-ifelse(cohort.ncc$id%in%ncc$id,cohort.ncc$x,NA)
  cohort.super$x<-ifelse(cohort.super$id%in%ncc$id,cohort.super$x,NA)
  cohort.super.ext$x<-ifelse(cohort.super.ext$id%in%ncc$id,
                             cohort.super.ext$x,NA)
  # ------------------------------
  # save NCC data sets
  # NCC sample
  saveRDS(ncc,file=paste0(filepath,"ncc",j,".rds")) 
  
  # NCC within full cohort
  saveRDS(cohort.ncc,file=paste0(filepath,"cohort.ncc",j,".rds")) 
  
  # NCC within superset ncc
  saveRDS(cohort.super,file=paste0(filepath,"cohort.super.ncc",j,".rds")) 
  
  # NCC within superset extended ncc
  saveRDS(cohort.super.ext,
          file=paste0(filepath,"cohort.super.ext.ncc",j,".rds"))
  
  # ------------------------------
  # GENERATE CASE-COHORT DATA
  # ------------------------------
  cohort.caco=cohort
  # ------------------------------
  # Generate subcohort
  n.subco=250
  cohort.caco$subco<-c(rep(1,n.subco),rep(0,n-n.subco))
  
  n.subco.super = 750
  cohort.caco$subco.super <- c(rep(1,n.subco.super),rep(0,n-n.subco.super))
  
  n.subco.super.ext = 1750
  cohort.caco$subco.super.ext <- c(rep(1,n.subco.super.ext),
                                   rep(0,n-n.subco.super.ext))
  # ------------------------------
  #make x1 missing in those outside the case-cohort sample
  cohort.caco$x<-ifelse(cohort.caco$subco==1|cohort.caco$d==1,
                        cohort.caco$x,NA)
  # ------------------------------
  # Generate data-set which is just the case-cohort substudy
  caco=cohort.caco[cohort.caco$subco==1|cohort.caco$d==1,]
  
  # Generate data-set which is the case-cohort supersets 
  cohort.super.caco=cohort.caco[cohort.caco$subco.super==1|
                                  cohort.caco$d==1,]
  
  cohort.super.ext.caco=cohort.caco[cohort.caco$subco.super.ext==1|
                                      cohort.caco$d==1,]
  
  cohort.super.caco$entertime=ifelse(cohort.super.caco$d==1 & 
                                     cohort.super.caco$subco.super==0,
                                     cohort.super.caco$t-0.001,0) 
  
  cohort.super.ext.caco$entertime=ifelse(cohort.super.ext.caco$d==1&
                                    cohort.super.ext.caco$subco.super.ext==0,
                                    cohort.super.ext.caco$t-0.001,0) 
  
  # ------------------------------
  # save case-cohort data sets
  saveRDS(caco,file=paste0(filepath,"caco",j,".rds"))  
  
  saveRDS(cohort.caco,file=paste0(filepath,"cohort.caco",j,".rds")) 
  
  saveRDS(cohort.super.caco,
          file=paste0(filepath,"cohort.super.caco",j,".rds")) 
  
  saveRDS(cohort.super.ext.caco,
          file=paste0(filepath,"cohort.super.ext.caco",j,".rds")) 
}
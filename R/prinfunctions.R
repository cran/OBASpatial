#setwd("C:/Users/alejo/Documents/Alejo/Doutorado/Tesis/oba package")
#source("auxfunctionst.R")
#source("auxfunctionsnorm.R")
################################################
#########densidade priori cjta (phi,nu)#########
################################################

##############
#### a data frame with the coordinates and the covariates (optional) is needed


dtsrprioroba=function(x,trend="cte",prior="reference",coords.col=1:2,kappa=0.5,cov.model="exponential",data){

####validation of arguments#####################################3

if(!is.numeric(x))  stop("x must be a numeric vector of length 2")
if(x[2]<=4) stop ("The degrees of freedom must be greater than 4")
if(x[1]<0) stop("Range parameter must be a real number in [0,Inf)")
if(length(x)!=2) stop("x must be a numeric vector of length 2")

if ( !is.numeric(coords.col)) stop("2D coordinates must be specified")
if ( length(coords.col)!=2) stop("2D coordinates must be specified")

if(!is.numeric(kappa))stop("kappa must be a real number in [0,Inf)")
if(kappa<0) stop("kappa must be a real number in [0,Inf)")

if (!(cov.model %in% c("matern","pow.exp","gaussian", "spherical","cauchy","exponential")))  stop('Valid covariance structures are (matern), (exponential), (gaussian), (spherical),
(pow.exp), (cauchy)')

if (!is.data.frame(data)) stop("data must be a data.frame")


coords=data[,coords.col]
H=as.matrix(dist(coords,upper=T,diag=T))

tau2=0
sigma=1

if(trend==0) if(prior!="jef.ind") stop("trend=0 is just valid for the jef.ind prior")
if(trend!=0) if(prior=="jef.ind") stop("The argument trend is not necessary when prior=jef.ind")

if(prior=="jef.ind"){
priorcand=jefpriortind2(x=x,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2)
return(priorcand)
}
else{

if(class(trend)!="formula") if (!(trend %in% c("cte","1st","2nd"))) stop('trend must be (cte), (1st),(2nd) or a valid formula')

if (!(prior %in% c("reference","jef.rul","jef.ind"))) stop("Accepted priors: reference, jef.rul, jef.ind")

c1=coords[,1]
c2=coords[,2]

if(trend=="cte") xmat=as.matrix(rep(1,dim(data)[1]))
if(trend=="1st") xmat=model.matrix(~c1+c2)
if(trend=="2nd") xmat= model.matrix(~ c1 + c2 + I(c1^2) + I(c2^2) + I(c1*c2))
if(class(trend)=="formula") xmat=model.matrix(trend,data=data)


#############Looking for NA#####################
if ((sum(is.na(xmat))+sum(is.na(coords)))>0) stop ("NAs not allowed")

if (det(t(xmat)%*%xmat)==0) stop("the columns of x must be linearly independent")



if(prior=="reference") priorcand=prioriref2(x=x,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
if(prior=="jef.rul") priorcand=jefpriort2(x=x,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat)
#if(prior=="jef.ind") priorcand=jefpriortind2(x=x,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2)
#if(prior=="vague") priorcand=vagprior(x=x,sigma=sigma,anu=anu,bnu=bnu,aphi=aphi,bphi=bphi)
return(priorcand)
}
}




################################################
#########densidade posteriori cjta (phi,nu)#########
################################################

##############
#### a data frame with the coordinates and the covariates (optional) is needed

dtsrposoba=function(x,formula,prior="reference",coords.col=1:2,kappa=0.5,cov.model="exponential",data,asigma=2.1,intphi="default",intnu="default"){

####validation of arguments#####################################3

if(!is.numeric(x))  stop("x must be a numeric vector of length 2")
if(x[2]<=4) stop ("The degrees of freedom must be greater than 4")
if(x[1]<0) stop("Range parameter must be a real number in [0,Inf)")
if(prior=="vague") if(length(x)!=3) stop("x must be a numeric vector of length 3")
if(prior %in% c("reference","jef.rul","jef.ind")) if(length(x)!=2) stop("x must be a numeric vector of length 2")

if (class(formula)!="formula")  stop('The argument formula must be a valid formula')

if (!(prior %in% c("reference","jef.rul","jef.ind","vague"))) stop("Accepted priors: reference, jef.rul, jef.ind")

if ( !is.numeric(coords.col)) stop("2D coordinates must be specified")
if ( length(coords.col)!=2) stop("2D coordinates must be specified")

if(!is.numeric(kappa))stop("kappa must be a real number in [0,Inf)")
if(kappa<0) stop("kappa must be a real number in [0,Inf)")

if (!(cov.model %in% c("matern","pow.exp","gaussian", "spherical","cauchy","exponential")))  stop('Valid covariance structures are (matern), (exponential), (gaussian), (spherical),
(pow.exp), (cauchy)')

if (!is.data.frame(data)) stop("data must be a data.frame")
  coords=data[,coords.col]
  H=as.matrix(dist(coords,upper=T,diag=T))

if(prior=="vague"){
  if(asigma<=2) stop("asigma must be a real number greater than 2 (proper posterior)")
  if(!is.numeric(asigma)) stop("asigma must be numeric")
if(!is.numeric(intphi)) stop("intphi must be a numeric 2 dimension vector")
if(!is.vector(intphi)) stop("intphi must be a numeric 2 dimension vector")
if ( length(intphi)!=2) stop("A valid interval for the range parameter must be provided")
if(intphi[1]>intphi[2]) stop("Invalid interval for the range parameter")
if(intphi[1]<=0 |intphi[2]<=0) stop("Range parameter is greater than 0")
  if(!is.numeric(intnu)) stop("intnu must be a numeric 2 dimension vector")
  if(!is.vector(intnu)) stop("intnu must be a numeric 2 dimension vector")
  if ( length(intnu)!=2) stop("A valid interval for the degrees of freedom must be provided")
  if(intnu[1]>intnu[2]) stop("Invalid interval for the degrees of freedom")
  if(intnu[1]<4 |intphi[2]<4) stop("degrees of reedom must be in (4,Inf)")
  anu=intnu[1]
  bnu=intnu[2]
aphi=intphi[1]
bphi=intphi[2]
  }







xmat=model.matrix(formula,data=data)
#if(prior=="jef.ind") if(max(xmat)==1 & min(xmat)==1) stop("This prior generates an improper posterior")
if(prior=="jef.ind") xmat=as.matrix(xmat[,-1])

y <-data[,all.vars(formula)[1]]

#############Looking for NA#####################
if ((sum(is.na(y))+sum(is.na(xmat))+sum(is.na(coords)))>0) stop ("NAs not allowed")

if (det(t(xmat)%*%xmat)==0) stop("the columns of x must be linearly independent")


tau2=0
sigma=1
cons=1
if(prior=="reference") posterior=posteriorref(y,cons=cons,sigma=sigma,x=x,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
if(prior=="jef.rul") posterior=posteriorjef(x=x,y=y,cons=cons,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
if(prior=="jef.ind") posterior=posteriorjefind(x=x,y=y,cons=cons,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
if(prior=="vague") posterior=posteriorvag(x=x,y=y,cons=cons,sigma=sigma,anu=anu,bnu=bnu,aphi=aphi,bphi=bphi,asigma=asigma,xmat=xmat,cov.model=cov.model,H=H,tau2=tau2,kappa=kappa)

#if(prior=="jef.ind") warning("The intercept was not considered for the X matrix (improper posterior)")
return(posterior)

}




################################################
#########Normal: densidade priori cjta (phi,nu)#########
################################################

##############
#### a data frame with the coordinates and the covariates (optional) is needed


dnsrprioroba=function(x,trend="cte",prior="reference",coords.col=1:2,kappa=0.5,cov.model="exponential",data){

####validation of arguments#####################################3

if(!is.numeric(x)) stop("Range parameter must be a real number in [0,Inf)")
if(x<0) stop("Range parameter must be a real number in [0,Inf)")
if(length(x)!=1) stop("Range parameter must be a real number in [0,Inf)")


if ( !is.numeric(coords.col)) stop("2D coordinates must be specified")
if ( length(coords.col)!=2) stop("2D coordinates must be specified")

if(!is.numeric(kappa))stop("kappa must be a real number in [0,Inf)")
if(kappa<0) stop("kappa must be a real number in [0,Inf)")

if (!(cov.model %in% c("matern","pow.exp","gaussian", "spherical","cauchy","exponential")))  stop('Valid covariance structures are (matern), (exponential), (gaussian), (spherical),
(pow.exp), (cauchy)')

if (!is.data.frame(data)) stop("data must be a data.frame")
if (!(prior %in% c("reference","jef.rul","jef.ind"))) stop("Accepted priors: reference, jef.rul, jef.ind")
coords=data[,coords.col]
H=as.matrix(dist(coords,upper=T,diag=T))
tau2=0
sigma=1
if(trend==0) if(prior!="jef.ind") stop("trend=0 is just valid for the jef.ind prior")
if(trend!=0) if(prior=="jef.ind") stop("The argument trend is not necessary when prior=jef.ind")

if(prior=="jef.ind"){
priorcand=priorijefindnorm(x=x,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2)
return(priorcand)
}
else{
if(class(trend)!="formula") if (!(trend %in% c("cte","1st","2nd"))) stop('trend must be (cte), (1st),(2nd) or a valid formula')
c1=coords[,1]
c2=coords[,2]
if(trend=="cte") xmat=as.matrix(rep(1,dim(data)[1]))
if(trend=="1st") xmat=model.matrix(~c1+c2)
if(trend=="2nd") xmat= model.matrix(~ c1 + c2 + I(c1^2) + I(c2^2) + I(c1*c2))
if(class(trend)=="formula") xmat=model.matrix(trend,data=data)



#############Looking for NA#####################
if ((sum(is.na(xmat))+sum(is.na(coords)))>0) stop ("NAs not allowed")

if (det(t(xmat)%*%xmat)==0) stop("the columns of x must be linearly independent")

if(prior=="reference") priorcand=priorirefnorm(x=x,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
if(prior=="jef.rul") priorcand=priorijefnorm(x=x,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
if(prior=="jef.ind") priorcand=priorijefindnorm(x=x,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2)
#if(prior=="vague") priorcand=vagprior(x=x,sigma=sigma,anu=anu,bnu=bnu,aphi=aphi,bphi=bphi)
return(priorcand)
}
}





################################################
#########densidade posteriori cjta (phi,nu)#####
################################################

##############
#### a data frame with the coordinates and the covariates (optional) is needed

dnsrposoba=function(x,formula,prior="reference",coords.col=1:2,kappa=0.5,cov.model="exponential",data,asigma=2.1,intphi){

####validation of arguments#####################################3

if(!is.numeric(x)) stop("Range parameter must be a real number in [0,Inf)")
if(x<0) stop("Range parameter must be a real number in [0,Inf)")
if(length(x)!=1) stop("Range parameter must be a real number in [0,Inf)")

if (class(formula)!="formula")  stop('The argument formula must be a valid formula')

if (!(prior %in% c("reference","jef.rul","jef.ind","vague"))) stop("Accepted priors: reference, jef.rul, jef.ind")

if ( !is.numeric(coords.col)) stop("2D coordinates must be specified")
if ( length(coords.col)!=2) stop("2D coordinates must be specified")

if(!is.numeric(kappa))stop("kappa must be a real number in [0,Inf)")
if(kappa<0) stop("kappa must be a real number in [0,Inf)")

if (!(cov.model %in% c("matern","pow.exp","gaussian", "spherical","cauchy","exponential")))  stop('Valid covariance structures are (matern), (exponential), (gaussian), (spherical),
(pow.exp), (cauchy)')

if (!is.data.frame(data)) stop("data must be a data.frame")

if(prior=="vague"){
if(!is.numeric(intphi)) stop("intphi must be a numeric 2 dimension vector")
if(!is.vector(intphi)) stop("intphi must be a numeric 2 dimension vector")
if(asigma<=2) stop("asigma must be a real number greater than 2 (proper posterior)")
if ( length(intphi)!=2) stop("A valid interval for the range parameter must be provided")
if(intphi[1]>intphi[2]) stop("Invalid interval for the range parameter")
if(intphi[1]<=0 |intphi[2]<=0) stop("Range parameter is greater than 0")
aphi=intphi[1]
bphi=intphi[2]
}

coords=data[,coords.col]
H=as.matrix(dist(coords,upper=T,diag=T))

xmat=model.matrix(formula,data=data)
if(prior=="jef.ind") if(max(xmat)==1 & min(xmat)==1) stop("This prior generates an improper posterior")
if(prior=="jef.ind") xmat=as.matrix(xmat)

y <-data[,all.vars(formula)[1]]

#############Looking for NA#####################
if ((sum(is.na(y))+sum(is.na(xmat))+sum(is.na(coords)))>0) stop ("NAs not allowed")

if (det(t(xmat)%*%xmat)==0) stop("the columns of x must be linearly independent")


tau2=0
sigma=1
cons=1
if(prior=="reference") posterior=posteriorrefnorm(y,cons=cons,sigma=sigma,x=x,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
if(prior=="jef.rul") posterior=posteriorjefnorm(x=x,y=y,cons=cons,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
if(prior=="jef.ind") posterior=posteriorjefindnorm(x=x,y=y,cons=cons,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
if(prior=="vague") posterior=posteriorvagnorm(x=x,y=y,cons=cons,sigma=sigma,aphi=aphi,bphi=bphi,asigma=asigma,xmat=xmat,cov.model=cov.model,H=H,tau2=tau2,kappa=kappa)

#if(prior=="jef.ind") warning("The intercept was not considered for the X matrix (improper posterior)")
return(posterior)

}



###########################################################################
############Integral M caso T##############################################
###########################################################################
intmT=function(formula,prior="reference",coords.col=1:2,kappa=0.5,cov.model="exponential",data,asigma=2.1,intphi="default",intnu=c(4.1,Inf),maxEval){

  ####validation of arguments#####################################3

  if(!is.numeric(maxEval)) stop("maxEval must be a positive integer value")
  if(maxEval<0) stop("maxEval must be a positive integer value")
  if(maxEval%%1 !=0 ) stop("maxEval must be a positive integer value")
  if(length(maxEval)!=1) stop("maxEval must be a positive integer value")

  if (class(formula)!="formula")  stop('The argument formula must be a valid formula')

  if (!(prior %in% c("reference","jef.rul","jef.ind","vague"))) stop("Accepted priors: reference, jef.rul, jef.ind")

  if ( !is.numeric(coords.col)) stop("2D coordinates must be specified")
  if ( length(coords.col)!=2) stop("2D coordinates must be specified")

  if(!is.numeric(kappa))stop("kappa must be a real number in [0,Inf)")
  if(kappa<0) stop("kappa must be a real number in [0,Inf)")

  if (!(cov.model %in% c("matern","pow.exp","gaussian", "spherical","cauchy","exponential")))  stop('Valid covariance structures are (matern), (exponential), (gaussian), (spherical),
                                                                                                    (pow.exp), (cauchy)')

  if (!is.data.frame(data)) stop("data must be a data.frame")

  coords=data[,coords.col]
  H=as.matrix(dist(coords,upper=T,diag=T))

  if(prior=="vague"){
    if(class(intphi)=="character"){
      if(intphi!="default") stop("A character different of 'default' is not allowed")
      if(cov.model=="spherical") stop("intphi argument must be specified for the spherical covariance case")
      distancias=unique(as.vector(H))
      bphi=try(find.phi(d=max(distancias),kappa=kappa,range=c(1e-04,100000),cut=0.05,cov.model=cov.model),silent=T)
      aphi=try(find.phi(d=min(distancias[distancias!=0]),kappa=kappa,range=c(1e-04,100000),cut=0.05,cov.model=cov.model),silent=T)
      if (class(aphi)=="try-error" | class(bphi)=="try-error") stop("error in calculating intphi by default")
    }else{
      if(!is.numeric(intphi)) stop("intphi must be a numeric 2 dimension vector")
      if(!is.vector(intphi)) stop("intphi must be a numeric 2 dimension vector")
      if(!is.numeric(intnu)) stop("intnu must be a numeric 2 dimension vector")
      if(!is.vector(intnu)) stop("intnu must be a numeric 2 dimension vector")
      if(asigma<=2) stop("asigma must be a real number greater than 2 (proper posterior)")
      if ( length(intphi)!=2) stop("A valid interval for the range parameter must be provided")
      if ( length(intnu)!=2) stop("A valid interval for the degrees of freedom must be provided")
      if(intphi[1]>intphi[2]) stop("Invalid interval for the range parameter")
      if(intphi[1]<=0 |intphi[2]<=0) stop("Range parameter is greater than 0")
      if(intnu[1]>intnu[2]) stop("Invalid interval for the degrees of freedom")
      if(intnu[1]<4 |intnu[2]<4) stop("degrees of freedom must be in (4,Inf)")
      aphi=intphi[1]
      bphi=intphi[2]
    }

    anu=intnu[1]
    bnu=intnu[2]

  }


  xmat=model.matrix(formula,data=data)
  if(prior=="jef.ind") if(max(xmat)==1 & min(xmat)==1) stop("This prior generates an improper posterior")
  if(prior=="jef.ind") xmat=as.matrix(xmat[,-1])

  y <-data[,all.vars(formula)[1]]

  #############Looking for NA#####################
  if ((sum(is.na(y))+sum(is.na(xmat))+sum(is.na(coords)))>0) stop ("NAs not allowed")

  if (det(t(xmat)%*%xmat)==0) stop("the columns of x must be linearly independent")


  tau2=0
  sigma=1

  if(prior=="reference") intm=intm1tref(y=y,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
  if(prior=="jef.rul") intm=intm1tjef(y=y,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
  if(prior=="jef.ind") intm=intm1tjefind(y=y,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
  if(prior=="vague") intm=intm1tvag(y=y,sigma=sigma,anu=anu,bnu=bnu,aphi=aphi,bphi=bphi,asigma=asigma,xmat=xmat,cov.model=cov.model,H=H,tau2=tau2,kappa=kappa,maxEval=maxEval)

  #if(prior=="jef.ind") warning("The intercept was not considered for the X matrix (improper posterior)")
  return(intm)

}




###########################################################################
############Integral M caso normal##############################################
###########################################################################

intmnorm=function(formula,prior="reference",coords.col=1:2,kappa=0.5,cov.model="exponential",data,asigma=2.1,intphi,maxEval){

####validation of arguments#####################################3

if(!is.numeric(maxEval)) stop("maxEval must be a positive integer value")
if(maxEval<0) stop("maxEval must be a positive integer value")
if(maxEval%%1 !=0 ) stop("maxEval must be a positive integer value")
if(length(maxEval)!=1) stop("maxEval must be a positive integer value")

if (class(formula)!="formula")  stop('The argument formula must be a valid formula')

if (!(prior %in% c("reference","jef.rul","jef.ind","vague"))) stop("Accepted priors: reference, jef.rul, jef.ind")

if ( !is.numeric(coords.col)) stop("2D coordinates must be specified")
if ( length(coords.col)!=2) stop("2D coordinates must be specified")

if(!is.numeric(kappa))stop("kappa must be a real number in [0,Inf)")
if(kappa<0) stop("kappa must be a real number in [0,Inf)")

if (!(cov.model %in% c("matern","pow.exp","gaussian", "spherical","cauchy","exponential")))  stop('Valid covariance structures are (matern), (exponential), (gaussian), (spherical),
(pow.exp), (cauchy)')

if (!is.data.frame(data)) stop("data must be a data.frame")

if(prior=="vague"){
if(!is.numeric(intphi)) stop("intphi must be a numeric 2 dimension vector")
if(!is.vector(intphi)) stop("intphi must be a numeric 2 dimension vector")
if(asigma<=2) stop("asigma must be a real number greater than 2 (proper posterior)")
if ( length(intphi)!=2) stop("A valid interval for the range parameter must be provided")
if(intphi[1]>intphi[2]) stop("Invalid interval for the range parameter")
if(intphi[1]<=0 |intphi[2]<=0) stop("Range parameter is greater than 0")
aphi=intphi[1]
bphi=intphi[2]
}

coords=data[,coords.col]
H=as.matrix(dist(coords,upper=T,diag=T))

xmat=model.matrix(formula,data=data)
if(prior=="jef.ind") if(max(xmat)==1 & min(xmat)==1) stop("This prior generates an improper posterior")
#if(prior=="jef.ind") xmat=as.matrix(xmat[,-1])

y <-data[,all.vars(formula)[1]]

#############Looking for NA#####################
if ((sum(is.na(y))+sum(is.na(xmat))+sum(is.na(coords)))>0) stop ("NAs not allowed")

if (det(t(xmat)%*%xmat)==0) stop("the columns of x must be linearly independent")


tau2=0
sigma=1

if(prior=="reference") intm=intm1normref(y=y,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
if(prior=="jef.rul") intm=intm1normjef(y=y,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
if(prior=="jef.ind") intm=intm1normjefind(y=y,sigma=sigma,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
if(prior=="vague") intm=intm1normvag(y=y,sigma=sigma,aphi=aphi,bphi=bphi,asigma=asigma,xmat=xmat,cov.model=cov.model,H=H,tau2=tau2,kappa=kappa,maxEval=maxEval)

#if(prior=="jef.ind") warning("The intercept was not considered for the X matrix (improper posterior)")
return(intm)

}

############################################################################
#####Amostrando da posteriori Caso normal###################################
############################################################################


nsroba=function(formula,proposal="unif",method="median",candpar,prior="reference",coords.col=1:2,kappa=0.5,cov.model="matern",data,asigma=2.1,intphi="default",ini.pars,burn=500,iter=5000,thin=10){

#xmat=xmat,proposal=proposal,candpar=candpar,prior=prior,asigma=asigma,y=y,coords=coordstot,covini=cov.ini,nuini=nuini,tau2=0,kappa=kappa,cov.model=cov.model,aphi=aphi,bphi=bphi,burn=burn,iter=iter,thin=thin

####validation of arguments#####################################3
x=ini.pars
if(!is.vector(x)) stop("ini.pars must be a numeric vector of length 2")
if(length(x)!=2) stop("ini.pars must be a numeric vector of length 2")
if(x[1]<0) stop ("sigma2 must be must be a real number in (0,Inf)")
if(x[2]<0) stop("Range parameter must be a real number in (0,Inf)")
if (!(method %in% c("mean","median","mode"))) stop("Estimation methods: mean, median, mode")
if (!(proposal %in% c("unif","gamma"))) stop("Accepted proposals: unif, gamma")

if(burn%%1 !=0 ) stop("burn must be a positive integer value")
if(iter%%1 !=0 ) stop("iter must be a positive integer value")
if(thin%%1 !=0 ) stop("thin must be a positive integer value")

if(burn> iter) stop("burning cannot be greater than number of iterations")
if(thin> iter) stop("thin cannot be greater than number of iterations")

if(proposal=="gamma") if(!is.numeric(candpar)) stop("candpar must be a numeric hyperparameter")

if (class(formula)!="formula")  stop('The argument formula must be a valid formula')

if (!(prior %in% c("reference","jef.rul","jef.ind","vague"))) stop("Accepted priors: reference, jef.rul, jef.ind")

if ( !is.numeric(coords.col)) stop("2D coordinates must be specified")
if ( length(coords.col)!=2) stop("2D coordinates must be specified")

if(!is.numeric(kappa))stop("kappa must be a real number in [0,Inf)")
if(kappa<0) stop("kappa must be a real number in [0,Inf)")

if (!(cov.model %in% c("matern","pow.exp","gaussian", "spherical","cauchy","exponential")))  stop('Valid covariance structures are (matern), (exponential), (gaussian), (spherical),
(pow.exp), (cauchy)')

if (!is.data.frame(data)) stop("data must be a data.frame")

coords=data[,coords.col]
H=as.matrix(dist(coords,upper=T,diag=T))


if(proposal=="unif"){
 if(class(intphi)=="character"){
if(intphi!="default") stop("A character different of 'default' is not allowed")
if(cov.model=="spherical") stop("intphi argument must be specified for the spherical covariance case")
distancias=unique(as.vector(H))
bphi=try(find.phi(d=max(distancias),kappa=kappa,range=c(1e-04,100000),cut=0.05,cov.model=cov.model),silent=T)
aphi=try(find.phi(d=min(distancias[distancias!=0]),kappa=kappa,range=c(1e-04,100000),cut=0.05,cov.model=cov.model),silent=T)
if (class(aphi)=="try-error" | class(bphi)=="try-error") stop("error in calculating intphi by default")
}



if(class(intphi)!="character"){
if(!is.numeric(intphi)) stop("intphi must be a numeric 2 dimension vector")
if(!is.vector(intphi)) stop("intphi must be a numeric 2 dimension vector")
if ( length(intphi)!=2) stop("A valid interval for the range parameter must be provided")
if(intphi[1]>intphi[2]) stop("Invalid interval for the range parameter")
if(intphi[1]<=0 |intphi[2]<=0) stop("Range parameter is greater than 0")
aphi=intphi[1]
bphi=intphi[2]
}

}






if(prior=="vague"){
if(proposal=="gamma") stop("For the vague prior, intphi must be specified")
if(!is.numeric(asigma)) stop("asigma must be a real number greater than 2 (proper posterior)")
if(asigma<=2) stop("asigma must be a real number greater than 2 (proper posterior)")
}



xmat=model.matrix(formula,data=data)
if(prior=="jef.ind") if(max(xmat)==1 & min(xmat)==1) stop("This prior generates an improper posterior")
if(prior=="jef.ind") xmat=as.matrix(xmat)

y <-data[,all.vars(formula)[1]]

#############Looking for NA#####################
if ((sum(is.na(y))+sum(is.na(xmat))+sum(is.na(coords)))>0) stop ("NAs not allowed")

if (det(t(xmat)%*%xmat)==0) stop("the columns of x must be linearly independent")


tau2=0
#sigma=1


if(proposal=="gamma") res=suppressWarnings(baysparefnorm(xmat=xmat,method=method,proposal=proposal,candpar=candpar,prior=prior,asigma=asigma,y=y,coords=coords,covini=ini.pars,tau2=0,kappa=kappa,cov.model=cov.model,aphi=aphi,bphi=bphi,burn=burn,iter=iter,thin=thin))
if(proposal=="unif") res=suppressWarnings(baysparefnorm(xmat=xmat,method=method,proposal=proposal,prior=prior,asigma=asigma,y=y,coords=coords,covini=ini.pars,tau2=0,kappa=kappa,cov.model=cov.model,aphi=aphi,bphi=bphi,burn=burn,iter=iter,thin=thin))
res$formula=formula
res$prior=prior
res$proposal=proposal
res=res

#if(prior=="jef.ind") warning("The intercept was not considered for the X matrix (improper posterior)")
 class(res)="nsroba"
return(res)


}


############################################################################
#####Amostrando da posteriori Caso T###################################
############################################################################


tsroba=function(formula,proposal="unif",method="median",candpar,prior="reference",coords.col=1:2,kappa=0.5,cov.model="matern",data,asigma=2.1,intphi="default",intnu="default",ini.pars,burn=500,iter=5000,thin=10){

#xmat=xmat,proposal=proposal,candpar=candpar,prior=prior,asigma=asigma,y=y,coords=coordstot,covini=cov.ini,nuini=nuini,tau2=0,kappa=kappa,cov.model=cov.model,aphi=aphi,bphi=bphi,burn=burn,iter=iter,thin=thin

####validation of arguments#####################################3
x=ini.pars
if(!is.vector(x)) stop("ini.pars must be a numeric vector of length 3")
if(length(x)!=3) stop("ini.pars must be a numeric vector of length 3")
if(x[1]<0) stop ("sigma2 must be must be a real number in (0,Inf)")
if(x[2]<0) stop("Range parameter must be a real number in (0,Inf)")
if(x[3]<=4) stop("degrees of freedom must be greater than 4")




if (!(proposal %in% c("unif","gamma"))) stop("Accepted proposals: unif, gamma")

if(burn%%1 !=0 ) stop("burn must be a positive integer value")
if(iter%%1 !=0 ) stop("iter must be a positive integer value")
if(thin%%1 !=0 ) stop("thin must be a positive integer value")

if(burn> iter) stop("burning cannot be greater than number of iterations")
if(thin> iter) stop("thin cannot be greater than number of iterations")

if(proposal=="gamma") if(!is.numeric(candpar)) stop("candpar must be a numeric hyperparameter")

if (class(formula)!="formula")  stop('The argument formula must be a valid formula')

if (!(prior %in% c("reference","jef.rul","jef.ind","vague"))) stop("Accepted priors: reference, jef.rul, jef.ind")
if (!(method %in% c("mean","median","mode"))) stop("Estimation methods: mean, median, mode")
if ( !is.numeric(coords.col)) stop("2D coordinates must be specified")
if ( length(coords.col)!=2) stop("2D coordinates must be specified")

if(!is.numeric(kappa))stop("kappa must be a real number in [0,Inf)")
if(kappa<0) stop("kappa must be a real number in [0,Inf)")

if (!(cov.model %in% c("matern","pow.exp","gaussian", "spherical","cauchy","exponential")))  stop('Valid covariance structures are (matern), (exponential), (gaussian), (spherical),
(pow.exp), (cauchy)')

if (!is.data.frame(data)) stop("data must be a data.frame")

coords=data[,coords.col]
H=as.matrix(dist(coords,upper=T,diag=T))


if(proposal=="unif"){
 if(class(intphi)=="character"){
if(intphi!="default") stop("A character different of 'default' is not allowed")
if(cov.model=="spherical") stop("intphi argument must be specified for the spherical covariance case")
distancias=unique(as.vector(H))
bphi=try(find.phi(d=max(distancias),kappa=kappa,range=c(1e-04,100000),cut=0.05,cov.model=cov.model),silent=T)
aphi=try(find.phi(d=min(distancias[distancias!=0]),kappa=kappa,range=c(1e-04,100000),cut=0.05,cov.model=cov.model),silent=T)
if (class(aphi)=="try-error" | class(bphi)=="try-error") stop("error in calculating intphi by default")
}

if(class(intphi)!="character"){
if(!is.numeric(intphi)) stop("intphi must be a numeric 2 dimension vector")
if(!is.vector(intphi)) stop("intphi must be a numeric 2 dimension vector")
if ( length(intphi)!=2) stop("A valid interval for the range parameter must be provided")
if(intphi[1]>intphi[2]) stop("Invalid interval for the range parameter")
if(intphi[1]<=0 |intphi[2]<=0) stop("Range parameter is greater than 0")
aphi=intphi[1]
bphi=intphi[2]
}
}

if(class(intnu)=="character"){
if(intnu!="default") stop("A character different of 'default' is not allowed")
anu=4
bnu=30
}

if(class(intnu)!="character"){
if(!is.numeric(intnu)) stop("intnu must be a numeric 2 dimension vector")
if(!is.vector(intnu)) stop("intnu must be a numeric 2 dimension vector")
if ( length(intnu)!=2) stop("A valid interval for the range parameter must be provided")
if(intnu[1]>intnu[2]) stop("Invalid interval for the degrees of freedom")
if(intnu[1]<4 |intnu[2]<4) stop("degrees of freedom for the proposal must be in (4,Inf]")
anu=intnu[1]
bnu=intnu[2]
}

if(prior=="vague"){
if(proposal=="gamma") stop("For the vague prior, intphi must be specified (just unif)")
if(!is.numeric(asigma)) stop("asigma must be a real number greater than 2 (proper posterior)")
if(asigma<=2) stop("asigma must be a real number greater than 2 (proper posterior)")
}



xmat=model.matrix(formula,data=data)
if(prior=="jef.ind") if(max(xmat)==1 & min(xmat)==1) stop("This prior generates an improper posterior")
if(prior=="jef.ind") xmat=as.matrix(xmat)

y <-data[,all.vars(formula)[1]]

#############Looking for NA#####################
if ((sum(is.na(y))+sum(is.na(xmat))+sum(is.na(coords)))>0) stop ("NAs not allowed")

if (det(t(xmat)%*%xmat)==0) stop("the columns of x must be linearly independent")


tau2=0
#sigma=1

p=ncol(xmat)
anujef=max(4,p)
if(proposal=="unif") {
if(prior=="reference"| prior=="jef.ind")  res=suppressWarnings(bayspaestT1(xmat=xmat,method=method,proposal=proposal,prior=prior,y=y,coords=coords,covini=ini.pars[1:2],nuini=ini.pars[3],tau2=0,kappa=kappa,cov.model=cov.model,aphi=aphi,bphi=bphi,anu=anu,bnu=bnu,burn=burn,iter=iter,thin=thin))
if(prior=="jef.rul") res=suppressWarnings(bayspaestTjef(proposal=proposal,xmat=xmat,method=method,y=y,coords=coords,covini=ini.pars[1:2],nuini=ini.pars[3],tau2=0,kappa=kappa,cov.model=cov.model,aphi=aphi,bphi=bphi,anu=anujef,bnu=bnu,burn=burn,iter=iter,thin=thin))
if(prior=="vague")   res= suppressWarnings(bayspaestTvag(xmat=xmat,asigma=asigma,y=y,method=method,coords=coords,covini=ini.pars[1:2],nuini=ini.pars[3],tau2=0,kappa=kappa,cov.model=cov.model,aphi=aphi,bphi=bphi,anu=anu,bnu=bnu,burn=burn,iter=iter,thin=thin))
}

if(proposal=="gamma") {
if(prior=="reference"| prior=="jef.ind")  res=suppressWarnings(bayspaestT1(xmat=xmat,proposal=proposal,method=method,candpar=candpar,prior=prior,y=y,coords=coords,covini=ini.pars[1:2],nuini=ini.pars[3],tau2=0,kappa=kappa,cov.model=cov.model,anu=anu,bnu=bnu,burn=burn,iter=iter,thin=thin))
if(prior=="jef.rul") res=suppressWarnings(bayspaestTjef(candpar=candpar,proposal=proposal,xmat=xmat,y=y,coords=coords,method=method,covini=ini.pars[1:2],nuini=ini.pars[3],tau2=0,kappa=kappa,cov.model=cov.model,anu=anujef,bnu=bnu,burn=burn,iter=iter,thin=thin))
if(prior=="vague")   res= suppressWarnings(bayspaestTvag(xmat=xmat,asigma=asigma,y=y,coords=coords,method=method,covini=ini.pars[1:2],nuini=ini.pars[3],tau2=0,kappa=kappa,cov.model=cov.model,aphi=aphi,bphi=bphi,anu=anu,bnu=bnu,burn=burn,iter=iter,thin=thin))
}

#if(prior=="jef.ind") warning("The intercept was not considered for the X matrix (improper posterior)")
res$formula=formula
res$prior=prior
res$proposal=proposal
res=res
class(res)="tsroba"
return(res)
}


tsrobapred=function(obj,xpred,coordspred){
  if(class(obj)!="tsroba") stop ("Argument obj must be of class tsroba")
  if(!is.numeric(xpred) & !is.data.frame(xpred)) stop ("xpred must be a numeric matrix or data.frame")
  if(!is.numeric(coordspred) & !is.data.frame(coordspred)) stop ("coordspred must be a numeric matrix or data.frame")
  if (!is.matrix(xpred)) xpred=as.matrix(xpred)
  if (!is.matrix(coordspred)) as.matrix(coordspred)
  if(ncol(coordspred)!=2) stop("2D coordinates must be specified")
  if(nrow(xpred)!=nrow(coordspred)) stop("xpred does not have the same number of lines than coordspred")
  if(sum(is.na(xpred)) > 0) stop("There are some NA values in xpred")
  if(sum(is.na(coordspred)) > 0) stop("There are some NA values in coordspred")
  aux=apply(obj$dist,1,FUN=prediction1,xobs=obj$X,z=obj$y,xpred=xpred,coordspred=coordspred,coordsobs=obj$coords,cov.model=obj$type,tau2=0,kappa=obj$kappa)
  obj.out=apply(aux,1,median)
  return(obj.out)
}


nsrobapred1=function(xpred,coordspred,obj){
  if(!is.numeric(xpred) & !is.data.frame(xpred)) stop ("xpred must be a numeric matrix or data.frame")
  if(!is.numeric(coordspred) & !is.data.frame(coordspred)) stop ("coordspred must be a numeric matrix or data.frame")
  if (!is.matrix(xpred)) xpred=as.matrix(xpred)
  if (!is.matrix(coordspred)) as.matrix(coordspred)
  if(ncol(coordspred)!=2) stop("2D coordinates must be specified")
  if(nrow(xpred)!=nrow(coordspred)) stop("xpred does not have the same number of lines than coordspred")
  if(class(obj)!="nsroba") stop("an object of the class nsroba must be provided")

  if(sum(is.na(xpred)) > 0) stop("There are some NA values in xpred")
  if(sum(is.na(coordspred)) > 0) stop("There are some NA values in coordspred")

  aux=apply(obj$dist,1,FUN=predictionnorm,xobs=obj$X,z=obj$y,xpred=xpred,coordspred=coordspred,coordsobs=obj$coords,cov.model=obj$type,tau2=0,kappa=obj$kappa)
  obj.out=apply(aux,1,median)
  return(obj.out)

}


summary.nsroba=function(object, ...){
  if(class(object)!="nsroba")  stop("an object of class nsroba  must be provided")
    cat('\n')
    call <- match.call()
    cat("Call:\n")
    print(call)
    cat('\n')
    cat('\n\n')
    cat('---------------------------------------------------\n')
    cat('  Objective Bayesian Analysis for NSR model  \n')
    cat('---------------------------------------------------\n')
    cat('\n')
    cat("*Formula:")
    cat('\n')
    print(object$formula)
    cat('\n')
    cat("*Prior:",object$prior)
    cat('\n')
    cat("*Proposal:",object$proposal)
    cat('\n')
    cat("*Covariance structure:",object$type)
    cat('\n')
    cat('---------\n')
    cat('Estimates\n')
    cat('---------\n')
    cat('\n')
    trends=object$X
    l = ncol(trends)
    lab = numeric(l+2)
    for (i in 1:l) lab[i] = paste('beta ',i-1,sep='')
    lab[l+1] = 'sigma2'
    lab[l+2] ='phi'
    tab = t(round(object$theta,4))
    colnames(tab)=t(lab)
    #colnames(tab)="Estimated"
    print(tab)
    invisible(list(mean.str=object$theta[1:l],var.str=object$theta[(l+1):(l+2)],betaF=object$betaF,sigmaF=object$sigmaF,phiF=object$phiF))
  }



  summary.tsroba=function(object, ...){
    if(class(object)!="tsroba")  stop("an object of class tsroba  must be provided")
  #Running the algorithm
  cat('\n')
  call <- match.call()
  cat("Call:\n")
  print(call)
  cat('\n')
  cat('\n\n')
  cat('---------------------------------------------------\n')
  cat('  Objective Bayesian Analysis for TSR model  \n')
  cat('---------------------------------------------------\n')
  cat('\n')
  cat("*Formula:")
  cat('\n')
  print(object$formula)
  cat('\n')
  cat("*Prior:",object$prior)
  cat('\n')
  cat("*Proposal:",object$proposal)
  cat('\n')
  cat("*Covariance structure:",object$type)
  cat('\n')
  cat('---------\n')
  cat('Estimates\n')
  cat('---------\n')
  cat('\n')
  trends=object$X
  l = ncol(trends)
  lab = numeric(l+2)
  for (i in 1:l) lab[i] = paste('beta ',i-1,sep='')
  lab[l+1] = 'sigma2'
  lab[l+2] ='phi'
  lab[l+3] ='nu'
  tab = t(round(object$theta,4))
  colnames(tab)=t(lab)
  #colnames(tab)="Estimated"
  print(tab)

  invisible(list(mean.str=object$theta[1:l],var.str=object$theta[(l+1):(l+3)],betaF=object$betaF,sigmaF=object$sigmaF,phiF=object$phiF,nuF=object$nuF))
}

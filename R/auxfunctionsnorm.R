
###############################################################
#########Caso Normal: densidade ref priori phi#####################
###############################################################

priorirefnorm=function(sigma,x,H,kappa,cov.model,tau2,xmat){
covpars=c(sigma,x)
p=ncol(xmat)
ndata=nrow(xmat)
nu1=ndata-p
cand=covpars
sigma=covpars[1]
phi=covpars[2]
phicand=phi
covacand=varcov.spatial(H=H,cov.model=cov.model,cov.pars=cand,nugget=tau2,kappa=kappa)
covacand=covacand+(10e-10*diag(1,ndata))

derivcovacand=sigma*derivcormatrix1(H,phi=phicand,kappa=kappa,type=cov.model)$dev1
covacandinv=solve(covacand,tol=10e-22)
Pcand=diag(ndata)-(xmat%*%solve(t(xmat)%*%covacandinv%*%xmat)%*%t(xmat)%*%covacandinv)
Wcand=derivcovacand%*%covacandinv%*%Pcand

A=Wcand%*%Wcand

rescand=sum(diag(A))- (1/nu1)*(sum(diag(Wcand))^2)

priorcand=sqrt(rescand)


return(priorcand)
}



###############################################################
#########Caso Normal: densidade ref posteriori phi#####################
###############################################################

posteriorrefnorm=function(y,cons,sigma,x,H,kappa,cov.model,tau2,xmat){
ycomp=y
n=nrow(xmat)
p=ncol(xmat)
v=n-p
last=c(sigma,x)
Psilastaux=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
#Psilast=Psilastaux
Psilast=Psilastaux

Rlast=(Psilast/sigma)-((tau2/sigma)*diag(n)) +(10e-10*diag(1,n))
Rlastinv=solve(Rlast,tol=10e-22)
vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
     vbetalast<-(vbetalast+t(vbetalast))/2

mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)



reflast=priorirefnorm(sigma=sigma,x=x,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
detrlast=det(Rlast)
if(detrlast<1e-323){
detrlast=1e-323
}

posterior=(((S2last)^(-v/2))*sqrt(det(vbetalast))/sqrt(detrlast))*reflast
return(posterior/cons)
}



########################################################################
#### Caso normal: (densidade priori Jeffreys)############################
########################################################################

priorijefnorm=function(sigma,x,H,kappa,cov.model,tau2,xmat){
covpars=c(sigma,x)
ndata=nrow(xmat)
n=ndata
cand=covpars
sigma=covpars[1]
phi=covpars[2]
phicand=phi
covacand=varcov.spatial(H=H,cov.model=cov.model,cov.pars=cand,nugget=tau2,kappa=kappa)
covacand=covacand+(10e-10*diag(1,ndata))

derivcovacand=sigma*derivcormatrix1(H,phi=phicand,kappa=kappa,type=cov.model)$dev1
covacandinv=solve(covacand,tol=10e-22)
Wcand=derivcovacand%*%covacandinv
compbeta=(t(xmat)%*%covacandinv%*%xmat)

detcompbeta=det(compbeta)
A=Wcand%*%Wcand

rescand=sum(diag(A))- (1/n)*(sum(diag(Wcand))^2)

priorcand=sqrt(rescand*detcompbeta)
return(priorcand)
}


######################################################
####Caso normal: Posterior Jeffreys###################
######################################################

posteriorjefnorm=function(y,cons,sigma,x,H,kappa,cov.model,tau2,xmat){
ycomp=y
n=nrow(xmat)
p=ncol(xmat)
v=n-p
last=c(sigma,x)
Psilastaux=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
#Psilast=Psilastaux
Psilast=Psilastaux

Rlast=(Psilast/sigma)-((tau2/sigma)*diag(n)) +(10e-4*diag(1,n))
Rlastinv=solve(Rlast,tol=10e-22)
vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
     vbetalast<-(vbetalast+t(vbetalast))/2

mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)



reflast=priorijefnorm(sigma=sigma,x=x,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
detrlast=det(Rlast)
if(detrlast<1e-323){
detrlast=1e-323
}

a= (p/2) + 1

posterior=(((S2last)^(-v/2-a+1))*sqrt(det(vbetalast))/sqrt(detrlast))*reflast
return(posterior/cons)

#posterior=(((S2last)^(-v/2))*sqrt(det(vbetalast))/sqrt(detrlast))*reflast
return(posterior/cons)
}





########################################################################
#### Caso normal: (densidade priori indep Jeffreys)############################
########################################################################

priorijefindnorm=function(sigma,x,H,kappa,cov.model,tau2){
covpars=c(sigma,x)
ndata=nrow(H)
n=ndata
cand=covpars
sigma=covpars[1]
phi=covpars[2]
phicand=phi
covacand=varcov.spatial(H=H,cov.model=cov.model,cov.pars=cand,nugget=tau2,kappa=kappa)
covacand=covacand+(10e-10*diag(1,ndata))

derivcovacand=sigma*derivcormatrix1(H,phi=phicand,kappa=kappa,type=cov.model)$dev1
covacandinv=solve(covacand,tol=10e-22)
Wcand=derivcovacand%*%covacandinv
A=Wcand%*%Wcand

rescand=sum(diag(A))- (1/n)*(sum(diag(Wcand))^2)

priorcand=sqrt(rescand)
return(priorcand)
}


######################################################
####Caso normal: Posterior ind Jeffreys###################
######################################################

posteriorjefindnorm=function(y,cons,sigma,x,H,kappa,cov.model,tau2,xmat){
ycomp=y
n=nrow(xmat)
p=ncol(xmat)
v=n-p
last=c(sigma,x)
Psilastaux=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
#Psilast=Psilastaux
Psilast=Psilastaux

Rlast=(Psilast/sigma)-((tau2/sigma)*diag(n)) +(10e-4*diag(1,n))
Rlastinv=solve(Rlast,tol=10e-22)
vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
     vbetalast<-(vbetalast+t(vbetalast))/2

mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)



reflast=priorijefindnorm(sigma=sigma,x=x,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2)
detrlast=det(Rlast)
if(detrlast<1e-323){
detrlast=1e-323
}

posterior=(((S2last)^(-v/2))*sqrt(det(vbetalast))/sqrt(detrlast))*reflast

return(posterior/cons)
}





######################################################
####Caso normal: Posterior vaga###################
######################################################

posteriorvagnorm=function(y,cons,sigma,x,H,kappa,cov.model,tau2,xmat,aphi,bphi,asigma){
ycomp=y
n=nrow(xmat)
p=ncol(xmat)
v=n-p
last=c(sigma,x)
Psilastaux=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
#Psilast=Psilastaux
Psilast=Psilastaux

Rlast=(Psilast/sigma)-((tau2/sigma)*diag(n)) +(10e-4*diag(1,n))
Rlastinv=solve(Rlast,tol=10e-22)
vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
     vbetalast<-(vbetalast+t(vbetalast))/2

mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)



reflast=dunif(x=x,aphi,bphi)
detrlast=det(Rlast)
if(detrlast<1e-323){
detrlast=1e-323
}

a= asigma
posterior=(((S2last)^(-v/2-a+1))*sqrt(det(vbetalast))/sqrt(detrlast))*reflast
return(posterior/cons)
}



#######################################
###Integral M caso Normal###################
#######################################



intm1normref=function(kappa,H,sigma,cov.model,tau2,xmat,y,maxEval){
consm1mat2=hcubature(priorirefnorm, lowerLimit=0, upperLimit=Inf,H=H,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
consm1posmat2=hcubature(posteriorrefnorm, lowerLimit=0, upperLimit=Inf,H=H,y=y,cons=consm1mat2$integral,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
ms=c(consm1posmat2$integral)
return(ms)
}


intm1normjef=function(kappa,H,sigma,cov.model,tau2,xmat,y,maxEval){
consm1mat2=hcubature(priorijefnorm, lowerLimit=0, upperLimit=Inf,H=H,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
consm1posmat2=hcubature(posteriorjefnorm, lowerLimit=0, upperLimit=Inf,H=H,y=y,cons=consm1mat2$integral,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
ms=c(consm1posmat2$integral)
return(ms)
}

intm1normjefind=function(kappa,H,sigma,cov.model,tau2,xmat,y,maxEval){
consm1mat2=hcubature(priorijefindnorm, lowerLimit=0, upperLimit=Inf,H=H,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,maxEval=maxEval)
consm1posmat2=hcubature(posteriorjefindnorm, lowerLimit=0, upperLimit=Inf,H=H,y=y,cons=consm1mat2$integral,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
ms=c(consm1posmat2$integral)
return(ms)
}

intm1normvag=function(kappa,H,sigma,cov.model,tau2,xmat,y,maxEval,aphi,bphi,asigma){
#consm1mat2=hcubature(priorirefnorm, lowerLimit=0, upperLimit=Inf,H=H,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
consm1posmat2=hcubature(posteriorvagnorm, lowerLimit=aphi, upperLimit=bphi,H=H,y=y,cons=1,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval,aphi=aphi,bphi=bphi,asigma=asigma)
ms=c(consm1posmat2$integral)
return(ms)
}





#############################################################################
##############amostrando da posteriori normal################################
#############################################################################

########usando a minha priori######

baysparefnorm=function(xmat,proposal,method,prior,y,coords,covini,nuini,tau2,kappa,cov.model,aphi,bphi,anu,bnu,burn,iter,thin,asigma){
p=ncol(xmat)
n=nrow(xmat)
H=as.matrix(dist(coords,upper=T,diag=T))
Psi1=varcov.spatial(H=H,cov.model=cov.model,cov.pars=covini,nugget=tau2,kappa=kappa)
#Psi1=varcov.spatial(coords,cov.model=cov.model,cov.pars=covini,nugget=tau2,kappa=kappa)$varcov
Psi1=(Psi1+t(Psi1))/2
beta1=solve(t(xmat)%*%solve(Psi1)%*%xmat)%*%t(xmat)%*%solve(Psi1)%*%y

betaF=matrix(0,iter+1,p)
betaF[1,]=beta1
sigmaF=c(covini[1],rep(0,iter))
phiF=c(covini[2],rep(0,iter))
count=0
countiter=0

#if(proposal=="unif"){
for(j in 2:iter){
ycomp=y

n=nrow(xmat)
p=ncol(xmat)
v=n-p

#phicand=rtrunc(1, spec="gamma",shape=phiF[j-1],scale=candpar,a=aphi,b=Inf)
phicand=runif(1,aphi,bphi)

cand=c(sigmaF[j-1],phicand)
Psicand=varcov.spatial(H=H,cov.model=cov.model,cov.pars=cand,nugget=tau2,kappa=kappa)

Rcand=(Psicand/sigmaF[j-1])-((tau2/sigmaF[j-1])*diag(n)) +(10e-4*diag(1,n))

Rcandinv=solve(Rcand)
vbetacand<- solve(t(xmat)%*%Rcandinv%*%xmat)
     vbetacand<-(vbetacand+t(vbetacand))/2

mubetacand<-   solve(t(xmat)%*%Rcandinv%*%xmat)%*%t(xmat)%*%Rcandinv%*%ycomp
S2cand=t(ycomp-xmat%*%mubetacand)%*%Rcandinv%*%(ycomp-xmat%*%mubetacand)/(n-p)

detrcand=det(Rcand)
if(detrcand<1e-323){
detrcand=1e-323
}

last=c(sigmaF[j-1],phiF[j-1])
Psilast=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)

Rlast=(Psilast/sigmaF[j-1])-((tau2/sigmaF[j-1])*diag(n)) +(10e-4*diag(1,n))
Rlastinv=solve(Rlast)
vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
     vbetalast<-(vbetalast+t(vbetalast))/2

mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)

philast=phiF[j-1]

detrlast=det(Rlast)
if(detrlast<1e-323){
detrlast=1e-323
}

if(prior=="reference"){
reflast=priorirefnorm(H=H,x=philast,sigma=1,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
refcand=priorirefnorm(H=H,x=phicand,sigma=1,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
a=1
}

if(prior=="jef.rul"){
reflast=priorijefnorm(H=H,x=philast,sigma=1,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
refcand=priorijefnorm(H=H,x=phicand,sigma=1,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
a=1+ (p/2)
}

if(prior=="jef.ind"){
reflast=priorijefindnorm(H=H,x=philast,sigma=1,kappa=kappa,cov.model=cov.model,tau2=tau2)
refcand=priorijefindnorm(H=H,x=phicand,sigma=1,kappa=kappa,cov.model=cov.model,tau2=tau2)
a=1
}

if(prior=="vague"){
reflast=dunif(philast,aphi,bphi)
refcand= dunif(phicand,aphi,bphi)
a=asigma
}


priorphicand=(((S2cand)^(-v/2-a+1))*sqrt(det(vbetacand))/sqrt(detrcand))*refcand
#priorphicand=(((S2cand)^(-v/2-a+1))*sqrt(det(vbetacand))/sqrt(detrcand))*refcand*dtrunc(nulast,spec="exp",rate=1/nulast,a=anu,b=bnu)*dtrunc(philast, spec="gamma",shape=phiF[j-1],scale=candpar,a=aphi,b=Inf)


priorphilast=(((S2last)^(-v/2-a+1))*sqrt(det(vbetalast))/sqrt(detrlast))*reflast
#priorphilast=(((S2last)^(-v/2))*sqrt(det(vbetalast))/sqrt(detrlast))*reflast*dtrunc(nucand,spec="exp",rate=1/nulast,a=anu,b=bnu)*dtrunc(phicand, spec="gamma",shape=phiF[j-1],scale=candpar,a=aphi,b=Inf)


Num=priorphicand
Den=priorphilast

unif=runif(1)
  if(unif<Num/Den)
  {
next.pi2F <- cand
mubeta=mubetacand
vbeta=vbetacand
S2=S2cand
count=count+1
}else
  {
    next.pi2F <- last
  mubeta=mubetalast
vbeta=vbetalast
S2=S2last
  }



betaF[j,]=rmvt(1,mu=t(mubeta),S=as.numeric(S2)*vbeta,df=v)
#betaF[j,]=rmvnorm(1,mean=mubeta,sigma=as.numeric(S2)*vbeta)

#sigmaaux2=rinvchisq(1,next.nuF,scale=S2)
#sigmaF[j]=rinvchisq(1,v,scale=S2)
sigmaF[j]=rinvchisq(1,(v+ 2*(a-1)),scale=4*S2)
phiF[j]=next.pi2F[2]
covinii=c(sigmaF[j],phiF[j])
Psi1=varcov.spatial(H=H,cov.model=cov.model,cov.pars=covinii,nugget=tau2,kappa=kappa)
Psi1=(Psi1+t(Psi1))/2

#print(c(betaF[j,],sigmaF[j],phiF[j],j))
countiter=countiter+1
cat("Iteration ",countiter," of ",iter,"\r")
}
#}

betaburn=as.matrix(betaF[(burn+1):iter,])
betaval=as.matrix(betaburn[seq((burn+1),iter-burn,thin),])
phiburn=phiF[(burn+1):iter]
phival=phiburn[seq((burn+1),iter-burn,thin)]
sigmaburn=sigmaF[(burn+1):iter]
sigmaval=sigmaburn[seq((burn+1),iter-burn,thin)]



if(method=="mode"){
  modebeta=apply(betaval,2,mlv)
  modephi=mlv(phival)
  modesigma=mlv(sigmaval)
  theta=c(modebeta,modesigma,modephi)
}


if(method=="mean"){
  modebeta=apply(betaval,2,mean)
  modephi=mean(phival)
  modesigma=mean(sigmaval)
  theta=c(modebeta,modesigma,modephi)
}


if(method=="median"){
  modebeta=apply(betaval,2,median)
  modephi=median(phival)
  modesigma=median(sigmaval)
  theta=c(modebeta,modesigma,modephi)
}

dist=cbind(betaval,sigmaval,phival)
return(list(prob=count/iter,dist=dist,betaF=betaval,sigmaF=sigmaval,phiF=phival,coords=coords,kappa=kappa,X=xmat,type=cov.model,theta=c(modebeta,modesigma,modephi),y=ycomp))
}





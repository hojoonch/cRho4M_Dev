function [misfit,sig,prho,ro]=pupd(data,J,par,yky,tri,es,akel,tev,k1,indx,V1,prho,npar,dd,so,p,C,lambda,Rd)
%Updating model parameters

%Damping factor
while lambda<0.01
    lambda=0.01;
end
%smoothness constrained least squares
%smoothness constrain is a second order laplacian

b=(J'*Rd'*Rd*dd-lambda*C*((1./prho(:))));
A=(J'*Rd'*Rd*J+lambda*C);
dp=A\b;
parg=1./((1./prho(:)).*exp(dp));
rhoort=exp(sum(log(parg)./length(parg)));
sigtmp(1:es)=1./(rhoort);
for s=1:npar
    sigtmp(par(s).ucg)=1./parg(s);
end
% Test the updated model
[J,rog]=forward(yky,tri,es,sigtmp,so,data.nel,akel,0,tev,k1,indx,V1,data,prho,npar,par,p);
misfitg=sqrt((Rd*dd)'*(Rd*dd)/data.nd)*100;

misfit=misfitg;
ro=rog;
sig=sigtmp;
prho=parg;


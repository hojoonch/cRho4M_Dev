function [prho,nu]=pupd_ip(J,par,es,prho,npar,dd,C,lambda,Rd)

while lambda<0.15
    lambda=0.15;
end
b=(J'*Rd'*Rd*dd-lambda*C*((prho(:))));
A=(J'*Rd'*Rd*J+lambda*C);
dp=A\b;
parg=prho(:)+(dp);

%log mu deðil mi
nuort=(sum((parg)./length(parg)));
% sigtmp(1:es)=1./rhoort;
% for s=1:npar
%     sigtmp(par(s).ucg)=1./parg(s);
% end
% Test the updated model
% [J,JIP,rog]=forward(yky,tri,es,sigtmp,so,data.nel,akel,0,tev,k1,indx,V1,data,prho,npar,par,p);


% [J,JIP,ro]=forward(yky,t,es,ssig,so,data.nel,akel,1,tev,k1,indx,V1,data,prhoR,npar,par,p);
%         ssig=sig.*(1-nug);
%         pprho=1./(1./prho(:)).*(1-pma(:)/1000);
%         [J1,JIP1,ron]=forward(yky,t,es,ssig,so,data.nel,akel,1,tev,k1,indx,V1,data,prho,npar,par,p);
%         mac=1000*(ron-roR)./ron;

% misfitg=sqrt((Rd*dd)'*(Rd*dd)/data.nd)*100;

% misfit=misfitg;
% ro=rog;
% sig=sigtmp;


homma=nuort;
nu(1:es)=homma;
for s=1:npar
    nu(par(s).ucg)=parg(s);
end

prho=abs(parg);


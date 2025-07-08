function [sig,es,ds,akel,V1,k1,prho,so,x,pma,nu]=initial(t,p,data,yky,npar)
% Initialize model and calculate the mesh response for 1ohm-m homogenous
% space
ds=length(p);
es=length(t);
I=2*pi;

[delta,b1,c1,b2,c2,b3,c3]=pdetrgm(p,t);
try
    [akel]=knnsearch(p',[data.xelek(:) data.zelek(:)]);
catch
    [tmp1,akel,tmp2]=    intersect(p',[data.xelek(:) data.zelek(:)],'rows');
end
so=spalloc(ds,data.nel,data.nel);
for ael=1:data.nel
    so(akel(ael),ael)=I;
end

if data.ip
    homma=sum(abs(data.ma))/data.nd;
    pma(1:npar)=homma;
    nu(1:es)=homma;
else
    homma=[];
    pma=[];
    nu=[];
end
sig(1:es)=(1./data.homro);


prho(1:npar)=data.homro;
[V1,k1,x]=k1f(t,es,delta,b1,b2,b3,c1,c2,c3,so,data.nel,akel,yky);


function [J,ro]=forward(yky,t,es,sig,so,nel,akel,opt,tev,k1,x,V1,data,prho,npar,par,p)
% Forward calculation routine
% This function calculates and returns the response of a 2D resistivity
% model and the Jacobian matrix J.
tt=zeros(nel,nel,npar);
xx2=zeros(nel,nel);
[delta,b1,c1,b2,c2,b3,c3]=pdetrgm(p,t);
p1_1=sig./(4.*delta);
el=1:es;
%% Finite Elements
for nky=1:length(yky)
    a=sig.*(1/6)*yky(nky).^2.*delta;
    b=a/2;
    k11(el)=p1_1.*((b1(el).^2+c1(el).^2))+a;
    k11(el)=p1_1.*((b1(el).^2+c1(el).^2))+a;
    k12(el)=p1_1.*((b1(el).*b2(el)+c1(el).*c2(el)))+b;
    k13(el)=p1_1.*((b1(el).*b3(el)+c1(el).*c3(el)))+b;
    k21(el)=k12(el);
    k22(el)=p1_1.*((b2(el).^2+c2(el).^2))+a;
    k23(el)=p1_1.*(b2(el).*b3(el)+c2(el).*c3(el))+b;
    k31(el)=k13(el);
    k32(el)=k23(el);
    k33(el)=p1_1.*(b3(el).^2+c3(el).^2)+a;
    k11=sparse([k11,k12,k13,k21,k22,k23,k31,k32,k33]);
    K1=accumarray(x,k11);
    fi1=K1\so;
    for ael=1:nel
        yuzpot1=fi1(akel,ael);
        xx1(:,ael)=full(yuzpot1');
    end
end
V0=xx1*yky(1)/pi;
%% Calculate apparent resistivities
[ro]=go2d(V0,V1,data.pat,data.eldiz);
%% Calculate Jacobian Matrix if it is requested
if opt==1
    for pno=1:npar
        ccd=sparse(length(fi1),length(fi1));
        uc=par(pno).ucg;
        for m=1:length(uc)
            CCD=sparse(length(fi1),length(fi1));
            d=uc(m):es:length(t)*9;
            cevap=full(k1(d));
            for id=1:9
                CCD(x(d(id),1),x(d(id),2))=cevap(id);
            end
            ccd=ccd+CCD;
        end
        uc=[];
        s1=-ccd*fi1;
        fitur=K1\s1;
        tt(:,:,pno)=(fitur(akel,:))';
    end
end

if opt==1
    VVT=tt.*yky(1)/pi;
    J=jacob(VVT,V1,data,ro,prho,npar);
else
    J=0;
end


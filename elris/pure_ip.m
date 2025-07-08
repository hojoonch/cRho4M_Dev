function [chamod,misfit,chaobs,iternum]=pure_ip(data,roR,sigR,JR,prhoR,C,es,akel,V1,k1,so,indx,pma,nu,tev,par,p,t,npar,Rd)
tic

itmax = 100;%
yky=1/data.zmax;%
lambda=std(log(abs(data.ma)));
sd=abs(data.ma).^-.0025;%sd=sd./max(sd);%
% Rd=diag((sd));
Rd=diag((sd));
pma=pma/1000;
nu=nu/1000;
for iter=1:itmax
    % Forward operator
    ssig=sigR.*(1-nu);
    pprho=prhoR(:).*(1-pma(:));
    [~,ron]=forward(yky,t,es,ssig,so,data.nel,akel,0,tev,k1,indx,V1,data,pprho,npar,par,p);
    mac=(ron-roR)./ron;
    
    for ii=1:length(ron)
        for jj=1:npar
            J(ii,jj)=(-pprho(jj))*(roR(ii)./ron(ii).^2)*JR(ii,jj);

        end
    end
    
    dd=(data.ma(:)/1000)-(mac(:));
    misfit1=sqrt((Rd*dd)'*(Rd*dd)/data.nd)*100;
    

%     misfit1=100*sqrt(sum((data.ma(:)/1000-(mac(:))).^2)/numel(mac))
    % Parameter update
    [pma,nu]=pupd_ip(J,par,es,pma,npar,dd,lambda*C,lambda,Rd);
    pprho=1./(1./prhoR(:)).*(1-abs(pma(:)));
    [~,ron]=forward(yky,t,es,ssig,so,data.nel,akel,0,tev,k1,indx,V1,data,pprho,npar,par,p);
    mac=(ron-roR)./ron;
    mfit(iter)=misfit1;
    
    
    cha(iter,:)=pma;
    ocha(:,iter)=mac;
    %     cizro=pma';
    

    [misfit,yer]=min(mfit);
    chamod=cha(yer,:)';
    chaobs=ocha(:,yer);
    iternum=yer;
    if iter>1
        if mfit(iter)>mfit(iter-1)            
            break            
        end
        
        if iter>=2
            lambda=lambda*.65;
        end
    end
end
toc
iter



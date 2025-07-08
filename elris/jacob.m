function J=jacob(VVT,V1,data,ro,prho,npar)
for i=1:npar
    M=VVT(:,:,i);
    J(:,i)=go2d(M,V1,data.pat,data.eldiz);
end

for ii=1:length(ro)
    for jj=1:npar
        J(ii,jj)=J(ii,jj)/prho(jj)/(ro(ii));
    end
end



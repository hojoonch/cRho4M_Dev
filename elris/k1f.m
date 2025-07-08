function [V1,k1,x]=k1f(t,es,delta,b1,b2,b3,c1,c2,c3,so,nel,akel,yky)
x1=t(1,:)';x2=t(2,:)';x3=t(3,:)';
x=[x1,x1;x1,x2;x1,x3;x2,x1;x2,x2;x2,x3;x3,x1;x3,x2;x3,x3];
p1_1=1./(4.*delta);
el=1:es;
for nky=1:length(yky)
    % p2_1=;
    a=(1/6)*yky(nky).^2.*delta;%p2_1;
    b=a/2;
    k11(el)=p1_1.*((b1(el).^2+c1(el).^2))+a;
    k11(el)=p1_1.*((b1(el).^2+c1(el).^2))+a;
    k12(el)=p1_1.*((b1(el).*b2(el)+c1(el).*c2(el)))+b;
    k13(el)=p1_1.*((b1(el).*b3(el)+c1(el).*c3(el)))+b;
    k21(el)=k12(el);%p1*(b(1)*b(2)+c(1)*c(2))+(1/12)*p2;
    k22(el)=p1_1.*((b2(el).^2+c2(el).^2))+a;
    k23(el)=p1_1.*(b2(el).*b3(el)+c2(el).*c3(el))+b;
    k31(el)=k13(el);%p1*(b(1)*b(3)+c(1)*c(3))+(1/12)*p2;
    k32(el)=k23(el);%p1*(b(2)*b(3)+c(2)*c(3))+(1/12)*p2;
    k33(el)=p1_1.*(b3(el).^2+c3(el).^2)+a;
    k1=sparse([k11,k12,k13,k21,k22,k23,k31,k32,k33]);
    K1=accumarray(x,k1);
    fi1=K1\so;
    for ael=1:nel
        yuzpot1=fi1(akel,ael);
        xx1(:,ael)=full(yuzpot1');
    end
    
end

V1=xx1*yky(1)/pi/1.5;




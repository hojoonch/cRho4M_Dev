function [ro]=go2d(V0,V1,IND,eldiz)

% pause
ro=zeros(1,length(IND.ind1));
switch eldiz
    case 1
        dv0=V0(IND.ind1)-V0(IND.ind2)-V0(IND.ind3)+V0(IND.ind4);
        dv1=V1(IND.ind1)-V1(IND.ind2)-V1(IND.ind3)+V1(IND.ind4);
        ro=dv0./dv1;
    case 2
        dv0=V0(IND.ind1);
        dv1=V1(IND.ind1);
        ro=dv0./dv1;
        
    case 3
        dv0=V0(IND.ind1)-V0(IND.ind2)-V0(IND.ind3)+V0(IND.ind4);
        dv1=V1(IND.ind1)-V1(IND.ind2)-V1(IND.ind3)+V1(IND.ind4);
        ro=dv0./dv1;       
    case 6
        dv0=V0(IND.ind1)-V0(IND.ind2);
        dv1=V1(IND.ind1)-V1(IND.ind2);
        ro=dv0./dv1;        
    case 7
        dv0=V0(IND.ind1)-V0(IND.ind2)-V0(IND.ind3)+V0(IND.ind4);
        dv1=V1(IND.ind1)-V1(IND.ind2)-V1(IND.ind3)+V1(IND.ind4);
        ro=dv0./dv1;
    case 11
         dv0=V0(IND.ind1)-V0(IND.ind2);
        dv1=V1(IND.ind1)-V1(IND.ind2);
        ro=dv0./dv1;
        
end


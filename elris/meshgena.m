function [p,t,nlay,tev,par,npar,z,xel,nxg,nzg]=meshgena(data)
% Mesh generator. 'fine' mode is selected. Model and outer space is constructed bey
% dividing rectangular blocks into two triangles. 
dz(1)=data.dz1/2;
for i=2:1000
    if sum(dz)<=data.zmax*1.1
        zk=dz(i-1)*1.1;
        dz(i)=zk;
        nlay=i;
    end
end

z=cumsum(dz);
oran=max(z)/data.zmax;
z=z/oran;
z=-[0 z];
xel=zeros(1,2*data.nel-1);
xelek=data.xelek;
fark=diff(xelek);
xel(1:2:end)=xelek;
xel(2:2:end)=xelek(1:end-1)+fark/2;
x1=xel(1)-[6*data.ela 3*data.ela data.ela];
x2=xel(end)+[ data.ela  3*data.ela 6*data.ela];
x_u= [x1 xel x2];
z_u=[z 1.25*z(end) 1.875*z(end) 2.8125*z(end)];  
nz=length(z_u);
[xa,za]=meshgrid(x_u,z_u);
tri(1,:)=[(nz+1) 1 2];
tri(2,:)=[(nz+1) 2 nz+2];
for i=3:2*(nz-1)
    tri(i,:)=tri(i-2,:)+1;
end
for j=2:(2*data.nel-2)+6
    sb=(j-1)*2*(nz-1)+1;
    ss=sb+2*(nz-1)-1;
    ib=sb-2*(nz-1);
    is=ib+2*(nz-1)-1;
    tri(sb:ss,:)=tri(ib:is,:)+(nz);
end
p=[xa(:)'; za(:)'];
nx=(2*data.nel-2);
npar=nlay*nx;
sira=reshape(1:npar,nlay,2*(data.nel-1));
clear tev
basx=1;
say=1;
bas=3*(nz-1)*2;
for k=1:nx
    basy=1;
    for j=1:nlay
        tev(say).par=say;
        yer=bas+(k-1)*(nz-1)*2+(j-1)*2+1;
        par(say).ucg=yer:yer+1;
        say=say+1;
    end
end
t=tri';
nxg=nx;
nzg=nlay;



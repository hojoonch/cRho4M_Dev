function data=read_data(varargin)
% Data input function for ELRIS
% Last revision 21.06.2014
% Able to read general array format
global pathname
if nargin==0
    if isempty(pathname)
        pathname=pwd;
    end

    fid=fopen([pathname filename]);

else
    filename=cell2mat(varargin);
    fid=fopen(filename);
end

if fid<0
    printf('Specified file cannot be found : %s \n', filename)
    data=[];
    return
elseif ischar(fid)
    data=[];
    pathname=pwd;
    return
end

data.prfadi=fgetl(fid);
data.ela=fscanf(fid,'%f',1);
data.eldiz=fscanf(fid,'%d',1);%1: Wenner, 2:Pole-pole 3: Dipole-dipole 6: Pole-dipole 7: Wenner-Schlumberger 11: Mixed
if data.eldiz==11
    data.subeldiz=fscanf(fid,'%d',1);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
end
data.nd=fscanf(fid,'%d',1);
data.mp=fscanf(fid,'%d',1);%MIDPOINT
data.ip=fscanf(fid,'%d',1);%IP

data=oku(data,fid); %Call data reader function
topog=fscanf(fid,'%d',1);
data.topog=topog;
if topog
    topsay=fscanf(fid,'%d',1);
    data.topo = fscanf(fid, '%g %g', [2 topsay])';    % Read topography data
else
    fscanf(fid,'%d',5) ;
    fscanf(fid,'%s',2);
    data.sd=fscanf(fid, '%g ', [1 data.nd])';
end
switch data.eldiz
    case 1 %wenner
        data.eldizc='Wenner';
        switch data.mp
            case 1
                tmp=data.nlev*data.ela;
                xc1=data.xd-1.5*tmp;
                xc2=data.xd+1.5*tmp;
                xp1=xc1+tmp;
                xp2=xp1+tmp;
                xelek=unique([xc1;xc2;xp1;xp2]);
                data.pat.kon=[xc1 xc2 xp1 xp2];
                data.xd=(xc1+xc2)/2;

            case 0
                tmp=data.nlev*data.ela;
                xc1=data.xd;
                xp1=xc1+tmp;
                xp2=xp1+tmp;
                xc2=xp2+tmp;
                xelek=unique([xc1;xc2;xp1;xp2]);
                data.pat.kon=[xc1 xc2 xp1 xp2];
                data.xd=(xc1+xc2)/2;
        end
        for j=1:data.nd
            data.pat.esira(j,:)=[find(xelek==xc1(j)) find(xelek==xc2(j)) find(xelek==xp1(j)) find(xelek==xp2(j))];%Elektrotlar�n s�ra numaralar�
        end
        data.nel=length(xelek);
        data.xelek=xelek;
        data.zmax=data.ela*(.519*max(data.nlev));
        data.psd=-data.ela*(.519*(data.nlev));

        data.pat.ind1=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,3));
        data.pat.ind2=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,4));
        data.pat.ind3=sub2ind([data.nel data.nel],data.pat.esira(:,2),data.pat.esira(:,3));
        data.pat.ind4=sub2ind([data.nel data.nel],data.pat.esira(:,2),data.pat.esira(:,4));
        data.dz1=0.519*data.ela;

    case 2 % Pol-Pol
%---------------------------------
         data.eldizc='Pol-Pol';
        switch data.mp
            case 1
                tmp=data.nlev*data.ela;
                xc1=data.xd-.5*tmp;
                xp1=xc1+tmp;

                xelek=unique([xc1;xp1]);
                data.pat.kon=[xc1 xp1];
            case 0
                tmp=data.nlev*data.ela;
                xc1=data.xd;
                xp1=xc1+tmp;
                xp2=xp1+tmp;
                xc2=xp2+tmp;
                xelek=unique([xc1;xc2;xp1;xp2]);
                data.pat.kon=[xc1 xc2 xp1 xp2];
                data.xd=(xc1+xc2)/2;
        end
        for j=1:data.nd
            data.pat.esira(j,:)=[find(xelek==xc1(j)) find(xelek==xp1(j))];%Electrode index
        end
        data.nel=length(xelek);
        data.xelek=xelek;
        data.zmax=data.ela*(.867*max(data.nlev));
        data.psd=-data.ela*(.867*(data.nlev));

        data.pat.ind1=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,2));
        data.dz1=0.867*data.ela;






%---------------------------------
    case 3 % Dipole dipole
        data.eldizc='Dipole-Dipole';
        xc1=data.xd-(data.nlev.*data.ela)/2;
        xc2=xc1-data.ela;
        xp1=xc1+data.nlev.*data.ela;
        xp2=xp1+data.mn;
        xelek=unique([xc1;xc2;xp1;xp2]);
        data.pat.kon=[xc1 xc2 xp1 xp2];
        for j=1:data.nd
            data.pat.esira(j,:)=[find(xelek==xc1(j)) find(xelek==xc2(j)) find(xelek==xp1(j)) find(xelek==xp2(j))];%Elektrotlar�n s�ra numaralar�
        end
        data.nel=length(xelek);
        data.pat.ind1=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,3));
        data.pat.ind2=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,4));
        data.pat.ind3=sub2ind([data.nel data.nel],data.pat.esira(:,2),data.pat.esira(:,3));
        data.pat.ind4=sub2ind([data.nel data.nel],data.pat.esira(:,2),data.pat.esira(:,4));
        data.zmax=data.ela*(.2587*max(data.nlev)+.1759);
        data.xelek=xelek;
        data.psd=-data.ela*(.2587*(data.nlev)+.1759);
        data.dz1=0.42*data.ela;

    case 6 %Pole Dipole
        data.eldizc='Pole-Dipole';
        switch data.mp
            case 1
                xc1=data.xd-((data.nlev+1).*data.ela)/2;
                xp1=xc1+data.nlev.*data.ela;
                xp2=xp1+data.mn;
                xelek=unique([xc1;xp1;xp2]);
                data.pat.kon=[xc1 xp1 xp2];%�l��m s�ras�nda elektrotlar�n
                %         konumland�r�ld�klar� noktan�n x koordinatlar�
            case 0
                xc1=data.xd;
                xp1=xc1+data.nlev.*data.ela;
                xp2=xp1+data.mn;
                xelek=unique([xc1;xp1;xp2]);
                data.pat.kon=[xc1 xp1 xp2];%�l��m s�ras�nda elektrotlar�n
                data.xd=(xc1+xp2)/2;
        end

        for j=1:data.nd
            data.pat.esira(j,:)=[find(xelek==xc1(j)) find(xelek==xp1(j)) find(xelek==xp2(j))];%Elektrotlar�n s�ra numaralar�
        end
        data.nel=length(xelek);
        data.xelek=xelek;
        data.pat.ind1=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,2));
        data.pat.ind2=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,3));
        data.zmax=data.ela*(.3874*max(data.nlev)+.1489);
        data.psd=-data.ela*(.3874*(data.nlev)+.1489);
        data.dz1=0.41*data.ela;
        data.pat.k=2*pi*data.nlev(data.nlev+1)*data.ela;


    case 7 %Wenner Schlumberger
        data.eldizc='Schlumberger';
        switch data.mp
            case 1
                xc1=data.xd-((data.nlev+.5).*data.ela);
                xc2=data.xd+((data.nlev+.5).*data.ela);
                xp1=xc1+data.nlev.*data.ela;
                xp2=xp1+data.mn;
                xelek=unique([xc1;xc2;xp1;xp2]);
                data.pat.kon=[xc1 xc2 xp1 xp2];
                data.xd=(xc1+xc2)/2;

            case 0
                xc1=data.xd;
                xp1=xc1+data.nlev.*data.ela;
                xp2=xp1+data.mn;
                xc2=xp2+data.nlev.*data.ela;
                xelek=unique([xc1;xc2;xp1;xp2]);
                data.pat.kon=[xc1 xc2 xp1 xp2];
                data.xd=(xc1+xc2)/2;
        end

        for j=1:data.nd
            data.pat.esira(j,:)=[find(xelek==xc1(j)) find(xelek==xc2(j)) find(xelek==xp1(j)) find(xelek==xp2(j))];%Elektrotlar�n s�ra numaralar�
        end
        data.nel=length(xelek);
        data.xelek=xelek;
        data.zmax=data.ela*(.3874*max(data.nlev)+.1489);
        data.psd=-data.ela*(.3874*(data.nlev)+.1489);
        data.pat.ind1=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,3));
        data.pat.ind2=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,4));
        data.pat.ind3=sub2ind([data.nel data.nel],data.pat.esira(:,2),data.pat.esira(:,3));
        data.pat.ind4=sub2ind([data.nel data.nel],data.pat.esira(:,2),data.pat.esira(:,4));
        data.dz1=0.5*data.ela;
    case 11
        clc
    otherwise
        errordlg('Not a supported array type!','modal')
end
data.homro=exp(mean(log(data.roa)));%exp(sum(log(data.roa)/length(data.roa)));%%
if data.topog
    data.zelek=pchip(data.topo(:,1),data.topo(:,2),data.xelek);
else
    data.zelek=zeros(size(data.xelek));
end
data.filename=filename;



function data=oku(data,fid)
if data.eldiz==6||data.eldiz==3||data.eldiz==7
    if data.ip
        fscanf(fid,'%s',1);
        fscanf(fid,'%s',1);

        fscanf(fid,'%f%c%f',[1 3]);
        veri = fscanf(fid, '%g %g %g %g %g', [5 data.nd])';    % Read apparent resistivity & chargeability
        data.xd=veri(:,1);
        data.mn=veri(:,2);
        data.nlev=veri(:,3);
        data.roa=veri(:,4);
        data.ma=veri(:,5);
    else
        veri = fscanf(fid, '%g %g %g %g %g', [4 data.nd])';    % Read apparent resistivity
        data.xd=veri(:,1);
        data.mn=veri(:,2);
        data.nlev=veri(:,3);
        data.roa=veri(:,4);
    end
elseif data.eldiz==1
    if data.ip
        fscanf(fid,'%s',1);
        fscanf(fid,'%s',1);

        fscanf(fid,'%f%c%f',[1 3]);
        veri = fscanf(fid, '%g %g %g %g', [4 data.nd])';     % Read apparent resistivity & chargeability
        data.xd=veri(:,1);
        data.mn=veri(:,2);
        data.nlev=veri(:,2)./data.ela;
        data.roa=veri(:,3);
        data.ma=veri(:,4);
    else
        veri = fscanf(fid, '%g %g %g %g %g', [3 data.nd])';    % Read apparent resistivity
        data.xd=veri(:,1);
        data.mn=veri(:,2)*data.ela;
        data.nlev=veri(:,2);
        data.roa=veri(:,3);
    end
elseif data.eldiz==2

    if data.ip
        fscanf(fid,'%s',1);
        fscanf(fid,'%s',1);

        fscanf(fid,'%f%c%f',[1 3]);
        veri = fscanf(fid, '%g %g %g %g', [4 data.nd])';     % Read apparent resistivity & chargeability
        data.xd=veri(:,1);
        data.mn=veri(:,2);
        data.nlev=veri(:,2)./data.ela;
        data.roa=veri(:,3);
        data.ma=veri(:,4);
    else
        veri = fscanf(fid, '%g %g %g', [3 data.nd])';    % Read apparent resistivity
        data.xd=veri(:,1);
        data.mn(1:data.nd)=9999;%veri(:,2)*data.ela;
        data.nlev=veri(:,2);
        data.roa=veri(:,3);
    end
elseif data.eldiz==11
    for k=1:data.nd
        ne=fscanf(fid,'%d',1);
        xc1(k)=fscanf(fid,'%f',1);
        zc1(k)=fscanf(fid,'%f',1);
        if ne==4
            xc2(k)=fscanf(fid,'%f',1);
            zc2(k)=fscanf(fid,'%f',1);
        else
            xc2(k)=NaN;
            zc2(k)=NaN;
        end

        xp1(k)=fscanf(fid,'%f',1);
        zp1(k)=fscanf(fid,'%f',1);
        if ne>2
            xp2(k)=fscanf(fid,'%f',1);
            zp2(k)=fscanf(fid,'%f',1);
        else
            xp2(k)=NaN;
            zp2(k)=NaN;
        end
        roa(k)=fscanf(fid,'%f',1);
    end
    switch data.subeldiz
        case 3
            data.eldizc='Dipol-Dipol';
            xelek=unique([xc1;xc2;xp1;xp2]);
            data.pat.kon=[xc1 xc2 xp1 xp2];
            for j=1:data.nd
                data.pat.esira(j,:)=[find(xelek==xc1(j)) find(xelek==xc2(j)) find(xelek==xp1(j)) find(xelek==xp2(j))];%Elektrotlar�n s�ra numaralar�
            end
            a2=(xp1-xc2)/data.ela;
            data.nlev=a2;%./a1;
            data.nel=length(xelek);
            data.pat.ind1=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,3));
            data.pat.ind2=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,4));
            data.pat.ind3=sub2ind([data.nel data.nel],data.pat.esira(:,2),data.pat.esira(:,3));
            data.pat.ind4=sub2ind([data.nel data.nel],data.pat.esira(:,2),data.pat.esira(:,4));
            data.zmax=data.ela*(.2587*max(data.nlev)+.1759);
            data.xelek=xelek;
            data.psd=-data.ela*(.2587*(data.nlev)+.1759);
            data.dz1=0.42*data.ela;
            data.roa=roa;
            data.xd=(xc1+xp2)/2;

        case 6
            xelek=unique([xc1;xp1;xp2]);
            data.pat.kon=[xc1 xp1 xp2];
            for j=1:data.nd
                data.pat.esira(j,:)=[find(xelek==xc1(j)) find(xelek==xc2(j)) find(xelek==xp1(j)) find(xelek==xp2(j))];%Elektrotlar�n s�ra numaralar�
            end
            a2=(xp1-xc1)/data.ela;
            data.nlev=a2;%./a1;
            data.nel=length(xelek);
            data.xelek=xelek;
            data.pat.ind1=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,2));
            data.pat.ind2=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,3));
            data.zmax=data.ela*(.3874*max(data.nlev)+.1489);
            data.psd=-data.ela*(.3874*(data.nlev)+.1489);
            data.dz1=0.41*data.ela;
            data.roa=roa;
            data.xd=(xc1+xp1)/2;
        case 7
            xelek=unique([xc1;xc2;xp1;xp2]);
            data.pat.kon=[xc1 xc2 xp1 xp2];
            a2=(xp1-xc1)/data.ela;
            data.nlev=a2;%./a1;
            for j=1:data.nd
                data.pat.esira(j,:)=[find(xelek==xc1(j)) find(xelek==xc2(j)) find(xelek==xp1(j)) find(xelek==xp2(j))];%Elektrotlar�n s�ra numaralar�
            end
            data.nel=length(xelek);
            data.zmax=data.ela*(.3874*max(data.nlev)+.1489);
            data.psd=-data.ela*(.3874*(data.nlev)+.1489);
            data.pat.ind1=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,3));
            data.pat.ind2=sub2ind([data.nel data.nel],data.pat.esira(:,1),data.pat.esira(:,4));
            data.pat.ind3=sub2ind([data.nel data.nel],data.pat.esira(:,2),data.pat.esira(:,3));
            data.pat.ind4=sub2ind([data.nel data.nel],data.pat.esira(:,2),data.pat.esira(:,4));
            data.dz1=0.5*data.ela;
            data.roa=roa;
            data.xd=(xc1+xc2)/2;


    end
end
sd=1./(data.roa).^.25;
sd=sd./max(sd);
Rd=diag(sd);
data.Rd=Rd;
% w=1./data.roa.^.25;
% data.Rd=diag(w);
% fclose (fid);

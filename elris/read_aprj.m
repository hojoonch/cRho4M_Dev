function data=read_aprj(varargin)
% Data input function for ELRIS
% Last revision 21.06.2014
% Able to read general array format
global pathname
if nargin==0
    data = [];
    return
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


% read json
%dataJ = load_json(fileread(filename));
%disp('read json using jsondecode')
dataJ = jsondecode(fileread(filename));

% for test
%dataJ.arrayType = 4000;

data.prfadi= [dataJ.area ':' dataJ.line];
data.ela=dataJ.elspacing;

at_conv = [4000 3000 0 -1 -1 1000 5000 -1 -1 6000 9000];
data.eldiz=find(at_conv == dataJ.arrayType);
data.nd= cal_ndata(data.eldiz, dataJ.head03_1, dataJ.head03_2);
data.mp=1;%MIDPOINT
data.ip=0;%IP

data=okuJ(data,dataJ,fid); %Call data reader function

topog = false; % skip reading topo data
data.topog = topog;

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
            data.pat.esira(j,:)=[find(xelek==xc1(j)) find(xelek==xc2(j)) find(xelek==xp1(j)) find(xelek==xp2(j))];
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
            data.pat.esira(j,:)=[find(xelek==xc1(j)) find(xelek==xc2(j)) find(xelek==xp1(j)) find(xelek==xp2(j))];
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
                data.pat.kon=[xc1 xp1 xp2];

            case 0
                xc1=data.xd;
                xp1=xc1+data.nlev.*data.ela;
                xp2=xp1+data.mn;
                xelek=unique([xc1;xp1;xp2]);
                data.pat.kon=[xc1 xp1 xp2];
                data.xd=(xc1+xp2)/2;
        end

        for j=1:data.nd
            data.pat.esira(j,:)=[find(xelek==xc1(j)) find(xelek==xp1(j)) find(xelek==xp2(j))];
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
            data.pat.esira(j,:)=[find(xelek==xc1(j)) find(xelek==xc2(j)) find(xelek==xp1(j)) find(xelek==xp2(j))];
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

     case 10 % modified pole pole
        data.eldizc='Modified Pole Pole';
        xc1=data.xd-(data.nlev.*data.ela)/2;
        %xc2=xc1-data.ela; %
        xc2 = xc1*0.0 + dataJ.sta0;
        xp1=xc1+data.nlev.*data.ela;
        %xp2=xp1+data.mn; %
        xp2 = xc2 + (dataJ.head03_2 + 2)*dataJ.elspacing;
        xelek=unique([xc1;xc2;xp1;xp2]);
        data.pat.kon=[xc1 xc2 xp1 xp2];
        for j=1:data.nd
            data.pat.esira(j,:)=[find(xelek==xc1(j)) find(xelek==xc2(j)) find(xelek==xp1(j)) find(xelek==xp2(j))];
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

    case 11
        clc
        errordlg('Mixed array : is not supported yet!','modal')
    otherwise
        errordlg('Not a supported array type!','modal')
end

% for test
% [xc1 xc2 xp1 xp2]

% remove negative roa
indP = find(data.roa > 0);
data.indP = indP;
data.nd = length(indP);
data.xd = data.xd(indP);
data.mn = data.mn(indP);
data.nlev = data.nlev(indP);
data.roaori=data.roa;
data.roa = data.roa(indP);

data.homro=exp(mean(log(data.roa)));%exp(sum(log(data.roa)/length(data.roa)));%%
if data.topog
    data.zelek=pchip(data.topo(:,1),data.topo(:,2),data.xelek);
else
    data.zelek=zeros(size(data.xelek));
end
data.filename=filename;


% 1: Wenner, 2:Pole-pole 3: Dipole-dipole 6: Pole-dipole 7: Wenner-Schlumberger 11: Mixed
function data=okuJ(data, dataJ, fid)
  veri = makeR2dinvData(dataJ);

if data.eldiz==6||data.eldiz==3||data.eldiz==7 || data.eldiz == 10 %P-Dp, Dp-Dp, WSch
    if data.ip
        %veri = fscanf(fid, '%g %g %g %g %g', [5 data.nd])';    % Read apparent resistivity & chargeability
        data.xd=veri(:,1);
        data.mn=veri(:,2);
        data.nlev=veri(:,3);
        data.roa=veri(:,4);
        data.ma=veri(:,5);
    else
        %veri = fscanf(fid, '%g %g %g %g %g', [4 data.nd])';    % Read apparent resistivity
        data.xd=veri(:,1);
        data.mn=veri(:,2);
        data.nlev=veri(:,3);
        data.roa=veri(:,4);
    end
elseif data.eldiz==1 % Wenner
    if data.ip
        data.xd=veri(:,1);
        data.mn=veri(:,2);
        data.nlev=veri(:,2)./data.ela;
        data.roa=veri(:,3);
        data.ma=veri(:,4);
    else
        data.xd=veri(:,1);
        data.mn=veri(:,2);
        data.nlev=veri(:,3);
        data.roa=veri(:,4);
    end
elseif data.eldiz==2 % pole pole
    if data.ip
        data.xd=veri(:,1);
        data.mn=veri(:,2);
        data.nlev=veri(:,2);
        data.roa=veri(:,3);
        data.ma=veri(:,4);
    else
        data.xd=veri(:,1);
        data.mn(1:data.nd)=9999;%veri(:,2)*data.ela;
        data.nlev=veri(:,3);
        data.roa=veri(:,4);
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

function ndata = cal_ndata(type, N, Nd1)
  %1: Wenner, 2:Pole-pole 3: Dipole-dipole 6: Pole-dipole 7: Wenner-Schlumberger 11: Mixed
  switch type
    case 1  % Wenner
      ndata = Nd1 * N - 3 * N * (N-1) / 2;
    case 2 % Pole-pole
      d1 = Nd1;
      d2 = d1 - N + 1;
      ndata = d1*(d1+1)/2 - (d2 - 1) * d2/2 ;
    case 3 % Dipole-dipole
      d1 = Nd1;
      d2 = d1 - N + 1;
      ndata = d1*(d1+1)/2 - (d2 - 1) * d2/2;
    case 6 % Pole-dipole
      d1 = Nd1;
      d2 = d1 - N + 1;
      ndata = d1*(d1+1)/2 - (d2 - 1) * d2/2 ;
    case 7 % Wenner-Schlimberger
      #ndata = Nd1 * N -  N * (N-1)
      ndata = (Nd1 - N + 1) * N;
   case 10 % modified pole pole
      d1 = Nd1;
      d2 = d1 - N + 1;
      ndata = d1*(d1+1)/2 - (d2 - 1) * d2/2 ;
    case 11 % Mixed
      d1 = Nd1;
      d2 = d1 - N + 1;
      ndata = d1*(d1+1)/2 - (d2 - 1) * d2/2 ;
  endswitch
  return

function veri = makeR2dinvData(dataJ)
  Rhoa01 = dataJ.Rhoa01;
  Rhoa11 = dataJ.Rhoa11;
  veri=[];
  if dataJ.arrayType == 5000  % schlumberger
    nf = 2;
  elseif dataJ.arrayType == 4000 % wenner
    nf = 3;
  else  % dpdp(0), p-p(3000), p-dp(1000)
    nf = 1;
  endif

  ind = 1;

  for i = 1 : double(dataJ.head03_1)
    for j = 1 : double(dataJ.head03_2) - (i-1)*nf
      if nf == 1
        c2 = j; % p-p, p-dp % C1=inf
        if dataJ.arrayType == 0 || dataJ.arrayType == 10 % dp-dp , mpp
          c2 = c2 + 1;
        endif
        p1 = c2 + i;

        mp = (c2+p1-2.0)*double(dataJ.elspacing)/2.0 + double(dataJ.sta0);
        mn = double(dataJ.elspacing);
        nn = i;

      elseif nf == 2 % schlumberger
        mp = (2*(j+i)-1) * double(dataJ.elspacing)/2.0 + double(dataJ.sta0);
        mn = dataJ.elspacing;
        nn = i;
      elseif nf == 3 % wenner
        mp = (2*(j-1) + 3*i) * double(dataJ.elspacing)/2.0 + double(dataJ.sta0);
        mn = i * dataJ.elspacing;
        nn = i;
      endif

      veri(ind,1) = mp;
      veri(ind,2) = mn;
      veri(ind,3) = nn;
      veri(ind,4) = Rhoa01(i,j);
      ind=ind+1;
    endfor
    #printf('\n')
  endfor

  return




function invert()%(varargin)
        %
        % rpi command line run
        %
        %

        args=argv();
        if length(args) ==3
           exepath = strtrim(args{1});
          datapath = strtrim(args{2});
          filename = strtrim(args{3});
        else
          exepath = pwd;
          datapath = pwd;
          filename = 'DC2023-09-21-15-03.aprj';
        endif

        % remove this part from 24. 9. 20
        % instead use builtin function ( >= v7.0) or octave package( < v7.0)
        if false
          cname = computer();
          rapidjson_dir = "";
          if strfind(cname, 'w64')
            % windows
            rapidjson_dir = strcat(exepath, '\rapidjson_win');
          elseif strfind(cname, 'arm')
            % raspberry pi
            rapidjson_dir = strcat(exepath , '/rapidjson_rpi');
				  elseif strfind(cname, 'pc-linux')
					  % ubuntu
					  rapidjson_dir = strcat(exepath , '/rapidjson_ubuntu');
				  elseif strfind(cname, 'apple')
					  % MacOS
					  rapidjson_dir = strcat(exepath , '/rapidjson_macos');
          end
          %rapidjson_dir
          addpath(rapidjson_dir);
        end
        % new codes for json

        octave_ver = cellfun(@str2double, strsplit(version(),"."));

        if octave_ver(1) >= 7
        else
          pkg load json;
        end

        fname = [datapath filesep filename];
        data = getResisitivityData(fname);

				itmax = 20;

        %if itmax==11
        %    itmax=15;
        %end

        mtype = 1;
        switch mtype
            case 0
                xa=1; za=1;
            case 1
                xa=2; za=1; % divides each cell into half
        end
        if ~isempty(data)
            alfax=1;
            alfaz=1;
            yky=1/data.zmax;%
            %             yky=1/((data.nel-1)*data.ela);
            %             yky=1/data.zmax;
            lambda=std(log(data.roa));
            % Mesh generator
            switch mtype
                case 0 %Fine mode selected
                    [p,t,nlay,tev,par,npar,z,xel,nx,nz]=meshgena(data);
                    parc=1:npar;
                    parc=reshape(parc,nlay,2*(data.nel-1));
                    parc=[parc;zeros(1,size(parc,2))];
                    parc=[zeros(size(parc,1),1),parc,zeros(size(parc,1),1)];
                    C=full(delsq(parc));
                    say=1;
                    for k=1:nx
                        for m=1:nz
                            yx1=(k-1)*xa+1;yx2=(k-1)*xa+xa+1;
                            yy1=(m-1)*za+1;yy2=(m-1)*za+za+1;
                            xp(say,:)=[xel(yx1) xel(yx2) xel(yx2) xel(yx1)];
                            zp(say,:)=[z(yy1) z(yy1) z(yy2) z(yy2)];
                            say=say+1;
                        end
                    end
                case 1 % Normal mode selected
                    [p,t,nlay,tev,par,npar,z]=meshgen(data,exepath,datapath );
                    parc=1:npar;
                    say=1;
                    for k=1:data.nel-1
                        for m=1:length(z)-1
                            xp(say,:)=[data.xelek(k) data.xelek(k+1) data.xelek(k+1) data.xelek(k)];
                            zp(say,:)=[z(m) z(m) z(m+1) z(m+1)];
                            say=say+1;
                        end

                    end
                    parc=reshape(parc,nlay,data.nel-1);
                    parc=[parc;zeros(1,size(parc,2))];
                    parc=[zeros(size(parc,1),1),parc,zeros(size(parc,1),1)];
                    C=full(delsq(parc));
            end

            [sig,es,ds,akel,V1,k1,prho,so,indx,pma,nu]=initial(t,p,data,yky,npar) ;

            sd=1./data.roa.^.025;
            %
            Rd=diag(sd);

            %tic
            disp('==== Inversion Start ====')
            for iter=1:itmax
                % Forward operator
                [J,ro]=forward(yky,t,es,sig,so,data.nel,akel,1,tev,k1,indx,V1,data,prho,npar,par,p);
                %dd=log(data.roa(:))-log(ro(:));
								ro = ro(data.indP);
								J = J(data.indP,:);
								dd=log(data.roa(:))-log(ro(:));
                misfit=sqrt((Rd*dd)'*(Rd*dd)/data.nd)*100;
                % Parameter update
                [misfit,sig,prho,ro]=pupd(data,J,par,yky,t,es,akel,tev,k1,indx,V1,prho',npar,dd,so,p,C,lambda,Rd);
                mfit(iter)=misfit;

								misfit

                g3 =1;
                switch g3
                    case 0
                        cizro=prho';
                    case 1
                        cizro=log10(prho');
                end

                if iter==1
                    alp=sum(abs(J),1);
                    alp=alp/max(alp);
                    alp1=repmat(alp,4,1);
                    alp1=(alp1(:));
                    alp1=alp1+(.91-min(alp1));
                    alp1(alp1>1)=1;
                else

                end

                oran=.2*(1/itmax)*iter;

                if iter>1
                    farkm=abs(mfit(iter)-mfit(iter-1))./mfit(iter);
                    if farkm<.025
                        break
                    end
                end
                if iter>=2
                    %lambda=lambda*.55; % original value
										lambda=lambda*.7;
                end
            end
            disp('==== Inversion Done ====')
            if data.ip
                [pma,misfit_ip,mac,iterx]=pure_ip(data,ro,sig,J,prho,C,es,akel,V1,k1,so,indx,pma,nu,tev,par,p,t,npar,Rd);
            end

            ndot = strfind(fname,'.');
            dadi = [fname(1:ndot) 'mat'];
						dadj = [fname(1:ndot) 'invj'];

            if data.ip==0
                %save (dadi,'data','xp','zp','prho','misfit','iter','ro','alp1','-mat')
								%roaori = data.roaori;
								%indP = data.indP;
								%save (dadi,'roaori','indP','xp','zp','prho','mfit','ro','-mat')
								outD = [data.xd data.mn data.nlev data.roa ro(data.indP)];
								model = [(xp(:,1)+xp(:,3))/2 xp(:,3)-xp(:,1) -(zp(:,1)+zp(:,3))/2 -(zp(:,3)-zp(:,1)) prho];
								%save (dadi,'model','outD','mfit','-mat')
								%jstr = save_json(struct("InvModel", model,"InOutData",outD,"misfit",mfit));
                %disp('save json using jsonencode')
                jstr = jsonencode(struct("InvModel", model,"InOutData",outD,"misfit",mfit));
								fid = fopen(dadj,'w+');
								fprintf(fid,jstr);
								fclose(fid);
            else
                save (dadi,'data','xp','zp','prho','misfit','iter','ro','alp1','pma','mac','misfit_ip','-mat')
            end

        end
end

%--------------------------------------------------------------------------
% getResisitivityData
%   This reads in all supported data files in the current directory
%--------------------------------------------------------------------------
    function record = getResisitivityData(filename)
      if length(strfind(filename,'aprj'))
        [record] = read_aprj(filename);
      else
        [record]=read_data(filename);
      end
    end


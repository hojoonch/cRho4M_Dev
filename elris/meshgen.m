function [p,t,nlay,tev,par,npar,zi]=meshgen(data,pathname,datapath)

% Mesh generator. 'normal' mode is selected. Model space is constructed by
% dividing rectangular blocks into two triangles. Outer part of the mesh is
% unstructured. Unstructured mesh is produced by Triangle program

useOriginal = false;


if useOriginal
	% adjust data.zmax
	data.zmax = data.ela * 5;
	data.dz1 = 2.5;
	dz(1)=data.dz1;
	%-----------------

	%dz(1)=data.dz1/1.5; % original
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
	zi=-[0 z];
else
	data.dz1 = 0.5*data.ela;
	dz = [0.5 0.5 1 1.5 2.5]*data.ela;
	nlay = 5;
	z=cumsum(dz);
	data.zmax=max(z);
	oran=max(z)/data.zmax;
	z=z/oran;
	zi=-[0 z];
end

ek=10;

x1=data.xelek(1)-ek*data.ela;
x2=data.xelek(1)*ones(size(zi));
x3=data.xelek(end)*ones(size(zi));
x4=data.xelek(end)+ek*data.ela;
zdoi=zi(end);

zcerceve=[0  zi ones(1,data.nel-2)*zi(end) fliplr(zi)  0 -6*data.zmax -6*data.zmax]';
xcerceve=[x1 x2 data.xelek(2:end-1)' x3 x4 x4 x1];
[x,y]=meshgrid(data.xelek(2:end-1),zi(1:end-1));

%preparing input file for Triangle

%pfix=data.filename(1:end-4);
pfix=data.filename(1:strfind(data.filename,'.')-1);
fidp=fopen([pfix,'.poly'],'w');
np=length(xcerceve);
fprintf(fidp,'%d 2 1 0\n',np);
for k=1:np
    fprintf(fidp,'%d %8.4f %8.4f 1\n',k,xcerceve(k),zcerceve(k));
end
fprintf(fidp,'%d 0\n',np);

for k=1:np-1
    fprintf(fidp,'%d %8.4f %8.4f 1\n',k,k,k+1);
end
fprintf(fidp,'%d %8.4f %8.4f 1\n1\n',np,np,1);
fprintf(fidp,'1 %8.4f %8.4f',-data.zmax/2,data.xelek(fix((2*data.nel-1)/2)));
fclose(fidp);

%eval(['!triangle -Q -q ',pfix,'.poly']) %matlab
%system(['! /opt/local/bin/triangle -Q -q ',pfix,'.poly']) %gnu octave
%system(['! /home/pi/cRho/triangle/triangle -Q -q ',pfix,'.poly']) %gnu octave
cname = computer();
if strfind(cname, "w64" )
  % windows
  system([pathname,'\triangle.exe -Q -q ',pfix,'.poly']);
elseif strfind(cname, "arm" ) || strfind(cname, "aarch64" )
  % Raspberry Pi
  [dum, os_bit] = system('getconf LONG_BIT');
  if str2num(os_bit) == 32
    disp('run triangle_rpi32')
    system([pathname,'/triangle_rpi32 -Q -q ',pfix,'.poly']);
  else
    disp('run triangle_rpi64')
    system([pathname,'/triangle_rpi64 -Q -q ',pfix,'.poly']);
  end

elseif strfind(cname, "pc-linux" )
  % Ubuntu
  system([pathname,'/triangle_ubuntu -Q -q ',pfix,'.poly']);
elseif strfind(cname, "apple" )
  % MacOs
  system([pathname,'/triangle_macos -Q -q ',pfix,'.poly']);
end

% p=load([pfix,'.1.node'])
[ node_num, marker ] = node_header_read ( [pfix,'.1.node']);
[ node_xy, node_marker ] = node_data_read ( [pfix,'.1.node'], node_num );
[ element_order, element_num ] = element_header_read ( [pfix,'.1.ele'] );
element_node = element_data_read ( [pfix,'.1.ele'], element_order, ...
    element_num );
t_dis=element_node;
p_dis=node_xy;

% check whether windows or not
icom = strfind(cname, 'w64');
if ~isempty(icom)
  % windows
  delete(strcat(datapath, '\*.poly'));
  delete(strcat(datapath,'\*.ele'));
  delete(strcat(datapath,'\*.node'));

else %cname(1:3) == 'arm'
  % Linux
  delete(strcat(datapath, '/*.poly'));
  delete(strcat(datapath,'/*.ele'));
  delete(strcat(datapath,'/*.node'));
end

nz=length(zi);
p_ic=[x(:)';y(:)'];

boy=2*nz+data.nel-2;
liste_dis=2:boy+1;
ic=reshape(length(p_dis)+1:length(p_dis)+numel(x),nz-1,data.nel-2);
alt=nz+2:liste_dis(end-nz);
T=[[2:nz+1]' [ic;alt] liste_dis(end:-1:end-nz+1)'] ;
tri=[];
for k=1:size(T,2)-1
    for m=1:size(T,1)-1
        t1=[T(m,k+1) T(m,k) T(m+1,k)];
        t2=[T(m,k+1) T(m+1,k) T(m+1,k+1)];
        tri=[tri;t1;t2];
    end
end
p=[p_dis p_ic];
t=[t_dis tri'];


%linking triangles with parameters
sira=[1:nlay]';
tev=sira;
for k=2:data.nel-1
    sira=sira+nlay;
    tev=[tev; sira];
end
npar=length(tev);
bas=length(t_dis)+1;
npar=length(tev);
for ip=1:npar
    par(ip).ucg=[bas bas+1];
    bas=bas+2;
end



end
function [ node_xy, node_marker ] = node_data_read ( node_filename, ...
    node_num )

%*****************************************************************************80
%
%% NODE_DATA_READ reads the data of a node file.
%  Licensing:
%    This code is distributed under the GNU LGPL license.
%  Modified:
%    01 November 2010
%  Author:
%    John Burkardt
%  Reference:
%    Jonathan Shewchuk,
%    Triangle: Engineering a 2D Quality Mesh Generator and
%    Delaunay Triangulator,
%    in Applied Computational Geometry: Towards Geometric Engineering,
%    edited by Ming Lin, Dinesh Manocha,
%    Lecture Notes in Computer Science, Volume 1148,
%    Springer, 1996,
%    ISBN: 354061785X,
%    LC: QA448.D38.A635.
%
%  Parameters:
%    Input, string NODE_FILENAME, the name of the file.
%    Input, integer NODE_NUM, the number of nodes.
%    Output, real NODE_COORD(2,NODE_NUM), ...
%    Output, integer NODE_MARKER(NODE_NUM), ...

node_xy = zeros ( 2, node_num );
node_marker = zeros ( node_num, 1 );

input = fopen ( node_filename, 'rt' );
text_num = 0;
read_header = 0;
node = 0;

while ( 1 )
    text = fgetl ( input );
    if ( text == -1 )
        break
    end
    text_num = text_num + 1;
    if ( text(1) == '#' )
        continue
    end

    if ( s_len_trim ( text ) <= 0 )
        continue
    end
    %
    %  Read (but ignore) header line.
    %
    if ( ~read_header )

        [ a, count ] = sscanf ( text, '%d', 4 );

        if ( count ~= 4 )
            fprintf ( 1, '\n' );
            fprintf ( 1, 'NODE_DATA_READ - Fatal error!\n' );
            fprintf ( 1, '  Could not read 4 integers from node header line.\n' );
            fprintf ( 1, '  File is "%s".\n', node_filename );
            fprintf ( 1, '  Line number = %d.\n', text_num );
            error ( 'NODE_DATA_READ - Fatal error!' );
        end

        read_header = 1;
        %
        %  Read data for next node.
        %
    else
        [ a, count ] = sscanf ( text, '%d  %f  %f %d', 4 );
        if ( count ~= 4 )
            fprintf ( 1, '\n' );
            fprintf ( 1, 'NODE_DATA_READ - Fatal error!\n' );
            fprintf ( 1, '  Could not read 4 values for node %d.\n',  node + 1 );
            fprintf ( 1, '  File is "%s".\n', node_filename );
            fprintf ( 1, '  Line number = %d.\n', text_num );
            error ( 'NODE_DATA_READ - Fatal error!' );
        end

        node = node + 1;
        node_xy(1,node) = a(2);
        node_xy(2,node) = a(3);
        node_marker(node) = a(4);

        if ( node_num <= node )
            break
        end

    end

end

fclose ( input );

return
end
function [ node_num, marker ] = node_header_read ( node_filename )

%*****************************************************************************80
%
%% NODE_HEADER_READ reads the header of a node file.
%
%  Discussion:
%
%    The header is a single line, of the form:
%
%      node_num  dim  att_num  marker
%
%    NODE_NUM is the number of nodes;
%    DIM is the spatial dimension, which should be 2;
%    ATT_NUM is the number of attributes;
%    MARKER is 1 to indicate that the last column of the node data
%    is used to identify boundary vertices and constrained vertices.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    01 November 2010
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Jonathan Shewchuk,
%    Triangle: Engineering a 2D Quality Mesh Generator and
%    Delaunay Triangulator,
%    in Applied Computational Geometry: Towards Geometric Engineering,
%    edited by Ming Lin, Dinesh Manocha,
%    Lecture Notes in Computer Science, Volume 1148,
%    Springer, 1996,
%    ISBN: 354061785X,
%    LC: QA448.D38.A635.
%
%  Parameters:
%
%    Input, string NODE_FILENAME, the name of the file.
%
%    Output, integer NODE_NUM, the number of nodes.
%
%    Output, integer MARKER, is 1 if boundary nodes are marked by the
%    value 1.
%
input = fopen ( node_filename, 'rt' );
text_num = 0;
while ( 1 )
    text = fgetl ( input );

    if ( text == -1 )
        break
    end
    text_num = text_num + 1;

    if ( text(1) == '#' )
        continue
    end
    if ( s_len_trim ( text ) <= 0 )
        continue
    end
    [ a, count ] = sscanf ( text, '%d', 4 );
    if ( count ~= 4 )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'NODE_HEADER_READ - Fatal error!\n' );
        fprintf ( 1, '  Could not read 4 integers from node header line.\n' );
        fprintf ( 1, '  File is "%s".\n', node_filename );
        fprintf ( 1, '  Line number = %d.\n', text_num );
        error ( 'NODE_HEADER_READ - Fatal error!' );
    end
    node_num = a(1);
    dim = a(2);
    att_num = a(3);
    marker = a(4);
    break
end
fclose ( input );
return
end
function element_node = element_data_read ( element_filename, element_order, ...
    element_num )

%*****************************************************************************80
%% ELEMENT_DATA_READ reads the data of an element file.
%  Discussion:
%  Licensing:
%    This code is distributed under the GNU LGPL license.
%  Modified:
%    01 November 2010
%  Author:
%    John Burkardt
%  Reference:
%    Jonathan Shewchuk,
%    Triangle: Engineering a 2D Quality Mesh Generator and
%    Delaunay Triangulator,
%    in Applied Computational Geometry: Towards Geometric Engineering,
%    edited by Ming Lin, Dinesh Manocha,
%    Lecture Notes in Computer Science, Volume 1148,
%    Springer, 1996,
%    ISBN: 354061785X,
%    LC: QA448.D38.A635.
%
%  Parameters:
%
%    Input, string ELEMENT_FILENAME, the name of the file.
%    Input, integer ELEMENT_ORDER, the order of the elements.
%    Input, integer ELEMENT_NUM, the number of elements.
%    Output, real ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), ...
%
element_node = zeros ( element_order, element_num );

input = fopen ( element_filename, 'rt' );
text_num = 0;
read_header = 0;
element = 0;

while ( 1 )
    text = fgetl ( input );
    if ( text == -1 )
        break
    end
    text_num = text_num + 1;
    if ( text(1) == '#' )
        continue
    end
    if ( s_len_trim ( text ) <= 0 )
        continue
    end
    %
    %  Read (but ignore) header line.
    if ( ~read_header )
        [ a, count ] = sscanf ( text, '%d', 3 );
        if ( count ~= 3 )
            fprintf ( 1, '\n' );
            fprintf ( 1, 'ELEMENT_DATA_READ - Fatal error!\n' );
            fprintf ( 1, '  Could not read 3 integers from element header line.\n' );
            fprintf ( 1, '  File is "%s".\n', element_filename );
            fprintf ( 1, '  Line number = %d.\n', text_num );
            error ( 'ELEMENT_DATA_READ - Fatal error!' );
        end
        read_header = 1;
        %  Read data for next element.
    else
        [ a, count ] = sscanf ( text, '%d', element_order + 1 );
        if ( count ~= 4 )
            fprintf ( 1, '\n' );
            fprintf ( 1, 'ELEMENT_DATA_READ - Fatal error!\n' );
            fprintf ( 1, '  Could not read %d values for element %d.\n',  ...
                element_order + 1, element + 1 );
            fprintf ( 1, '  File is "%s".\n', element_filename );
            fprintf ( 1, '  Line number = %d.\n', text_num );
            error ( 'ELEMENT_DATA_READ - Fatal error!' );
        end
        element = element + 1;
        element_node(1:element_order,element) = a(2:element_order+1);
        if ( element_num <= element )
            break
        end
    end
end
fclose ( input );
return
end
function len = s_len_trim ( s )

%*****************************************************************************80
%
%% S_LEN_TRIM returns the length of a character string to the last nonblank.
%
%  Licensing:
%    This code is distributed under the GNU LGPL license.
%  Modified:
%    14 June 2003
%  Author:
%  John Burkardt
%  Parameters:
%    Input, string S, the string to be measured.
%    Output, integer LEN, the length of the string up to the last nonblank.
len = length ( s );
while ( 0 < len )
    if ( s(len) ~= ' ' )
        return
    end
    len = len - 1;
end
return
end
function [ element_order, element_num ] = element_header_read ( ...
    element_filename )

%*****************************************************************************80
%
%% ELEMENT_HEADER_READ reads the header of an element file.%
%  Discussion:
%    The header is a single line, of the form:
%    element_num  element_order  att_num
%    ELEMENT_NUM is the number of elements;
%    ELEMENT_ORDER is the element order (number of nodes per element);
%    ATT_NUM is the number of attributes;
%
%  Licensing:
%  This code is distributed under the GNU LGPL license.
%  Modified:
%  01 November 2010
%
%  Author:
%  John Burkardt
%  Reference:
%
%    Jonathan Shewchuk,
%    Triangle: Engineering a 2D Quality Mesh Generator and
%    Delaunay Triangulator,
%    in Applied Computational Geometry: Towards Geometric Engineering,
%    edited by Ming Lin, Dinesh Manocha,
%    Lecture Notes in Computer Science, Volume 1148,
%    Springer, 1996,
%    ISBN: 354061785X,
%    LC: QA448.D38.A635.
%
%  Parameters:
%    Input, string ELEMENT_FILENAME, the name of the file.
%    Output, integer ELEMENT_ORDER, the element order.
%    Output, integer ELEMENT_NUM, the number of elements.

input = fopen ( element_filename, 'rt' );
text_num = 0;

while ( 1 )
    text = fgetl ( input );
    if ( text == -1 )
        break
    end
    text_num = text_num + 1;
    if ( text(1) == '#' )
        continue
    end
    if ( s_len_trim ( text ) <= 0 )
        continue
    end
    [ a, count ] = sscanf ( text, '%d', 3 );
    if ( count ~= 3 )
        fprintf ( 1, '\n' );
        fprintf ( 1, 'ELEMENT_HEADER_READ - Fatal error!\n' );
        fprintf ( 1, '  Could not read 3 integers from element header line.\n' );
        fprintf ( 1, '  File is "%s".\n', element_filename );
        fprintf ( 1, '  Line number = %d.\n', text_num );
        error ( 'ELEMENT_HEADER_READ - Fatal error!' );
    end

    element_num = a(1);
    element_order = a(2);
    att_num = a(3);
    break
end
fclose (input );
return
end

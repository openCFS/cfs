function [ file, volume, dimension ] = generateShearedCrossExact(point,filepath,nx,ny)
% GENERATECROSS  -  Generates a sheared cross.
%
% @param:
%       point(1)  thickness of horizontal bar in [0,1]
%       point(2)  thickness of vertical bar in [0,1]
%       point(3)  shearing angle in [0,1] (will be mapped to [-pi/2,pi/2])
%       filepath  path to generated mesh file (optional)
%       nx        resolution of mesh (optional)
% 
% 
%                                    o  angle o o
%                                    o     xxxxoxxxx                
%                                    o    xxxxoxxxx                
%                                    o   xxxxoxxxx                
%                                    o  xxxxoxxxx                
%                                    o xxxxoxxxx                
%                                    oxxxxoxxxx                
%                                    oxxxoxxxx                
%                                   xoxxoxxxx                
%                  xxxxxxxxxxxxxxxxxxoxoxxxxxxxxxxxxxxxxxxxx
%                 xxxxxxxxxxxxxxxxxxxooxxxxxxxxxxxxxxxxxxxx
%            s1  xxxxxxxxxxxxxxxxxxxxoxxxxxxxxxxxxxxxxxxxx
%               xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%              xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%                             xxxxxxxxx                
%                            xxxxxxxxx                
%                           xxxxxxxxx                
%                          xxxxxxxxx                
%                         xxxxxxxxx                
%                        xxxxxxxxx                
%                       xxxxxxxxx                
%                      xxxxxxxxx                
%                         s2
% 

if nargin < 3
    nx = 128;
    ny = 128;
end
if nargin < 2
    filepath = '.';
end

if length(point) ~= 3
    throw( MException( 'generateFramedCross:wrongParameterCount', 'point does not contain three parameters.' ));
end

if point(1) < 0 || point(1) > 1
    throw( MException( 'generateFramedCross:parameterOutOfRange', sprintf('point(1) has to be in [0;1], but is %f.', point(1)) ));
end
if point(2) < 0 || point(2) > 1
    throw( MException( 'generateFramedCross:parameterOutOfRange', sprintf('point(2) has to be in [0;1], but is %f.', point(3)) ));
end
if point(3) <= 0 || point(3) >= 1
    throw( MException( 'generateFramedCross:parameterOutOfRange', sprintf('point(3) has to be in (0;1), but is %f.', point(3)) ));
end

s1 = point(1);
s2 = point(2);
s3 = (.5-point(3))*pi;

volume = s1 + s2 - s1*s2;

dimension = 2;

% Initial number of elements
nx = min(max(nx,ceil(128/500/min(s1,s2*3))),4096);
ny = min(max(ny,ceil(128/500/min(s1*3,s2))),4096);

% Nodes
% Nodes are separated into three parts:
% 0 to remainder / remainder to remainder + s / remainder + s to 1
% Thus the second part has excactly the width s.

yremainder = (1-s1)/2;
xremainder = (1-s2)/2;

% Points in each part
px1 = ceil(xremainder*nx)+1;
px2 = max(ceil(s2*nx)+1,3);
px3 = px1;
py1 = ceil(yremainder*ny)+1;
py2 = max(ceil(s1*ny)+1,3);
py3 = py1;

% Coordinates of points in each part
x1 = linspace(0,xremainder,px1);
x2 = linspace(x1(end),x1(end)+s2,px2);
x3 = linspace(x2(end),1,px3);
x = [x1,x2(2:end),x3(2:end)];

y1 = linspace(0,yremainder,py1);
y2 = linspace(y1(end),y1(end)+s1,py2);
y3 = linspace(y2(end),1,py3);
y = [y1,y2(2:end),y3(2:end)];

% New number of elements
npx = size(x,2);
npy = size(y,2);
nx = npx - 1;
ny = npy - 1;
numNodes = npx*npy;


% Density matrix
A = zeros(ny,nx);
vertBar = numel(x1):numel(x1)+numel(x2)-2;
horzBar = numel(y1):numel(y1)+numel(y2)-2;
A(horzBar,:) = 1;
A(:,vertBar) = 1;


% Node coordinates
ycoords = reshape(repmat(y,npx,1),1,numNodes);
if s3 == 0 % no shearing -> just repeat
    xcoords = repmat(x,1,npy);
else % shearing -> add offset
    xcoords = zeros(1,numNodes);
    for i=1:npy
        idx = (1:npx)+(i-1)*npx;
        xcoords(idx) = x + y(i) * tan(s3);
    end
end
nodes = [1:numNodes; xcoords; ycoords; zeros(1,numNodes)];

% Boundary nodes
bottom = 1:npx;
top = numNodes-npx+1:numNodes;
left = 1:npx:numNodes;
right = npx:npx:numNodes;
numNodeBC = numel(bottom) + numel(top) + numel(left) + numel(right);

% Elements
num2DElements = nx*ny;
elems = zeros(4,num2DElements);

for i=1:ny
    z = 1:npx;
    a = z+(i-1)*npx;
    b = z+i*npx;

    % Repeat inner node numbers
    x = [a(1),reshape([a(2:end-1);a(2:end-1)],1,2*(nx-1)),a(end)];
    y = [b(1),reshape([b(2:end-1);b(2:end-1)],1,2*(nx-1)),b(end)];

    % Get element indeces
    idx = (1:nx)+(i-1)*nx;
    
    % Get nodes of elements
    elems(:,idx) = [reshape(x,2,[]);flipud(reshape(y,2,[]))];
end
elems = [1:num2DElements;elems];


% Write mesh and density file
[~,filename] = fileparts(tempname);

% Get absolute path of mesh file
oldpath=pwd;
cd(filepath)
fullpath=pwd;
cd(oldpath)
filename = fullfile(fullpath,filename);

file = [filename,'.dens'];
meshfile = [filename,'.mesh'];

% Write density file
Homogenization.matrixToDensity(A,file);

% Write mesh file
fid = fopen(meshfile,'wt');

fprintf(fid,'[Info]\n');
fprintf(fid,'Version 1\n');
fprintf(fid,'Dimension 2\n');
fprintf(fid,'NumNodes %d\n',numNodes);
fprintf(fid,'Num3DElements 0\n');
fprintf(fid,'Num2DElements %d\n',num2DElements);
fprintf(fid,'Num1DElements 0\n');
fprintf(fid,'NumNodeBC %d\n',numNodeBC);
fprintf(fid,'NumSaveNodes 0\n');
fprintf(fid,'NumSaveElements 0\n');
fprintf(fid,'Num 2d-line      : 0\n');
fprintf(fid,'Num 2d-line,quad : 0\n');
fprintf(fid,'Num 3d-line      : 0\n');
fprintf(fid,'Num 3d-line,quad : 0\n');
fprintf(fid,'Num triangle     : 0\n');
fprintf(fid,'Num triangle,quad: 0\n');
fprintf(fid,'Num quadr        : %d\n',num2DElements);
fprintf(fid,'Num quadr,quad   : 0\n');
fprintf(fid,'Num tetra        : 0\n');
fprintf(fid,'Num tetra,quad   : 0\n');
fprintf(fid,'Num brick        : 0\n');
fprintf(fid,'Num brick,quad   : 0\n');
fprintf(fid,'Num pyramid      : 0\n');
fprintf(fid,'Num pyramid,quad : 0\n');
fprintf(fid,'Num wedge        : 0\n');
fprintf(fid,'Num wedge,quad   : 0\n');
fprintf(fid,'\n');

fprintf(fid,'[Nodes]\n');
fprintf(fid,'#NodeNr x-coord y-coord z-coord\n');
if ~isempty(nodes)
    fprintf(fid,'% 8d %f %f %f\n',nodes);
end
fprintf(fid,'\n');

fprintf(fid,'[1D Elements]\n');
fprintf(fid,'#ElemNr  ElemType  NrOfNodes  Level\n');
fprintf(fid,'#Node1 Node2 ... NodeNrOfNodes\n');
fprintf(fid,'\n');

fprintf(fid,'[2D Elements]\n');
fprintf(fid,'#ElemNr  ElemType  NrOfNodes  Level\n');
fprintf(fid,'#Node1 Node2 ... NodeNrOfNodes\n');
if ~isempty(elems)
    fprintf(fid,'%d 6 4 mech\n%d %d %d %d\n',elems);
end
fprintf(fid,'\n');

fprintf(fid,'[3D Elements]\n');
fprintf(fid,'#ElemNr  ElemType  NrOfNodes  Level\n');
fprintf(fid,'#Node1 Node2 ... NodeNrOfNodes\n');
fprintf(fid,'\n');

fprintf(fid,'[Node BC]\n');
fprintf(fid,'#NodeNr Level\n');
if ~isempty(bottom)
    fprintf(fid,'% 8d bottom\n',int32(bottom));
    fprintf(fid,'% 8d top\n',int32(top));
end
if ~isempty(left)
    fprintf(fid,'% 8d left\n',int32(left));
    fprintf(fid,'% 8d right\n',int32(right));
end
fprintf(fid,'\n');

fprintf(fid,'[Save Nodes]\n');
fprintf(fid,'#NodeNr Level\n');
fprintf(fid,'\n');

fprintf(fid,'[Save Elements]\n');
fprintf(fid,'#ElemNr Level\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

fclose(fid);

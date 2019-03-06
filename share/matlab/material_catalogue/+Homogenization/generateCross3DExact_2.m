function [ file, volume, dimension ] = generateCross3DExact(point,filepath,nx,ny,nz)
% GENERATECROSS  -  Generates a sheared cross.
%
% @param:
%       point(1)  thickness of bar in x direction in [0,1]
%       point(2)  thickness of bar in y direction in [0,1]
%       point(3)  thickness of bar in z direction in [0,1]
%       filepath  path to generated mesh file (optional)
%       nx        resolution of mesh in x direction (optional)
%       ny        resolution of mesh in y direction (optional)
%       nz        resolution of mesh in z direction (optional)
% 
% 
%                      xxxxxxxxx                
%                      xxxxxxxxx                
%                      xxxxxxxxx                
%                      xxxxxxxxx                
%                      xxxxxxxxx                
%                      xxxxxxxxx                
%                      xxxxxxxxx                
%                      xxxxxxxxx                
%      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%                      xxxxxxxxx                
%                      xxxxxxxxx                
%                      xxxxxxxxx                
%                      xxxxxxxxx                
%                      xxxxxxxxx                
%                      xxxxxxxxx                
%                      xxxxxxxxx                
%                      xxxxxxxxx                
% 

if nargin < 3
    nx = 16;
    ny = 16;
    nz = 16;
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
if point(3) < 0 || point(3) > 1
    throw( MException( 'generateFramedCross:parameterOutOfRange', sprintf('point(3) has to be in [0;1], but is %f.', point(3)) ));
end

s1 = point(1);
s2 = point(2);
s3 = point(3);

volume = CalcCross3DVolume(s1, s2, s3);

dimension = 3;

% Initial number of elements
nx = min(max(nx,ceil(128/500/min(s1,s2*3,s3*3))),4096);
ny = min(max(ny,ceil(128/500/min(s1*3,s2,s3*3))),4096);
nz = min(max(nz,ceil(128/500/min(s1*3,s2*3,s3))),4096);

% Nodes
% Nodes are eparated into three parts:
% 0 to remainder / remainder to remainder + s / remainder + s to 1
% Thus the second part has excactly the width s.

xremainder = (1-s1)/2;
yremainder = (1-s2)/2;
zremainder = (1-s3)/2;

% Points in each part
px1 = max(ceil(xremainder*nx)+1,3);
px2 = max(ceil(s1*nx)+1,3);
px3 = max(ceil(xremainder*nx)+1,3);
py1 = max(ceil(yremainder*ny)+1,3);
py2 = max(ceil(s2*ny)+1,3);
py3 = max(ceil(yremainder*ny)+1,3);
pz1 = max(ceil(zremainder*nz)+1,3);
pz2 = max(ceil(s3*nz)+1,3);
pz3 = max(ceil(zremainder*nz)+1,3);

% Coordinates of points in each part
x1 = linspace(0,xremainder,px1);
x2 = linspace(x1(end),x1(end)+s1,px2);
x3 = linspace(x2(end),1,px3);
x = [x1,x2(2:end),x3(2:end)];
y1 = linspace(0,yremainder,py1);
y2 = linspace(y1(end),y1(end)+s2,py2);
y3 = linspace(y2(end),1,py3);
y = [y1,y2(2:end),y3(2:end)];
z1 = linspace(0,zremainder,pz1);
z2 = linspace(z1(end),z1(end)+s3,pz2);
z3 = linspace(z2(end),1,pz3);
z = [z1,z2(2:end),z3(2:end)];

% New number of elements
npx = size(x,2);
npy = size(y,2);
npz = size(z,2);
nx = npx - 1;
ny = npy - 1;
nz = npz - 1;
numNodes = npx*npy*npz;


% Density matrix
A = zeros(nx,ny,nz);
xBar = numel(x1):numel(x1)+numel(x2)-2;
yBar = numel(y1):numel(y1)+numel(y2)-2;
zBar = numel(z1):numel(z1)+numel(z2)-2;
A(xBar,:,:) = 1;
A(:,yBar,:) = 1;
A(:,:,zBar) = 1;


% Node coordinates
xcoords = reshape(repmat(x,1,npy*npz),1,numNodes);
ycoords = repmat(reshape(repmat(y,npx,1),1,npx*npy),1,npz);
zcoords = reshape(repmat(z,npx*npy,1),1,numNodes);
nodes = [1:numNodes; xcoords; ycoords; zcoords];

% Boundary nodes
left = 1:npx:numNodes;
right = npx:npx:numNodes;
bottom = zeros(1,npx*npz);
for i = 1:npz
    for j = 1:npx
        bottom( (i-1)*npx + j ) = (i-1) * npy * npx + j;
    end
end
top = zeros(1,npx*npz);
for i = 1:npz
    for j = 1:npx
        top( (i-1)*npx + j ) = (i * npy - 1) * npx + j;
    end
end
back = 1:npx*npy;
front = (npz-1)*npx*npy+1:numNodes;

numNodeBC = numel(bottom) + numel(top) + numel(left) + numel(right) + numel(back) + numel(front);

% Elements
num3DElements = nx*ny*nz;
elems = zeros(8,num3DElements);

idx = 0;
for z=0:nz-1
    for y=0:ny-1
        for x=0:nx-1
            % lower-left-back of current element
            ll = npx*npy*z + npx*y + x;

            % Get nodes of element
            elemnodes = [ll+npx, ll+1+npx, ll+1+(npx*npy)+npx, ll+(npx*npy)+npx, ...
                ll, ll+1, ll+1+(npx*npy), ll+(npx*npy)] + 1;
                
            % Get element indeces
            idx = idx + 1;

            elems(:,idx) = elemnodes';
        end
    end
end
elems = [1:num3DElements;elems];


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
fprintf(fid,'Dimension 3\n');
fprintf(fid,'NumNodes %d\n',numNodes);
fprintf(fid,'Num3DElements %d\n',num3DElements);
fprintf(fid,'Num2DElements 0\n');
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
fprintf(fid,'Num quadr        : 0\n');
fprintf(fid,'Num quadr,quad   : 0\n');
fprintf(fid,'Num tetra        : 0\n');
fprintf(fid,'Num tetra,quad   : 0\n');
fprintf(fid,'Num brick        : %d\n',num3DElements);
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
fprintf(fid,'\n');

fprintf(fid,'[3D Elements]\n');
fprintf(fid,'#ElemNr  ElemType  NrOfNodes  Level\n');
fprintf(fid,'#Node1 Node2 ... NodeNrOfNodes\n');
if ~isempty(elems)
    fprintf(fid,'%d 10 8 mech\n%d %d %d %d %d %d %d %d\n',elems);
end
fprintf(fid,'\n');

fprintf(fid,'[Node BC]\n');
fprintf(fid,'#NodeNr Level\n');
if ~isempty(bottom)
    fprintf(fid,'% 8d nodes3\n',int32(bottom));
    fprintf(fid,'% 8d nodes4\n',int32(top));
end
if ~isempty(left)
    fprintf(fid,'% 8d nodes5\n',int32(left));
    fprintf(fid,'% 8d nodes6\n',int32(right));
end
if ~isempty(back)
    fprintf(fid,'% 8d nodes7\n',int32(back));
    fprintf(fid,'% 8d nodes8\n',int32(front));
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

end % function


function [vol] = CalcCross3DVolume(stiff1, stiff2, stiff3)
  
if (stiff1 >= stiff2 && stiff1 >= stiff3)
    vol = stiff1*stiff1 + stiff2*stiff2 + stiff3*stiff3 - stiff1*stiff3*stiff3 - stiff1*stiff2*stiff2;
elseif (stiff2 >= stiff1 && stiff2 >= stiff3)
    vol = stiff1*stiff1 + stiff2*stiff2 + stiff3*stiff3 - stiff2*stiff3*stiff3 - stiff2*stiff1*stiff1;
else % if (stiff3 >= stiff1 && stiff3 >= stiff2) 
    vol = stiff1*stiff1 + stiff2*stiff2 + stiff3*stiff3 - stiff3*stiff2*stiff2 - stiff3*stiff1*stiff1;
end

end
function [ meshfile, volume ] = generateShearedCross(point,filepath)
%GENERATE_MESH generates an ANSYS mesh for a sheared cross.
%
%   meshfile = generate_mesh(filepath,nx,ny,a,b,phi) generates a .mesh
%   file in ANSYS notation. The file is stored at filepath. The
%   resolution of the micromesh is nx x ny and the sheared cross is given
%   by the width of the vertical beam a in [0,1], the width of the
%   horizontal beam b in [0,1] and the shearing angle phi in (0,1). No
%   shearing applies, when phi = 0.5.
%
%   The width a of the vertical beam is the normal onto its long side, so
%   the actual width ah of the sheared beam may differ.
%   
%         <--ah->
%         xxxxxxx              ^
%          xxxxxxx             |
%           xxxxxxx            |
%    ^ xxxxxxxxxxxxxxxxxxx     |
%    |  xxxxxxxxxxxxxxxxxxx    |
%    b   xxxxxxxxxoxxxxxxxxx   ny
%    |    xxxxxxxxxxxxxxxxxxx  |
%    v     xxxxxxxxxxxxxxxxxxx |
%                 xxxxxxx      |
%                  xxxxxxx     |
%                   xxxxxxx    v
%      <----------nx--------->
%
%   Example:
%   Homogenization.generate_mesh('.',32,32,0.6,0.4,0.7)
%

nx = 128;
ny = nx;

a = point(1);
b = point(2);
phi = point(3);

dx = 1/nx;
dy = 1/ny;

% Shift phi and convert it to radians
% phi = phi/180*pi;
phi = (phi - .5)*pi;



% The given a is the actual width of the vertical beam. For the mesh we need
% the horizontal width of this beam.
% ah = a/cos(phi);
ah = a;

% if ah > 1
%     disp('generate_mesh: phi too big.');
%     meshfile = -1;
%     return;
% end

% To ensure periodic boundary conditions meshes are not allowed to have
% only one void element per row or column (which would be located at the 
% boundary). So if ah or b exceeds a given threshold fill complete row or
% column with material.
if round(nx*ah) >= nx*(1-dx)
    xCoordMin = 0;
    xCoordMax = 1;
else
    xCoordMin = round(nx*(1-ah)/2)*dx;
    xCoordMax = round(nx*(xCoordMin+ah))*dx;
end
if round(ny*b) >= ny*(1-dy)
    yCoordMin = 0;
    yCoordMax = 1;
else
    yCoordMin = round(ny*(1-b)/2)*dy;
    yCoordMax = round(ny*(yCoordMin+b))*dy;
end

xCoord = xCoordMin:dx:xCoordMax;
yCoord = yCoordMin:dy:yCoordMax;
numXCoord = numel(xCoord);
numYCoord = numel(yCoord);

% If a beam exists only of one column/row of nodes and no elements, this
% beam does not exist.
if numXCoord <= 1
    numXCoord = 0;
    xCoordMax = xCoordMax - dx;
end
if numYCoord <= 1
    numYCoord = 0;
end

% number of nodes and elements
numNodesVertBeam = (nx+1)*numXCoord;
numNodesHorzBeam = (ny+1-numXCoord)*numYCoord;
numNodes = numNodesVertBeam + numNodesHorzBeam;
if numYCoord ~= 0
    num2DElementsHorzBeamPerRow1 = round(nx*xCoordMin);
    num2DElementsHorzBeamPerRow2 = nx-num2DElementsHorzBeamPerRow1-(numXCoord-1);
    if numNodesVertBeam == 0
        num2DElementsHorzBeamPerRow2 = num2DElementsHorzBeamPerRow2 - 1;
    end
else
    num2DElementsHorzBeamPerRow1 = 0;
    num2DElementsHorzBeamPerRow2 = 0;
end
% number of 2D elements = elements in vert beam + elements in horz beam
num2DElements = ny*max(numXCoord-1,0) + (nx-max(numXCoord-1,0))*max(numYCoord-1,0);

% Get boundary nodes
bottom = 1:numXCoord;
top = (1:numXCoord) +ny*numXCoord;
% If horizontal beam has width 1 (i.e. full material), add nodes to bottom 
% and top which are members of horizontal beam but not (yet) of vertical beam
if num2DElementsHorzBeamPerRow1 ~= 0 && numYCoord == ny+1
    tmp1 = 1:num2DElementsHorzBeamPerRow1;
    tmp2 = 1:num2DElementsHorzBeamPerRow2;
    offset = (ny+1)*(numXCoord+num2DElementsHorzBeamPerRow1);
    bottom = [bottom,tmp1+(ny+1)*numXCoord,tmp2+offset];
    top = [top,tmp1+offset-num2DElementsHorzBeamPerRow1,tmp2+numNodes-num2DElementsHorzBeamPerRow2];
end
if num2DElementsHorzBeamPerRow1==0
    if numXCoord == nx+1
        left = 1:nx+1:(nx+1)*(ny+1);
        right = left + nx;
    else
        left = [];
        right = [];
    end
else
    left = (1:num2DElementsHorzBeamPerRow1:num2DElementsHorzBeamPerRow1*numYCoord) +numNodesVertBeam;
    right = left(end)+num2DElementsHorzBeamPerRow1+num2DElementsHorzBeamPerRow2-1:num2DElementsHorzBeamPerRow2:numNodes;
end
numNodeBC = (numel(bottom)+numel(left)) * 2;


% Get all nodes
nodeNum = 1;
nodes = zeros(4,numNodes);
% Nodes in vertical beam
if numNodesVertBeam > 0
    nodeNum = 1:numXCoord;
    for i=1:ny+1
        nodes(:,nodeNum) = [nodeNum;xCoord + (i-1)*dy*tan(phi);(i-1)*dy*ones(1,numXCoord);zeros(1,numXCoord)];
        nodeNum = nodeNum + numXCoord;
    end
end

% Nodes in horizontal beam
if numNodesHorzBeam > 0
    indeces = 1:num2DElementsHorzBeamPerRow1;
    nodeNum = min(nodeNum)-1 + indeces;
    for i=1:numYCoord
        nodes(:,nodeNum) = [nodeNum;(indeces-1)*dx + yCoord(i)*tan(phi);yCoord(i)*ones(1,numel(indeces));zeros(1,numel(indeces))];
        nodeNum = nodeNum + numel(indeces);
    end
    indeces = round(nx*xCoordMax+1):nx;
    nodeNum = min(nodeNum)-1 + (1:numel(indeces));
    for i=1:numYCoord
        nodes(:,nodeNum) = [nodeNum;indeces*dx + yCoord(i)*tan(phi);yCoord(i)*ones(1,numel(indeces));zeros(1,numel(indeces))];
        nodeNum = nodeNum + numel(indeces);
    end
end


% Get all elements
count = 1;
elems = zeros(4,num2DElements);

% Elements in vertical beam
if numNodesVertBeam > 0
    if numXCoord > 0
        tmp = 1:numXCoord-1;
        for i=1:ny
            A = [tmp;tmp+1;tmp+1+numXCoord;tmp+numXCoord] + (i-1)*numXCoord;
            elems(:,(1:numXCoord-1) + (i-1)*(numXCoord-1)) = A;
            count = count+numXCoord-1;
        end
    end
end

% Elements in horizontal beam
if numNodesHorzBeam > 0
    for i=1:numYCoord-1
        tmp = (1:num2DElementsHorzBeamPerRow1) + (i-1)*num2DElementsHorzBeamPerRow1 + numNodesVertBeam;
        if numNodesVertBeam > 0
            A = [tmp, int32(ny*yCoordMin*numXCoord+1+(i-1)*numXCoord);
                tmp + num2DElementsHorzBeamPerRow1, int32(ny*yCoordMin*numXCoord+1+i*numXCoord)];
        else
            A = [tmp, int32(num2DElementsHorzBeamPerRow1*numYCoord+1+(i-1)*(num2DElementsHorzBeamPerRow2+1));
                tmp + num2DElementsHorzBeamPerRow1, int32(num2DElementsHorzBeamPerRow1*numYCoord+1+i*(num2DElementsHorzBeamPerRow2+1))];
        end
        for j=1:num2DElementsHorzBeamPerRow1
            elems(:,count) = [A(1,j);A(1,j+1);A(2,j+1);A(2,j)];
            count = count + 1;
        end
        if numNodesVertBeam > 0
            tmp = (1:num2DElementsHorzBeamPerRow2) + (i-1)*num2DElementsHorzBeamPerRow2 + numNodesVertBeam + num2DElementsHorzBeamPerRow1*numYCoord;
            A = [int32(ny*yCoordMin*numXCoord+1+i*numXCoord-1), tmp;
                int32(ny*yCoordMin*numXCoord+1+(i+1)*numXCoord-1), tmp + num2DElementsHorzBeamPerRow2];
        else
            tmp = (1:num2DElementsHorzBeamPerRow2+1) + (i-1)*(num2DElementsHorzBeamPerRow2+1) + numNodesVertBeam + num2DElementsHorzBeamPerRow1*numYCoord;
            A = [tmp;tmp + num2DElementsHorzBeamPerRow2+1];
        end
        for j=1:num2DElementsHorzBeamPerRow2
            elems(:,count) = [A(1,j);A(1,j+1);A(2,j+1);A(2,j)];
            count = count + 1;
        end
    end
end

%%%%% Correct ? %%%%%
% If their are no nodes at the boundary we have to add one, so CFS can apply periodic BCs
if isempty(left)
    numNodes = numNodes + 2;
    numNodeBC = numNodeBC + 2;
end
if isempty(bottom)
    numNodes = numNodes + 2;
    numNodeBC = numNodeBC +2;
end
%%%%% Correct ? %%%%%


[~,temp] = fileparts(tempname);
meshName = [temp,'.mesh'];

% Get absolute path of mesh file
oldpath=pwd;
cd(filepath)
fullpath=pwd;
cd(oldpath)
meshfile = fullfile(fullpath,meshName);

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
fprintf(fid,'% 8d %f %f %f\n',nodes);

% If their are no nodes at the boundary we have to add one, so CFS can apply periodic BCs
nodeNum = min(nodeNum);
if isempty(left)
    fprintf(fid,'% 8d %f %f %f\n',nodeNum,.5*tan(phi),.5,0);
    left = nodeNum;
    nodeNum = nodeNum + 1;
    fprintf(fid,'% 8d %f %f %f\n',nodeNum,.5*tan(phi)+1,.5,0);
    right = nodeNum;
    nodeNum = nodeNum + 1;
end
if isempty(bottom)
    fprintf(fid,'% 8d %f %f %f\n',nodeNum,.5,0,0);
    bottom = nodeNum;
    nodeNum = nodeNum + 1;
    fprintf(fid,'% 8d %f %f %f\n',nodeNum,.5+tan(phi),1,0);
    top = nodeNum;
    nodeNum = nodeNum + 1;
end
fprintf(fid,'\n');


fprintf(fid,'[1D Elements]\n');
fprintf(fid,'#ElemNr  ElemType  NrOfNodes  Level\n');
fprintf(fid,'#Node1 Node2 ... NodeNrOfNodes\n');
fprintf(fid,'\n');

fprintf(fid,'[2D Elements]\n');
fprintf(fid,'#ElemNr  ElemType  NrOfNodes  Level\n');
fprintf(fid,'#Node1 Node2 ... NodeNrOfNodes\n');
fprintf(fid,'%d 6 4 mech\n%d %d %d %d\n',[1:num2DElements;elems]);

fprintf(fid,'\n');

fprintf(fid,'[3D Elements]\n');
fprintf(fid,'#ElemNr  ElemType  NrOfNodes  Level\n');
fprintf(fid,'#Node1 Node2 ... NodeNrOfNodes\n');
fprintf(fid,'\n');

fprintf(fid,'[Node BC]\n');
fprintf(fid,'#NodeNr Level\n');
fprintf(fid,'% 8d nodes3\n',int32(bottom));
fprintf(fid,'% 8d nodes4\n',int32(top));
fprintf(fid,'% 8d nodes5\n',int32(left));
fprintf(fid,'% 8d nodes6\n',int32(right));
fprintf(fid,'\n');

fprintf(fid,'[Save Nodes]\n');
fprintf(fid,'#NodeNr Level\n');
fprintf(fid,'\n');

fprintf(fid,'[Save Elements]\n');
fprintf(fid,'#ElemNr Level\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

fclose(fid);

volume = num2DElements/nx/ny;

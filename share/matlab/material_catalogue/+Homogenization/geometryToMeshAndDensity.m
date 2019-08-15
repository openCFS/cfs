function [ densityfile ] = geometryToMeshAndDensity(geom,meshfile,holemarker)
% GEOMETRYTOMESH  -  Generates a (sparse) ANSYS mesh file out of a geometry.
%
% @param:
%       geom            Decomposed Constructive Solid Geometry
%       meshfile        name of generated .mesh file
%                       (optional, default = 'mesh.mesh')
%

if nargin < 2
    meshfile = 'mesh.mesh';
end
if nargin < 3
    holemarker = 2;
end

% Generate a triangular mesh (Delaunay)
[p, e, t] = initmesh(geom, 'Hmax', 0.02, 'Jiggle', 'on');

numNodes = size(p,2);
num2DElements = size(t,2);

nodes = [1:numNodes; p; zeros(1,numNodes)];
elems = [1:num2DElements; t(1:3,:)];

% Boundary nodes
boundarynodes = unique([e(1,:),e(2,:)]);
bp = nodes(:,boundarynodes)';
bottom = bp( bp(:,3) == 0, 1 );
top = bp( bp(:,3) == 1, 1 );
left = bp( bp(:,2) == 0, 1 ); % only works for no shearing
right = bp( bp(:,2) == 1, 1 ); % only works for no shearing
numNodeBC = numel(left) + numel(right) + numel(bottom) + numel(top);

% Density
density = ones(1,size(t,2));
density(t(4,:)==holemarker) = 0; % holes

% Write density file
densityfile = [meshfile(1:end-4),'dens'];
Homogenization.writeDensityFile(density,densityfile);

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
fprintf(fid,'Num triangle     : %d\n', num2DElements);
fprintf(fid,'Num triangle,quad: 0\n');
fprintf(fid,'Num quadr        : 0\n');
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
    fprintf(fid,'%d 4 3 mech\n%d %d %d\n',elems);
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

function matrixToMesh(A,meshfile,shearingAngle)
% MATRIXTOMESH  -  Generates a (sparse) ANSYS mesh file out of a given density matrix.
%
% @param:
%       A               density matrix (entry has to be 0 for no material)
%       meshfile        name of generated .mesh file
%                       (optional, default = 'mesh.mesh')
%       shearingAngle   shearing angle of mesh in (-pi/2,pi/2)
%                       (optional, default = 0)
%

if nargin < 3
    shearingAngle = 0.0;
end

if nargin < 2
    meshfile = 'mesh.mesh';
end

A = flipud(A);
[m,n] = size(A);
hy = 1/m;
hx = 1/n;

% Get nodes
nodes = zeros((m+1)*(n+1),5);
nodeIdx = zeros(m+1,n+1);
nodeNum = 1;
for i = 1:m+1
    for j = 1:n+1
        % Check if a surrounding element has material
        if A(min(i,m),min(j,n)) || A(min(i,m),max(j-1,1)) || A(max(i-1,1),min(j,n)) || A(max(i-1,1),max(j-1,1))
            nodes(nodeNum,1:3) = [nodeNum,(j-1)*hx,(i-1)*hy];
            % Shearing
            nodes(nodeNum,2) = nodes(nodeNum,2) + (i-1) / m * tan(shearingAngle);
            % Set boundary flags
            if j == 1
                % left
                nodes(nodeNum,5) = 5;
            elseif j == n+1
                % right
                nodes(nodeNum,5) = 6;
            end
            nodeIdx(i,j) = nodeNum;
            nodeNum = nodeNum + 1;
        end
    end
end
numNodes = nodeNum - 1;
nodes = nodes( nodes(:,1) ~= 0, : );

% Boundary nodes
bottom = nodes( nodes(:,3) == 0, 1 );
top = nodes( nodes(:,3) == 1, 1 );
left = nodes( nodes(:,5) == 5, 1 );
right = nodes( nodes(:,5) == 6, 1 );
numNodeBC = numel(left) + numel(right) + numel(bottom) + numel(top);

% Get (quadratic) elements
elems = zeros(m*n,4);
c = 1;
for i = 1:m
    for j = 1:n
        if A(i,j)
            elems(c,:) = [nodeIdx(i,j),nodeIdx(i,j+1),nodeIdx(i+1,j+1),nodeIdx(i+1,j)];
            c = c + 1;
        end
    end
end
num2DElements = c - 1;
elems = elems( elems(:,1) ~= 0, : );

nodes = nodes(:,1:4)';
elems = elems';
elems = [1:num2DElements;elems];

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

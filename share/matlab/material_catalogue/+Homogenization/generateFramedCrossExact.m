function [ file, volume, dimension ] = generateFramedCrossExact(point,filepath,~)
% GENERATEFRAMEDCROSS  -  Generates a quadratic frame overlayed with an
% orthogonal cross, which is rotated by 45 degrees.
%
% @param:
%       point(1)  thickness of horizontal frame parts in [0,1]
%       point(2)  thickness of vertical frame parts in [0,1]
%       point(3)  thickness of the bar from the upper left corner to the
%                 lower right corner in normal direction in [0,1]
%       point(4)  thickness of the bar from the upper right corner to the
%                 lower left corner in normal direction in [0,1]
%       filepath  path to generated mesh file (optional)
%       nx        resolution of mesh (optional)
% 
% 
%       s2 
%      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%   s1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%      xxxxxx                           xxxxxxxx
%      xxxx xxx                       xxxxxxxxxx
%      xxxx   xxx                   xxxxxxx xxxx
%      xxxx     xxx               xxxxxxx   xxxx
%      xxxx        s3             s4        xxxx
%      xxxx         xxx       xxxxxxx       xxxx
%      xxxx           xxx   xxxxxxx         xxxx
%      xxxx             xxxxxxxxx           xxxx
%      xxxx             xxxoxxx             xxxx
%      xxxx           xxxxxxxxx             xxxx
%      xxxx         xxxxxxx   xxx           xxxx
%      xxxx       xxxxxxx       xxx         xxxx
%      xxxx     xxxxxxx           xxx       xxxx
%      xxxx   xxxxxxx               xxx     xxxx
%      xxxx xxxxxxx                   xxx   xxxx
%      xxxxxxxxxx                       xxx xxxx
%      xxxxxxxx                           xxxxxx
%      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% 

if nargin < 2
    filepath = '.';
end

if length(point) ~= 4
    throw( MException( 'generateFramedCross:wrongParameterCount', 'point does not contain four parameters.' ));
end

s1 = point(1);
s2 = point(2);
s3 = point(3);
s4 = point(4);

dimension = 2;

% Describe shape with Constructive Solid Geometry
% basicshape = [id; #edges; x-coords; y-coords; fill]
% Lower horizontal bar
rect11 = [3; 4; 0; 1; 1; 0; 0; 0; s1/2; s1/2; zeros(4,1)];
% Upper horizontal bar
rect12 = [3; 4; 0; 1; 1; 0; 1-s1/2; 1-s1/2; 1; 1; zeros(4,1)];
% Left vertical bar
rect21 = [3; 4; 0; s2/2; s2/2; 0; 0; 0; 1; 1; zeros(4,1)];
% Right vertical bar
rect22 = [3; 4; 1-s2/2; 1; 1; 1-s2/2; 0; 0; 1; 1; zeros(4,1)];
% Diagonal bar from upper left to lower right corner
diag11 = [2; 6; 0; 0; 1-s3/2; 1; 1; s3/2; 1; 1-s3/2; 0; 0; s3/2; 1];
diag12 = [2; 3; 0; 0; s3/2; s3/2; 0; 0; zeros(6,1)];
diag13 = [2; 3; 1; 1; 1-s3/2; 1-s3/2; 1; 1; zeros(6,1)];
% Diagonal bar from lower left to upper right corner
diag21 = [2; 6; 0; s4/2; 1; 1; 1-s4/2; 0; 0; 0; 1-s4/2; 1; 1; s4/2];
diag22 = [2; 3; 1-s4/2; 1; 1; 0; 0; s4/2; zeros(6,1)];
diag23 = [2; 3; s4/2; 0; 0; 1; 1; 1-s4/2; zeros(6,1)];

rect0 = [3; 4; 0; 1; 1; 0; 0; 0; 1; 1; zeros(4,1)];

% Combine shapes
gd = [rect11, rect12, rect21, rect22, diag11, diag12, diag13, diag21, diag22, diag23];

% Delete empty shapes
idx = csgchk(gd);
gd = gd(:,idx==0);

assert(~isempty(gd),'Assertion failed: Geometry is empty.')

% Get decomposed CSG
[geom,bt] = decsg(gd);
gd = [gd, rect0];
[fullgeom,fullbt] = decsg(gd);
% Remove subdomain boundaries
geom = csgdel(geom,bt);
fullgeom = csgdel(fullgeom,fullbt);

% Get volume
% We calculate the volume of each triangle by taking the crossproduct of
% two defining vectors.
if nargout == 2
    [p, ~, t] = initmesh(geom);
    v1 = p(:,t(2,:)) - p(:,t(1,:));
    v2 = p(:,t(3,:)) - p(:,t(1,:));
    volume = sum( abs(v1(1,:).*v2(2,:)-v1(2,:).*v2(1,:)) / 2 );
end

% Mark holes as subdomain
holemarker = max(max(geom(6:7,:))) + 1;
borders = ~(ismember(geom(1:end-2,:)',fullgeom(1:end-2,:)','rows')...
          | ismember(geom([1,3,2,5,4],:)',fullgeom(1:end-2,:)','rows'));
A = fullgeom(6:7,:);
A(A == 1) = holemarker;
fullgeom(6:7,:) = A;
A = geom(6:7,borders);
A(A == 0) = holemarker;
geom(6:7,borders) = A;
boundaries = ~(ismember(fullgeom(1:end-2,:)',geom(1:end-2,:)','rows')...
          | ismember(fullgeom(1:end-2,:)',geom([1,3,2,5,4],:)','rows'));
geom = [geom,fullgeom(:,boundaries)];


% Write mesh (and possible density) file
[~,filename] = fileparts(tempname);

% Get absolute path of mesh file
oldpath=pwd;
cd(filepath)
fullpath=pwd;
cd(oldpath)
filename = fullfile(fullpath,filename);
meshfile = [filename,'.mesh'];

file = Homogenization.geometryToMeshAndDensity(geom,meshfile,holemarker);

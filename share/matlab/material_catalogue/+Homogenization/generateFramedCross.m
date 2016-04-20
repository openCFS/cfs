function [ file, volume, dimension ] = generateFramedCross(point,filepath,nx)
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

if nargin < 3
    nx = 128;
end

if nargin < 2
    filepath = '.';
end

if length(point) ~= 4
    throw( MException( 'generateFramedCross:wrongParameterCount', 'point does not contain four parameters.' ));
end

density = zeros(nx);

% For periodic boundaries we generate the structure for a smaller grid
% first and copy the first row and column later.
nx = nx - 1;

s3 = point(3)*nx;
s4 = point(4)*nx;

dimension = 2;

% Cross
for i = 1:nx
    bar1 = mod( (i-round(s3/2) : i+round(s3/2)-1) - 1, nx ) + 1;
    density(i,bar1) = 1;
    bar2 = mod( (i-round(s4/2) : i+round(s4/2)-1) - 1, nx ) + 1;
    density(i,nx-bar2+1) = 1;
end

nx = nx + 1;

% Copy first row and column to the opposite boundary -> periodicity!
density(nx,:) = density(1,:);
density(:,nx) = density(:,1);

% Frame
s1 = point(1)*nx;
s2 = point(2)*nx;
horzBar = [1:round(s1/2),nx-round(s1/2)+1:nx];
vertBar = [1:round(s2/2),nx-round(s2/2)+1:nx];
density(horzBar,:) = 1;
density(:,vertBar) = 1;


volume = sum(sum(density))/nx^2;

% Write mesh (and possible density) file
[~,filename] = fileparts(tempname);

% Get absolute path of mesh file
oldpath=pwd;
cd(filepath)
fullpath=pwd;
cd(oldpath)
filename = fullfile(fullpath,filename);

file = Homogenization.matrixToMeshAndDensity(density,filename);

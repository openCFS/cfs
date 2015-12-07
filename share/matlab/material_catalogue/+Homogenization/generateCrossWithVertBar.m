function [ file, volume ] = generateCrossWithVertBar(point,filepath,nx)
% GENERATECROSSWITHVERTBAR  -  Generates an orthogonal cross, which is rotated
% by 45 degrees overlayed with vertical bars at the left and right edge.
%
% @param:
%       point(1)  thickness of horizontal frame parts
%       point(2)  thickness of vertical frame parts
%       point(3)  thickness of the bar from the upper left corner to the
%                 lower right corner in normal direction
%       point(4)  thickness of the bar from the upper right corner to the
%                 lower left corner in normal direction
%       filepath  path to generated mesh file (optional)
%       nx        resolution of mesh (optional)
% 
% 
%       s2 
%      xxxx                                 xxxx
%      xxxx                               xxxxxx
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
%      xxxx                                 xxxx
%      xxxx                                 xxxx
% 

if nargin < 3
    nx = 128;
end

if nargin < 2
    filepath = '.';
end

if length(point) ~= 3
    throw( MException( 'generateFramedCross:wrongParameterCount', 'point does not contain three parameters.' ));
end

density = zeros(nx);

s2 = point(1)*nx;
s3 = point(2)*nx;
s4 = point(3)*nx;

% Frame
% horzBar = [1:round(s1/2),nx-round(s1/2)+1:nx];
vertBar = [1:round(s2/2),nx-round(s2/2)+1:nx];
density(:,vertBar) = 1;

% Cross
for i = 1:nx
    bar1 = mod( (i-round(s3/2)+1 : i+round(s3/2)-1) - 1, nx ) + 1;
    density(i,bar1) = 1;
    bar2 = mod( (i-round(s4/2)+1 : i+round(s4/2)-1) - 1, nx ) + 1;
    density(i,nx-bar2+1) = 1;
end

volume = sum(sum(density))/nx^2;

% Write mesh (and possible density) file
[~,filename] = fileparts(tempname);

% Get absolute path of mesh file
oldpath=pwd;
cd(filepath)
fullpath=pwd;
cd(oldpath)
filename = fullfile(fullpath,filename);

file = Homogenization.matrixToMeshOrDensity(density,filename);

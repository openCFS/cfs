function [ file, volume ] = generateCross(point,filepath,nx)
% GENERATECROSS  -  Generates a vertical cross.
%
% @param:
%       point(1)  thickness of horizontal bar in [0,1]
%       point(2)  thickness of vertical bar in [0,1]
%       filepath  path to generated mesh file (optional)
%       nx        resolution of mesh (optional)
% 
% 
%                         s2
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
%   s1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
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
    nx = 128;
end

if nargin < 2
    filepath = '.';
end

if length(point) ~= 2
    throw( MException( 'generateFramedCross:wrongParameterCount', 'point does not contain two parameters.' ));
end

density = zeros(nx);

s1 = round(point(1)*nx);
s2 = round(point(2)*nx);

% To ensure periodic boundary conditions meshes are not allowed to have
% one and only one void row or column (which would be located at the 
% boundary). So if s1 or s2 would leave one and only one row or column
% empty, we increase the value to its maximum.
if s1 == nx - 1
    s1 = nx;
end
if s2 == nx - 1
    s2 = nx;
end

nxhalf = (nx+1)/2;

% Cross
horzBar = round( nxhalf - s1/2 + (1:s1) ) - 1;
vertBar = round( nxhalf - s2/2 + (1:s2) ) - 1;
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

file = Homogenization.matrixToMeshOrDensity(density,filename);

function [ file, volume, dimension ] = generateFrame(point,filepath,nx)
% GENERATEFRAME  -  Generates a quadratic frame.
%
% @param:
%       point(1)  thickness of horizontal frame parts in [0,1]
%       point(2)  thickness of vertical frame parts in [0,1]
%       filepath  path to generated mesh file (optional)
%       nx        resolution of mesh (optional)
% 
% 
%       s2 
%      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%   s1 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxx                                 xxxx
%      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%      xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
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

s1 = point(1)*nx;
s2 = point(2)*nx;

dimension = 2;

% Frame
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

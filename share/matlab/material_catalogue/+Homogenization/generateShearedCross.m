function [ file, volume, dimension ] = generateShearedCross(point,filepath,nx,ny)
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

% We first generate an unsheared cross like in generateCross.m and at the
% end pass the shearing angle to the mesh generation function.

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

while point(1) > 1 - 2/ny && point(1) ~= 1
    ny = ny*2;
end
while point(2) > 1 - 2/nx && point(2) ~= 1
    nx = nx*2;
end

density = zeros(ny,nx);


s1 = round(point(1)*ny);
s2 = round(point(2)*nx);
angle = (.5 - point(3))*pi;

% To ensure periodic boundary conditions meshes are not allowed to have
% one and only one void row or column (which would be located at the 
% boundary). So if s1 or s2 would leave one and only one row or column
% empty, we increase the value to its maximum.
% if s1 == nx - 1
%     s1 = nx;
% end
% if s2 == nx - 1
%     s2 = nx;
% end

nxhalf = (nx+1)/2;
nyhalf = (ny+1)/2;

% Cross
horzBar = round( nyhalf - s1/2 + (1:s1) ) - 1;
vertBar = round( nxhalf - s2/2 + (1:s2) ) - 1;
density(horzBar,:) = 1;
density(:,vertBar) = 1;

volume = sum(sum(density))/nx/ny;

dimension = 2;

% Write mesh (and possible density) file
[~,filename] = fileparts(tempname);

% Get absolute path of mesh file
oldpath=pwd;
cd(filepath)
fullpath=pwd;
cd(oldpath)
filename = fullfile(fullpath,filename);

file = Homogenization.matrixToMeshAndDensity(density,filename,angle);

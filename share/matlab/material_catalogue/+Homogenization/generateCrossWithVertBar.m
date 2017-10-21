function [ file, volume, dimension ] = generateCrossWithVertBar(point,filepath,nx)
% GENERATECROSSWITHVERTBAR  -  Generates an orthogonal cross, which is rotated
% by 45 degrees overlayed with a vertical bar.
%
% @param:
%       point(1)  thickness of vertical bar
%       point(2)  thickness of the bar from the upper left corner to the
%                 lower right corner in normal direction
%       point(3)  thickness of the bar from the upper right corner to the
%                 lower left corner in normal direction
%       filepath  path to generated mesh file (optional)
%       nx        resolution of mesh (optional)
% 
% 
%                          s2 
%      xx                xxxxx              xxxx
%       xxx              xxxxx            xxxxxx
%         xxx            xxxxx          xxxxxxx
%           xxx          xxxxx        xxxxxxx
%             xxx        xxxxx      xxxxxxx 
%               xxx      xxxxx    xxxxxxx   
%                  s3    xxxxx    s4        
%                   xxx  xxxxxxxxxxxx       
%                     xxxxxxxxxxxxx         
%                       xxxxxxxxx           
%                       xxxoxxx             
%                     xxxxxxxxx             
%                   xxxxxxxxxxxxx           
%                 xxxxxxxxxxxx  xxx         
%               xxxxxxx  xxxxx    xxx       
%             xxxxxxx    xxxxx      xxx     
%           xxxxxxx      xxxxx        xxx   
%         xxxxxxx        xxxxx          xxx 
%       xxxxxxx          xxxxx            xxx
%      xxxxxx            xxxxx              xxx
%      xxxx              xxxxx                xx
% 

if nargin < 3
    nx = 128;
end

if nargin < 2
    filepath = '.';
end

if length(point) ~= 3
    disp(point)
    throw( MException( 'generateFramedCross:wrongParameterCount', 'point does not contain three parameters.' ));
end

density = zeros(nx);

% For periodic boundaries we generate the structure for a smaller grid
% first and copy the first row and column later.
nx = nx - 1;

s2 = point(1)*nx;
s3 = point(2)*nx;
s4 = point(3)*nx;

dimension = 2;

% Frame
vertBar = round((nx-s2)/2)+1:round((nx+s2)/2);
density(:,vertBar) = 1;

% Cross
for i = 1:nx
    bar1 = mod( (i-round(s3/2)+1 : i+round(s3/2)-1) - 1, nx ) + 1;
    density(i,bar1) = 1;
    bar2 = mod( (i-round(s4/2)+1 : i+round(s4/2)-1) - 1, nx ) + 1;
    density(i,nx-bar2+1) = 1;
end

nx = nx + 1;

% Copy first row and column to the opposite boundary -> periodicity!
density(nx,:) = density(1,:);
density(:,nx) = density(:,1);

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

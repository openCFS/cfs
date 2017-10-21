function validateCFS()


% Sheared cross
points = zeros(9,3);

% Full material without shearing
points(1,:) = [1,1,.5];
% Full material with shearing
points(2,:) = [1,1,.25];
points(3,:) = [1,1,.95];
% Horizontal bar with width 0.5 without shearing
points(4,:) = [.5,0,.5];
% Horizontal bar with width 0.5 with shearing
points(5,:) = [.5,0,.25];
% Vertical bar with width 0.5 without shearing
points(6,:) = [0,.5,.5];
% Vertical bar with width 0.5 with shearing (has to be rotated for comparison)
points(7,:) = [0,.5,.25];
% Horizontal bar with width 0.001 with shearing
points(8,:) = [0.001,0,.25];
% Horizontal bar with width 0.999 with shearing
points(9,:) = [0.999,0,.25];

% % Framed Cross
% points = zeros(10,4);
% 
% % Full material 1
% points(1,:) = [1,0,0,0];
% % Full material 2
% points(2,:) = [0,1,0,0];
% % Full material 3
% points(3,:) = [0,0,1,0];
% % Full material 4
% points(4,:) = [0,0,0,1];
% % Horizontal bar with width 0.5
% points(5,:) = [.5,0,0,0];
% % Vertical bar with width 0.5
% points(6,:) = [0,.5,0,0];
% % Diagonal bar 1 with width 0.5
% points(7,:) = [0,0,.5,0];
% % Diagonal bar 1 with width 0.5
% points(8,:) = [0,0,0,.5];
% % Diagonal bar 1 with width 0.5
% points(9,:) = [0.01,0.01,.5,0];
% % Diagonal bar 1 with width 0.5
% points(10,:) = [0.01,0.01,0,.5];

% points(:,:) = [-0.500000,0.000000,-0.500000,0.500000];

% points = (points+1)/2;

sz = size(points);
% givenByLevelAndIndex, dimension, level, #points
points = [0, sz(2), 0, sz(1), zeros(1,sz(2)-4); points, zeros(sz(1), 4-sz(2))];

Homogenization.getInterpolationValues(points, 'tests')

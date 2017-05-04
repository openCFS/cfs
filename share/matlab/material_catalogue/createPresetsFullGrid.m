dim  = 2;
n = 3;
mini = .01;
maxi = .99;

npoints = 2^n-1;
sz = [ npoints^dim, dim ];

fprintf('number of grid points:  %d\n', sz(1));

pointsPerDim = 1:npoints;

if dim == 3
    x = reshape(repmat(pointsPerDim,npoints^2,1),[],1);
    y = repmat(reshape(repmat(pointsPerDim,npoints,1),[],1),npoints,1);
    z = repmat(pointsPerDim',npoints^2,1);
    points = [x,y,z];
elseif dim == 2
    x = reshape(repmat(pointsPerDim,npoints,1),[],1);
    y = repmat(pointsPerDim',npoints,1);
    points = [x,y];
end

points = points/(npoints+1)*(maxi-mini)+mini;

% Write coordinates
if sz(2) > 3
    data = [dim, npoints, sz(1), zeros(1,sz(2)-3); points];
elseif sz(2) == 3
    data = [dim, npoints, sz(1); points];
else
    data = [dim, npoints, sz(1); points, zeros(sz(1), 3-sz(2))];
end
% csvwrite('presets',data);
if ~exist('presets','dir')
    mkdir('presets')
end
dlmwrite( sprintf('presets/byCoords/presets%dD_%d',dim,npoints), data, 'delimiter', ',', 'precision', '%.10f' );

% Clean up
clear

dim  = 3;
npoints = 2^8-1; % 2^n-1
sz = [ npoints^dim, dim ];

% npoints = 6;
% sz = [ 6^2, dim ];

fprintf('number of grid points:  %d, 4096\n', sz(1));

pointsPerDim = 1:npoints;

% points = zeros(sz);
% idx = 1;
% for i=pointsPerDim
%     for j=pointsPerDim
%         for k=pointsPerDim
% %             points(idx,:) = [ 2/128*2-1, 2/128*2-1, i/(npoints+1)*2-1, j/(npoints+1)*2-1 ];
% %             points(idx,:) = [ 2/128*2-1, 2/128*2-1, i*2-1, j*2-1 ];
%             points(idx,:) = [ i, j, k ];
%             idx = idx + 1;
%         end
%     end
% end

x = reshape(repmat(pointsPerDim,npoints^2,1),[],1);
y = repmat(reshape(repmat(pointsPerDim,npoints,1),[],1),npoints,1);
z = repmat(pointsPerDim',npoints^2,1);

points = [x,y,z];
points = points/(npoints+1)*2-1;

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
dlmwrite( sprintf('presets/presets%dD_%d',dim,npoints), data, 'delimiter', ',', 'precision', '%f' );

% Clean up
clear

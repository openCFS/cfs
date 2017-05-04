dim  = 3;
level = 7;
blevel = 3;

% Get current path
path = fileparts(which('createPresetsSparseGrid.m'));

% Generate points with python script
tmp = pwd;
cd('generate_points/');
[status,result] = system( sprintf('python generate_points.py %d %d',dim,level) );
disp(result);
if status ~= 0
    disp('Fehler beim Aufruf von generate_points.py');
    return;
end
cd(tmp);

% Read levels and indices
format = ['[ ',repmat('%d, %d, ',1,dim-1),'%d, %d ]\n'];
fid = fopen([path,'/generate_points/grid_points.csv']);
A = fscanf(fid,format);
fclose(fid);

A = reshape(A,2*dim,size(A,1)/2/dim)';

grid.Xl = uint64(A(:,1:2:2*dim));
grid.Xi = uint64(A(:,2:2:2*dim));

% Calculate coordinates
grid.X = pow2(-double(grid.Xl)) .* double(grid.Xi) * 2 - 1;

% Map to intervall from 0 to 1
grid.X = (grid.X + 1) / 2;
sz = size(grid.X);

% Write coordinates
coords = [0, dim, level, sz(1), zeros(1,sz(2)-4); grid.X, zeros(sz(1), 4-sz(2))];
A = [1, dim, level, sz(1), zeros(1,2*dim-4); A, zeros(sz(1), 4-2*dim)];
% csvwrite('presets',data);
if ~exist('presets','dir')
    mkdir('presets')
end
if ~exist('presets/byCoords','dir')
    mkdir('presets/byCoords')
end
if ~exist('presets/byLevelAndIndex','dir')
    mkdir('presets/byLevelAndIndex')
end
dlmwrite( sprintf('presets/byCoords/presets%dD_L%d_b%d',dim,level,blevel), coords, 'delimiter', ',', 'precision', '%.10f' );
dlmwrite( sprintf('presets/byLevelAndIndex/presets%dD_L%d_b%d',dim,level,blevel), A, 'delimiter', ',', 'precision', '%d' );

% Clean up
clear tmp format status result fid ans A sz
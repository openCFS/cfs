dim  = 3;
level = 3;

% Generate points with python script
tmp = pwd;
cd('generate_points/');
[status,result] = system( sprintf('python generate_points.py %d %d',dim,level) );
if status ~= 0
    disp('Fehler beim Aufruf von generate_points.py');
    disp(result);
    return;
end
cd(tmp);

% Read levels and indices
format = ['[ ',repmat('%d, %d, ',1,dim-1),'%d, %d ]\n'];
fid = fopen('/home/daniel/Masterarbeit/generate_points/grid_points.csv');
A = fscanf(fid,format);
fclose(fid);

A = reshape(A,2*dim,[]);

grid.Xl = uint64(A(1:2:2*dim,:)');
grid.Xi = uint64(A(2:2:2*dim,:)');

% Calculate coordinates
grid.X = pow2(-double(grid.Xl)) .* double(grid.Xi) * 2 - 1;
sz = size(grid.X);

% Write coordinates
data = [dim, level, sz(1), zeros(1,sz(2)-3); grid.X, zeros(sz(1), 3-sz(2))];
% csvwrite('presets',data);
if ~exist('presets','dir')
    mkdir('presets')
end
dlmwrite( sprintf('presets/presets%dD_L%d',dim,level), data, 'delimiter', ',', 'precision', '%f' );

% Clean up
clear tmp format fid A B C Proj

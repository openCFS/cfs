clc

% MICRO MESH
lxmicro = 1;
lymicro = 1;

dim = 4;

cfsWorkingDirectory = '/home/daniel/Masterarbeit/Matlab/+Homogenization/CFS_Working_Directory';
        
for level = 4:4
    fprintf('\n --------------------------------------------------\n');
    fprintf('Dimension: %d, Level: %d\n',dim,level);

    % MICRO MESH
    nelxmicro = 2^max(7,level);
    nelymicro = 2^max(7,level);
    % GET INTERPOLANTS FOR ELASTICITY TENSOR
    try
        tic;
        file = sprintf('%s/grid_dim_%d_level_%d',cfsWorkingDirectory,dim,level);
        time = Homogenization.computeInterpolationValues(nelxmicro,nelymicro,file);
        toc;
    catch ME
        fprintf('line %d: %s\n', ME.stack.line, ME.message);
    end
end

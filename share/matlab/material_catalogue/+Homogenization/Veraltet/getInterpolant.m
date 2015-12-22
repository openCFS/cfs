addpath inviso/data
addpath inviso/fem
addpath inviso/snopt
tic;
% MICRO MESH
nelxmicro = 20;
nelymicro = 20;
lxmicro = 1;
lymicro = 1;

fhandle = @createCross;

% GET INTERPOLANTS FOR ELASTICITY TENSOR
Eh = computeInterpolant(fhandle,nelxmicro,nelymicro,lxmicro,lymicro);

str = strcat('Eh',func2str(fhandle),num2str(nelxmicro));

save(str,'Eh')
toc;
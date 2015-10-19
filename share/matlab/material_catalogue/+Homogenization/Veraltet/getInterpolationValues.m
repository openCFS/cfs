clc

% MICRO MESH
lxmicro = 1;
lymicro = 1;

dim = 3;

for level = 2:2
    fprintf('\n --------------------------------------------------\n');
    fprintf('Dimension: %d, Level: %d\n',dim,level);

    % MICRO MESH
    nelxmicro = 2^max(7,level);
    nelymicro = 2^max(7,level);
%     presetsa = 1/2^level:1/2^level:1;
%     presetsb = presetsa; 
    presetsa = .1:.1:1;
    presetsb = presetsa;
    % presetsphi = 1/2^level:1/2^level:1-(1/2^level);
    presetsphi = .1:.1:.9;
    
    % GET INTERPOLANTS FOR ELASTICITY TENSOR
    try
        tic;
        time = Homogenization.computeInterpolationValues(nelxmicro,nelymicro,presetsa,presetsb,presetsphi);
        toc;
    catch ME
        disp(ME.message);
    end
end


clc

cfsWorkingDirectory = '/home/daniel/Masterarbeit/Matlab/+Homogenization/CFS_Working_Directory';

referenceRadius = .5;

for level = 6
    fprintf('\n --------------------------------------------------\n');
    fprintf('Dimension: 2, Level: %d\n',level);

    % GET INTERPOLANTS FOR ELASTICITY TENSOR
    try
        tic;
        file = sprintf('presets8D_%d',level);
        time = Homogenization.computeInterpolationValuesForHoleblock(file, referenceRadius);
        toc;
    catch ME
        fprintf('line %d: %s\n', ME.stack.line, ME.message);
    end
end

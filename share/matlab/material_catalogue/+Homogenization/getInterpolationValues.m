function getInterpolationValues(gridfile, threadID)

% meshgenerationfunc = @Homogenization.generateFrame;
% meshgenerationfunc = @Homogenization.generateCross;
% meshgenerationfunc = @Homogenization.generateShearedCross;
meshgenerationfunc = @Homogenization.generateShearedCrossExact;
% meshgenerationfunc = @Homogenization.generateFramedCross;
% meshgenerationfunc = @Homogenization.generateFramedCrossExact;
    
cfsworkingdirectory = '/home/daniel/Masterarbeit/Matlab/+Homogenization/CFS_Working_Directory';

% GET INTERPOLANTS FOR ELASTICITY TENSOR
try
    tic;
    createcataloguetime = Homogenization.computeInterpolationValues(gridfile, meshgenerationfunc, cfsworkingdirectory, threadID);
    toc;
catch ME
    fprintf('line %d: %s\n', ME.stack.line, ME.message);
end

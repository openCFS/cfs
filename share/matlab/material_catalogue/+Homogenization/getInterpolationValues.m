function getInterpolationValues(gridfile, threadID)

% meshgenerationfunc = @Homogenization.generateFrame;
% meshgenerationfunc = @Homogenization.generateCross;
meshgenerationfunc = @Homogenization.generateShearedCross;
% meshgenerationfunc = @Homogenization.generateShearedCrossExact;
% meshgenerationfunc = @Homogenization.generateFramedCross;
% meshgenerationfunc = @Homogenization.generateFramedCrossExact;

% Get current path
path = fileparts(which('+Homogenization/getInterpolationValues.m'));

cfsworkingdirectory = [path,'/CFS_Working_Directory'];
% GET INTERPOLANTS FOR ELASTICITY TENSOR
try
    tic;
    createcataloguetime = Homogenization.computeInterpolationValues(gridfile, meshgenerationfunc, cfsworkingdirectory, threadID);
    toc;
catch ME
    fprintf('%s\n', ME.getReport());
end

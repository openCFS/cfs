function [createcataloguetime] = computeInterpolationValues(points, meshgenerationfunc, cfsworkingdirectory, threadID)
% COMPUTEINTERPOLATIONVALUES Computes values of homogenized elasticity
%                            tensor using CFS in VOIGT notation.
% 
%    computeInterpolationValues(points, meshgenerationfunc, cfsworkingdirectory, threadID)
%    computes the homogenized elasticity tensors for the data matrix points.
%    Each row of points contains one data point. meshgenerationfunc is a
%    handle to a function for the generation of the mesh of the micro cell.
%    cfsworkingdirectory is a directory containing a mat.xml and
%    inv_tensor.xml file where CFS++ will work. The resulting tensors are
%    stored in 'catalogues/detailed_stats_threadID'.
% 
%    computeInterpolationValues(file, ...) computes the tensors for the
%    data points given in file. The values in the file are ranged in [-1;1]
%    and will be mapped to [0;1]! Data files can be created with the
%    scripts createPresetsFullGrid and createPresetsSparseGrid.
%
%    computeInterpolationValues(data, ...) computes the tensors for the
%    data points given in the matrix data. The values in the matrix are
%    ranged in [0;1], no mapping applies!
% 
% 
% Example:
% 
%    points = rand(5,3);
%    meshgenerationfunc = @Homogenization.generateShearedCrossExact;
%    cfsworkingdirectory = '/home/daniel/code/cfs/share/matlab/material_catalogue/+Homogenization/CFS_Working_Directory';
%    threadID = 1;
%    
%    Homogenization.computeInterpolationValues(points, meshgenerationfunc, cfsworkingdirectory, threadID)
% 

% At the moment the catalogue generation only works on unix systems. For
% other platforms use C = computer; if (C == ... and modify
% Homogenization.getElasticityTensorOfMicroCell

if ~isunix
    warning('computeInterpolationValues:NoUnixSystem',...
        'computeInterpolationValues using CFS++ only works on UNIX systems!');
    return;
end

if ischar(points)
    % Read gridpoints from file
    file = points;
    gridPoints = load(file);
else
    % Gridpoints are given as input matrix
    gridPoints = points;
end

% First line contains dimension and level (and maybe number of points)
if gridPoints(1,1) > 1
    dim = gridPoints(1,1);
    level = gridPoints(1,2);
    gridPoints(1,:) = [];
else
    dim = size(gridPoints,2);
    level = 1;
end

if size(gridPoints,2) > dim
    gridPoints = gridPoints(:,1:dim);
end

% Map to intervall from 0 to 1
if ischar(points)
    gridPoints = (gridPoints + 1) / 2;
end

numPoints = size(gridPoints,1);

% Check if catalogue fits to meshgenerationfunc by calling the function once
try
    meshfile = meshgenerationfunc(gridPoints(1,:), '.');
    [meshfilepath, meshfilename] = fileparts(meshfile);
    if exist( sprintf('%s/%s.dens', meshfilepath, meshfilename), 'file' )
        delete( sprintf('%s/%s.dens', meshfilepath, meshfilename) );
    end
    delete( sprintf('%s/%s.mesh', meshfilepath, meshfilename) );
catch ME
    fprintf('line %d: %s\n', ME.stack.line, ME.message);
    createcataloguetime = -1;
    return
end

Efull = Homogenization.getElasticityTensorOfMaterial( sprintf('%s/inv_tensor.xml',cfsworkingdirectory) );

%if ~isempty(level)
%    id = strcat(num2str(dim),'D_level_',num2str(level));
%else
%    id = strcat(num2str(dim),'D_numpoints_',num2str(numPoints));
%end

% Create folder 'catalogues' if necessary
if ~exist('catalogues','dir')
    mkdir('catalogues');
end    

% Start timer
createcataloguetime = tic;

% Compute homogenized elasticity tensors
Tensors = zeros(numPoints,7);
cache = [];
if ischar(threadID)
    filename = strcat('catalogues/detailed_stats_', threadID);
else
    filename = strcat('catalogues/detailed_stats_', num2str(threadID));
end
fid = fopen(filename,'wt');
% fprintf(fid,'%dD\tL%d\t%d\t%s\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n',dim,level,numPoints,'voigt',zeros(1,11-dim));
for i=1:numPoints
    point = gridPoints(i,:);
%     Take advantage of symmetry
%     [ism, idx] = ismemberf(point([1 3 2]),gridPoints(1:i-1,:),'rows');
%     disp(ism)
%     if ism && ~ismemberf(point([2 1]),cache,'rows');
%         Tensor = zeros(3,3);
%         vec = Tensors(idx,:);
%         Tensor(1,2:3) = vec(2:3);
%         Tensor(2,3) = vec(5);
%         Tensor = Tensor + Tensor';
%         Tensor(1,1) = vec(1);
%         Tensor(2,2) = vec(4);
%         Tensor(3,3) = vec(6);
%         % rotate tensor
%         Eh = rotatevoigt(Tensor, pi/2);
%         volume = vec(end);
%     else
        [Eh, volume, meshfilename] = Homogenization.getElasticityTensorOfMicroCell(point, meshgenerationfunc, Efull, cfsworkingdirectory);
%     end
    if ~isempty(Eh)
        Tensors(i,:) = [Eh(1,1) Eh(1,2) Eh(1,3) Eh(2,2) Eh(2,3) Eh(3,3) volume];
        fprintf(fid,'%.10f\t',point);
        fprintf(fid,'%e\t%e\t%e\t',Eh(1,1),Eh(1,2),Eh(1,3));
        fprintf(fid,'%e\t%e\t%e\t',Eh(2,2),Eh(2,3),Eh(3,3));
        fprintf(fid,'%e\n',volume);
    else
        cache = [cache;point];
        warning( strcat('No tensor available for point', sprintf('\t%.6f',point)) );
    end

    showProgress(i,numPoints);
    if abs(volume-1) >= 1e-14
        delete( sprintf('%s/inv_tensor_%s.info.xml', cfsworkingdirectory, meshfilename) );
    end
end

% Write failed points to file
if ~isempty(cache)
    if ischar(points)
        fprintf(fid,'%.10f,%.10f,%.10f,%.10f\n',cache'*2-1);
    else
        fprintf(fid,'%.10f,%.10f,%.10f,%.10f\n',cache');
    end
end
fclose(fid);

% End timer
createcataloguetime = toc(createcataloguetime);

fprintf('material catalogue file:  %s\n', filename);
fprintf('time for creating catalogue:     %f s\n', createcataloguetime);

end



function showProgress(counter,complete)
progress = round(counter/complete*100);
% clc;
fprintf('Progress: %d %%\n',progress);
end

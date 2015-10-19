function [cfsTime] = computeInterpolationValues(nelx,nely,arg1,arg2,arg3)
% COMPUTEINTERPOLATIONVALUES Computes values of homogenized elasticity
%                            tensor using CFS.
%     cfsTime = computeInterpolationValues(nelx,nely,file)
%     computes the tensors for the data points given in file on a micromesh
%     with resolution nelx x nely.
%
%     cfsTime = computeInterpolationValues(nelx,nely,presetsa,presetsb)
%     computes the tensors for the data points given by the presetsa and 
%     presetsb on a micromesh with resolution nelx x nely. Phi is
%     automatically assigned the value 0.5.
%
%     cfsTime = computeInterpolationValues(nelx,nely,presetsa,presetsb,presetsphi)
%     computes the tensors for the data points given by the presets on a
%     micromesh with resolution nelx x nely.
%
%     See also Homogenization.generate_mesh
%

if ~isunix
    warning('computeInterpolationValues:NoUnixSystem',...
        'computeInterpolationValues using CFS only works on UNIX systems!');
    return;
end

if nargin < 5
    arg3 = 0.5;
end

cfsWorkingDirectory = '/home/daniel/Masterarbeit/Matlab/+Homogenization/CFS_Working_Directory';

% Remove in further versions
dim = 4;

if nargin < 4
    if ischar(arg1)
        % Read gridpoints from file
        file = arg1;
        data = load(file);
    else
        % Gridpoints are given as input matrix
        data = arg1;
    end
    
    % If only two columns are given we have to add the third
    if size(data,2) == 2
        data = [ data,arg3*ones(size(data,1),1) ];
    end
    
    % First line contains dimension and level (and maybe number of points)
    if data(1,1) > 1
        level = data(1,2);
        data(1,:) = [];
    else
        level = getLevel(data);
    end
else
    level = getLevel(arg1,arg2);
    % Presets of points are given
    [tmp1,tmp2,tmp3] = meshgrid(arg1,arg2,arg3);
    numPoints = numel(tmp1);
    data = zeros(numPoints,3);
    for i = 1:numPoints
        data(i,:) = [tmp1(i),tmp2(i),tmp3(i)];
    end
end

numPoints = size(data,1);

% Set accumulated time
preprocessing = 0;
cfs = 0;
postprocessing = 0;
readFromXML = 0;
writeStats = 0;

Ehf = Homogenization.getElasticityTensor(sprintf('%s/inv_tensor.xml',cfsWorkingDirectory));

% Compute homogenized elasticity tensors
Tensors = cell(numPoints,1);
parfor i=1:numPoints
    vec = data(i,:);
    a = vec(1);
    b = vec(2);
    phi = vec(3);
            
    % For full material without shearing Eh equals the elasticity tensor.
    if a / cos((phi-.5)*pi) >= 1
        Tensors{i} = Ehf;
        continue
    end

    tpreprocessing = tic;
    % Generate sparse mesh
    meshFile = Homogenization.generate_mesh(cfsWorkingDirectory,nelx,nely,a,b,phi);
    [~,meshFileName] = fileparts(meshFile);
    % Call CFS++ to calculate homogenized tensor
    path = pwd;
    cd(cfsWorkingDirectory);
    invFileName = sprintf('inv_tensor_%s',meshFileName);
    invFile = strcat(invFileName,'.xml');
    [status,message] = copyfile('inv_tensor.xml',invFile);
    if status ~= 1
        disp('Fehler beim Kopieren von inv_tensor.xml');
        disp(message);
    end
    preprocessing = preprocessing + toc(tpreprocessing);
    tCfs = tic;
    [status,result] = unix( sprintf('LC_ALL=C cfs.rel -m %s %s',meshFile,invFileName) );
    cfs = cfs + toc(tCfs);
    cd(path);
    if status ~= 0
        disp('Fehler in CFS');
        disp(result);
    else
        tReadFromXML = tic;
        % Read homogenized tensor from xml file
        Eh = Homogenization.read_matrix_from_xml(...
            sprintf('%s/%s.info.xml',cfsWorkingDirectory,invFileName) );
        Tensors{i} = Eh;
        readFromXML = readFromXML + toc(tReadFromXML);
    end
    tpostprocessing = tic;
    delete(meshFile);
    delete( sprintf('%s/%s',cfsWorkingDirectory,invFile) );
    postprocessing = postprocessing + toc(tpostprocessing);
%    showProgress(counter,complete);
end
delete( sprintf('%s/inv_tensor_*',cfsWorkingDirectory) );

cfsTime = cfs;

if ~isempty(level)
    id = strcat(num2str(dim),'D_level_',num2str(level));
else
    id = strcat(num2str(dim),'D_numpoints_',num2str(numPoints));
end

% Save homogenized tensors in a mat file
% save( strcat('Tensors_',id), 'Tensors' );

% Write tensors to a file
twriteStats = tic;
filename = strcat('detailed_stats_',id);
fid = fopen(filename,'wt');
fprintf(fid,'%dD\tL%d\t%d\t%s\t%e\t%e\t%e\t%e\t%e\n',dim,level,numPoints,'voigt',0,0,0,0,0);
for i=1:numPoints
    fprintf(fid,'%.10f\t%.10f\t%.10f\t',data(i,1),data(i,2),data(i,3));
    Eh = Tensors{i};
    try
        fprintf(fid,'%e\t%e\t%e\t',Eh(1,1),Eh(1,2),Eh(1,3));
        fprintf(fid,'%e\t%e\t%e\n',Eh(2,2),Eh(2,3),Eh(3,3));
    catch ME
        disp(ME.message);
    end
end
fclose(fid);
writeStats = writeStats + toc(twriteStats);

fprintf('time for preprocessing:   %f s\n',preprocessing);
fprintf('time spent in cfs:        %f s\n',cfs);
fprintf('time for postprocessing:  %f s\n',postprocessing);
fprintf('time for reading xml:     %f s\n',readFromXML);
fprintf('time for writing stats:   %f s\n',writeStats);

end



function showProgress(counter,complete)
progress = round(counter/complete*100);
clc;
fprintf('Progress: %d %%',progress);
end



function level = getLevel(arg1,arg2)
% Helper method to get level of data discretization
level = [];

if nargin == 2
    n1 = numel(arg1);
    n2 = numel(arg2);
    % Check if number of points is equal in both directions
    if n1 ~= n2
        return
    end
    
    % Check if points are equidistant
    if any( abs( sort(arg1)-linspace(min(arg1),max(arg1),numel(arg1) ) ) > eps )
        return
    end
    if any( abs( sort(arg2)-linspace(min(arg2),max(arg2),numel(arg2) ) ) > eps )
        return
    end
    
    % Get level
    [f,e] = log2(n1);
    if abs(f-.5)<eps
        level = e-1;
    end    
else
    % Extract slice of data with most points (usually mid slice)
    slices = unique(arg1(:,3));
    pointsperslice = zeros(1,numel(slices));
    for i=1:numel(slices)
        pointsperslice(i) = numel( find( arg1(:,3)==slices(i) ) );
    end
    numPoints = max(pointsperslice);
    % Check if number of points in slice equals number of points in a two
    % dimensional sparse grid
    k = 1;
    tmp = 0;
    while tmp < numPoints
        tmp = 2^k*(k-1)+1;
        k = k+1;
    end
    if tmp == numPoints
        level = k-1;
    end
end

end

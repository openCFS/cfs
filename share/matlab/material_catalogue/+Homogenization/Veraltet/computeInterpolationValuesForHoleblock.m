function [cfsTime] = computeInterpolationValuesForHoleblock(arg1,referenceRadius)
% COMPUTEINTERPOLATIONVALUES Computes values of homogenized elasticity
%                            tensor using CFS.
%     cfsTime = computeInterpolationValuesForHoleblock(file)
%     computes the tensors for the data points given in file.
%
%     cfsTime = computeInterpolationValuesForHoleblock(data)
%     computes the tensors for the data points given in the matrix data.
% 

if ~isunix
    warning('computeInterpolationValues:NoUnixSystem',...
        'computeInterpolationValues using CFS only works on UNIX systems!');
    return;
end

cfsWorkingDirectory = '/home/daniel/Masterarbeit/Matlab/+Homogenization/CFS_Working_Directory';

addpath(genpath('/home/daniel/examples/shape/shapeopt'));

if ischar(arg1)
    % Read gridpoints from file
    file = arg1;
    data = load(file);
else
    % Gridpoints are given as input matrix
    data = arg1;
end

% First line contains dimension and level (and maybe number of points)
if data(1,1) > 1
    level = data(1,2);
    data(1,:) = [];
else
    level = 1;
end
    
numPoints = size(data,1);
dim = size(data,2);

phi = 2 * pi * 2/dim * (0:dim/2-1);
reference = reshape( referenceRadius * [cos(phi-pi/4);sin(phi-pi/4)], 1, dim);

% Set accumulated time
preprocessing = 0;
cfs = 0;
postprocessing = 0;
readFromXML = 0;
writeStats = 0;

Ehf = Homogenization.getElasticityTensor(sprintf('%s/inv_tensor.xml',cfsWorkingDirectory));

% Compute homogenized elasticity tensors
Tensors = cell(numPoints,1);
cerr = zeros(1,numPoints);
for i=1:numPoints
    cp = data(i,:) + reference;
    
    % For full material without shearing Eh equals the elasticity tensor.
    if (numel(unique(cp(1:2:end))) == 1) && (numel(unique(cp(2:2:end))) == 1)
        Tensors{i} = Ehf;
        continue
    end

    tpreprocessing = tic;
    % Generate sparse mesh
    try
        meshFile = holeblockmesh(cfsWorkingDirectory,cp,i);
    catch ME
        fprintf('Fehler beim Erzeugen des Mesh: %s\n',ME.message);
%         fprintf('In line %d of function %s.\n',ME.stack(1).line,ME.stack(1).file);
        Tensors{i} = zeros(3);
        cerr(i) = 1;
        continue
    end
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
        disp('Fehler in CFS:');
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
    delete( meshFile );
    delete( sprintf('%s/%s',cfsWorkingDirectory,invFile) );
    delete( sprintf('%s/inv_tensor_%s.las',cfsWorkingDirectory,meshFileName) );
    delete( sprintf('%s/inv_tensor_%s.info',cfsWorkingDirectory,meshFileName) );
    delete( sprintf('%s/inv_tensor_%s.plot.dat',cfsWorkingDirectory,meshFileName) );
    postprocessing = postprocessing + toc(tpostprocessing);
    showProgress(i,numPoints);
end
clear('showProgress')
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
fprintf(fid,'%dD\tL%d\t%d\t%s\t',dim,level,numPoints,'voigt');
fprintf(fid,'%e\t',zeros(1,size(data,2)+6-4));
fprintf(fid,'\n');
for i=1:numPoints
    fprintf(fid,'%.10f\t',data(i,:));
    Eh = Tensors{i};
    try
        fprintf(fid,'%e\t%e\t%e\t',Eh(1,1),Eh(1,2),Eh(1,3));
        fprintf(fid,'%e\t%e\t%e\n',Eh(2,2),Eh(2,3),Eh(3,3));
    catch ME
        fprintf('Fehler beim Schreiben der Tensoren in den Katalog: %s\n',ME.message);
        fprintf('In line %d of function %s.\n',ME.stack(1).line,ME.stack(1).file);
    end
end
fclose(fid);
writeStats = writeStats + toc(twriteStats);

fprintf('Errors occurred during mesh generation: %d out of %d data points\n',sum(cerr),numPoints);
fprintf('time for preprocessing:   %f s\n',preprocessing);
fprintf('time spent in cfs:        %f s\n',cfs);
fprintf('time for postprocessing:  %f s\n',postprocessing);
fprintf('time for reading xml:     %f s\n',readFromXML);
fprintf('time for writing stats:   %f s\n',writeStats);

end

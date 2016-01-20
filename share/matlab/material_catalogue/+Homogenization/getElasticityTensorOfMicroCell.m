function [Tensor, volume, meshfilename] = getElasticityTensorOfMicroCell(point, meshgenerationfunc, Efull, cfsworkingdirectory)
% GETELASTICITYTENSOROFMICROCELL Computes a homogenized elasticity tensor
%                                using CFS in VOIGT notation.
% 
%    getElasticityTensorOfMicroCell(point, meshgenerationfunc, Efull, cfsworkingdirectory)
%    computes the homogenized elasticity tensor for the data vector point.
%    meshgenerationfunc is a handle to a function for the generation of the
%    mesh of the micro cell. Efull is the tensor of full material which can
%    be obtained by Homogenization.getElasticityTensorOfMaterial.
%    cfsworkingdirectory is a directory containing a mat.xml and
%    inv_tensor.xml file where CFS++ will work.
% 
%    [Tensor, volume, meshfilename] = getElasticityTensorOfMicroCell(...)
%    returns the homogenized elasticity tensor, the relative volume of the
%    material in the micro cell and the name of the generated meshfile
%    (without extension).
% 
% Example:
% 
%    point = rand(1,3);
%    meshgenerationfunc = @Homogenization.generateShearedCrossExact;
%    cfsworkingdirectory = '/home/daniel/code/cfs/share/matlab/material_catalogue/+Homogenization/CFS_Working_Directory';
%    Efull = Homogenization.getElasticityTensorOfMaterial(sprintf('%s/inv_tensor.xml',cfsworkingdirectory));
%    
%    getElasticityTensorOfMicroCell(point, meshgenerationfunc, Efull, cfsworkingdirectory)
% 

% Generate sparse mesh
[meshfile, volume] = meshgenerationfunc(point, cfsworkingdirectory);
%close;Homogenization.plotmesh(meshfile);
[meshfilepath, meshfilename, ext] = fileparts(meshfile);
% For an empty grid Eh equals the all zero tensor.
if volume < 1e-14
    Tensor = zeros(3,3);
    delete( sprintf('%s/%s.dens', meshfilepath, meshfilename) );
    delete( meshfile );
    return
end
% For full material Eh equals the elasticity tensor.
if 1-volume < 1e-14
    Tensor = Efull;
    delete( meshfile );
    return
end
% Homogenization.plotmesh(meshfile);
% Call CFS++ to calculate homogenized tensor
path = pwd;
cd(cfsworkingdirectory);
invfilename = sprintf('inv_tensor_%s', meshfilename);
invfile = strcat(invfilename, '.xml');
[status,message] = copyfile('inv_tensor.xml', invfile);
if status ~= 1
    disp('Fehler beim Kopieren von inv_tensor.xml');
    disp(message);
end
if strcmp(ext,'.mesh')
    [status,result] = system( sprintf('LC_ALL=C cfs.rel -m %s %s', meshfile, invfilename) );
else
    [status,result] = system( sprintf('LC_ALL=C cfs.rel -m %s.mesh -x %s %s', meshfilename, meshfile, invfilename) );
end
cd(path);

if status ~= 0
    disp('Fehler in CFS');
    disp(point)
    disp(result);
    Tensor = [];
else
    % Read homogenized tensor from xml file
    try
        Eh = Homogenization.read_matrix_from_xml( sprintf('%s/%s.info.xml', cfsworkingdirectory, invfilename) );
        Tensor = Eh;
    catch ME
        Tensor = [];
        disp(ME.message);
    end
end
densfile = sprintf('%s/%s.dens', meshfilepath, meshfilename);
infofile = sprintf('%s/%s.info', cfsworkingdirectory, invfilename);
if exist(densfile, 'file')
    delete(densfile);
end
delete( sprintf('%s/%s.mesh', meshfilepath, meshfilename) );
delete( sprintf('%s/%s', cfsworkingdirectory, invfile) );
delete( sprintf('%s/%s.las', cfsworkingdirectory, invfilename) );
if exist(infofile, 'file')
    delete(infofile);
end
delete( sprintf('%s/%s.plot.dat', cfsworkingdirectory, invfilename) );
delete( sprintf('%s/results_hdf5/%s.h5', cfsworkingdirectory, invfilename) );

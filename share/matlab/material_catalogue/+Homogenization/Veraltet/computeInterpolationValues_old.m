function [cfsTime] = computeInterpolationValues_old(nelx,nely,presetsa,presetsb,presetsphi)
% Computes interpolation values of Eh using CFS.

if ~isunix
    warning('computeInterpolationValues:NoUnixSystem',...
        'computeInterpolationValues using CFS only works on UNIX systems!');
    return;
end

cfsWorkingDirectory = '/home/daniel/Masterarbeit/Matlab/+Homogenization/CFS_Working_Directory';

preprocessing = 0;
cfs = 0;
postprocessing = 0;
readFromXML = 0;
writeStats = 0;

npresetsa = numel(presetsa);
npresetsb = numel(presetsb);
npresetsphi = numel(presetsphi);

% Compute some homogenized elasticity tensors
Tensors = cell(npresetsa,npresetsb,npresetsphi);
for loopa = 1:npresetsa
    a = presetsa(loopa);
    for loopb = 1:npresetsb
        b = presetsb(loopb);
        parfor loopphi = 1:npresetsphi;
            
            % SPEEDHACK
            % at this time only works for 99line material
            if a == 1 || b == 1
                Eh = [ 1, .3, 0; .3, 1, 0; 0, 0, (1-.3)/2 ]./(1-.3^2);
                Tensors{loopa,loopb,loopphi} = Eh;
                continue
            end
            
            phi = presetsphi(loopphi);
            tpreprocessing = tic;
            % Generate sparse mesh
            meshFile = Homogenization.generate_mesh(cfsWorkingDirectory,nelx,nely,a,b,phi)
            [~,meshFileName] = fileparts(meshFile);
            % Call CFS to calculate homogenized tensor
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
                Eh = Homogenization.read_matrix_from_xml( sprintf('%s/%s.info.xml',cfsWorkingDirectory,invFileName) );
                Tensors{loopa,loopb,loopphi} = Eh;
                readFromXML = readFromXML + toc(tReadFromXML);
            end
            tpostprocessing = tic;
            delete(meshFile);
            postprocessing = postprocessing + toc(tpostprocessing);
        end
        delete( sprintf('%s/inv_tensor_*',cfsWorkingDirectory) );
        delete( sprintf('%s/results_hdf5/inv_Tensor_*',cfsWorkingDirectory) );
%         showProgress(counter,complete);
    end
end

cfsTime = cfs;

% Save homogenized tensors in a file
twriteStats = tic;
if numel(presetsphi) > 1
    dim = '3D';
else
    dim = '2D';
end
filename = strcat('Interpolation/detailed_stats_',num2str(npresetsa),'_',dim);
fid = fopen(filename,'wt');
fprintf(fid,'%d\t%d\t%d\t%s\t%e\t%e\t%e\t%e\t%e\n',npresetsa,npresetsb,npresetsphi,'voigt',0,0,0,0,0);
for loopa = 1:npresetsa
    for loopb = 1:npresetsb
        for loopphi = 1:npresetsphi;
            fprintf(fid,'%f\t%f\t%f\t',presetsa(loopa),presetsb(loopb),presetsphi(loopphi));
            Eh = Tensors{loopa,loopb,loopphi};
            try
                fprintf(fid,'%e\t%e\t%e\t',Eh(1,1),Eh(1,2),Eh(1,3));
                fprintf(fid,'%e\t%e\t%e\n',Eh(2,2),Eh(2,3),Eh(3,3));
            catch ME
                disp(ME.message);
            end
        end
    end
end
fclose(fid);
writeStats = writeStats + toc(twriteStats);

preprocessing
cfs
postprocessing

readFromXML
writeStats

end



function showProgress(counter,complete)
progress = round(counter/complete*100);
clc;
display(sprintf('Progress: %d %%',progress));
end

%    FftHdf5Core
%
%      Time data filtering tool for CFS++ HDF5 data.
%
%      This script reads transient results from HDF5 files and performs a FFT
%      on them. Afterwards either the harmonic data or the filtered,
%      retransformed transient data is written to a (different) HDF5 file.
%
%      First the script tries to determine the size of the transient data in
%      megabytes: size_in_mb = num_items * numsteps * 8 / 1024 / 1024;
%      Then it determines how many iterations are needed to perform the FFT
%      with the given buffer size (bufsize) and how many items
%      (items_per_iter=nodes) can be processed in one iteration. Next the
%      script allocates a matrix of size numsteps*items_per_iter and performs
%      a FFT on it. Then entries outside of the given frequency range are set
%      to zero. In case of mode=1 the harmonic data are written to temporary
%      files in real-imag format. In case of mode=2 an inverse FFT is applied
%      and transient data are written to temporary files. At the end all
%      temporary files are combined into one file.
%    
% Input Parameters
%   * mode      - 1: save harmonic data; 2: save filtered transient data
%   * infile    - path of input HDF5 file
%   * outfile   - path of output HDF5 file
%   * quantity  - which quantity to convert
%   * region    - region the quantity is defined on
%   * lowfreq   - lowest frequency to be stored (0 for unlimited)
%   * highfreq  - highest frequency to be stored (0 for unlimited)
%   * bufsize   - maximum memory consumption (in megabytes)
%
% Return Value
%   None
%
% About
%   * Created:  Jan 2006
%   * Authors:  Max Escobar, Simon Triebenbacher, Jens Grabinger
%   * Revision: $Id$


function [] =  FftHdf5Core(mode, infile, outfile, quantity, region, lowfreq, highfreq, bufsize)


if (mode < 1) || (mode > 2)
  error('Unknown mode: %d', mode);
end

if exist(infile, 'file') ~= 2
  error('File not found: %s', infile);
end

% add path to HDF5Tools from Matlab Central
% http://www.mathworks.com/matlabcentral/fileexchange/17172-hdf5tools
thisfile = mfilename('fullpath');
[here, ~, ~] = fileparts(thisfile);
addpath([here '/hdf5tools'])

multistep = 1;
step = 1;
fileinfo = hdf5info(infile);
toplevel = fileinfo.GroupHierarchy;

% check that requested result is really present in input file
[found resgroup restype msgroup datafile] = FindPathHDF5(toplevel, multistep, step, quantity, region);
if found < 8
  errorstr = 'Cannot find requested dataset.';
  if found >= 3
    errorstr = sprintf('%s Search aborted at %s', errorstr, resgroup.Name);
  end
  error(errorstr);
  return
end

% store path of multistep and of result dataset
basepath = msgroup.Name;

% calculate no of time and frequency steps
numsteps = size(msgroup.Groups,2)-1;
numharmsteps = floor(numsteps / 2)+1;

% store result type of quantity
switch restype
case 1 % nodes
  restypestr = 'Nodes';
case 4 % elements
  restypestr = 'Elements';
end

% make sure we don't overwrite results in output file
if exist(outfile, 'file') == 2
  outfileinfo = hdf5info(outfile);
  toplevelout = outfileinfo.GroupHierarchy;

  [found2, resgroup2, ~, ~, datafile2] = FindPathHDF5(toplevelout, multistep, step, quantity, region);

  if found2 == 8
    error('Quantity %s already present in file %s under path: %s.', quantity, datafile2.Filename, resgroup2.Name);
  end
end

% get delta_t of the transient simulation
time_step1 = h5attget(infile, [basepath '/Step_1'], 'StepValue');
time_step2 = h5attget(infile, [basepath '/Step_2'], 'StepValue');
dt = time_step2 - time_step1;

% adapt no. of frequency steps according to maximum frequency
if highfreq > 0
  harmstepbnd = floor(highfreq*numsteps*dt)+1;
  if (harmstepbnd > 0) && (harmstepbnd < numharmsteps)
    numharmsteps = harmstepbnd;
  end
end

% adapt no. of frequency steps according to minimum frequency
firstharmstep = 0;
if lowfreq > 0
  firstharmstep = ceil(lowfreq*numsteps*dt);
  if (firstharmstep > 0)
    if (firstharmstep < numharmsteps)
      numharmsteps = numharmsteps - firstharmstep;
    else
      firstharmstep = 0;
      warning('Lower frequency boundary out of range. Ignoring it.'); %#ok<*WNTAG>
    end
  end
end

% display info on no. of steps
fprintf('No. of time steps:      %d\n', numsteps)
if mode == 1
  fprintf('No. of frequency steps: %d\n', numharmsteps)
end
disp(' ')

% read first time step
dataset = sprintf('%s/Real', resgroup.Name);
ds = h5varget(datafile.Filename, dataset);

% Number of scalars in dataset
num_items = length(ds);
clear ds

% Size in megabytes of whole transient dataset
size_in_mb = num_items*numsteps*8 / 1024 / 1024;

% make sure that buffer size makes sense
if bufsize < 1
  warning('Invalid buffer size. Reset to 256 M.');
  bufsize = 256;
end

% Number of iterations required to perform the FFT on the whole dataset
numiter = ceil(size_in_mb / bufsize * 6); % scale by 6, because we need
                                          % 6*bufsize for fft

% Delete temp files from last run.
[~, name, ~] = fileparts(outfile);
tmpdir = tempname();
exec(sprintf('rm -rf %s', tmpdir));
exec(sprintf('mkdir %s', tmpdir));
                                         
% Number of scalars treated in one iteration (= chunk size)
items_per_iter = ceil(num_items / numiter);

% initialize chunk counters
if items_per_iter < num_items
  item_end = items_per_iter;
else
  item_end = num_items;
end

% initialize chunk counters while takeing into account
% a possible restart (right now hardcoded as its only prepared)

%flag if we really wanna do this...
resumefile =0; 

if resumefile == 1;
  %specify the chunk number to START!! with
  %NOT the number of chunks to skip...
  startChunkNum = 3;
  itercount = startChunkNum;
  item_start = ((startChunkNum-1)*items_per_iter) + 1;
  item_end = (startChunkNum)*items_per_iter;
  if item_end > num_items
    item_end = num_items;
  end 
else
  item_start = 1;
  itercount = 1;
  exec(sprintf('rm -rf %s', tmpdir));
  exec(sprintf('mkdir %s', tmpdir));
end



% create buffer with chunk size
mat = zeros(numsteps, items_per_iter);

% read in transient data divided into chunks
for iter=itercount:numiter

  for i=1:numsteps
  
    fprintf('Reading step %d of %d, chunk %d of %d\n', i, numsteps, iter, numiter)

    [found, resgroup, restype, ~, datafile] = FindPathHDF5(toplevel, multistep, i, quantity, region);
    if found < 8
      errorstr = sprintf('Cannot find dataset of time step %d.', i);
      if found >= 3
        errorstr = sprintf('%s Search aborted at %s', errorstr, resgroup.Name);
      end
      error(errorstr);
      return
    end

    dataset = sprintf('%s/Real', resgroup.Name);
    ds = h5varget(datafile.Filename, dataset);

    mat(i,1:(item_end-item_start+1)) = ds(item_start:item_end);

    clear ds
  end

  % move counters to next chunk
  item_start = item_end + 1;
  item_end = item_end + items_per_iter;
  if item_end > num_items
    item_end = num_items;
  end

  % perform FFT
  fprintf('\nPerforming FFT\n')
  MAT = fft(mat);
  
  % mode 1: write out harmonic data
  outfile_iter = sprintf('%s/%s_%d.h5', tmpdir, name, iter);
  if mode == 1
    % split result into real and imaginary part
    real_MAT = 2.*real(MAT)./numsteps;
    imag_MAT = 2.*imag(MAT)./numsteps;
    clear MAT

    % write back harmonic datasets, one file for each chunk
    if numiter > 1
      fprintf('Buffering result of chunk %d\n', iter)
    end

    for i=1:numharmsteps

      outpath = sprintf('%s/Step_%d/%s/%s/%s', basepath, i, quantity, region, restypestr);
      dataset = sprintf('%s/Real', outpath);

      ds = real_MAT(i+firstharmstep,:);
      if exist(outfile_iter, 'file') ~= 2
        hdf5write(outfile_iter, dataset, ds, 'WriteMode', 'overwrite');
      else
        hdf5write(outfile_iter, dataset, ds, 'WriteMode', 'append');
      end

      dataset = sprintf('%s/Imag', outpath);

      ds = imag_MAT(i+firstharmstep,:);
      hdf5write(outfile_iter, dataset, ds, 'WriteMode', 'append');

    end

    clear real_MAT imag_MAT ds

  else % mode 2: filter out some frequencies and transfrom back
    if firstharmstep > 0
      MAT(1:firstharmstep,:) = 0;
    end
    if numharmsteps < size(MAT,1) - firstharmstep
      MAT(firstharmstep+numharmsteps+1:size(MAT,1),:) = 0;
    end

    fprintf('Performing inverse FFT\n\n')

    mat = real(ifft(MAT));

    clear MAT

    % write back transient datasets, one file for each chunk
    if numiter > 1
      fprintf('Buffering result of chunk %d\n', iter)
    end

    for i=1:numsteps
      outpath = sprintf('%s/Step_%d/%s/%s/%s/Real', basepath, i, quantity, region, restypestr);

      if exist(outfile_iter, 'file') ~= 2
        hdf5write(outfile_iter, outpath, mat(i,:), 'WriteMode', 'overwrite');
      else
        hdf5write(outfile_iter, outpath, mat(i,:), 'WriteMode', 'append');
      end
    end

    clear mat
  end % if mode 1/2
end

% create buffer for whole dataset of real part
ds_real = zeros(1, num_items);

% number of temporary files for mode 2
numfiles = numsteps;

if mode == 1
  % create buffer for whole dataset of imaginary part
  ds_imag = zeros(1, num_items);

  % number of temporary files for mode 1
  numfiles = numharmsteps;

  % calculate frequency steps
  f = double((firstharmstep:numharmsteps-1+firstharmstep)/(numsteps*dt));
else % mode == 2
  % calculate time steps
  t = double((1:numsteps)*dt);
end

% Write back harmonic datasets
for i=1:numfiles

  steppath = sprintf('%s/Step_%d', basepath, i);
  outpath = sprintf('%s/%s/%s/%s', steppath, quantity, region, restypestr);

  % store chunks into one dataset
  for iter=1:numiter
    outfile_iter = sprintf('%s/%s_%d.h5', tmpdir, name, iter);

    idx = (iter-1)*items_per_iter;
    idxend = idx+items_per_iter;
  
    dataset = sprintf('%s/Real', outpath);
    ds = h5varget(outfile_iter, dataset);
    ds_real(idx+1:idxend) = ds;

    if mode == 1
      dataset = sprintf('%s/Imag', outpath);
      ds = h5varget(outfile_iter, dataset);
      ds_imag(idx+1:idxend) = ds;
    end

  end

  clear ds

  fprintf('Writing step %d of %d\n', i, numfiles)

  % write to final output file
  dataset = sprintf('%s/Real', outpath);
  if exist(outfile, 'file') ~= 2
    hdf5write(outfile, dataset, ds_real(1:num_items), 'WriteMode', 'overwrite');
  else
    hdf5write(outfile, dataset, ds_real(1:num_items), 'WriteMode', 'append');
  end

  if mode == 1
    dataset = sprintf('%s/Imag', outpath);
    hdf5write(outfile, dataset, ds_imag(1:num_items), 'WriteMode', 'append');
  end

  % store current frequency in step attribute
  attr_details.AttachedTo = steppath;
  attr_details.AttachType = 'group';
  attr_details.Name = 'StepValue';
  if mode == 1
    step_value = f(i);
  else
    step_value = t(i);
  end
  try
    fattr=h5attget(outfile, steppath, 'StepValue');
    if fattr ~= step_value
      warning('Attribute %s/StepValue already exists in output file and has a different value', ...
              steppath);
    end
  catch %#ok<CTCH>
    h5attput(outfile, steppath, 'StepValue', step_value);
  end

end

clear ds_real ds_imag

fprintf('\nFinalizing output file\n')

% set required attributes of results
try
  oldlaststep = h5attget(outfile, basepath, 'LastStepNum');
catch %#ok<CTCH>
  oldlaststep = 0;
end
if i > oldlaststep
  h5attput(outfile, basepath, 'LastStepNum', uint32(i));
if mode == 1
    h5attput(outfile, basepath, 'LastStepValue', f(i));
else
    h5attput(outfile, basepath, 'LastStepValue', t(i));
  end
end
if mode == 1
  analtype = 'harmonic';
else
  analtype = 'transient';
end
try
  oldatype = h5attget(outfile, basepath, 'AnalysisType');
catch %#ok<CTCH>
  oldatype = {analtype};
end
if oldatype{1} ~= analtype
  error('Output file has "%s" set as analysis type, but you want to write %s data', ...
      oldatype{1}, analtype)
end
h5WriteVLStrAtt(outfile, basepath, 'AnalysisType', analtype, true);
try
  oldextfiles = h5attget(outfile, '/Results/Mesh', 'ExternalFiles');
catch %#ok<CTCH>
  oldextfiles = 0;
end
if oldextfiles ~= 0
  error('Output file uses external step files, but this is not supported here')
end
h5attput(outfile, '/Results/Mesh', 'ExternalFiles', uint32(0));

% write ResultDescription of quantity
rdpath = [basepath '/ResultDescription/'  quantity];
% create first dataset explicitly, so that groups get created as well
try
  h5datacreate(outfile, [rdpath '/DefinedOn'], 'type', 'uint32', 'size', 1);
  h5varput(outfile, [rdpath '/DefinedOn'], uint32(restype));
catch %#ok<CTCH>
end
try
  h5datacreate(outfile, [rdpath '/EntryType'], 'type', 'uint32', 'size', 1);
  h5varput(outfile, [rdpath '/EntryType'], uint32(1));
catch %#ok<CTCH>
end
try
  h5datacreate(outfile, [rdpath '/NumDOFs'], 'type', 'uint32', 'size', 1);
  h5varput(outfile, [rdpath '/NumDOFs'], uint32(1));
catch %#ok<CTCH>
end
try
  h5datacreate(outfile, [rdpath '/StepNumbers'], 'type', 'uint32', 'size', numfiles);
  h5varput(outfile, [rdpath '/StepNumbers'], uint32(1:numfiles));
catch %#ok<CTCH>
end
try
  h5datacreate(outfile, [rdpath '/StepValues'], 'type', 'double', 'size', numfiles);
if mode == 1
    h5varput(outfile, [rdpath '/StepValues'], f);
else
    h5varput(outfile, [rdpath '/StepValues'], t);
  end
catch %#ok<CTCH>
end

h5WriteVLStrDset(outfile, rdpath, 'DOFNames', '-', true);

try
  oldregnames = h5varget(outfile, [rdpath '/EntityNames']);
  regnames = [oldregnames; {region}];
catch %#ok<CTCH>
  regnames = {region};
end
h5WriteVLStrDset(outfile, rdpath, 'EntityNames', regnames, true);

h5WriteVLStrDset(outfile, rdpath, 'Unit', 'unknown', true);

% copy mesh to output file
try
  h5attget(outfile, '/Mesh', 'Dimension');
catch %#ok<CTCH>
  h5copy(infile, outfile, '/Mesh', '/Mesh');
end

% write FileInfo group to outfile
WriteFileInfoHDF5(outfile, [1 2]);

% Delete temp files
exec(sprintf('rm -rf %s', tmpdir));

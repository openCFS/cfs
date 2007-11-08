% ==============================================================
%
%    Time2FreqHDF5
%
%
%    GENERAL
%    Time data filtering tool for CFS++ HDF5 data.
%
%    This script reads transient results from our HDF5 files
%    and performs a FFT on them. Afterwards the harmonic data
%    is written either to a (different) HDF5 file.
%    First the script tries to determine the size of the
%    transient data in megabytes:
%    size_in_mb = num_items * numsteps * 8 / 1024 / 1024;
%    Then it determines how many iterations are needed to
%    perform the FFT with the given buffer size (bufsize)
%    and how many items (items_per_iter=nodes) can be processed
%    in one iteration. Next the script allocates a matrix of
%    size numsteps*items_per_iter performs a FFT on it and
%    writes the harmonic data to temporary files in amplitude-
%    phase format.
%    At the end all temporary files are combined into one file.
%    
%    PATH must include h5tool for this script to work correctly.
%
%
%    INPUT/S
%      infile    - path of input HDF5 file
%      outfile   - path of output HDF5 file
%      quantity  - which quantity to convert
%      region    - region the quantity is defined on
%      bufsize   - buffer size for transient data in megabytes
%      dt        - time step size in seconds (not used)
%        
%    OUTPUT/S
%
%      
%    ABOUT
%
%      -Created:     Jan 2006
%      -Last update: 07 Nov 2007
%      -Revision:    0.5
%      -Authors:     Max Escobar, Simon Triebenbacher, Jens Grabinger
%
% =======================================================================


function [amplitude] =  Time2FreqHDF5(infile, outfile, quantity, region, bufsize, dt)


multistep = 1
step = 1
fileinfo = hdf5info(infile);
toplevel = fileinfo.GroupHierarchy;

% check that requested result is really present in input file
[found resgroup restype msgroup datafile] =  FindPathHDF5(toplevel, multistep, step, quantity, region);
if found < 8
  errorstr = 'Cannot find requested dataset.';
  if found >= 3
    errorstr = sprintf('%s Search aborted at %s', errorstr, resgroup.Name);
  end
  error(errorstr);
  return;
end

% store path of multistep and of result dataset
basepath = msgroup.Name;
respath = resgroup.Name;

% calculate no of time and frequency steps
numsteps = size(msgroup.Groups,2)-1
numharmsteps = floor(numsteps / 2)

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

  [found2 resgroup2 restype2 msgroup2] = FindPathHDF5(toplevelout, multistep, step, quantity, region);

  if found2 == 8
    errorstr = sprintf('Quantity %s already present in file %s under path: %s.', quantity, outfile, respath);
    error(errorstr);
    return;
  end
end

% get delta_t of the transient simulation
time_step1 = hdf5read(infile, strcat(basepath, '/Step_1/StepValue'));
time_step2 = hdf5read(infile, strcat(basepath, '/Step_2/StepValue'));
dt = time_step2 - time_step1;

% read first time step
dataset = sprintf('%s/Real', respath);
ds = hdf5read(datafile.Filename, dataset);

% Number of scalars in dataset
num_items = length(ds);

% Size in megabytes of whole transient dataset
size_in_mb = num_items*numsteps*8 / 1024 / 1024;

% Number of iterations required to perform the FFT on the whole dataset
numiter = ceil(size_in_mb / bufsize);

% Number of scalars treated in one iteration (= chunk size)
items_per_iter = ceil(num_items / numiter);

% initialize chunk counters
item_start = 1;
if items_per_iter < num_items
  item_end = items_per_iter;
else
  item_end = num_items;
end

% create buffer with chunk size
mat = zeros(numsteps, items_per_iter);

% read in transient data divided into chunks
for iter=1:numiter

  for i=1:numsteps
  
    dataset = sprintf('%s/Real', respath);
    ds = hdf5read(datafile.Filename, dataset);
  
    mat(i,1:item_end-item_start+1) = ds(item_start:item_end);

    disp(sprintf('reading step %d of %d, chunk %d', i, numsteps, iter))

  end

  % move counters to next chunk
  item_start = item_end + 1;
  item_end = item_end + items_per_iter;
  if item_end > num_items
    item_end = num_items;
  end

  % perform FFT
  MAT = fft(mat);

  % split result into real and imaginary part
  real_MAT = real(MAT);
  imag_MAT = imag(MAT);

  % write back harmonic datasets, one file for each chunk
  outfile_iter = sprintf('%s_%d', outfile, iter);

  for i=1:numharmsteps
  
    outpath = sprintf('%s/Step_%d/%s/%s/%s', basepath, i, quantity, region, restypestr);
    dataset = sprintf('%s/Real', outpath);

    ds = real_MAT(i,:);
    if exist(outfile_iter, 'file') ~= 2
      hdf5write(outfile_iter, dataset, ds, 'WriteMode', 'overwrite');
    else
      hdf5write(outfile_iter, dataset, ds, 'WriteMode', 'append');
    end

    dataset = sprintf('%s/Imag', outpath);

    ds = imag_MAT(i,:);
    hdf5write(outfile_iter, dataset, ds, 'WriteMode', 'append');

  end

end

% create buffers for whole datasets
ds_real = zeros(1, num_items);
ds_imag = zeros(1, num_items);

% calculate frequency steps
f = double((0:numharmsteps-1)/(numsteps*dt));

% Write back harmonic datasets
for i=1:numharmsteps

  steppath = sprintf('%s/Step_%d', basepath, i);
  inoutpath = sprintf('%s/%s/%s/%s', steppath, quantity, region, restypestr);

  % store chunks into one dataset
  for iter=1:numiter

    outfile_iter = sprintf('%s_%d', outfile, iter);

    idx = (iter-1)*items_per_iter;
    idxend = idx+items_per_iter;
  
    dataset = sprintf('%s/Real', inoutpath);
    ds = hdf5read(outfile_iter, dataset);
    ds_real(idx+1:idxend) = ds;

    dataset = sprintf('%s/Imag', inoutpath);
    ds = hdf5read(outfile_iter, dataset);
    ds_imag(idx+1:idxend) = ds;

  end

  % write to final output file
  dataset = sprintf('%s/Real', inoutpath);
  if exist(outfile, 'file') ~= 2
    hdf5write(outfile, dataset, ds_real(1:num_items), 'WriteMode', 'overwrite');
  else
    hdf5write(outfile, dataset, ds_real(1:num_items), 'WriteMode', 'append');
  end

  dataset = sprintf('%s/Imag', inoutpath);
  hdf5write(outfile, dataset, ds_imag(1:num_items), 'WriteMode', 'append');

  % store current frequency in step attribute
  attr_details.AttachedTo = steppath;
  attr_details.AttachType = 'group';
  attr_details.Name = 'StepValue';
  hdf5write(outfile, attr_details, f(i), 'WriteMode', 'append');

  disp(sprintf('Writing step %d of %d', i, numharmsteps))

end

% set required attributes of results
attr_details.AttachedTo = basepath;
attr_details.AttachType = 'group';
attr_details.Name = 'LastStepNum';
hdf5write(outfile, attr_details, uint32(i), 'WriteMode', 'append');
attr_details.Name = 'LastStepValue';
hdf5write(outfile, attr_details, f(i), 'WriteMode', 'append');

% construct cfg filename for h5tool
tmpfile = strcat(infile, '.h5cfg');

% use h5tool to write analysis type
fid = fopen(tmpfile, 'w');
attr_cfg = fprintf(fid, 'attribute\ncreate\n%s\n%s\nstring\nAnalysisType\nappend\nharmonic\n', outfile, basepath);
fclose(fid);
exec(sprintf('h5tool < %s > /dev/null', tmpfile));

% declare not to use external files
attr_details.AttachedTo = '/Results/Mesh';
attr_details.AttachType = 'group';
attr_details.Name = 'ExternalFiles';
hdf5write(outfile, attr_details, uint32(0), 'WriteMode', 'append');

% write ResultDescription of quantity
rdpath = sprintf('%s/ResultDescription/%s', basepath, quantity);
ds_name = sprintf('%s/%s', rdpath, 'DefinedOn');
hdf5write(outfile, ds_name, uint32(restype), 'WriteMode', 'append');
ds_name = sprintf('%s/%s', rdpath, 'EntryType');
hdf5write(outfile, ds_name, uint32(1), 'WriteMode', 'append');
ds_name = sprintf('%s/%s', rdpath, 'NumDOFs');
hdf5write(outfile, ds_name, uint32(1), 'WriteMode', 'append');
ds_name = sprintf('%s/%s', rdpath, 'StepNumbers');
hdf5write(outfile, ds_name, uint32(1:numharmsteps), 'WriteMode', 'append');
ds_name = sprintf('%s/%s', rdpath, 'StepValues');
hdf5write(outfile, ds_name, f, 'WriteMode', 'append');

% use h5tool to write string datasets in ResultDescription
fid = fopen(tmpfile, 'w');
attr_cfg = fprintf(fid, 'dataset\ncreate\n%s\n%s/ResultDescription/%s\nstring\nDOFNames\nappend\n-\n', outfile, basepath, quantity);
fclose(fid);
if exec(sprintf('h5tool < %s > /dev/null', tmpfile))
  warning('h5tool failed, result file may be incomplete.');
end

fid = fopen(tmpfile, 'w');
attr_cfg = fprintf(fid, 'dataset\ncreate\n%s\n%s/ResultDescription/%s\nstring\nEntityNames\nappend\n%s\n', outfile, basepath, quantity, region);
fclose(fid);
if exec(sprintf('h5tool < %s > /dev/null', tmpfile))
  warning('h5tool failed, result file may be incomplete.');
end

fid = fopen(tmpfile, 'w');
attr_cfg = fprintf(fid, 'dataset\ncreate\n%s\n%s/ResultDescription/%s\nstring\nUnit\nappend\nunknown\n', outfile, basepath, quantity);
fclose(fid);
if exec(sprintf('h5tool < %s > /dev/null', tmpfile))
  warning('h5tool failed, result file may be incomplete.');
end

% copy mesh to output file
CopyTreeHDF5(infile, '/Mesh', outfile, '/Mesh');

% write FileInfo group to outfile
WriteFileInfoHDF5(outfile, [1 2]);

% Delete temp files
for iter=1:numiter
  outfile_iter = sprintf('%s_%d', outfile, iter);
  delete(outfile_iter);
end
delete(tmpfile);

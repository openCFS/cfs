% ==============================================================
%
%    GENERAL
%    Time data filtering tool for CFS++ HDF5 data.
%
%    This script reads transient results from our HDF5 files
%    and performs a FFT on them. Afterwards the harmonic data
%    is written either to same HDF5 file (outfile = infile) or
%    to a different file.
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
%      -Last update: Jun 2007
%      -Revision:    0.4
%      -Author:      Max Escobar, Simon Triebenbacher
% ==============================================================


function [amplitude] =  Time2FreqHDF5(infile, outfile, quantity, region, bufsize, dt)



progress_bar_position = 0;

% estimate for first iteration
time_for_this_iteration = 0.01;
% quantity = 'acouRHSval';
multistep = 1
fileinfo = hdf5info(infile);
toplevel = fileinfo.GroupHierarchy;
step = 1

[found stepname msgroup] =  FindPathHDF5(toplevel, multistep, step, quantity, region);

basepath = msgroup.Name;
numsteps = size(msgroup.Groups,2)-1
numharmsteps = numsteps / 2

sprintf('%d %s %s %s', found, basepath, stepname, msgroup.Name)


dataset = sprintf('%s/Step_1/%s/%s', basepath, quantity, region);
ds = hdf5read(infile, dataset);
% Number of scalars in dataset
num_items = length(ds)

% Size in megabytes of whole transient dataset
size_in_mb = num_items*numsteps*8 / 1024 / 1024

% Buffersize for FFT
% bufsize = 1;

% Number of iterations required to perform the FFT on the whole dataset
numiter = ceil(size_in_mb / bufsize)

% Number of scalars treated in one iteration
items_per_iter = ceil(num_items / numiter)

item_start = 1
if items_per_iter < num_items
  item_end = items_per_iter
else
  item_end = num_items
end

mat = zeros(numsteps, items_per_iter);

if exist(outfile, 'file') == 2
  outfileinfo = hdf5info(outfile);
  toplevelout = outfileinfo.GroupHierarchy;

  qtity = sprintf('%s_abs', quantity);
  [found stepname msgroup] =  FindPathHDF5(toplevelout, multistep, step, qtity, region);

  if found > 5
    errorstr = sprintf('Quantity %s already present in file %s under path: %s/%s.', qtity, outfile, basepath, stepname);
    error(errorstr);
  end
end

for iter=1:numiter

  for i=1:numsteps
    tic;
    % loading data
  
    dataset = sprintf('%s/Step_%d/%s/%s', basepath, i, quantity, region);

    echo off
    ds = hdf5read(infile, dataset);
  
    mat(i,1:item_end-item_start+1) = ds(item_start:item_end);

    %  percent = i / steps_per_percent;
    %  if check == 0
    %    percent = percent + 1
    %  end
    sprintf('In iter %d, reading step %d of %d', iter, i, numsteps)

  %  attr = 'Some random data';
  %  dataset = sprintf('/Data/Volume/Multistep_1/Step_%d',i);
  %  attr_details.AttachedTo = dataset;
  %  attr_details.AttachType = 'group';

  end

  item_start = item_end + 1
  item_end = item_end + items_per_iter
  if item_end > num_items
    item_end = num_items
  end

  MAT = fft(mat);

  abs_MAT = 2*abs(MAT) / numsteps;
  phase_MAT = unwrap(angle(MAT))*180/(pi);

  % Write back harmonic datasets
  outfile_iter = sprintf('%s_%d', outfile, iter);
  for i=1:numharmsteps
    tic;
  
    dataset = sprintf('%s/Step_%d/%s_abs/%s', basepath, i, quantity, region);

    echo off
    ds = abs_MAT(i,:);
    if exist(outfile_iter, 'file') == 0
      hdf5write(outfile_iter, dataset, ds, 'WriteMode', 'overwrite');
    else
      hdf5write(outfile_iter, dataset, ds, 'WriteMode', 'append');
    end

    dataset = sprintf('%s/Step_%d/%s_phase/%s', basepath, i, quantity, region);

    echo off
    ds = phase_MAT(i,:);
    hdf5write(outfile_iter, dataset, ds, 'WriteMode', 'append');
  
    sprintf('Writing step %d of %d', i, numharmsteps)

  end

end

ds_abs = zeros(1, num_items);
ds_phase = zeros(1, num_items);

% Write back harmonic datasets
for i=1:numharmsteps

  for iter=1:numiter

    outfile_iter = sprintf('%s_%d', outfile, iter);

    idx = (iter-1)*items_per_iter
    idxend = idx+items_per_iter
  
    dataset = sprintf('%s/Step_%d/%s_abs/%s', basepath, i, quantity, region);
    ds = hdf5read(outfile_iter, dataset);
    ds_abs(idx+1:idxend) = ds;

    dataset = sprintf('%s/Step_%d/%s_phase/%s', basepath, i, quantity, region);
    ds = hdf5read(outfile_iter, dataset);
    ds_phase(idx+1:idxend) = ds;
  end

  echo off
  dataset = sprintf('%s/Step_%d/%s_abs/%s', basepath, i, quantity, region);
  if exist(outfile, 'file') == 0
    hdf5write(outfile, dataset, ds_abs(1:num_items), 'WriteMode', 'overwrite');
  else
    hdf5write(outfile, dataset, ds_abs(1:num_items), 'WriteMode', 'append');
  end

  dataset = sprintf('%s/Step_%d/%s_phase/%s', basepath, i, quantity, region);

  hdf5write(outfile, dataset, ds_phase(1:num_items), 'WriteMode', 'append');
  
  sprintf('Writing step %d of %d', i, numsteps)

end

% Delete temp files
for iter=1:numiter
  outfile_iter = sprintf('%s_%d', outfile, iter);
  delete(outfile_iter);
end
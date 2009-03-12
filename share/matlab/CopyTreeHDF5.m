% ==============================================================
%
%    CopyTreeHDF5
%
%
%    GENERAL
%    Copies a given data tree from one HDF5 file to another.
%    Known Issues (due to restrictions of Matlab):
%    - Cannot copy string datasets
%    - Cannot copy attributes of datasets
%
%    INPUT/S
%      srcfile   - Top level node (/) of a HDF5 file
%      srcpath   - which multistep to search for
%      destfile  - which step to search for
%      destpath  - what quantity
%        
%    OUTPUT/S
%
%      
%    ABOUT
%
%      -Created:     30 Oct 2007
%      -Revision:    0.1
%      -Author:      Jens Grabinger
%
% ==============================================================

function [] = CopyTreeHDF5(srcfile, srcpath, destfile, destpath)

% source file must exist
if exist(srcfile, 'file') ~= 2
  errstr = sprintf('File not found: %s', srcfile);
  error(errstr);
  return;
end

% srcpath must not end in '/', but srcpath_ must
if strcmp(srcpath(length(srcpath)), '/')
  srcpath_ = srcpath;
  srcpath = srcpath(1:length(srcpath)-1);
else
  srcpath_ = strcat(srcpath, '/');
end
srcpath_len = length(srcpath_);

srcinfo = hdf5info(srcfile);
toplevel = srcinfo.GroupHierarchy;

% search for srcpath
i = 1;
cur_group = toplevel;
while ~strcmp(cur_group.Name, srcpath)
  if i > length(cur_group.Groups)
    break;
  end

  test_group = cur_group.Groups(i);
  cmpstr = strcat(test_group.Name, '/');
  
  len = min(length(cmpstr), srcpath_len);
  if strncmp(cmpstr, srcpath_, len)
    i = 1;
    cur_group = test_group;
  else
    i = i + 1;
  end
end

if ~strcmp(cur_group.Name, srcpath)
  errstr = sprintf('source path not found: %s', srcpath);
  error(errstr);
  return;
end

% copy tree recursively
CopyGroupHDF5(cur_group, length(srcpath), destfile, destpath);


function [] = CopyGroupHDF5(srcgroup, srcpathlen, destfile, destpath)

if ~strcmp(destpath(length(destpath)), '/')
  destpath = strcat(destpath, '/');
end

% first copy all datasets in current group
ndatasets = length(srcgroup.Datasets);
for dataset=1:ndatasets
  ds = srcgroup.Datasets(dataset);
  srcpath = ds.Name;
  ds = hdf5read(ds.Filename, srcpath);

  if isnumeric(ds)
    fullpath = strcat(destpath, srcpath(srcpathlen+2:length(srcpath)));

    hdf5write(destfile, fullpath, ds, 'WriteMode', 'append');
  end
end

% handle all subgroups
nsubgroups = length(srcgroup.Groups);
for subgroup=1:nsubgroups
  CopyGroupHDF5(srcgroup.Groups(subgroup), srcpathlen, destfile, destpath);
end

% copy attributes of current group
nattribs = length(srcgroup.Attributes);
for attrib=1:nattribs
  attr = srcgroup.Attributes(attrib);
  attr_name = attr.Name;
  srcpath = srcgroup.Name;

  fullpath = strcat(destpath, srcpath(srcpathlen+2:length(srcpath)));
  if strcmp(fullpath(length(fullpath)), '/')
    fullpath = fullpath(1:length(fullpath)-1);
  end
  details.AttachedTo = fullpath;
  details.AttachType = 'group';
  details.Name = attr_name(length(srcpath)+2:length(attr_name));

  hdf5write(destfile, details, attr.Value, 'WriteMode', 'append');
end

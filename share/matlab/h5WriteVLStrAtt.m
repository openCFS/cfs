% h5WriteVLStrAtt writes a string of variable length to an attribute in an
% HDF5 file.
%
% h5WriteVLStrAtt(FILENAME, PATH, NAME, VALUE, OVERWRITE) opens FILENAME,
% which may be a string containing the path to an HDF5 file or a file
% handle returned by the H5F.open function. Then it creates the attribute
% called NAME and attaches it to a group or dataset given by PATH. The
% attribute gets assigned VALUE, which may be either a string or a cell
% array of strings. If OVERWRITE is true, the atrribute will be overwritten
% if it already exists.

% $Id$


function [] = h5WriteVLStrAtt(filename, path, name, value, overwrite)

% make sure we have HDF5 1.8.0 at least
[maj, min, rel] = H5.get_libversion();
if (maj == 1) && (min < 8)
  error('Incompatible HDF5 library version (current: %d.%d.%d; required: >= 1.8.0)', ...
        maj, min, rel)
end

% check that arguments are valid
if ~ ischar(path)
  error('Attribute path must be a string parameter')
end
if ~ ischar(name)
  error('Attribute name must be a string parameter')
end
if iscellstr(value)
  dims = size(value);
  numdim = length(dims);
  attrval = value;
elseif ischar(value)
  dims = 1;
  numdim = 1;
  attrval = {value};
else
  error('Attribute value must be either a string or cell array of strings')
end
if isempty(overwrite)
  overwrite = false;
elseif ~ (islogical(overwrite) && isscalar(overwrite))
  error('Overwrite parameter must be a logical value')
end

if ischar(filename)
  % open HDF5 file
  if ~ H5F.is_hdf5(filename)
    error('"%s" is not a valid HDF5 file', filename)
  end
  h5file = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
elseif isa(filename, 'H5ML.id')
  % make sure it is actually is a file handle by re-opening
  h5file = H5F.reopen(filename);
  % known issue: writing doesn't work if the file hasn't been opened with r/w
  % access. One could provide an error message, if MATLAB had an interface to
  % the H5Fget_intent function.
else
  error('File parameter must be either a string or a H5F handle')
end

% open parent group or dataset
try
  % first try opening as a group
  parent = H5G.open(h5file, path);
catch %#ok<CTCH>
  % it isn't a group, so try opening as a datset
  slash = strfind(path, '/');
  if isempty(slash)
    H5F.close(h5file);
    error('%s is not a valid HDF5 path', path)
  end
  group_path = path(1:(slash(length(slash))-1));
  ds_path = path((slash(length(slash))+1):length(path));
  try
    group = H5G.open(h5file, group_path);
  catch %#ok<CTCH>
    H5F.close(h5file);
    error('Could not open group "%s"', group_path)
  end
  try
    parent = H5D.open(group, ds_path);
  catch %#ok<CTCH>
    H5F.close(h5file);
    error('Could not open dataset "%s"', path)
  end
end

attr_exists = true;
try
  attr = H5A.open_name(parent, name);
catch %#ok<CTCH>
  attr_exists = false;
end
if attr_exists
  H5A.close(attr);
  if overwrite
    H5A.delete(parent, name);
  else
    H5F.close(h5file);
    error('Attribute "%s" already exists', name)
  end
end

strtype = H5T.copy('H5T_C_S1');         % create a C-like string type
H5T.set_size(strtype, 'H5T_VARIABLE');  % set string size to "variable"
dspace = H5S.create_simple(numdim, dims, []); % create a dataspace

try
  attr = H5A.create(parent, name, strtype, dspace, 'H5P_DEFAULT');
  H5A.write(attr, strtype, attrval);
catch ex
  H5F.close(h5file);
  throw(ex)
end

H5A.close(attr);
H5S.close(dspace);
H5T.close(strtype);
try
  H5G.close(parent);
  H5D.close(parent);
catch %#ok<CTCH>
end
if exist('group', 'var')
  H5G.close(group);
end
H5F.close(h5file);

% h5WriteVLStrDset writes a string of variable length to a dataset in an
% HDF5 file.
%
% h5WriteVLStrDset(FILENAME, PATH, NAME, DATA, OVERWRITE) opens FILENAME,
% which may be a string containing the path to an HDF5 file or a file
% handle returned by the H5F.open function. Then it creates the dataset
% called NAME inside the group given by PATH. The contents of the dataset
% will be a 1x1 array of a variable-length string, if DATA is a string. If
% DATA is a cell array, the dataset will contain an array of equal
% dimensions as the cell array. If OVERWRITE is true, the dataset will be
% overwritten if it already exists.

% $Id$


function [] = h5WriteVLStrDset(filename, path, name, data, overwrite)

% make sure we have HDF5 1.8.0 at least
[maj, min, rel] = H5.get_libversion();
if (maj == 1) && (min < 8)
  error('Incompatible HDF5 library version (current: %d.%d.%d; required: >= 1.8.0)', ...
        maj, min, rel)
end

% check that arguments are valid
if ~ ischar(path)
  error('Dataset path must be a string parameter')
end
if ~ ischar(name)
  error('Dataset name must be a string parameter')
end
if iscellstr(data)
  dims = fliplr(size(data));
  numdim = length(dims);
  dsval = transpose(data);
elseif ischar(data)
  dims = 1;
  numdim = 1;
  dsval = {data};
else
  error('Data must be either a string or cell array of strings')
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

% open parent group
try
  parent = H5G.open(h5file, path);
catch %#ok<CTCH>
  H5F.close(h5file);
  error('Could not open group "%s"', path)
end

strtype = H5T.copy('H5T_C_S1');         % create a C-like string type
H5T.set_size(strtype, 'H5T_VARIABLE');  % set string size to "variable"
dspace = H5S.create_simple(numdim, dims, []); % create a dataspace

try
  ds = H5D.open(parent, name);
  ds_exists = true;
catch %#ok<CTCH>
  ds_exists = false;
end
if ds_exists
  if overwrite
    ds_delete = false;
    oldtype = H5D.get_type(ds);
    if H5T.equal(oldtype, strtype)
      oldspace = H5D.get_space(ds);
      if H5S.is_simple(oldspace)
        if H5S.get_simple_extent_ndims(oldspace) == numdim
          [~, olddims, ~] = H5S.get_simple_extent_dims(oldspace);
          if any(fliplr(olddims) ~= dims)
            ds_delete = true;
          end
        else
          ds_delete = true;
        end
      else
        ds_delete = true;
      end
      H5S.close(oldspace);
    else
      ds_delete = true;
    end
    H5T.close(oldtype);
    if ds_delete
      H5D.close(ds);
      H5L.delete(parent, name, []);
      ds_exists = false;
    end
  else
    H5F.close(h5file);
    error('Dataset "%s" already exists', name)
  end
end
if ~ds_exists
  dcpl = H5P.create('H5P_DATASET_CREATE');
  H5P.set_layout(dcpl, 'H5D_COMPACT');
  ds = H5D.create(parent, name, strtype, dspace, dcpl);
end

try
  H5D.write(ds, strtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', dsval);
catch ex
  H5F.close(h5file);
  throw(ex)
end

if exist('dcpl','var')
  H5P.close(dcpl);
end
H5D.close(ds);
H5S.close(dspace);
H5T.close(strtype);
H5G.close(parent);
H5F.close(h5file);

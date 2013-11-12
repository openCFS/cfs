function info = h5imfinfo(fileName,varname)
%H5IMFINFO returns information about HDF5 images.
%
%   INFO = H5IMFINFO(FILENAME) returns a structure whose fields
%   contain information about images in an HDF5 file.  FILENAME
%   is a string that specifies the name of the HDF5 file.
%
%   INFO = H5IMFINFO(FILENAME,VARNAME) returns a structure whose fields
%   contain information about only the specified image dataset.

%   Copyright 2007-2010 The MathWorks, Inc.


global imageVars;
imageVars = cell(0);

file_id = H5F.open (fileName, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

if nargin == 1

    % Iterate through all the datasets.
    H5G.iterate (file_id, '/',[] , @findHdf5Images);

else
    dset_id = H5D.open(file_id,    varname);
    if checkIfDatasetIsHdf5Image(dset_id)
        imageVars = { varname };
    end
    H5D.close(dset_id);
end

H5F.close (file_id);

if numel(imageVars) > 0
    info = getInfo(fileName,imageVars);
else
    error('HDF5TOOLS:h5imfinfo:notAnImageFile', ...
          'No images detected in %s.\n', fileName );
end

return


%--------------------------------------------------------------------------
function info = getInfo(hfile,vars)

hinfo = struct ( 'name', '', ...
               'class', '', ...
               'image_version', '', ...
               'image_subclass', '', ...
               'interlace_mode', '', ...
               'display_origin', [], ...
               'image_white_is_zero', [], ...
               'image_minmaxrange', [], ...
               'image_backgroundindex', [], ...
               'image_transparency', [], ...
               'image_aspectratio', [], ...
               'image_colormodel', [], ...
               'image_gammacorrection', [], ...
               'height', [], ...
               'width', [], ...
               'num_palettes', [] );
info = repmat(hinfo,numel(vars),1);
for j = 1:numel(vars)
    info(j) = getInfoSingleVar(hfile,vars{j});
end

return

%--------------------------------------------------------------------------
function hinfo = getInfoSingleVar(hfile,varname)

hinfo = struct ( 'name', '', ...
               'class', '', ...
               'image_version', '', ...
               'image_subclass', '', ...
               'interlace_mode', '', ...
               'display_origin', [], ...
               'image_white_is_zero', [], ...
               'image_minmaxrange', [], ...
               'image_backgroundindex', [], ...
               'image_transparency', [], ...
               'image_aspectratio', [], ...
               'image_colormodel', [], ...
               'image_gammacorrection', [], ...
               'height', [], ...
               'width', [], ...
               'num_palettes', [] );

hinfo.name = varname;

% Most of these attributes are optional
attributes = { 'CLASS', 'IMAGE_VERSION', ...
               'IMAGE_SUBCLASS', 'INTERLACE_MODE', 'DISPLAY_ORIGIN', ...
               'IMAGE_WHITE_IS_ZERO', 'IMAGE_MINMAXRANGE', ...
               'IMAGE_BACKGROUNDINDEX', 'IMAGE_TRANSPARENCY', ...
               'IMAGE_ASPECTRATIO', 'IMAGE_COLORMODEL', ...
               'IMAGE_GAMMACORRECTION' };
for j = 1:numel(attributes)
    try
        attval = h5attget(hfile,varname,attributes{j});

        % Some of these attributes come back with a null 
        % character tacked onto the end.
        if ischar(attval) && (attval(end) == char(0))
            attval = attval(1:end-1);
        end
        hinfo.(lower(attributes{j})) = attval;
    catch me
        switch ( version('-release') )
            case { '2008b', '2008a', '2007b', '2007a', '2006b' }
                % On these releases, if we get this error identifier,
                % then it just means that the attribute did not
                % exist.
                if ~strcmp(me.identifier,'MATLAB:H5ML_hdf5:invalidID')
                    rethrow(me);
                end
        end
    end
end

if ~isfield(hinfo,'class')
    error('MATLAB:h5imageinfo:missingClass', ...
          'The image %s is missing the required CLASS attribute.', ...
          varname );
end
if ~isfield(hinfo,'image_version')
    error('MATLAB:h5imageinfo:missingVersion', ...
          'The image %s is missing the required IMAGE_VERSION attribute.', ...
          varname );
end

hinfo = getWidthAndHeight(hfile,varname,hinfo);
if strcmp(hinfo.image_subclass,'IMAGE_INDEXED')
    hinfo.num_palettes = determineNumberOfPalettes(hfile,varname);
end

return

%---------------------------------------------------------------------------
function n = determineNumberOfPalettes(hfile,varname)
% How many pallettes are linked to the image?


file_id = H5F.open (hfile, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
dset_id = H5D.open(file_id,varname);
attr_id = H5A.open_name(dset_id,'PALETTE');

attr_space_id = H5A.get_space(attr_id);
dims = H5S.get_simple_extent_dims(attr_space_id);

if dims == 0
    n = 1;
else
    n = prod(dims);
end


H5S.close(attr_space_id);
H5A.close(attr_id);
H5D.close(dset_id);
H5F.close(file_id);
%---------------------------------------------------------------------------
function hinfo = getWidthAndHeight(hfile,varname,hinfo)

% Get the width and height.
flags = 'H5F_ACC_RDONLY';
plist_id = 'H5P_DEFAULT';
fileId = H5F.open ( hfile, flags, plist_id );
datasetId=H5D.open(fileId,varname);
spaceId = H5D.get_space(datasetId);
[numdims dims] = H5S.get_simple_extent_dims(spaceId); %#ok<ASGLU>
H5S.close(spaceId);
H5D.close(datasetId);
H5F.close(fileId);

switch hinfo.image_subclass
    case { 'IMAGE_INDEXED', 'IMAGE_GRAYSCALE' }
        hinfo.height = dims(1);
        hinfo.width = dims(2);

    case 'IMAGE_TRUECOLOR'
        if isfield(hinfo,'interlace_mode') && strcmp(hinfo.interlace_mode,'INTERLACE_PIXEL')
            % component values are contiguous
            hinfo.height = dims(1);
            hinfo.width = dims(2);
        elseif isfield(hinfo,'interlace_mode') ...
          && strcmp(hinfo.interlace_mode,'INTERLACE_PLANE')
            % component values are stored as a plane
            hinfo.height = dims(2);
            hinfo.width = dims(3);
        else
            error ( 'true color, but interlace mode ''%s''?', hinfo.interlace_mode);
        end

    otherwise
        error ( 'unhandled image class ''%s''', hinfo.image_subclass );
end
return

%--------------------------------------------------------------------------
function tf = correctImageDatatype(dsetId)

tf = false;

% The datatype must be floating point or integer.
dtype = H5D.get_type(dsetId);

if H5T.equal(dtype,'H5T_NATIVE_UCHAR')
    tf = true;
elseif H5T.equal(dtype,'H5T_NATIVE_USHORT')
    tf = true;
elseif H5T.equal(dtype,'H5T_NATIVE_ULONG')
    tf = true;
elseif H5T.equal(dtype,'H5T_NATIVE_FLOAT')
    tf = true;
elseif H5T.equal(dtype,'H5T_NATIVE_DOUBLE')
    tf = true;
end

H5T.close(dtype);

return

%--------------------------------------------------------------------------
function tf = correctAttributeClass(dsetId)
% The dataset for an image is distinguished from other datasets by
% giving it an attribute "CLASS=IMAGE".  In addition, the Image dataset
% may have an optional attribute "PALETTE" that is an array of object
% references for zero or more palettes.

tf = false;

attr_id = H5A.open_name(dsetId,'CLASS');
imageClass = H5A.read(attr_id,'H5ML_DEFAULT');
if ~ischar(imageClass)
    return
end

if ~( strcmp(imageClass','IMAGE') ...
      || strcmp(imageClass',['IMAGE' char(0)]))
    return
end

H5A.close(attr_id);

tf = true;

return

%--------------------------------------------------------------------------
function tf = checkIfDatasetIsHdf5Image(dset_id)

tf = correctImageDatatype(dset_id) ...
     && correctAttributeClass(dset_id);

return

%--------------------------------------------------------------------------
function tf = checkIfHdf5Image(loc_id,name)

dsetId = H5D.open(loc_id,name);

tf = checkIfDatasetIsHdf5Image(dsetId);

H5D.close(dsetId);

return

%--------------------------------------------------------------------------
function status=findHdf5Images (loc_id, name)

global imageVars;

statbuf=H5G.get_objinfo (loc_id, name, 0);

switch (statbuf.type)
    case H5ML.get_constant_value('H5G_GROUP')
        % Check this group out.
        H5G.iterate (loc_id, name,[] , @findHdf5Images);
        
        
    case H5ML.get_constant_value('H5G_DATASET')
        tf = checkIfHdf5Image(loc_id,name);
        if tf
            n = numel(imageVars);
            fullname = sprintf ( '%s/%s', H5I.get_name(loc_id), name );
            if (numel(fullname)>2) && ( strcmp(fullname(1:2),'//') ) 
                fullname = fullname(2:end);
            end
            imageVars{n+1} = fullname;
        end
        
    otherwise
        % Ignore anything else.
end

status=0;

return

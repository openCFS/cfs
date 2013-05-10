function [X,cmap] = h5imread(filename,varname,start,count,stride)
%H5IMREAD Read an image from an HDF5 file.
%
%   [X,MAP] = H5IMREAD(FILENAME) returns the first image it finds
%   in the HDF5 file specified by FILENAME.  X will be a 2-D uint8 
%   array if the specified data set contains an 8-bit image.  It 
%   will be an M-by-N-by-3 uint8 array if the specified data set 
%   contains a 24-bit image.  MAP contains the colormap if present; 
%   otherwise it is empty. 
%
%   [X,MAP] = H5IMREAD(FILENAME,DSETNAME) reads the specified
%   image from the HDF5 file specified by the string variable 
%   FILENAME.  
%
%   H5IMREAD assumes that the HDF5 file conforms to the HDF5 image
%   specification.  See http://www.hdfgroup.org/HDF5/Tutor/h5image.html
%   for details.
%
%   See also H5IMFINFO, IMREAD, IMWRITE, IMFINFO.

%   Copyright 2007-2010 The MathWorks, Inc.


switch(nargin)
	case 1
    	% restrict to the first image.
		hinfo = h5imfinfo(filename);
		hinfo = hinfo(1);
		varname = hinfo.name;
		hargs = { filename, varname };

	case 2
		hargs = { filename, varname };
		hinfo = h5imfinfo(filename,varname);

	case 4
		hargs = { filename, varname, start, count };
		hinfo = h5imfinfo(filename,varname);
	case 5
		hargs = { filename, varname, start, count, stride };
		hinfo = h5imfinfo(filename,varname);
end

X = h5varget(hargs{:});

X = permute(X, ndims(X):-1:1 );


if isfield(hinfo,'image_subclass') && strcmp(hinfo.image_subclass,'IMAGE_INDEXED')
	cmap = retrieveColormap(filename,varname);
end




%===============================================================================
function cmap = retrieveColormap(hfile,varname)

cmap = [];

flags = 'H5F_ACC_RDONLY';
plist_id = 'H5P_DEFAULT';
fid = H5F.open ( hfile, flags, plist_id );
datasetId=H5D.open(fid,varname);


try
	% If there is an IMAGE_SUBCLASS attribute witht the value 'IMAGE_INDEXED',
	% then we have an associated colormap.
    attId = H5A.open_name(datasetId,'IMAGE_SUBCLASS');
	attval = H5A.read(attId,'H5ML_DEFAULT');
	H5A.close(attId);

	% If a string, the attribute will be a column vector.  Compare it to
	% a regular row string.
	switch(attval')
		case ['IMAGE_INDEXED' char(0)]
			% YES!  Now retrieve the 'PALETTE' attribute.  That's a reference
			% that point to the associated colormap.
    		attId = H5A.open_name(datasetId,'PALETTE');
			ref = H5A.read(attId,'H5ML_DEFAULT');
			H5A.close(attId);

			% dereference the palette.  This had better point to a dataset.
			cmapId = H5R.dereference(datasetId,'H5R_OBJECT',ref);
			cmap = H5D.read(cmapId, 'H5ML_DEFAULT', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
			H5D.close(cmapId);

			% Transpose the colormap, normalize, convert to double precision.
			cmap = double(cmap')/256;

	end
catch me %#ok<NASGU>
	H5D.close(datasetId);
	H5F.close(fid);
	return
end

H5D.close(datasetId);
H5F.close(fid);








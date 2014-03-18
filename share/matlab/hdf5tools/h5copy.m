function h5copy(ifile,ofile,srcObjName,dstObjName)
%H5COPY recursively copies an object from one HDF5 file to another file.
%
%   h5copy(INPUTFILE,OUTPUTFILE,SRCOBJ,DSTOBJ) copies the source object 
%   name SRCOBJ in the input HDF5 file INPUTFILE to the destination 
%   identified by DSTOBJ in OUTPUTFILE.
%
%   The behavior of the copy is to recursively copy all objects below the 
%   group, keep soft links as they are, keep external links as they are, 
%   update only the values of object references, and to copy attributes 
%   along with the object.
%
%   This program was translated from the source code of "h5copy.c" written 
%   by Pedro Vicente Nunes and Quincey Koziol of the HDF Group

%   Copyright 2010 The MathWorks, Inc.

% open input file
ifid = H5F.open(ifile,'H5F_ACC_RDONLY','H5P_DEFAULT');

% open output file
if exist(ofile,'file')
	ofid = H5F.open(ofile,'H5F_ACC_RDWR','H5P_DEFAULT');
else
	fcpl = H5P.create('H5P_FILE_CREATE');
	ofid = H5F.create(ofile,'H5F_ACC_TRUNC',fcpl,'H5P_DEFAULT');
	H5P.close(fcpl);
end

% create property lists for copy
ocpl = H5P.create('H5P_OBJECT_COPY');

% create link creation property list
lcpl = H5P.create('H5P_LINK_CREATE');

% Set the intermediate creation property
H5P.set_create_intermediate_group(lcpl,true);

% do the copy
H5O.copy(ifid,srcObjName,ofid,dstObjName,ocpl,lcpl);

H5P.close(lcpl);
H5P.close(ocpl);
H5F.close(ifid);
H5F.close(ofid);

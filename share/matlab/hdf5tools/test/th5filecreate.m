function th5filecreate
fprintf('\tTesting H5FILECREATE...');
test_simple;
test_userblock;
test_userblock_not_power_of_two;
test_truncate;
test_truncate_not;
test_buffered;

fprintf('OK\n');


%--------------------------------------------------------------------------
function test_simple()
% test h5filecreate with no options
hfile = 'a.h5';
fid = h5filecreate(hfile);
H5F.close(fid);

% just being able to open the file is good enough
fid = H5F.open(hfile,'H5F_ACC_RDONLY', 'H5P_DEFAULT');
H5F.close(fid);


%--------------------------------------------------------------------------
function test_userblock()
% Simple usage of user block size
hfile = 'a.h5';
fid = h5filecreate(hfile,'userblock',512);
H5F.close(fid);

% 
fid = H5F.open(hfile,'H5F_ACC_RDONLY', 'H5P_DEFAULT');
fapl = H5F.get_create_plist(fid);
usize = H5P.get_userblock(fapl);
if ( usize ~= 512 )
    error('did not get the expected size of userblock');
end
H5F.close(fid);



%--------------------------------------------------------------------------
function test_userblock_not_power_of_two()
% Simple usage of user block size
hfile = 'a.h5';
fid = h5filecreate(hfile,'userblock',10191);
H5F.close(fid);

% 
fid = H5F.open(hfile,'H5F_ACC_RDONLY', 'H5P_DEFAULT');
fapl = H5F.get_create_plist(fid);
usize = H5P.get_userblock(fapl);
if ( usize ~= 16384 )
    error('did not get the expected size of userblock');
end
H5F.close(fid);



%--------------------------------------------------------------------------
function test_truncate()
% the truncate option should blow away any existing file
hfile = 'a.h5';
if exist(hfile,'file')
    delete(hfile);
end
fid = fopen(hfile,'w');
fwrite(fid,5);
fclose(fid);

fid = h5filecreate(hfile,'truncate',true);
H5F.close(fid);

% It's good enough to be able to open the file.
fid = H5F.open(hfile,'H5F_ACC_RDONLY', 'H5P_DEFAULT');
H5F.close(fid);


%--------------------------------------------------------------------------
function test_truncate_not()
% the truncate option should blow away any existing file
hfile = 'a.h5';
if exist(hfile,'file')
    delete(hfile);
end
fid = fopen(hfile,'w');
fwrite(fid,5);
fclose(fid);

try
    fid = h5filecreate(hfile,'truncate',false); %#ok<NASGU>
catch me
	switch(me.identifier)
		case { 'MATLAB:H5ML_hdf5:invalidID', ...
		       'MATLAB:hdf5lib2:H5Fcreate:fileCreationFailed', ...
			   'MATLAB:hdf5lib2:H5Fcreate:failure' }
			% 2007b, 
			% 2009a
			% 2009b and up
		otherwise
            rethrow(me)
    end
end

%--------------------------------------------------------------------------
function test_buffered()
% Simple usage of data buffering
hfile = 'a.h5';

% Can't specify buffer size without buffered.
try
    fid = h5filecreate(hfile,'buffer_size',1000); %#ok<NASGU>
catch me
    if ~strcmp(me.identifier,'HDF5TOOLS:bufferSizeWithoutBuffered')
        rethrow(me)
    end
end

fid = h5filecreate(hfile,'buffered',true,'buffer_size',100191);

fapl = H5F.get_access_plist(fid);
usize = H5P.get_sieve_buf_size(fapl);
if ( usize ~= 100191 )
    error('did not get the expected size of userblock');
end
H5F.close(fid);

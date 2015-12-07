function th5datacreate()
fprintf('\tTesting H5DATACREATE...  ');
%test_complex;
test_double;

fprintf('OK\n');

%--------------------------------------------------------------------------
function test_double()

hfile = 'a.h5';

fid = h5filecreate(hfile);
dsid = h5datacreate(fid,'/DS1','Size',[20 30]);
H5D.close(dsid);
H5F.close(fid);

fid = H5F.open(hfile,'H5F_ACC_RDONLY','H5P_DEFAULT');
dsetId = H5D.open(fid,'/DS1');
spaceId = H5D.get_space(dsetId);

[n,dims] = H5S.get_simple_extent_dims(spaceId); %#ok<ASGLU>

% make it vertical so that this can pass on 2007b.  What were they thinking?
dims = flipud(dims(:));
if ( dims(1) ~= 20 ) || ( dims(2) ~= 30 ) 
    error('Did not create a 20x30 double precision dataset');
end


%--------------------------------------------------------------------------
function test_complex()
% Verify that we create a complex dataset according to our convention.
hfile = 'a.h5';

fid = h5filecreate(hfile);
dsid = h5datacreate(fid,'DS1','Size',[20 30],'type','complex');
H5D.close(dsid);
H5F.close(fid);

fid = H5F.open(hfile,'H5F_ACC_RDONLY','H5P_DEFAULT');
dsetId = H5D.open(fid,'DS1');
spaceId = H5D.get_space(dsetId);

[n,dims] = H5S.get_simple_extent_dims(spaceId); %#ok<ASGLU>
dims = fliplr(dims);
if ( dims ~= 2 )
    error('The complex dataset does not consist of the required two references.');
end

% ok, read the references
dref = H5D.read(dsetId,'H5T_STD_REF_OBJ','H5S_ALL','H5S_ALL','H5P_DEFAULT');

% The first dereference should be to the real part.
dsetR = H5R.dereference(dsetId,'H5R_OBJECT',dref(:,1));
if ~strcmp(H5I.get_name(dsetR),'/_DS1_Real_')
    error('Did not dereference the real part correctly.');
end
dsetI = H5R.dereference(dsetId,'H5R_OBJECT',dref(:,2));
if ~strcmp(H5I.get_name(dsetI),'/_DS1_Imag_')
    error('Did not dereference the real part correctly.');
end

% check the comment
comment = H5O.get_comment(dsetR);
if ~strcmp(comment,'real part of complex dataset');
    error('failed');
end
comment = H5O.get_comment(dsetI);
if ~strcmp(comment,'imaginary part of complex dataset');
    error('failed');
end

spaceId = H5D.get_space(dsetR);
[n,dims] = H5S.get_simple_extent_dims(spaceId); %#ok<ASGLU>
dims = fliplr(dims);
if ( dims(1) ~= 20 ) || ( dims(2) ~= 30 ) 
    error('Did not create a 20x30 double precision dataset');
end

spaceId = H5D.get_space(dsetI);
[n,dims] = H5S.get_simple_extent_dims(spaceId); %#ok<ASGLU>
dims = fliplr(dims);
if ( dims(1) ~= 20 ) || ( dims(2) ~= 30 ) 
    error('Did not create a 20x30 double precision dataset');
end

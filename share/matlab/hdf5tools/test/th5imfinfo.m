function th5imfinfo()

%   Copyright 2008-2009 The MathWorks, Inc.
test_paletteImage;    
test_multipleImages;    
test_rgb;    

test_noImagesInFile;    

fprintf ( 1, 'All H5IMFINFO tests succeeded.\n\n' );

%--------------------------------------------------------------------------
function test_paletteImage()

actInfo = h5imfinfo('data/ex_image1.h5');

load('th5imfinfo.mat');
expInfo = expInfo.paletteImage;

if ~isequal(expInfo,actInfo)
    error('Palette info struct did not match');
end


fprintf ( 1, 'Retrieval of palette image info succeeded .\n' );

%--------------------------------------------------------------------------
function test_multipleImages()

actInfo = h5imfinfo('data/ex_image2.h5');

load('th5imfinfo.mat');
expInfo = expInfo.multipleImages;

if ~isequal(expInfo,actInfo)
    error('Multiple image info struct did not match');
end


fprintf ( 1, 'Retrieval of multiple image info structure succeeded .\n' );

%--------------------------------------------------------------------------
function test_rgb()

actInfo = h5imfinfo('data/ex_image2.h5','/image24bitpixel');

load('th5imfinfo.mat');
expInfo = expInfo.rgb;

if ~isequal(expInfo,actInfo)
    error('RGB image info struct did not match');
end


fprintf ( 1, 'Retrieval of RGB image info succeeded .\n' );

%--------------------------------------------------------------------------
function test_noImagesInFile()

switch ( version('-release') ) 
    case { '2007b', '2007a', '2006b', '2006a' }
		try
    		actInfo = h5imfinfo('data/h5ex_d_fillval.h5');
		catch 
			[a,b] = lasterr;
            if ~strcmp(b,'HDF5TOOLS:h5imfinfo:notAnImageFile')
                error(b,a);
            end
		end
	
	otherwise
		try
    		actInfo = h5imfinfo('data/h5ex_d_fillval.h5');
		catch me
            if ~strcmp(me.identifier,'HDF5TOOLS:h5imfinfo:notAnImageFile')
                rethrow(me)
            end
		end
	
end

fprintf ( 1, 'Negative test on non-image file succeeded.\n' );


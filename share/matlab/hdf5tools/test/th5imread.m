function th5imread()

%   Copyright 2008-2009 The MathWorks, Inc.
test_paletteImage;    
test_rgb;    

fprintf ( 1, 'All H5IMREAD tests succeeded.\n\n' );

%--------------------------------------------------------------------------
function test_paletteImage()

[x,map] = h5imread('data/ex_image1.h5');

image(x); colormap(map);
yn = input ('Does it look like a series of horizontal bars? (y/n)\n', 's' );
if strcmp(yn,'n')
	error('Stopping tests at %s', mfilename);
end
delete(gcf);


fprintf ( 1, 'Retrieval of palette image succeeded .\n' );

%--------------------------------------------------------------------------
function test_rgb()

x = h5imread('data/ex_image2.h5','/image24bitpixel');

image(x);
yn = input ('Does it look like a rose? (y/n)\n', 's' );
if strcmp(yn,'n')
	error('Stopping tests at %s', mfilename);
end
delete(gcf);

fprintf ( 1, 'Retrieval of RGB image succeeded .\n' );



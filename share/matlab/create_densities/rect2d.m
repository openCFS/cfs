function rect2d(sx, sy, ux, uy, filename)
# parameters:
# sx = steps in x direction
# sy = steps in y direction
# ux, uy = x and y coords of upper right corner
# assumes lower left corner = 0,0
# 
# we assume a lexicographic ordering
# 11 12 ...
# 6  7  8  9  10
# 1  2  3  4  5 

xml_file = fopen(filename, 'w');
fprintf(xml_file, '<?xml version="1.0"?>\n  <cfsErsatzMaterial>\n    <header>\n    <design constant="false" initial="0.5" lower="1e-3" name="density" region="mech" scale="false" upper="1"/>\n    <transferFunction application="mech" design="density" param="1" type="simp"/>\n  </header>\n\n');

###############################
fprintf(xml_file, '<set id="full">\n');
  for x = 1:sx*sy;
    # number of current element
    fprintf(xml_file, '  <element nr="%d" type="density" design="1.0"/>\n', x);
  endfor
fprintf(xml_file, '</set>\n\n');

#fprintf(xml_file, '</cfsErsatzMaterial>\n');
#fclose(xml_file);

#exit
###############################
# width of element
w = ux/sx;
# height of element
h = uy/sy;

# center coords
cenx = ux/2;
ceny = uy/2;

# maximal distance of any element to the center
# = distance of element 1 to the center
dmax = sqrt((cenx - 0.5*w)^2 + (ceny - 0.5*h)^2);
dmin = 0.001;

fprintf(xml_file, '<set id="center">\n');
'center'
for y = 1:sy;
    # center y of current element
    cy = (y - 1) * h + 0.5*h;
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

    # center x of current element
    cx = (x - 1) * w + 0.5*w;

    # distance of current center to center of rect
    dist = sqrt((cenx - cx)^2 + (ceny - cy)^2);

    v = dist/dmax;
    if v < dmin 
      v = dmin; 
    end
    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="center-inverted">\n');
'center-inverted'
for y = 1:sy;
    # center y of current element
    cy = (y - 1) * h + 0.5*h;
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

    # center x of current element
    cx = (x - 1) * w + 0.5*w;

    # distance of current center to center of rect
    dist = sqrt((cenx - cx)^2 + (ceny - cy)^2);

    v = 1 - dist/dmax;
    if v < dmin 
      v = dmin; 
    end
    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="vertical">\n');
'vertical'
for y = 1:sy;
  for x = 1:sx;
    # center x of current element
    v = (y - 1) * h + 0.5*w;
    if v < dmin 
      v = dmin;
		endif
    # number of current element
    num = (y - 1) * sx + x;

    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="diagonal">\n');
'diagonal'
for y = 1:sy;
    # center y of current element
    cy = (y - 1) * h + 0.5*h;
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

    # center x of current element
    cx = (x - 1) * w + 0.5*w;

    v = min(cx, cy);
    if v < dmin
      v = dmin;
		endif

    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="checkerboard">\n');
'checkerboard'
needoffset = 0;
if(mod(sx, 2) == 0)
	needoffset = 1;
endif
offset = 0;
for y = 1:sy;
	if(needoffset == 1)
		if(mod(y, 2) == 0)
			offset = 1;
		else
			offset = 0;
		endif
	endif
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

		if(x == 1 || x == sx || y == 1 || y == sy)
			v = 1;
		else
			if(mod(x + offset, 2))
				v = dmin;
			else 
				v = 1;
			endif
		endif

    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="vbars">\n');
'vbars'
for y = 1:sy;
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

		if(x == 1 || x == sx || y == 1 || y == sy)
			v = 1;
		else
			if(mod(x, 2))
				v = dmin;
			else 
				v = 1;
			endif
		endif

    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="hbars">\n');
'hbars'
for y = 1:sy;
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

		if(x == 1 || x == sx || y == 1 || y == sy)
			v = 1;
		else
			if(mod(y, 2))
				v = dmin;
			else 
				v = 1;
			endif
		endif

    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="random">\n');
'random'
for y = 1:sy;
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

		if(x == 1 || x == sx || y == 1 || y == sy)
			v = 1;
		else
			v = rand();
      # do not use dmin here, material might get too weak
			if(v < 0.1)
				v = dmin;
			endif
		endif

    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="random2">\n');
'random2'
for y = 1:sy;
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

		if(x == 1 || x == sx || y == 1 || y == sy)
			v = 1;
		else
			v = rand();
      # do not use dmin here, material might get too weak
			if(v < 0.1)
				v = dmin;
			endif
		endif

    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="circle">\n');
'circle'
rad = sqrt((1.0-0.38)/pi);
for y = 1:sy;
    # center y of current element
    cy = (y - 1) * h + 0.5*h;
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

    # center x of current element
    cx = (x - 1) * w + 0.5*w;

    # distance of current center to center of rect
    dist = sqrt((cenx - cx)^2 + (ceny - cy)^2);

    if dist > rad
      v = 1.0;
    else
      v = dmin; 
    end
    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="donut">\n');
'donut'
rad1 = 0.5;
rad2 = 0.4;
for y = 1:sy;
    # center y of current element
    cy = (y - 1) * h + 0.5*h;
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

    # center x of current element
    cx = (x - 1) * w + 0.5*w;

    # distance of current center to center of rect
    dist = sqrt((cenx - cx)^2 + (ceny - cy)^2);

    if dist > rad2 && dist < rad1
      v = 1.0;
    else
      v = dmin; 
    end
    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="frame">\n');
'frame'
# frame should be 10 percent of width
xframe = 0.1 * sx;
yframe = 0.1 * sy;

for y = 1:sy;
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

		if(x < xframe || x > sx-xframe || y < yframe || y > sy - yframe)
			v = 1;
		else
		  v = dmin;
		end

    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="orthox">\n');
'orthox'
for y = 1:sy;
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

    if(y > sy/2 - sy/10 && y < sy/2 + sy/10)
      v = 1;
    else 
      v = dmin;
		endif

    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="orthoy">\n');
'orthoy'
for y = 1:sy;
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

    if(x > sx/2 - sx/10 && x < sx/2 + sx/10)
      v = 1;
    else 
      v = dmin;
		endif

    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '<set id="cross">\n');
'cross'
for y = 1:sy;
  for x = 1:sx;
    # number of current element
    num = (y - 1) * sx + x;

    if((y > sy/2 - sy/10 && y < sy/2 + sy/10) || (x > sx/2 - sx/10 && x < sx/2 + sx/10) )
      v = 1;
    else 
      v = dmin;
		endif

    fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '</cfsErsatzMaterial>\n');
fclose(xml_file);

endfunction

function cube3d(sx, sy, sz, ux, uy, uz, vol, filename)
# parameters:
# sx = steps in x direction
# sy = steps in y direction
# sz = steps in z direction
# ux, uy, uz = x, y and z coords of upper right corner
# v = approximate toal volume
# assumes lower left corner = 0,0,0
# 
# we assume a lexicographic ordering
# 11 12 ...
# 6  7  8  9  10
# 1  2  3  4  5 

xml_file = fopen(filename, 'w');
fprintf(xml_file, '<?xml version="1.0"?>\n  <cfsErsatzMaterial>\n    <header>\n    <design constant="false" initial="0.5" lower="1e-3" name="density" region="mech" scale="false" upper="1"/>\n    <transferFunction application="mech" design="density" param="1" type="simp"/>\n  </header>\n\n');

# width of element
w = ux/sx;
# height of element
h = uy/sy;
# depth of element
d = uz/sz;

# center coords
cenx = ux/2;
ceny = uy/2;
cenz = uz/2;

# maximal distance of any element to the center
# = distance of element 1 to the center
dmax = sqrt((cenx - 0.5*w)^2 + (ceny - 0.5*h)^2 + (cenz - 0.5*d)^2);
dmin = 0.05;

fprintf(xml_file, '<set id="center">\n');
for z = 1:sz;
  # center z of current element
  cz = (z - 1) * d + 0.5*d;
  for y = 1:sy;
    # center y of current element
    cy = (y - 1) * h + 0.5*h;
    for x = 1:sx;
      # center x of current element
      cx = (x - 1) * w + 0.5*w;

      # number of current element
      num = (z - 1) * sx * sy + (y - 1) * sx + x;

      # distance of current center to center of rect
      dist = sqrt((cenx - cx)^2 + (ceny - cy)^2 + (cenz - cz)^2);

      v = dist/dmax;
      # assume 0.58 as normal volume
      v *= vol;
      if v < dmin 
        v = dmin; 
      end
      fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
    endfor
  endfor
endfor
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '</cfsErsatzMaterial>\n');
fclose(xml_file);

endfunction

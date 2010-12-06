function orthox(sx, sy, sz, ux, uy, uz, filename)
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
fprintf(xml_file, '<?xml version="1.0"?>\n  <cfsErsatzMaterial>\n    <header>\n    <mesh x="%d" y="%d" z="%d"/>\n    <design constant="false" initial="0.5" lower="1e-3" name="density" region="mech" scale="false" upper="1"/>\n    <transferFunction application="mech" design="density" param="1" type="simp"/>\n  </header>\n\n', sx, sy, sz);

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

fprintf(xml_file, '<set id="orthox">\n');
'orthox'
for z = 1:sz;
  for y = 1:sy;
    for x = 1:sx;
      v = 0.001;
      # number of current element
      num = (z - 1) * sx * sy + (y - 1) * sx + x;

      if (x == floor(sx/2))
        if (y == floor(sy/2)) || (z == floor(sz/2))
          v = 1.0;
        end
      end

      if y < (sy/2+3) && z < (sz/2+3)
        if y > (sy/2-3) && z > (sz/2-3)
          v = 1.0;
        end
      end
      fprintf(xml_file, '  <element nr="%d" type="density" design="%g"/>\n', num, v);
    end
  end
end
fprintf(xml_file, '</set>\n\n');

fprintf(xml_file, '</cfsErsatzMaterial>\n');
fclose(xml_file);

endfunction

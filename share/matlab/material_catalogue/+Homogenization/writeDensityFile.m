function writeDensityFile(density,densityfile,x,y)
% WRITEDENSITYFILE  -  Generates a density file out of a given density vector or matrix.
%
% @param:
%       density         density (entry has to be 0 for no material)
%       densityfile     name of generated .dens file
%

if nargin < 3
    [n,m] = size(density);
else
    n = x;
    m = y;
end

if islogical(density)
    density = double(density);
end

% Set density for weak material
density(density < 1e-10) = 1e-7;

% If density is column vector, transpose
if m == 1
    density = density';
else
    density = flipud(density);
end

% Write density file
fid = fopen(densityfile,'wt');

fprintf(fid,'<?xml version="1.0"?>\n');
fprintf(fid,'<cfsErsatzMaterial>\n');
fprintf(fid,'  <header>\n');
fprintf(fid,'    <mesh x="%d" y="%d" z="1"/>\n',n,m);
fprintf(fid,'    <design initial="0.5" lower="1e-3" name="density" region="mech" upper="1"/>\n');
fprintf(fid,'    <transferFunction application="mech" design="density" param="1" type="simp"/>\n');
fprintf(fid,'  </header>\n');
fprintf(fid,'  <set id="4">\n');
fprintf(fid,'    <element nr="%d" type="density" design="%e"/>\n',[1:n*m;density(:)']);
fprintf(fid,'  </set>\n');
fprintf(fid,'</cfsErsatzMaterial>\n');

fclose(fid);

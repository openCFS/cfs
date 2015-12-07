function writeDensityFile(density,densityfile)
% WRITEDENSITYFILE  -  Generates a density file out of a given density vector.
%
% @param:
%       density         density vector (entry has to be 0 for no material)
%       densityfile     name of generated .dens file
%

if islogical(density)
    density = double(density);
end

% Set density for weak material
density(density < 1e-10) = 1e-7;

% If density is column vector, transpose
if size(density,2) == 1
    density = density';
end

[n,m] = size(density);

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
fprintf(fid,'    <element nr="%d" type="density" design="%e"/>\n',[1:size(density,2);density]);
fprintf(fid,'  </set>\n');
fprintf(fid,'</cfsErsatzMaterial>\n');

fclose(fid);

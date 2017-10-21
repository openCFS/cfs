function matrixToDensity(A,densityfile)
% MATRIXTODENSITY  -  Generates a density file out of a given density matrix.
%
% @param:
%       A               density matrix (entry has to be 0 for no material)
%       densityfile     name of generated .dens file
%

% Set density for weak material
A(A < 1e-10) = 1e-7;

if ismatrix(A)

    % Rotate matrix to fit element mesh numbering (lower left to upper right, line by line)
    % to matlab linear indexing (uper left to lower right, column by column)
    A = flipud(A)';

    [m,n] = size(A);

    % Write density file
    fid = fopen(densityfile,'wt');

    fprintf(fid,'<?xml version="1.0"?>\n');
    fprintf(fid,'<cfsErsatzMaterial>\n');
    fprintf(fid,'  <header>\n');
    fprintf(fid,'    <mesh x="%d" y="%d" z="1"/>\n',m,n);
    fprintf(fid,'    <design initial="0.5" lower="1e-3" name="density" region="mech" upper="1"/>\n');
    fprintf(fid,'    <transferFunction application="mech" design="density" param="1" type="simp"/>\n');
    fprintf(fid,'  </header>\n');
    fprintf(fid,'  <set id="4">\n');

    for i = 1:m*n
        fprintf(fid,'    <element nr="%d" type="density" design="%e"/>\n',i,A(i));
    end

    fprintf(fid,'  </set>\n');
    fprintf(fid,'</cfsErsatzMaterial>\n');

    fclose(fid);

else
    
    [m,n,o] = size(A);

    % Write density file
    fid = fopen(densityfile,'wt');

    fprintf(fid,'<?xml version="1.0"?>\n');
    fprintf(fid,'<cfsErsatzMaterial>\n');
    fprintf(fid,'  <header>\n');
    fprintf(fid,'    <mesh x="%d" y="%d" z="%d"/>\n',m,n,o);
    fprintf(fid,'    <design initial="0.5" lower="1e-3" name="density" region="mech" upper="1"/>\n');
    fprintf(fid,'    <transferFunction application="mech" design="density" param="1" type="simp"/>\n');
    fprintf(fid,'  </header>\n');
    fprintf(fid,'  <set id="4">\n');

    for i = 1:m*n*o
        fprintf(fid,'    <element nr="%d" type="density" design="%e"/>\n',i,A(i));
    end

    fprintf(fid,'  </set>\n');
    fprintf(fid,'</cfsErsatzMaterial>\n');

    fclose(fid);
end
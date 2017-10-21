function array_to_densfile(arr)

[m,n] = size(arr);
filename = 'output.dens.xml';
fid = fopen(filename,'wt'); 
fprintf(fid,'<?xml version="1.0"?>\n<cfsErsatzMaterial>\n<header>\n<mesh x="%d" y="%d" z="1"/>',m,n);
fprintf(fid,'\n<design initial="0.5" lower="1e-9" name="density" region="mech" upper="1"/>');
fprintf(fid,'\n<transferFunction application="mech" design="density" param="1" type="simp"/>');
fprintf(fid,'\n</header>\n<set id="set">\n');
count = 1;
for i=1:m
    for j=1:n
        if arr(i,j)<1e-9
            arr(i,j) = 1e-9;
        end
        fprintf(fid,'<element nr="%d" type="density" design="%g"/>\n',count,arr(i,j));
        count = count +1;
    end
end
fprintf(fid,'</set>\n</cfsErsatzMaterial>');
fclose(fid);

end

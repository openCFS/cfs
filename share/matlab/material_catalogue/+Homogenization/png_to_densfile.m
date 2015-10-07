function png_to_densfile(name)

arr = imread(name);
% Rescale to values between 0 and 1
arr = 1/255 * arr;
% Flip black and white
arr = mod(1,arr);
[m,n] = size(arr);
filename = 'output.dens.xml';
fid = fopen(filename,'wt'); 
fprintf(fid,'<?xml version="1.0"?>\n <cfsErsatzMaterial>\n<header>\n<mesh x="%d" y="%d" z="1"/>',m,n);
fprintf(fid,'\n<design initial="0.5" lower="1e-6" name="density" region="mech" upper="1"/>');
fprintf(fid,'\n <transferFunction application="mech" design="density" param="1" type="simp"/>');
fprintf(fid,'\n</header>\n<set id="set">\n');
count = 1;
for i=1:m
    for j=1:n
        fprintf(fid,'<element nr="%d" type="density" design="%f"/>\n',count,arr(i,j));
        count = count +1;
    end
end
fprintf(fid,'</set>\n</cfsErsatzMaterial>');
fclose(fid);

end

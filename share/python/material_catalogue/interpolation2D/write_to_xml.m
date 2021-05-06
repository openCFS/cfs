function write_to_xml(outputfile, params, Coeff11, Coeff22, Coeff33, Coeff23, Coeff13, Coeff12, VolCoeff)
% Erzeugt .xml file mit namen file fuer CFS

fid = fopen(outputfile,'wt');

fprintf(fid, '<homRectC1 notation="voigt">\n');
for j=1:length(params)
    m = length(params{j});
    if m > 1
        fprintf(fid, '<param%d>\n\t<matrix dim1="%d" dim2="1">\n\t\t<real>\n', j, m);
        for i=1:m
            fprintf(fid, '%.16f ', params{j}(i));
        end
        fprintf(fid,'\n\t\t</real>\n\t</matrix>\n</param%d>\n',j);
    end
end

ndatapoints = prod(max(1,cellfun('length', params)-1));
ncoeffs = 4^length(params);

fprintf(fid,'<coeff11>\n\t<matrix dim1="%d" dim2="%d">\n\t\t<real>\n', ndatapoints, ncoeffs);
for i=1:ndatapoints
    fprintf(fid,'%.16f ', Coeff11(i,:));
    fprintf(fid,'\n');
end
fprintf(fid,'\t\t</real>\n\t</matrix>\n</coeff11>\n');

fprintf(fid,'<coeff22>\n\t<matrix dim1="%d" dim2="%d">\n\t\t<real>\n', ndatapoints, ncoeffs);
for i=1:ndatapoints
    fprintf(fid,'%.16f ', Coeff22(i,:));
    fprintf(fid,'\n');
end
fprintf(fid,'\t\t</real>\n\t</matrix>\n</coeff22>\n');

fprintf(fid,'<coeff33>\n\t<matrix dim1="%d" dim2="%d">\n\t\t<real>\n', ndatapoints, ncoeffs);
for i=1:ndatapoints
    fprintf(fid,'%.16f ', Coeff33(i,:));
    fprintf(fid,'\n');
end
fprintf(fid,'\t\t</real>\n\t</matrix>\n</coeff33>\n');

fprintf(fid,'<coeff23>\n\t<matrix dim1="%d" dim2="%d">\n\t\t<real>\n', ndatapoints, ncoeffs);
for i=1:ndatapoints
    fprintf(fid,'%.16f ', Coeff23(i,:));
    fprintf(fid,'\n');
end
fprintf(fid,'\t\t</real>\n\t</matrix>\n</coeff23>\n');

if numel(Coeff13) > 0
    fprintf(fid, '<coeff13>\n\t<matrix dim1="%d" dim2="%d">\n\t\t<real>\n', ndatapoints, ncoeffs);
    for i=1:ndatapoints
        fprintf(fid,'%.16f ', Coeff13(i,:));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\t\t</real>\n\t</matrix>\n</coeff13>\n');
end

if numel(Coeff12) > 0
    fprintf(fid,'<coeff12>\n\t<matrix dim1="%d" dim2="%d">\n\t\t<real>\n', ndatapoints, ncoeffs);
    for i=1:ndatapoints
        fprintf(fid,'%.16f ', Coeff12(i,:));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\t\t</real>\n\t</matrix>\n</coeff12>\n');
end

if numel(VolCoeff) > 0
    fprintf(fid,'<coeffvol>\n\t<matrix dim1="%d" dim2="%d">\n\t\t<real>\n', ndatapoints, ncoeffs);
    for i=1:ndatapoints
        fprintf(fid,'%.16f ', VolCoeff(i,:));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\t\t</real>\n\t</matrix>\n</coeffvol>\n');
end

fprintf(fid,'</homRectC1>');

fclose(fid);
end

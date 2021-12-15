function write_to_xml(outputfile, params, coeffs)
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
ncoeffs = 4^sum(cellfun('length',params)>0);

names = fieldnames(coeffs);
for k=1:length(names)
    fprintf(fid,'<%s>\n\t<matrix dim1="%d" dim2="%d">\n\t\t<real>\n', names{k}, ndatapoints, ncoeffs);
    for i=1:ndatapoints
        fprintf(fid,'%.16f ', coeffs.(names{k})(i,:));
        fprintf(fid,'\n');
    end
    fprintf(fid,'\t\t</real>\n\t</matrix>\n</%s>\n', names{k});
end

fprintf(fid,'</homRectC1>');

fclose(fid);
end

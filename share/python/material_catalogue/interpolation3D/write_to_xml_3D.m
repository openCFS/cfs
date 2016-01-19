function write_to_xml_3D(file,m,n,o,a,b,c,Coeff11,Coeff12,Coeff22,Coeff33)
filename = sprintf('%s/coeff.xml',file);
fid = fopen(filename,'wt'); 
fprintf(fid,'<homRectC1>\n <a>\n<matrix dim1="%d" dim2="1">\n<real>\n',m+1);
for i=1:m+1
fprintf(fid,'%f ',a(i));
end
fprintf(fid,'\n </real>\n </matrix>\n </a>\n');
fprintf(fid,'<b>\n<matrix dim1="%d" dim2="1">\n<real>\n',n+1);
for i=1:n+1
fprintf(fid,'%f ',b(i));
end
fprintf(fid,'\n </real>\n</matrix>\n</b>');
fprintf(fid,'<c>\n<matrix dim1="%d" dim2="1">\n<real>\n',o+1);
for i=1:o+1
fprintf(fid,'%f ',c(i));
end
fprintf(fid,'\n </real>\n</matrix>\n</c>');
ende = m*n*o;
fprintf(fid,'<coeff11>\n<matrix dim1="%d" dim2="64">\n <real>\n',ende);
for i=1:ende
    for j=1:64
        fprintf(fid,'%f ', Coeff11(i,j));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'</real>\n </matrix>\n </coeff11>\n');
fprintf(fid,'<coeff12>\n<matrix dim1="%d" dim2="64">\n <real>\n',ende);
for i=1:ende
    for j=1:64
        fprintf(fid,'%f ', Coeff12(i,j));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'</real>\n </matrix>\n </coeff12>\n');
fprintf(fid,'<coeff22>\n<matrix dim1="%d" dim2="64">\n <real>\n',ende);
for i=1:ende
    for j=1:64
        fprintf(fid,'%f ', Coeff22(i,j));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'</real>\n </matrix>\n </coeff22>\n');

fprintf(fid,'<coeff33>\n<matrix dim1="%d" dim2="64">\n <real>\n',ende);
for i=1:ende
    for j=1:64
        fprintf(fid,'%f ', Coeff33(i,j));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'</real>\n </matrix>\n </coeff33>\n');
fprintf(fid,'</homRectC1>');
fclose(fid);
end
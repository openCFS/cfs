function write_to_xml_msfem(file,m,n,a,b,Coeff)
% Erzeugt .xml file mit namen file fuer CFS
filename = file;
fid = fopen(filename,'wt'); 
fprintf(fid,'<MSFEMC1>\n <a>\n<matrix dim1="1" dim2="%d">\n<row id = "1">\n',m+1);
for i=1:m+1
fprintf(fid,'<col id ="%d" data="%.16f"/>\n',i,a(i));
end
fprintf(fid,'</row>\n </matrix>\n </a>\n');
fprintf(fid,'<b>\n<matrix dim1="1" dim2="%d">\n<row id = "1">\n',n+1);
for i=1:n+1
fprintf(fid,'<col id ="%d" data="%.16f"/>\n',i,b(i));
end
fprintf(fid,'</row>\n </matrix>\n </b>\n');
index = {'11','12','13','14','15','16','17','18','22','23','24','25','26','27','28'...
    ,'33','34','35','36','37','38','44','45','46','47','48','55','56','57','58','66','67','68'...
    ,'77','78','88'};
ende = m*n;
for k = 1:36
    fprintf(fid,'<coeff%s>\n<matrix dim1="%d" dim2="16">\n',index{k},ende);
    for i=1:ende
        fprintf(fid,'<row id="%d">\n',i);
        for j=1:16
            fprintf(fid,'<col id = "%d" data="%.16f"/>\n',j, Coeff(k,i,j));
        end
        fprintf(fid,'</row>\n');
    end
    fprintf(fid,'</matrix>\n </coeff%s>\n',index{k});
end
fprintf(fid,'</MSFEMC1>');
fclose(fid);
end


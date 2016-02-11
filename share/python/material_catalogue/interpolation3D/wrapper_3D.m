function wrapper_3D(file)
load detailed_stats
list= detailed_stats;
m = list(1,1)-1;
n = list(1,2)-1;
o = list(1,3)-1;
da = 1/m;
db = 1/n;
dc = 1/o;
a = [0:da:1];
b = [0:db:1];
c = [0:dc:1];

E11 = zeros(m+1,n+1);
for i=2:size(list,1)
    E11(list(i,1)+1,list(i,2)+1) = list(i,3);
end
E12 = zeros(m+1,n+1);
for i=2:size(list,1)
    E12(list(i,1)+1,list(i,2)+1) = list(i,4);
end

E22 = zeros(m+1,n+1);
for i=2:size(list,1)
    E22(list(i,1)+1,list(i,2)+1) = list(i,5);
end
E33 = zeros(m+1,n+1);
for i=2:size(list,1)
    E33(list(i,1)+1,list(i,2)+1) = list(i,6);
end
% Coefficients for bicubic interpolation polynomial
[Coeff11] = tricubic_offline(a,b,c,E11);
[Coeff12] = tricubic_offline(a,b,c,E12);
[Coeff22] = tricubic_offline(a,b,c,E22);
[Coeff33] = tricubic_offline(a,b,c,E33);
%test(a,b,E11);


write_to_xml(file,m,n,o,a,b,c,Coeff11,Coeff12,Coeff22,Coeff33);

%Save interpolation coefficients for E11
% filename = sprintf('%s/coeff11.txt',file); 
% fid = fopen(filename,'wt');
% fprintf(fid,'%d %d\n',m,n);
% ende = m*n;
% for i=1:ende
%     for j=1:16
%         fprintf(fid,'%f ', Coeff11(i,j));
%     end
%     fprintf(fid,'\n');
% end
% fclose(fid);
% %Save interpolation coefficients for E12
% filename = sprintf('%s/coeff12.txt',file); 
% fid = fopen(filename,'wt');
% fprintf(fid,'%d %d\n',m,n);
% ende = m*n;
% for i=1:ende
%     for j=1:16
%         fprintf(fid,'%f ', Coeff12(i,j));
%     end
%     fprintf(fid,'\n');
% end
% fclose(fid); 
% %Save interpolation coefficients for E22
% filename = sprintf('%s/coeff22.txt',file); 
% fid = fopen(filename,'wt');
% fprintf(fid,'%d %d\n',m,n);
% ende = m*n;
% for i=1:ende
%     for j=1:16
%         fprintf(fid,'%f ', Coeff22(i,j));
%     end
%     fprintf(fid,'\n');
% end
% fclose(fid); 
% %Save interpolation coefficients for E33
% filename = sprintf('%s/coeff33.txt',file); 
% fid = fopen(filename,'wt');
% fprintf(fid,'%d %d\n',m,n);
% ende = m*n;
% for i=1:ende
%     for j=1:16
%         fprintf(fid,'%f ', Coeff33(i,j));
%     end
%     fprintf(fid,'\n');
% end
%fclose(fid); 

end

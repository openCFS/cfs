function [Coeff11,Coeff12,Coeff22,Coeff33,a,b,E11,E12,E22,E33] = main(file)
%Einlesen_der Materialtensoren aus dem Materialkatalog
load detailed_stats_5_cross
% Rotation des Materialkatalogs
%angle = 0;%pi/4;
%list = detailed_stats_10;
%list = rotate_list(detailed_stats_10,angle);
list = detailed_stats_5_cross;
m = list(1,1);
n = list(1,2);
da = 1/m;
db = 1/n;
a = [0:da:1];
b = [0:db:1];
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
[Coeff11] = bicubic_offline(a,b,E11,1);
[Coeff12] = bicubic_offline(a,b,E12,3);
[Coeff22] = bicubic_offline(a,b,E22,2);
[Coeff33] = bicubic_offline(a,b,E33,3);
%test(a,b,E11);


write_to_xml(file,m,n,a,b,Coeff11,Coeff12,Coeff22,Coeff33);

end

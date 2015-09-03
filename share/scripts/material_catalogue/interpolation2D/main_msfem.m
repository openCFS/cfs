function [Coeff,a,b,A] = main_msfem(inputfile,outputfile,opt)

if nargin < 3
% option if only 1 entry and not the full material catalogue was calculated
    opt = 0; 
end
%read material catalogue
list = load(inputfile);
m = 1/list(1,1);
n = 1/list(1,2);
da = 1/m;
db = 1/n;
a = [0:da:1];
b = [0:db:1];
A = zeros(36,m+1,n+1);
for i=3:38
    for j=2:size(list,1)
        if opt == 1
            A(i-2,list(2,1)+1,list(2,2)+1) = list(2,i);
        else
            A(i-2,list(j,1)+1,list(j,2)+1) = list(j,i); 
        end
    end
end
Coeff = zeros(36,(m)*(n),16);
% Coefficients for bicubic interpolation polynomial
E = zeros(m+1,n+1);
for i =1:36
    E(:,:) = A(i,:,:);
    Coeff(i,:,:) = bicubic_offline(a,b,E,[],[],[],[]);
end
write_to_xml_msfem(outputfile,m,n,a,b,Coeff);
end

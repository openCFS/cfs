function [res] = bicubic_int(Coeff,a,b,x1,x2)
% Auswertung des Interpolationspolynoms an der Stelle x1,x2,
% diskretisierte Dickenvektoren a, b
% Interpolationskoeffizienten Coeff
m = length(a);
n = length(b);
for i=1:m-1
    if a(i) <= x1 && x1 < a(i+1)
        al = a(i);
        au = a(i+1);
        j=i;
    elseif x1 == a(m)
        al = a(m-1);
        au = a(m);
        j=m-1;
    elseif x1 > a(m)
       display(sprintf('x1 out of bounds')); 
       break;
    end
end
for i=1:(n-1)
    if b(i) <= x2 && x2 < b(i+1)
        bl = b(i);
        bu = b(i+1);
        k = i;
    elseif x2 == b(n)
        bl = b(n-1);
        bu = b(n);
        k=n-1;
    elseif x2 > b(n)
       display(sprintf('x2 out of bounds')); 
       break;
    end
end
ct = Coeff((n-1)*(j-1)+k,:);
da = au-al;
db = bu-bl;
t = (x1-al)/da;
u = (x2-bl)/db;
res = 0;
for i=0:3
    for j=0:3
        res = res + ct(i + 4*j + 1) * t^i * u^j;
    end
end
end

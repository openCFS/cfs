function [res] = tricubic_int(E,a,b,c,x1,x2,x3)
% Tricubic interpolation within a grid square. 
m = length(a);
n = length(b);
o = length(c);
j=-1;
k=-1;
l = -1;
for i=1:m-1
    if a(i) <= x1 && x1 < a(i+1)
       j=i;
       al = j;
       au = j+1;
    elseif x1 == a(m)
        j=m-1;
        al = j;
       au = j+1;
    elseif x1 > a(m)
       display(sprintf('x1 out of bounds')); 
       break;
    end
end
for i=1:(n-1)
    if b(i) <= x2 && x2 < b(i+1)
        k = i;
        bl=k;
        bu = k+1;
    elseif x2 == b(n)
        k=n-1;
        bl=k;
        bu = k+1;
    elseif x2 > b(n)
       display(sprintf('x2 out of bounds')); 
       break;
    end
end
for i=1:(o-1)
    if c(i) <= x3 && x3 < c(i+1)
        l = i;
        cl=l;
        cu = l+1;
    elseif x3 == c(n)
        l=o-1;
        cl=l;
        cu = l+1;
    elseif x3 > c(n)
       display(sprintf('x3 out of bounds')); 
       break;
    end
end
da = au-al;
db = bu-bl;
dc = cu-cl;
[Eint Eda Edb Edc Edadb Edadc Edbdc Edadbdc] = tricubic_parderiv(j,k,l,E,da,db,dc); 
ct = tricubic_coeff(Eint, Eda, Edb, Edc, Edadb,Edadc,Edbdc,Edadbdc, da, db,dc);
%ct = Coeff((n-1)*(o-1)*(j-1)+(o-1)*(k-1)+l,:);
res = 0;
 for i=0:3
    for j=0:3
      for k=0:3 
          res = res + ct(1+i+4*j+16*k)*x1^i*x2^j*x3^k;
      end
    end
 end
end
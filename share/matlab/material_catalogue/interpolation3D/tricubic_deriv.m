function [deriv1,deriv2, deriv3] = tricubic_deriv(Coeff,a,b,c,x1,x2,x3)
% Tricubic interpolation within a grid square.
m = length(a);
n = length(b);
o = length(c);
% Choose the right intervall al, au and bl and bu
j=-1;
k=-1;
l = -1;
for i=1:m-1 
    if a(i) <= x1 && x1 < a(i+1)
       j=i;
    elseif x1 == a(m)
        j=m-1;
    elseif x1 > a(m)
       display(sprintf('x1 out of bounds')); 
       break;
    end
end
for i=1:n-1
    if b(i) <= x2 && x2 < b(i+1)
        k = i;
    elseif x2 == b(n)
        k=n-1;
    elseif x2 > b(n)
       display(sprintf('x2 out of bounds')); 
       break;
    end
end
for i=1:o-1
    if c(i) <= x3 && x3 < c(i+1)
        l = i;
    elseif x3 == c(o)
        l=o-1;
    elseif x3 > c(o)
       display(sprintf('x3 out of bounds')); 
       break;
    end
end
ct = Coeff((o-1)*(n-1)*(j-1)+(o-1)*(k-1)+l,:);
deriv1 = 0;
deriv2 = 0;
deriv3 = 0;
 for i=0:3
    for j=0:3
      for k=0:3
        deriv1 = deriv1 + ct(1+i+4*j+16*k)*i*x^(i-1)*y^j*z^k;
        deriv2 = deriv2 + ct(1+i+4*j+16*k)*x^i*j*y^(j-1)*z^k;
        deriv3 = deriv3 + ct(1+i+4*j+16*k)*x^i*y^j*k*z^(k-1);
      end
    end
 end
end
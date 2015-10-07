function [deriv1,deriv2] = bicubic_deriv(Coeff,a,b,x1,x2)
% Auswertung des Interpolationspolynoms an der Stelle x1 und x2 für die zwei Ableitungen

m = length(a);
n = length(b);
% Choose the right intervall al, au and bl and bu
j=-1;
k=-1;
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
for i=1:n-1
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
da=au-al;
db=bu-bl;
ct = Coeff((n-1)*(j-1)+k,:);
c = zeros(4,4);
for i=1:4
    for j=1:4
       c(i,j) = ct((i-1)*4+j);
    end 
end
t=(x1-al)/da;
u=(x2-bl)/db;
% Derivative with respect to x1
deriv1 = 0;
for i = 4:-1:2
deriv1=t*(deriv1)+((i-1)/da)*(((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1));
end
% Derivative with respect to x2
deriv2 = 0;
for i = 4:-1:1
deriv2=t*(deriv2)+(1/db)*((c(i,4)*3*u+2*c(i,3))*u+c(i,2));
end
end

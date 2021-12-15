function [E_out,elem_vol] = tricubic_int(a,b,c,x1,x2,x3,coeff11,coeff12,coeff22,coeff13,coeff23,coeff33,coeff44,coeff55,coeff66,coeffvol)
% Tricubic interpolation of 3D 6x6 orthotropic stiffness tensor within a
% grid square + interpolated scalar volume of the element
m = length(a);
n = length(b);
o = length(c);
j=-1;
k=-1;
l = -1;
for i=1:m-1
    if a(i) <= x1 && x1 < a(i+1)
       j=i;
       break;
    elseif x1 == a(m)
        j=m-1;
    elseif x1 > a(m)
        j = m-1;
        if x1 > 1.01
            display(sprintf('x1 out of bounds, SNOPT violates bounds')); 
        end
       x1 = 1;
       break;
    elseif x1 <0
        j = 1;
        if x1 < -0.01
            display(sprintf('x1 out of bounds, SNOPT violates bounds')); 
        end
        x1 = 0;
        break;
    end
end
for i=1:(n-1)
    if b(i) <= x2 && x2 < b(i+1)
        k = i;
    elseif x2 == b(n)
        k=n-1;
    elseif x2 > b(n)
        k = n-1;
        if x2 > 1.01
            display(sprintf('x2 out of bounds, SNOPT violates bounds')); 
        end
       x2 = 1;
       break;
    elseif x2 <0
        k = 1;
        if x2 < -0.01
            display(sprintf('x2 out of bounds, SNOPT violates bounds')); 
        end
        x2 = 0;
        break;
    end
end
for i=1:(o-1)
    if c(i) <= x3 && x3 < c(i+1)
        l = i;
    elseif x3 == c(n)
        l=o-1;
    elseif x3 > c(n)
        l = o-1;
        if x3 > 1.01
            display(sprintf('x3 out of bounds, SNOPT violates bounds')); 
        end
        x3 = 1;
       break;
    elseif x3 <0
        l = 1;
        if x3 < -0.01
            display(sprintf('x3 out of bounds, SNOPT violates bounds')); 
        end
        x3 = 0;
        break;
    end
end
u = (x1 - a(j))/(a(j)-a(j-1));
t = (x2 - b(k))/(b(k)-b(k-1));
v = (x3 - c(l))/(c(l)-c(l-1));
E_out = zeros(6,6);
elem_vol = 0;
 for ii=0:3
    for jj=0:3
      for kk=0:3 
          E_out(1,1) = E_out(1,1) + coeff11((n-1)*(o-1)*j+(o-1)*k+l,ii+4*jj+16*kk+1)*u^ii*t^jj*v^kk;
          E_out(1,2) = E_out(1,2) + coeff12((n-1)*(o-1)*j+(o-1)*k+l,ii+4*jj+16*kk+1)*u^ii*t^jj*v^kk;
          E_out(1,3) = E_out(1,3) + coeff13((n-1)*(o-1)*j+(o-1)*k+l,ii+4*jj+16*kk+1)*u^ii*t^jj*v^kk;
          E_out(2,2) = E_out(2,2) + coeff22((n-1)*(o-1)*j+(o-1)*k+l,ii+4*jj+16*kk+1)*u^ii*t^jj*v^kk;
          E_out(2,3) = E_out(2,3) + coeff23((n-1)*(o-1)*j+(o-1)*k+l,ii+4*jj+16*kk+1)*u^ii*t^jj*v^kk;
          E_out(3,3) = E_out(3,3) + coeff33((n-1)*(o-1)*j+(o-1)*k+l,ii+4*jj+16*kk+1)*u^ii*t^jj*v^kk;
          E_out(4,4) = E_out(4,4) + coeff44((n-1)*(o-1)*j+(o-1)*k+l,ii+4*jj+16*kk+1)*u^ii*t^jj*v^kk;
          E_out(5,5) = E_out(5,5) + coeff55((n-1)*(o-1)*j+(o-1)*k+l,ii+4*jj+16*kk+1)*u^ii*t^jj*v^kk;
          E_out(6,6) = E_out(6,6) + coeff66((n-1)*(o-1)*j+(o-1)*k+l,ii+4*jj+16*kk+1)*u^ii*t^jj*v^kk;
          elem_vol = elem_vol + coeffvol((o-1)*(n-1)*j+(o-1)*k+l,ii+4*jj+16*kk+1)*u^ii*t^jj*v^kk;
      end
    end
 end
 E_out(2,1) = E_out(1,2);
 E_out(3,1) = E_out(1,3);
 E_out(3,2) = E_out(2,3);
end
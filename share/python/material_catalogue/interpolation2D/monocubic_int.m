function [res] = monocubic_int(Coeff, a, x)
% Auswertung des Interpolationspolynoms an der Stelle x
% diskretisierte Dickenvektoren a
% Interpolationskoeffizienten Coeff
m = length(a);
k = -1;
for i=1:m-1
    if a(i) <= x && x < a(i+1)
        al = a(i);
        au = a(i+1);
        k=i;
    elseif x == a(m)
        al = a(m-1);
        au = a(m);
        k=m-1;
    elseif x > a(m)
        fprintf('x out of bounds\n');
        break;
    end
end

c = Coeff(k,:);
t = (x-al)/(au-al);
res = 0;

for i=0:3
    res = res + c(i+1) * t^i;
end

end

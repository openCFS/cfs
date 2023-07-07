function [res] = monocubic_int(Coeff, a, x, deriv)
% Auswertung des Interpolationspolynoms an der Stelle x
% diskretisierte Dickenvektoren a
% Interpolationskoeffizienten Coeff

if isrow(a)
    a = a';
end
if isrow(x)
    x = x';
end

if nargin < 4
    deriv = 0;
end

switch 'fast'
    case 'fast'
        if any(x < a(1) | x > a(end))
            outOfBoundEntries = [find(x < a(1) | x > a(end)), x(x < a(1) | x > a(end))]
            error('x out of bounds');
        end
        [~,k] = max(x' < a);
        k = k-1; k(k==0) = length(a)-1;
        al = a(k);
        au = a(k+1);

    case 'slow'
        if length(x) > 1
            res = zeros(size(x));
            for i=1:length(x)
                y = x(i);
                res(i) = monocubic_int(Coeff, a, y, deriv);
            end
            return
        else
            m = length(a);
            k = -1;
            for i=1:m-1
                if a(i) <= x && x < a(i+1)
                    al = a(i);
                    au = a(i+1);
                    k=i;
                    break
                elseif x == a(m)
                    al = a(m-1);
                    au = a(m);
                    k=m-1;
                    break
                elseif x < a(1) || x > a(m)
                    x
                    a(1)
                    a(m)
                    error('x out of bounds');
            %         al = a(m-1);
            %         au = a(m);
            %         k=m-1;
                    break
                end
            end
        end
end

c = Coeff(k,:);
t = (x-al)./(au-al);
res = 0;

if deriv == 0
    for i=0:3
        res = res + c(:,i+1) .* t.^i;
    end
elseif deriv == 1
    for i=1:3
        res = res + c(:,i+1) .* t.^(i-1) * i ./ (au-al);
    end
elseif deriv == 2
    for i=2:3
        res = res + c(:,i+1) .* t.^(i-2) * i * (i-1) ./ (au-al) ./ (au-al);
    end
elseif deriv == 3
    res = res + c(:,4) * 6;
end

end

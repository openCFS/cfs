function polarplottensor(T, xoff, yoff, c, showlowortho)

if nargin < 3
    xoff = 0;
    yoff = 0;
end
if nargin < 4
    c = 1;
end
if nargin < 5
    showlowortho = 0;
end

x = 0:0.01:2*pi;
pl = T(1,1)*cos(x).^4 + T(2,2)*sin(x).^4 + (4*T(3,3)+T(1,2)+T(2,1))*cos(x).^2.*sin(x).^2 - 2*(T(1,3)+T(3,1))*cos(x).^3.*sin(x) - 2*(T(2,3)+T(3,2))*cos(x).*sin(x).^3;
[x,y] = pol2cart(x, pl);
x = x + xoff;
y = y + yoff;
%polar(x, pl);
patch(x,y,c*ones(size(x)));

if showlowortho
    t = 0:0.001:2*pi;
    cost = cos(t);
    sint = sin(t);
    t31 = (T(1,1)-T(1,2))*cost.^3.*sint + (T(2,1)-T(2,2))*cost.*sint.^3 + (T(3,2)-T(3,1))*2*cost.^2.*sint.^2 + (T(1,3)*cost.^2+T(2,3)*sint.^2).*(cost.^2-sint.^2) - 2*cost.*sint.*(cost.^2-sint.^2);
    t32 = (T(1,1)-T(2,1))*cost.*sint.^3 + (T(1,2)-T(2,2))*cost.^3.*sint + (T(1,3)-T(2,3))*2*cost.^2.*sint.^2 + (T(3,1)*sint.^2+T(3,2)*cost.^2).*(cost.^2-sint.^2) + 2*cost.*sint.*(cost.^2-sint.^2);
    
    r = t31.^2 + t32.^2;
    
    [x,y] = pol2cart(t(r < 1e-3*trace(T)), 1);
    set(line([xoff - x; xoff + x], [yoff - y; yoff + y]), 'Color', 'b');
end



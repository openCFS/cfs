function polarplotdata(D,nx,ny,showlowortho)

% example:
% D = load('textfile');           containing one line for each tensor entry
% for each element
% polarplotdata(D, 20, 40, 1);

R = reshape(D, 6, nx, ny);

for x = 1:nx
    for y = 1:ny
        h = R(:,y,x);
        T = [h(1), h(2), h(4); h(2), h(3), h(5); h(4), h(5), h(6)];
        polarplottensor(T, 40 - 2*y, 2*x, trace(T), showlowortho);
    end
end
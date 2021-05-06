function monoplot_interp(inputEhom, inputXML)

list = load(inputEhom);
coeff = xml2struct(inputXML);
Coeff11 = str2num(coeff.homRectC1.coeff11.matrix.real.Text);
Coeff22 = str2num(coeff.homRectC1.coeff22.matrix.real.Text);
Coeff33 = str2num(coeff.homRectC1.coeff33.matrix.real.Text);
Coeff12 = str2num(coeff.homRectC1.coeff12.matrix.real.Text);


m = list(1,1)-1;
a = list(2:end,1);
E11_grid = zeros(m+1,1);
E12_grid = zeros(m+1,1);
E22_grid = zeros(m+1,1);
E33_grid = zeros(m+1,1);
for i=2:size(list,1)
    E11_grid(i-1) = list(i,2);
    E22_grid(i-1) = list(i,3);
    E33_grid(i-1) = list(i,4);
    E12_grid(i-1) = list(i,7);
end

figure;
XX = 0:.01:1;
ZZ = zeros(101,1);
for ii=1:101
    ZZ(ii) = monocubic_int(Coeff11, a, XX(ii));
end
plot(XX,ZZ);
hold on;
plot(a, E11_grid, '*');
hold off;

figure;
ZZ = zeros(101,1);
for ii=1:101
    ZZ(ii) = monocubic_int(Coeff22, a, XX(ii));
end
plot(XX,ZZ);
hold on;
plot(a, E22_grid, '*');
hold off;

figure;
ZZ = zeros(101,1);
for ii=1:101
    ZZ(ii) = monocubic_int(Coeff33, a, XX(ii));
end
plot(XX,ZZ);
hold on;
plot(a, E33_grid, '*');
hold off;

figure;
ZZ = zeros(101,1);
for ii=1:101
    ZZ(ii) = monocubic_int(Coeff12, a, XX(ii));
end
plot(XX,ZZ);
hold on;
plot(a, E12_grid, '*');
hold off;

clf;
exact = readmatrix('Jexact_xz.ref.txt','NumHeaderLines',1);
model = readmatrix('Jmodel_xz.ref.txt','NumHeaderLines',1);

angles = 0:pi/180:pi;
nangles = size(angles,2);
valuesModel = reshape(model(:,4), [nangles nangles]);
valuesExact = reshape(exact(:,4), [nangles nangles]);

[xangles,zangles] = meshgrid(angles,angles);

surf(xangles,zangles,valuesExact,'DisplayName','exact','LineStyle','none')
hold on
surf(xangles,zangles,valuesModel,'DisplayName','model')

xlabel('angle about z-axis')
ylabel('angle about x-axis')
zlabel('compliance')

legend

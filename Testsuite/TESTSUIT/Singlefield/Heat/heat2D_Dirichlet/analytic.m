clear all;
close all;

L = 1;
x = 0:0.001:1;
y = 0:0.5:L;
[X,Y] = meshgrid(x,y);

% truncation of eigenvalue expansion
N = 50;
M = 50;

% coefficients
a_m = @(m) m*pi/L;
b_n = @(n) (2*n-1)*pi/2;
bm = @(m) 4*sin(a_m(m)*L)/((2*L-sin(2*a_m(m)*L))*a_m(m));
bmn = @(m,n) bm(m)*(cos(b_n(n))-1)/((0.5-sin(2*b_n(n))/(4*b_n(n)))*b_n(n));
b0n = @(n) (cos(b_n(n))-1)/((0.5-sin(2*b_n(n))/(4*b_n(n)))*b_n(n));

% time step for evaluation
t = 0.1;

%% compute section
T1=0; 
for n = 1 : N
    T1 = T1 + b0n(n)*exp(-(b_n(n)^2*t)).*sin(b_n(n).*X);
end
T2 = 0;
for m = 1 : N
    T3 = bm(m);
    for n = 1 : N
        T3 = T3 + (b_n(n)^2*exp(-(b_n(n)^2+a_m(m)^2)*t)+a_m(m)^2)/(b_n(n)^2+a_m(m)^2)*bmn(m,n)*sin(b_n(n).*X);
    end
    T2 = T2 + T3.*cos(a_m(m).*Y);
end

T = 1 + T1 + T2;

surf(X,Y,T);

axis([0 1 0 L 0 1])
xlabel('x');
ylabel('y');
zlabel('T');
drawnow;
colormap gray
colorbar
shading interp
view(2)
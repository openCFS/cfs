function h = monoplot_interp(inputEhom, inputXML)

list = load(inputEhom);
coeffs = readXML(inputXML);
Coeff11 = coeffs.coeff11;
Coeff22 = coeffs.coeff22;
Coeff33 = coeffs.coeff33;
Coeff12 = coeffs.coeff12;
Coeff66 = coeffs.coeff66;


m = list(1,1)-1;
a = list(2:end,1);
E11_grid = zeros(m+1,1);
E12_grid = zeros(m+1,1);
E22_grid = zeros(m+1,1);
E33_grid = zeros(m+1,1);
E66_grid = zeros(m+1,1);
for i=2:size(list,1)
    E11_grid(i-1) = list(i,2);
    E22_grid(i-1) = list(i,3);
    E33_grid(i-1) = list(i,4);
    E12_grid(i-1) = list(i,7);
    E12_grid(i-1) = list(i,3);
    E66_grid(i-1) = list(i,end);
end

% figure;
XX = linspace(min(a),max(a),101);
ZZ = zeros(101,1);
for ii=1:101
    ZZ(ii) = monocubic_int(Coeff11, sort(a), XX(ii));
end
h = plot(XX,ZZ,'LineWidth',2);
hold on;
plot(a, E11_grid, '+', 'LineWidth',2,'MarkerSize',10);
% plot(XX,XX.^3*list(end,2))
% plot(XX,XX.^2*list(end,2))
% q = 2.8;
% plot(XX,XX./(1+q*(1-XX))*list(end,2))
% legend('int','dat', 'rho^2','rho^3','RAMP')
hold off;

set(gca,'FontSize',18)

% 
% figure;
% ZZ = zeros(101,1);
% for ii=1:101
%     ZZ(ii) = monocubic_int(Coeff22, sort(a), XX(ii));
% end
% plot(XX,ZZ);
% hold on;
% plot(a, E22_grid, '*');
% hold off;
% 
% figure;
% ZZ = zeros(101,1);
% for ii=1:101
%     ZZ(ii) = monocubic_int(Coeff33, sort(a), XX(ii));
% end
% plot(XX,ZZ);
% hold on;
% plot(a, E33_grid, '*');
% hold off;
% 
% figure;
% ZZ = zeros(101,1);
% for ii=1:101
%     ZZ(ii) = monocubic_int(Coeff12, sort(a), XX(ii));
% end
% plot(XX,ZZ);
% hold on;
% plot(a, E12_grid, '*');
% hold off;

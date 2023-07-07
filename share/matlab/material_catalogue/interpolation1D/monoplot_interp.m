function h = monoplot_interp(inputEhom, inputXML)

[coeffs, params] = readXML(inputXML);
is2d = norm(coeffs.coeff66)/norm(coeffs.coeff11) < 1e-4;

if(inputEhom ~= "")
    list = load(inputEhom);
    m = list(1,1);
    a = list(2:end,1);
else
    m = length(params.param1);
    a = params.param1;
end

sel = a<=1;

Coeff11 = coeffs.coeff11;
Coeff22 = coeffs.coeff22;
Coeff33 = coeffs.coeff33;
Coeff12 = coeffs.coeff12;
Coeff66 = coeffs.coeff66;
CoeffMLF = coeffs.microloadfactor;

if(inputEhom ~= "")
    E11_grid = zeros(m,1);
    E12_grid = zeros(m,1);
    E22_grid = zeros(m,1);
    E33_grid = zeros(m,1);
    E66_grid = zeros(m,1);
    MLF_grid = zeros(m,1);
    for i=2:size(list,1)
        if size(list,2) > 2
            if is2d
                E11_grid(i-1) = list(i,2);
                E22_grid(i-1) = list(i,3);
                E33_grid(i-1) = list(i,4);
                E12_grid(i-1) = list(i,7);
        %         E12_grid(i-1) = list(i,3);
                E66_grid(i-1) = list(i,end);
            else
                E11_grid(i-1) = list(i,2);
                E22_grid(i-1) = list(i,8);
                E33_grid(i-1) = list(i,13);
                E12_grid(i-1) = list(i,3);
        %         E12_grid(i-1) = list(i,3);
                E66_grid(i-1) = list(i,22);
                MLF_grid(i-1) = list(i,end);
            end
        else
            MLF_grid(i-1) = list(i,2);
        end
    end
end

figure;
%XX = linspace(min(a(sel)),max(a(sel)),101);
x0 = 0.6;
XX = linspace(min(a(sel)),x0,101);
ZZ = zeros(101,1);
for ii=1:101
    ZZ(ii) = monocubic_int(CoeffMLF, sort(a(sel)), XX(ii));
end
ev0 = monocubic_int(CoeffMLF, sort(a(sel)), x0);
dev0 = monocubic_int(CoeffMLF, sort(a(sel)), x0, 1);
ddev0 = monocubic_int(CoeffMLF, sort(a(sel)), x0, 2);
XX2 = linspace(max(a(sel)),1,101);
XX2 = linspace(x0,1,101);
ZZ2 = ddev0/2 * (XX2-x0).^2 + dev0 * (XX2-x0) + ev0;

XX = [XX,XX2]; ZZ = [ZZ;ZZ2'];

% monocubic_int(CoeffMLF, sort(a), 1)
plot(XX,ZZ,'LineWidth',2);
if(inputEhom ~= "")
    hold on;
    plot(a(sel), MLF_grid(sel), '+', 'LineWidth',2,'MarkerSize',10);
    sel(end) = 1;
%     plot(a(~sel), MLF_grid(~sel), 'o', 'LineWidth',2,'MarkerSize',10);
    hold off;
end
title('microloadfactor')
% plot(XX,XX.^3*list(end,2))
% plot(XX,XX.^2*list(end,2))
% q = 2.8;
% plot(XX,XX./(1+q*(1-XX))*list(end,2))
% legend('int','dat', 'rho^2','rho^3','RAMP')
set(gca,'FontSize',18)


% figure;
% XX = linspace(min(a(sel)),max(a(sel)),101);
% ZZ = zeros(101,1);
% for ii=1:101
%     ZZ(ii) = monocubic_int(Coeff11, sort(a), XX(ii));
% end
% plot(XX,ZZ,'LineWidth',2);
% if(inputEhom ~= "")
%     hold on;
%     plot(a, E11_grid, '*');
%     hold off;
% end
% title('coeff11')

% figure;
% ZZ = zeros(101,1);
% for ii=1:101
%     ZZ(ii) = monocubic_int(Coeff22, sort(a), XX(ii));
% end
% plot(XX,ZZ);
% hold on;
% plot(a, E22_grid, '*');
% hold off;
% title('coeff22')
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
% title('coeff33')
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
% title('coeff12')

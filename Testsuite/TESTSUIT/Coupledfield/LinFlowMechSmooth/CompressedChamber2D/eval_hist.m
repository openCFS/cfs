%% Settings
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')

sf = 1;

fig_pos_x = 2;
fig_pos_y = 2;
fig_width = 16*sf;
fig_height = 9*sf;
real_fs = 12*sf;
lw = 1*sf;

ftnr = 1;

print_bit = 0;

rho0 = 1.2;
c0 = 340;
p0 = 1e5;
kappa = 1.4;
V0 = 1e-3;

%% read in simulation data

f_path_S1 = 'history/CompressedChamber2D-fluidMechPressure-node-1-S1.hist';
f_path_S1_smooth = 'history/CompressedChamber2D-smoothDisplacement-node-1-S1.hist';

% read files
fid = fopen(f_path_S1,'r');
S1_pres = [];
tline = fgets(fid);
while ischar(tline)
    parts = textscan(tline, '%f %f;');
    if numel(parts{1}) > 0
        S1_pres = [ S1_pres ; parts{:} ];
    end
    tline = fgets(fid);
end
fclose(fid);

fid = fopen(f_path_S1_smooth,'r');
S1_smooth = [];
tline = fgets(fid);
while ischar(tline)
    parts = textscan(tline, '%f %f;');
    if numel(parts{1}) > 0
        S1_smooth = [ S1_smooth ; parts{:} ];
    end
    tline = fgets(fid);
end
fclose(fid);

%% eval

dp = zeros(length(S1_pres(:,1)),1);
for ii=1:length(S1_pres(:,1))-1
    dp(ii+1) = S1_pres(ii+1,2)-S1_pres(ii,2);
end

time_vec = S1_pres(:,1);

dV_sim = S1_smooth(:,2); % smoothDisplacement
V_sim = (V0-dV_sim);

%% reference

steps = length(S1_pres);
deltaT = time_vec(2)-time_vec(1);
tau = 25e-6;
%p_exc = 10.*(1-exp(-time_vec./tau));

V_ref = V0*nthroot(p0./(p0+S1_pres(:,2)),kappa);
p_ref = (S1_pres(1,2)+p0)*(1./(1-dV_sim/V0)).^kappa;

%plot(S1_pres(:,1),S1_pres(:,2))

%% plot

fig(1) = figure('units','centimeters','PaperUnits','centimeters','position',...
    [fig_pos_x fig_pos_y fig_width fig_height],'PaperPosition',...
    [0 0 fig_width fig_height],'PaperSize',[fig_width fig_height]); 
plot(V_sim,p_ref,'LineWidth',lw)
hold on
plot(V_sim,S1_pres(:,2)+p0,'LineWidth',lw)
hold off
xlabel('Volume expressed as smoothDisplacment in mm') 
ylabel('fluidMechPressure in Pa') 
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','off')
set(gca, 'LabelFontSizeMultiplier',1.0)
set(gca,'fontsize',real_fs)
grid on
legend('calculated Reference','Simulation results LinFlow ALE')

p_inter = interp1(V_sim,p_ref,V_sim,'linear','extrap');
p_err = abs(S1_pres(:,2)+p0-p_inter)./p_inter;

fig(2) = figure('units','centimeters','PaperUnits','centimeters','position',...
    [fig_pos_x fig_pos_y fig_width fig_height],'PaperPosition',...
    [0 0 fig_width fig_height],'PaperSize',[fig_width fig_height]); 
plot((1-dV_sim)*100,p_err*100,'LineWidth',lw)
xlabel('Rest volume in %') 
ylabel('Error in %') 
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','off')
set(gca, 'LabelFontSizeMultiplier',1.0)
set(gca,'fontsize',real_fs)
grid on
legend('Reference','LinFlow ALE')


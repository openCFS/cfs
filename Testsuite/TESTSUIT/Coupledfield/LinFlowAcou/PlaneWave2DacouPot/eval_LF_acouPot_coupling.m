%% Settings
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')

fig_pos_x = 2;
fig_pos_y = 2;
fig_width = 16;
fig_height = 9;
real_fs = 12;
lw = 1;

% from mat file
rho0 = 1.225;
c0 = sqrt(1.4261e5/rho0);

% simulation setup
S1_x = 400e-3; % distance of the mic-position
dt = 5e-5;
steps = 50;

%% Read in data

f_path = ['history' filesep 'LF_Acou_transient-acouPotentialD1-node-103-S1.hist']; % use filesep for compatibility
data_S1 = dlmread(f_path,'',3,0);

t = data_S1(:,1);
p = data_S1(:,2)*rho0;

%% Process

ta = 0:dt:(steps-1)*dt;
travel_time = S1_x/c0+dt; %we have to add "dt" since CFS already starts one step early (step 0 is at dt and not at 0)

% simple 1D solution (p_\mathrm{a} = \rho_0 c_0 v_\mathrm{a}, for the
% acouPotentialD1 we multiply by \rho_0
Pa = sin(2*pi*1000*(t)).*(1-exp(-(t)/(3e-3)))*c0*rho0;

figure('units','centimeters','PaperUnits','centimeters','position',...
    [fig_pos_x fig_pos_y fig_width fig_height],'PaperPosition',...
    [0 0 fig_width fig_height],'PaperSize',[fig_width fig_height]); 
plot(1e3*t,p,'LineWidth',lw)
hold on
plot(1e3*(ta+travel_time),Pa,'--','LineWidth',lw)
hold off
xlabel('Time in ms') 
ylabel('Sound pressure in Pa') 
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','off')
set(gca, 'LabelFontSizeMultiplier',1.0)
set(gca,'fontsize',real_fs)
grid on
xlim(1e3*[t(1) t(end)])
legend('Simulation','1D reference solution','Location','NorthWest')


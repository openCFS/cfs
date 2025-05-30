%% Settings
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')

real_fs = 10;
save_quality = '-r600';

fig_pos_x = 1;
fig_pos_y = 2;
fig_width = 16;
fig_heigth = 9;


%% Evaluation

eval_max_pos = 1; % evaluate the third local maximum for the wave number

f_path_S1 = 'history/PlaneWave2DGridVelC2-fluidMechPressure-node-201-S1.hist';
f_path_S2 = 'history/PlaneWave2DGridVelC2-fluidMechPressure-node-101-S2.hist';

% read files
fid = fopen(f_path_S1,'r');
S1 = [];
tline = fgets(fid);
while ischar(tline)
    parts = textscan(tline, '%f %f;');
    if numel(parts{1}) > 0
        S1 = [ S1 ; parts{:} ];
    end
    tline = fgets(fid);

end
fclose(fid);

fid = fopen(f_path_S2,'r');
S2 = [];
tline = fgets(fid);
while ischar(tline)
    parts = textscan(tline, '%f %f;');
    if numel(parts{1}) > 0
        S2 = [ S2 ; parts{:} ];
    end
    tline = fgets(fid);

end
fclose(fid);

% since we evaluate at certain nodes at a certain time, we won't get the
% peak exactly. Therefore, we evaluate the peak, search for the next
% zero-crossing and linearly interpolate to get the interpolated time for
% the zero crossing. This is more accurate than evaluating just the nodes.

TF_S1 = islocalmax(S1(:,2),'MinProminence',0.05);
TF_S2 = islocalmax(S2(:,2),'MinProminence',0.05);

TF_S1_ind = find(TF_S1==1,eval_max_pos);
TF_S1_ind = TF_S1_ind(end);

zc_S1_ind = find(S1(TF_S1_ind:end,2)<0);
zc_S1_ind_vec = [zc_S1_ind(1)-1 zc_S1_ind(1)];
zc_S1_t = interp1(S1(TF_S1_ind+zc_S1_ind_vec-1,2),S1(TF_S1_ind+zc_S1_ind_vec-1,1),0);

TF_S2_ind = find(TF_S2==1,eval_max_pos);
TF_S2_ind = TF_S2_ind(end);

zc_S2_ind = find(S2(TF_S2_ind:end,2)<0);
zc_S2_ind_vec = [zc_S2_ind(1)-1 zc_S2_ind(1)];
zc_S2_t = interp1(S2(TF_S2_ind+zc_S2_ind_vec-1,2),S1(TF_S2_ind+zc_S2_ind_vec-1,1),0);

% plot
figure('units','centimeters','PaperUnits','centimeters','position',...
    [fig_pos_x fig_pos_y fig_width fig_heigth],'PaperPosition',...
    [0 0 fig_width fig_heigth])
plot(10^3*S1(:,1),S1(:,2),'LineWidth',1)
hold on
plot(10^3*S2(:,1),S2(:,2),'LineWidth',1)
plot(10^3*S1(TF_S1_ind,1),S1(TF_S1_ind,2),'b*')
plot(10^3*zc_S1_t,0,'b*')
plot(10^3*S2(TF_S2_ind,1),S2(TF_S2_ind,2),'r*')
plot(10^3*zc_S2_t,0,'r*')
hold off
xlabel('Time in ms') 
ylabel('Sound pressure in Pa') 
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','off')
set(gca, 'LabelFontSizeMultiplier',1.0)
set(gca,'fontsize',real_fs)
grid on


%% calculate the wave number

ds = 0.5; % 0.5 m distance between the evaluation points
f = 1000;
omega = 2*pi*f;
v_g = 50; % prescribed grid velocity
K = 1.42709e5; % used on simulation
rho0 = 1.225; % used on simulation
c0 = sqrt(K/rho0);

dt = zc_S2_t-zc_S1_t;
c_sim = ds/dt;
k_sim = omega/c_sim;

k_analytic = omega*(v_g+sqrt(v_g^2+4*c0^2))/(2*c0^2);

disp(['k (Simulation):' num2str(k_sim)])
disp(['k (Analytic):' num2str(k_analytic)])

c_analytic = omega/k_analytic;

disp(['c (Simulation):' num2str(c_sim)])
disp(['c (Analytic):' num2str(c_analytic)])

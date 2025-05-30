%% Settings
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'tex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')

fig_pos_x = 2;
fig_pos_y = 2;
fig_width = 16;
fig_height = 9;
real_fs = 12;
lw = 1;

scale_height = 1;

save_quality = '-r600';
print_bit = 0;

d = 1e-2; % channel width
% verified with sf 25 (simulation parameter)


fig_test = figure;
color_mat = get(gca, 'ColorOrder');
close(fig_test)


%% define paths

normIntensPresExc_path = 'history/IntensityAndPower-fluidMechNormalIntensityPressureOnly-surfElement-874-Exc_LF.hist';
normIntensExc_path = 'history/IntensityAndPower-fluidMechNormalIntensity-surfElement-874-Exc_LF.hist';

normIntensPresIF_path = 'history/IntensityAndPower-fluidMechNormalIntensityPressureOnly-surfElement-2112-IF.hist';
normIntensIF_path = 'history/IntensityAndPower-fluidMechNormalIntensity-surfElement-2112-IF.hist';


powerPresExc_path = 'history/IntensityAndPower-fluidMechPowerPressureOnly-surfRegion-Exc_LF.hist';
powerExc_path = 'history/IntensityAndPower-fluidMechPower-surfRegion-Exc_LF.hist';

powerPresIF_path = 'history/IntensityAndPower-fluidMechPowerPressureOnly-surfRegion-IF.hist';
powerIF_path = 'history/IntensityAndPower-fluidMechPower-surfRegion-IF.hist';


intensPresExc_path = 'history/IntensityAndPower-fluidMechSurfaceIntensityPressureOnly-surfElement-874-Exc_LF.hist';
intensExc_path = 'history/IntensityAndPower-fluidMechSurfaceIntensity-surfElement-874-Exc_LF.hist';

intensPresIF_path = 'history/IntensityAndPower-fluidMechSurfaceIntensityPressureOnly-surfElement-2112-IF.hist';
intensIF_path = 'history/IntensityAndPower-fluidMechSurfaceIntensity-surfElement-2112-IF.hist';


viscDissPowerCh1_path = 'history/IntensityAndPower-fluidMechViscousDissipation-region-Channel_1.hist';
viscDissPowerCh2_path = 'history/IntensityAndPower-fluidMechViscousDissipation-region-Channel_2.hist';


%% read results


normIntensPresExc = readmatrix(normIntensPresExc_path,'FileType','text','Delimiter','  ');
normIntensExc = readmatrix(normIntensExc_path,'FileType','text','Delimiter','  ');

normIntensPresIF = readmatrix(normIntensPresIF_path,'FileType','text','Delimiter','  ');
normIntensIF = readmatrix(normIntensIF_path,'FileType','text','Delimiter','  ');


powerPresExc = readmatrix(powerPresExc_path,'FileType','text','Delimiter','  ');
powerExc = readmatrix(powerExc_path,'FileType','text','Delimiter','  ');

powerPresIF = readmatrix(powerPresIF_path,'FileType','text','Delimiter','  ');
powerIF = readmatrix(powerIF_path,'FileType','text','Delimiter','  ');


intensPresExc = readmatrix(intensPresExc_path,'FileType','text','Delimiter','  ');
intensExc = readmatrix(intensExc_path,'FileType','text','Delimiter','  ');

intensPresIF = readmatrix(intensPresIF_path,'FileType','text','Delimiter','  ');
intensIF = readmatrix(intensIF_path,'FileType','text','Delimiter','  ');


viscDissPowerCh1 = readmatrix(viscDissPowerCh1_path,'FileType','text','Delimiter','  ');
viscDissPowerCh2 = readmatrix(viscDissPowerCh2_path,'FileType','text','Delimiter','  ');


%% intensity

% intensity Exc
fig_plot1 = figure('units','centimeters','PaperUnits','centimeters','position',...
    [fig_pos_x fig_pos_y fig_width scale_height*fig_height],'PaperPosition',...
    [0 0 fig_width scale_height*fig_height],'PaperSize',[fig_width scale_height*fig_height]);
plot(1e3*normIntensPresExc(:,1),-normIntensPresExc(:,2),'--','Color',color_mat(1,:),'Linewidth',1)
hold on
plot(1e3*intensPresExc(:,1),intensPresExc(:,2),'-.','Color',color_mat(2,:),'Linewidth',1)
plot(1e3*intensPresExc(:,1),intensPresExc(:,3),':','Color',color_mat(2,:),'Linewidth',1)

plot(1e3*normIntensExc(:,1),-normIntensExc(:,2),'--','Color',color_mat(3,:),'Linewidth',1)
plot(1e3*intensExc(:,1),intensExc(:,2),'-.','Color',color_mat(4,:),'Linewidth',1)
plot(1e3*intensExc(:,1),intensExc(:,3),':','Color',color_mat(4,:),'Linewidth',1)
hold off
grid on
xlabel('Time in ms')
ylabel('Intensity in W/m$^2$')
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','off')
set(gca, 'LabelFontSizeMultiplier',1.0)
set(gca,'fontsize',real_fs)
title('Exc LF intensity comparison')
legend('normIntensPres','intensPresX','intensPresY','normIntens','intensX','intensY')


% intensity IF
fig_plot2 = figure('units','centimeters','PaperUnits','centimeters','position',...
    [fig_pos_x fig_pos_y fig_width scale_height*fig_height],'PaperPosition',...
    [0 0 fig_width scale_height*fig_height],'PaperSize',[fig_width scale_height*fig_height]);
plot(1e3*normIntensPresIF(:,1),normIntensPresIF(:,2),'--','Color',color_mat(1,:),'Linewidth',1)
hold on
plot(1e3*intensPresIF(:,1),intensPresIF(:,2),'-.','Color',color_mat(2,:),'Linewidth',1)
plot(1e3*intensPresIF(:,1),intensPresIF(:,3),':','Color',color_mat(2,:),'Linewidth',1)

plot(1e3*normIntensIF(:,1),normIntensIF(:,2),'--','Color',color_mat(3,:),'Linewidth',1)
plot(1e3*intensIF(:,1),intensIF(:,2),'-.','Color',color_mat(4,:),'Linewidth',1)
plot(1e3*intensIF(:,1),intensIF(:,3),':','Color',color_mat(4,:),'Linewidth',1)
hold off
grid on
xlabel('Time in ms')
ylabel('Intensity in W/m$^2$')
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','off')
set(gca, 'LabelFontSizeMultiplier',1.0)
set(gca,'fontsize',real_fs)
title('IF intensity comparison')
legend('normIntensPres','intensPresX','intensPresY','normIntens','intensX','intensY')


%% power


powerPresExc_man = normIntensPresExc(:,2)*d;
powerExc_man = normIntensExc(:,2)*d;

powerPresIF_man = normIntensPresIF(:,2)*d;
powerIF_man = normIntensIF(:,2)*d;


% intensity IF
fig_plot2 = figure('units','centimeters','PaperUnits','centimeters','position',...
    [fig_pos_x fig_pos_y fig_width scale_height*fig_height],'PaperPosition',...
    [0 0 fig_width scale_height*fig_height],'PaperSize',[fig_width scale_height*fig_height]);
plot(1e3*powerPresExc(:,1),powerPresExc(:,2),'--','Color',color_mat(1,:),'Linewidth',1)
hold on
plot(1e3*powerPresExc(:,1),powerPresExc_man,'-.','Color',color_mat(2,:),'Linewidth',1)

plot(1e3*powerExc(:,1),powerExc(:,2),'--','Color',color_mat(3,:),'Linewidth',1)
plot(1e3*powerExc(:,1),powerExc_man,'-.','Color',color_mat(4,:),'Linewidth',1)

plot(1e3*powerPresIF(:,1),powerPresIF(:,2),'--','Color',color_mat(5,:),'Linewidth',1)
plot(1e3*powerPresIF(:,1),powerPresIF_man,'-.','Color',color_mat(6,:),'Linewidth',1)

plot(1e3*powerIF(:,1),powerIF(:,2),'--','Color',color_mat(7,:),'Linewidth',1)
plot(1e3*powerIF(:,1),powerIF_man,'-.','Color','k','Linewidth',1)
hold off
grid on
xlabel('Time in ms')
ylabel('Power in W')
set(gca,'TickLabelInterpreter','latex')
set(gca,'XMinorTick','off')
set(gca, 'LabelFontSizeMultiplier',1.0)
set(gca,'fontsize',real_fs)
title('IF intensity comparison')
legend('powerPresExc','powerPresExcMan','powerExc','powerExcMan','powerPresIF','powerPresIFMan','powerIF','powerIFMan')


%% energy


viscDissEnergyCh1 = cumtrapz(viscDissPowerCh1(:,1),viscDissPowerCh1(:,2));
viscDissEnergyCh2 = cumtrapz(viscDissPowerCh2(:,1),viscDissPowerCh2(:,2));

energyIn = cumtrapz(powerExc(:,1),powerExc(:,2));
energyOut = cumtrapz(powerIF(:,1),powerIF(:,2));

energyInPres = cumtrapz(powerPresExc(:,1),powerPresExc(:,2));
energyOutPres = cumtrapz(powerPresIF(:,1),powerPresIF(:,2));


disp(['Energy difference (simplified): ' num2str(energyInPres(end)+energyOutPres(end))])
disp(['Energy difference: ' num2str(energyIn(end)+energyOut(end))])
disp(['Dissipated energy: ' num2str(viscDissEnergyCh1(end))])





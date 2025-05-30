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

K = 1.42709e5;
rho0 = 1.2;
c0 = 340;
kappa = 1.406;
p0 = K/kappa;
V0 = 1e-3;

%% read in simulation data

f_path_S1 = 'history/ElecForceCompressedChamber2D-fluidMechPressure-node-1-S1.hist';
f_path_S1_smooth = 'history/ElecForceCompressedChamber2D-smoothDisplacement-node-1-S1.hist';

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

time_vec = S1_pres(:,1);

% could be used to approximate the whole system
%p_ref = (p0*V0./(V0-x)-p0)*kappa+p0 = eps0*V^2/(2*(V0-x)^2)+p_mech

% here we go the simplified route and just check if we get the desired
% equilibrium

eps0 = 8.85419E-12;
s = 1e-3;
A = s*1;
V = 25000;
D = 1e-3-S1_smooth(end,2);

F_el = eps0*A*V^2/(2*D^2);
p_el = F_el/A;
T = 2e-5;
t_eval = S1_pres(end,1);

p_mech = 25000*(1-exp(-t_eval/T));

p_exc = p_mech+p_el;

disp(['Excitation pressure: ' num2str(p_exc/1000) ' kPa'])
disp(['Fluid pressure: ' num2str(S1_pres(end,2)/1000) ' kPa'])





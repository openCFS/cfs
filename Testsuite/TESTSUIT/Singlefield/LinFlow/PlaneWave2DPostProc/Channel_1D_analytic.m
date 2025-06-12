%% 1D solution of the LinFlow equations with non-zero viscosity

oversamp = 100; % oversample to get a finer resolution in time and space
f = 1000;
omega = 2*pi*f;
step = 28;
dt = 5e-5/oversamp;
h = 1e-0;
n_x = 100*oversamp;
dx = h/n_x;
x = 0:dx:h;

% mat. param.
mu = 10;
lambda = 5;
b_c = 2e-3; % width of the channel, used for force calculation

rho0 = 1.225;
c0 = 343;

plot_bit=0;


v_t_end = zeros(length(x),1);
for ii=1:length(x)
    v_t_end(ii) = vx(step*oversamp*dt,x(ii),omega,mu,lambda,rho0,c0);
end

v_x_start = zeros(step,1);
for ii=1:step*oversamp
    v_x_start(ii) = vx(dt*ii,0,omega,mu,lambda,rho0,c0);
end

v_t_end_diff = diff(real(v_t_end))/dx;

dvx_x_start = zeros(step*oversamp,1);
for ii=1:step*oversamp
    v_x1 = vx(dt*ii,0,omega,mu,lambda,rho0,c0);
    v_x2 = vx(dt*ii,dx,omega,mu,lambda,rho0,c0);
    dvx = (v_x2-v_x1)/dx;
    dvx_x_start(ii) = dvx;
end
p_a_ = cumtrapz(dt,-real(dvx_x_start)*rho0*c0^2);
p_a_offset = (max(p_a_)+min(p_a_))/2;
p_a = p_a_-p_a_offset;


if plot_bit==1
    figure
    plot(x,real(v_t_end))

    figure
    plot(x(1:end-1),v_t_end_diff)
end

%% calc stresses for x=0

strain = v_t_end_diff(1);
disp(['Strain at x = 0 m and step ' num2str(step) ': ' num2str(strain)])

stress = 2*mu*v_t_end_diff(1);
disp(['Stress at x = 0 m and step ' num2str(step) ': ' num2str(stress)])

comp_stress = (lambda-2/3*mu)*v_t_end_diff(1);
disp(['Compressible stress at x = 0 m and step ' num2str(step) ': ' num2str(comp_stress)])

visc_stress = stress+comp_stress;
disp(['Viscous stress at x = 0 m and step ' num2str(step) ': ' num2str(visc_stress)])

p_acou = p_a(end);
disp(['Pressure at x = 0 m and step ' num2str(step) ': ' num2str(p_acou)])

total_stress = -p_acou+visc_stress;
disp(['Total stress at x = 0 m and step ' num2str(step) ': ' num2str(total_stress)])

total_normal_stress = -total_stress; % note: this result differs slightly due to extremely steep gradient (change between 2 steps)
disp(['Total normal stress at x = 0 m and step ' num2str(step) ': ' num2str(total_normal_stress)])

force = total_normal_stress*b_c; % note: this result differs slightly due to extremely steep gradient (change between 2 steps)
disp(['Total normal force at x = 0 m and step ' num2str(step) ': ' num2str(force)])



%% analytical solution

function v=vx(t,x,omega,mu,lambda,rho0,c0)
    v_hat = -1i;

    k = 1i*sqrt(rho0*omega^2/(rho0*c0^2+(4/3*mu+lambda)*1i*omega));
    v = exp(1i*omega*t)*v_hat*exp(-k*x);
end

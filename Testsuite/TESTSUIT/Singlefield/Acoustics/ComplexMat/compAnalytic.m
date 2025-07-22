% example in book of Allard, 23ff)
d = 10e-2;
sigma = 10e3;

rho0 = 1.2;
K0 = 138720.0;
c0 = sqrt(K0/rho0);

f = 100:10:4e3;
omega = 2*pi*f;
X = (rho0*f)/sigma;
%impedance of porous material
Z2 = rho0*c0*( 1 + 0.057*X.^(-0.754) - i*0.087*X.^(-0.732) );
%complex wave number of porous material
k2 = (omega/c0).*( 1 + 0.0978*X.^(-0.7) - i*0.189*X.^(-0.595) );

%surface impedance
Zs = -i*Z2.*cot(k2.*d);

%effective parameters
K2   = Z2.*omega./k2;
rho2 = Z2.*k2./omega;

%refelction coefficient
Z0 = rho0*c0;
k0 = omega/c0;
R  = (Zs - Z0)./(Zs+Z0);

%absorption coefficient
alpha = 1 - (abs(R)).^2;

%transmission loss
TC = (2*i*Z2./(Z0 - 2*i*Z2.*cot(k2*d) + Z2.^2/Z0) ).*(sin(k2*d).^(-1)).*exp(-i*k0*d);
TL = -20*log(abs(TC));

plot(f,real(Zs),f,imag(Zs)),grid
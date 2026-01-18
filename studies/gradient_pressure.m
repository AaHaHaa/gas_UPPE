close all; clearvars;

Pin = 1; % input pressure (must be >Pout here)
Pout = 0; % output pressure

e = (Pin-Pout)/Pin;

L = 100; % fiber length
z = linspace(0,L,10000)'; % length sampling points
xi = z/L;

g = (1 - sqrt(1-2*e*(1-e/2)*xi))/e;

P = Pin - g*(Pin-Pout);

figure;
plot(z,P)

%% incompressible assumption in the literature
Psqrt = sqrt( Pin^2 + z/L*(Pout^2 - Pin^2) );
hold on;
plot(z,Psqrt);
hold off;

%%
h = 203.6;
R = 15e-6;
beta = R/L;
mu = 1.78e-5; % viscosity of nitrogen
M = 14; % atomic mass of nitrogen
rho0 = Pin*M/1.38e-23/(273.15+25);
kappa = R^3*rho0*(Pin-Pout)/L/mu^2;
f = @(xi,g) -2*h/R^2*(1-e*g)./(beta*e*h^2 - kappa*(1-e*g).^2);
[zz,gg] = ode45(f,[0,L],0);

P = Pin - gg*(Pin-Pout);

hold on;
plot(zz,P);
hold off;
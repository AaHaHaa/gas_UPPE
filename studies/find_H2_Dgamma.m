% This code finds the derivative of polarizability anisotropy at the
% equilibrium position, d(gamma)dR.
%
% The values are taken from "Polarizability of the Hydrogen Molecule" by W.
% Kolos, and L. Wolniewicz (2004)

close all; clearvars;

% 1 a.u. (length) = 5.29177210903(80)e-11 m
R_eq = 1.3984; % a.u.

R = [0.4,0.6,0.8,1,1.2,1.35,1.4,1.45,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4]';
% 1 a.u. (polarizability) = 1.64878e-41 Fm^2
alpha = [1.87653,2.35170,2.92314,3.58995,4.34482,4.96334,5.17862,5.39808,6.07864,7.02546,7.99537,8.95233,9.85682,10.66489,11.33265,11.82627,12.12397,12.22809,12.16123,11.95988,11.66968]';
gamma = [0.07862,0.20747,0.42182,0.74680,1.20255,1.63864,1.80281,1.97632,2.55312,3.44247,4.45353,5.52798,6.59998,7.58741,8.38281,8.92021,9.13487,9.01594,8.59397,7.93447,7.12298]';

% Change units
R_eq = R_eq*5.2917721090380e-11; % m
R = R*5.2917721090380e-11; % m
alpha = alpha*1.64878e-41; % Fm^2
gamma = gamma*1.64878e-41; % Fm^2

figure;
plot(R,[alpha,gamma]);
legend('\alpha','\gamma');
xlabel('R'); ylabel('polarizability');

R_range = 1*5.2917721090380e-11; % this can't be too small; otherwise, the number of data points is too small
around_eq = R<R_eq+R_range & R>R_eq-R_range;
n = 6;
[p_alpha,S_alpha,mu_alpha] = polyfit(R(around_eq),alpha(around_eq),n);
[p_gamma,S_gamma,mu_gamma] = polyfit(R(around_eq),gamma(around_eq),n);

m_H2 = 1/6.02e23/1e3; % kg, the reduced mass of H2
Dalpha = polyval(polyder(p_alpha)/mu_alpha(2),R_eq,[],mu_alpha)/sqrt(m_H2);
Dgamma = polyval(polyder(p_gamma)/mu_alpha(2),R_eq,[],mu_gamma)/sqrt(m_H2);
fprintf('Dalpha=%6.4g (Fm/sqrt(kg))\n',Dalpha); % Fm/sqrt(kg)
fprintf('Dgamma=%6.4g (Fm/sqrt(kg))\n',Dgamma); % Fm/sqrt(kg)

fit_alpha = polyval(p_alpha,R,[],mu_alpha);
fit_gamma = polyval(p_gamma,R,[],mu_gamma);

hold on;
plot(R,[fit_alpha,fit_gamma]);
xlim([min(R(around_eq)),max(R(around_eq))]);

%% Verification
permittivity0 = 8.85418782e-12; % F/m
fprintf('\nThe following value matches with the one in Table I in\n"Absolute measurement of the ultrafast nonlinear electronic and rovibrational response in H2 and D2"\nby Wahlstrand et al. (2015)\n');
fprintf(':Dalpha=%6.4g\n',polyval(polyder(p_alpha),R_eq,[],mu_alpha)/mu_alpha(2)/(4*pi*permittivity0)*1e20); % 1e-16 cm^2
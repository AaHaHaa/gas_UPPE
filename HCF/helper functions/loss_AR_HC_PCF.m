function confinement_loss = loss_AR_HC_PCF(wavelength,...
                                           num_tubes,r_tube,t_tube,delta_tube)
%LOSS_AR_HC_PCF Summary of this function goes here
%   
% References:
% 1. Vincetti et al., "Empirical formulas for calculating loss in hollow
%    core tube lattice fibers," Opt. Express 24(10), 10313-10325 (2016)
% 2. Vincetti et al., "A simple analytical model for confinement loss
%    estimation in hollow-core Tube Lattice Fibers," Opt. Express 27(4),
%    5230-5237 (2019)

k = 1 + delta_tube/(2*r_tube);
r_core = (k/sin(pi/num_tubes) - 1)*r_tube;

n_silica = calc_n_silica(wavelength*1e9,false);
n_silica = real(n_silica);

normalized_confinement_loss = 3e-4; % dB; from [2], but [1] employs 5e-4.
min_confinement_loss = real(normalized_confinement_loss/(r_core/r_tube)^4*(wavelength/r_tube).^4.5/(1 - t_tube/r_tube)^12.*sqrt(n_silica.^2-1)/t_tube.*exp(2*wavelength/r_tube./(real(n_silica).^2-1))); % dB/m

normf = 2*t_tube./wavelength.*sqrt(n_silica.^2 - 1);

u = 0:10;

A = 2e3*exp(-0.05*abs(u-1).^2.6);
tr = t_tube/r_tube;

p = 0;
for v = 1:20
    p = p + real(sum(A.*(Lorentzian(normf-cutoff_normf_HE(u,v,tr,n_silica)) + Lorentzian(normf-cutoff_normf_EH(u,v,tr,n_silica))),2));
end
confinement_loss = min_confinement_loss.*p; % dB/m
confinement_loss = (confinement_loss/10)*log(10)/2;
confinement_loss(n_silica<1.1) = 1e3; % just a large value; the ARHCF has high loss as n_silica is close to air's index; BTW, the loss of silica is high in this regime

end

function cutoff_normf = cutoff_normf_HE(u,v,tr,n)
%
% tr: t_tube/r_tube

if v == 1
    cutoff_normf = abs(0.21 + 0.175*u - 0.1./(u-0.35).^2).*tr.^(0.55+5e-3*sqrt(n.^4-1)) + 0.04*sqrt(u)*tr;
else
    cutoff_normf = 0.3./n.^0.3*(2/v)^1.2*abs(u-0.8).*tr + v - 1;
end

end

function cutoff_normf = cutoff_normf_EH(u,v,tr,n)
%
% tr: t_tube/r_tube

if v == 1
    cutoff_normf = (0.73 + 0.57*(u.^0.8+1.5)/4 - 0.04./(u-0.35)).*tr.^(0.5 - (n-1)/10./(u+0.5).^0.1);
else
    cutoff_normf = 11.5/v^1.2/(7.75-v)*(0.34+u/4.*(n/1.2).^1.15)./(u+0.2./n).^0.15.*tr.^(0.75+0.06./n.^1.15+0.1*sqrt(1.44./n)*(v-2)) + v - 1;
end

end

function output = Lorentzian(normf)

gamma = 0.003;
output = gamma^2./(gamma^2 + normf.^2);

end
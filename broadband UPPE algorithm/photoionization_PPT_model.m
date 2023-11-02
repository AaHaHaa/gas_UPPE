function [ne,DneDt] = photoionization_PPT_model(A_t, inverse_Aeff, ionization_energy, f0, dt, Ng,...
                                                erfi_x, erfi_y)
%PHOTOIONIZATION_PPT_MODEL This code computes the ionized free electrons
%with the Perelomov-Popov-Terent'ev (PPT) photoionization model. It's used
%as a helper function in the gas-UPPE model.
%   
% Input:
%
%   A_w: (Nt,1); the electric field under the frequency domain (in ifft order)
%   omega_pulse: (Nt,1); the instantaneous angular frequency of the field (2*pi*Hz=rad)
%   inverse_Aeff: scalar; 1/Aeff=SR value (Aeff: mode-field area) (1/m^2)
%   ionization_energy: scalar (J)
%   dt: scalar (ps)
%   Ng: the gas number density (1/m^3)
%   erfi_x, erfi_y: (Nt,1); the lookup table for the imaginary error function, erfi()
%
% Output:
%   ne: (Nt,1); the free electron number density (1/m^3)
%   DneDt: (Nt,1); the time derivative of ne (1/m^3/s)

[Nt,num_modes] = size(A_t);
if num_modes ~= 1
    error('Photoionization model works only for single mode.');
end

% Find instantaneous frequency of the pulse
pulse_phase = unwrap(angle(A_t));
pulse_phase = conv(pulse_phase,ones(floor(Nt/100),1)/floor(Nt/100),'same'); % smoothing is required; otherwise, it's too noisy such that an erroneous high-frequency number is calculated
omega_pulse = -(pulse_phase(3:end)-pulse_phase(1:end-2))/(2*dt)+2*pi*f0; % THz; I use "central difference" to calculate the slope here
omega_pulse = [omega_pulse(1);omega_pulse;omega_pulse(end)]*1e12; % Hz

me = 9.1093837e-31; % kg
e = 1.60217663e-19; % Coulomb
permittivity0 = 8.85418782e-12; % m^(-3)/kg*s^4*A^2
c = 299792458; % m/s
k = 4*pi*permittivity0;
h = 6.62607015e-34; % m^2*kg/s
hbar = h/2/pi;
a0 = k*hbar^2/me/e^2; % Bohr radius
U_H = e^2/k/a0/2; % hydrogen ionization energy = 13.6 eV

n = sqrt(U_H/ionization_energy); % effective principal quantum number

I = abs(A_t).^2*inverse_Aeff; % intensity; W/m^2

% This modification is for I too small to avoid spurious result in 1./I computation
I(I<max(I)/1e5) = max(I)/1e5;

ponderomotive_energy = e^2/2/me/permittivity0/c*I./omega_pulse.^2;
Keldysh_parameter = sqrt(ionization_energy/2./ponderomotive_energy);

kappa = 4*ionization_energy*sqrt(2*me*ionization_energy)/hbar/e;

v = ionization_energy/hbar./omega_pulse.*(1+1/2./Keldysh_parameter.^2);
n_v = ceil(v)-v + (0:10); % the minimum positive number of n-v+S, where S, an integer, adds n-v until it becomes a positive number

beta = 2*Keldysh_parameter./sqrt(1+Keldysh_parameter.^2);
g = 3/2./Keldysh_parameter.*((1+1/2./Keldysh_parameter.^2).*asinh(Keldysh_parameter) - sqrt(1+Keldysh_parameter.^2)/2./Keldysh_parameter);

% the PPT correction factor
A0 = sum(2/sqrt(3)*Keldysh_parameter.^2./(1+Keldysh_parameter.^2).*exp(-2*n_v.*asinh(Keldysh_parameter)).*interp1(erfi_x,erfi_y,sqrt(beta.*n_v)),2);
A0(Keldysh_parameter<0.8) = 1; % A0 should be close to 1 at small Keldysh parameter

W = 2^(2*n)/n/gamma(n)/gamma(n+1)*ionization_energy/hbar*sqrt(6/pi)*A0.*(sqrt(permittivity0*c/2./I)*kappa./sqrt(1+Keldysh_parameter.^2)).^(2*n-1.5).*exp(-sqrt(permittivity0*c/2./I)*kappa/3.*g);
W(I<max(I)/1e4) = 0; % avoid non-negligible W when there is no or weak field due to numerical precision error

ne = Ng*cumsum(W)*(dt*1e-12);
% Below is the more accurate version with integrating factor.
% It's important only when relative_ne approaches the gas number density, Ng.
%integrating_factor = exp(cumsum(W)*(dt*1e-12));
%ne = Ng./integrating_factor.*cumsum(W.*integrating_factor)*(dt*1e-12);
DneDt = W.*(Ng - ne);

end
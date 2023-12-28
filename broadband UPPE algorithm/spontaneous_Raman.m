function sponRS = spontaneous_Raman(Nt,dt,sim,fiber,gas,gas_eqn)
%SPONTANEOUS_RAMAN It computes the spontaneous Raman term.
%SPONTANEOUS_RAMAN It computes the spontaneous Raman term.
% Spontaneous Raman scattering is equivalent to scattering with a field
% with one photon per frequency band.
%
%   It creates the spontaneous-Raman counterpart of SR*|A|^2, A: the pulse field
%
%   I consider one-photon noise per frequency band.
%   spectral noise field = counterpart of sqrt( SR*|A|^2 )
%                        = sponRS_prefactor{1}.*randn(size(sponRS_prefactor{1})).*exp(1i*2*pi*rand(size(sponRS_prefactor{1})))
%   noise temporal intensity = abs( fft(spectral noise field) ).^2
%   Transform into the spectral domain for convolution = ifft( noise temporal intensity ).*sponRS_prefactor{2}
%   sponRS_prefactor{2} modifies the spectral noise according to the Bose-Einstein distribution and Stokes generation.
%   sponRS_Gamma = fft(haw.*sponRS) finishes the spontaneous-Raman convolution
%       such that
%   sponRS_Gamma*A is the spontaneous-Raman field in the GMMNLSE/UPPE.
%
% -------------------------------------------------------------------------
%   They're summarized and implemented as the following the stepping functions:
%
%      sponRS = ifft(abs(fft(sponRS_prefactor{1}.*randn(size(sponRS_prefactor{1})).*exp(1i*2*pi*rand(size(sponRS_prefactor{1}))))).^2).*sponRS_prefactor{2};
%      sponRS_Gamma = fft(haw.*sponRS);

h = 6.62607015e-34; % J*s
hbar = h/(2*pi); % J*s
k = 1.38064852e-23; % Boltzmann constant (MKS unit)

time_window = Nt*dt; % ps
f = ifftshift((-Nt/2:Nt/2-1)'/time_window,1); % THz
real_omegas = 2*pi*(f+sim.f0)*1e12; % Hz

% I use analytical-signal representation for solving GMMNLSE, so the field
% covers only the positive frequencies.
real_omegas(real_omegas<0) = 0; % no noise at negative frequencies

% The quantum noise is uniformly distributed over the cross section
num_spatial_modes = size(fiber.SR,1);
if sim.scalar
    num_modes = num_spatial_modes;
else
    num_modes = num_spatial_modes*2;
end
SR = zeros(1,num_modes);
for midx = 1:num_modes
    spatial_midx = ceil(int32(midx)/2);
    SR(midx) = fiber.SR(spatial_midx,spatial_midx,spatial_midx,spatial_midx);
end

% phonon number based on Bose-Einstein distribution
nth = 1./(exp(h*abs(f*1e12)./k/gas.temperature)-1);
nth(isinf(nth)) = 0; % if f=0
Heaviside = double(f<0); % Stokes wave

sponRS = {sqrt(SR*hbar.*real_omegas/(time_window*1e-12)),...
          nth+Heaviside};

% In gas simulations, the frequency window is larger than the input to avoid aliasing due to the large Raman frequency shift.
sponRS = {cat(1,sponRS{1}(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,     sponRS{1}(gas_eqn.n+1:end,:)),...
          cat(1,sponRS{2}(1:gas_eqn.n,:),gas_eqn.upsampling_zeros(:,1),sponRS{2}(gas_eqn.n+1:end,:))};
end
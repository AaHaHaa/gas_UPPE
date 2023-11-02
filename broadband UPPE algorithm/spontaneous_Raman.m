function sponRS = spontaneous_Raman(Nt,dt,sim,fiber,gas,gas_eqn)
%SPONTANEOUS_RAMAN It computes the spontaneous Raman term.

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
QR = zeros(1,num_modes);
for midx = 1:num_modes
    spatial_midx = ceil(int32(midx)/2);
    QR(midx) = fiber.SR(spatial_midx,spatial_midx,spatial_midx,spatial_midx)./sim.mode_profiles.norms(1).^4;
end

% phonon number based on Bose-Einstein distribution
nth = 1./(exp(h*abs(f*1e12)./k/gas.temperature)-1);
nth(isinf(nth)) = 0; % if f=0
Heaviside = double(f<0); % Stokes wave

sponRS = {sqrt(hbar*real_omegas.*QR/(time_window*1e-12)),...
          nth+Heaviside};

% In gas simulations, the frequency window is larger than the input to avoid aliasing due to the large Raman frequency shift.
sponRS = {cat(1,sponRS{1}(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,     sponRS{1}(gas_eqn.n+1:end,:)),...
          cat(1,sponRS{2}(1:gas_eqn.n,:),gas_eqn.upsampling_zeros(:,1),sponRS{2}(gas_eqn.n+1:end,:))};
end
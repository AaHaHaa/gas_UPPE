function sponRS = spontaneous_Raman(Nt,dt,sim,gas,gas_eqn)
%SPONTANEOUS_RAMAN It computes the spontaneous Raman term.
% Spontaneous Raman scattering is equivalent to scattering with a field
% with one photon per frequency band.

h = 6.62607015e-34; % J*s
hbar = h/(2*pi); % J*s
k = 1.38064852e-23; % Boltzmann constant (MKS unit)

time_window = Nt*dt; % ps
f = ifftshift((-Nt/2:Nt/2-1)'/time_window,1); % THz
real_omegas = 2*pi*(f+sim.f0)*1e12; % Hz

% I use analytical-signal representation for solving GMMNLSE, so the field
% covers only the positive frequencies.
real_omegas(real_omegas<0) = 0; % no noise at negative frequencies

% phonon number based on Bose-Einstein distribution
if gas.model == 0 % "gas.model=0" uses acyclic convolution which has an extended time window, so f needs to be re-defined.
    ratio = gas_eqn.acyclic_conv_stretch(Nt)/Nt;
    f = (-floor(gas_eqn.acyclic_conv_stretch(Nt)/2):ceil(gas_eqn.acyclic_conv_stretch(Nt)/2)-1)'/time_window/ratio; % THz
else
    f = fftshift(f,1);
end
nth = 1./(exp(h*abs(f*1e12)./k/gas.temperature)-1);
nth(isinf(nth)) = 0; % if f=0
Heaviside = double(f<0); % Stokes wave

sponRS = {sqrt(hbar.*real_omegas/(time_window*1e-12)),...
          nth+Heaviside};

% In gas simulations, the frequency window is larger than the input to avoid aliasing due to the large Raman frequency shift.
sponRS{1} = cat(1,sponRS{1}(1:gas_eqn.n,:),...
                  gas_eqn.upsampling_zeros(:,1),...
                  sponRS{1}(gas_eqn.n+1:end,:)...
               );
end
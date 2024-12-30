function At_noise = shot_noise(gas_Nt,gas_dt,sim,gas,gas_eqn,num_modes)
%SHOT_NOISE It computes the shot noise included in the governing equation

h = 6.62607015e-34; % J*s
k = 1.38064852e-23; % Boltzmann constant (MKS unit)

time_window = gas_Nt*gas_dt; % ps
f = ifftshift((-gas_Nt/2:gas_Nt/2-1)'/time_window,1); % THz
real_f = (f+sim.f0)*1e12; % Hz

% I use analytical-signal representation for solving GMMNLSE, so the field
% covers only the positive frequencies.
real_f(real_f<0) = 0; % no noise at negative frequencies

% phonon number based on Bose-Einstein distribution
if gas.model == 0 % "gas.model=0" uses acyclic convolution which has an extended time window, so f needs to be re-defined.
    ratio = gas_eqn.acyclic_conv_stretch(gas_Nt)/gas_Nt;
    f = (-floor(gas_eqn.acyclic_conv_stretch(gas_Nt)/2):ceil(gas_eqn.acyclic_conv_stretch(gas_Nt)/2)-1)'/time_window/ratio; % THz
else
    f = fftshift(f,1);
end

noise_amplitude = sqrt(h*real_f/(time_window*1e-12));

if isequal(sim.step_method,'RK4IP')
    At_noise = fft(noise_amplitude.*randn(gas_Nt,num_modes).*exp(1i*2*pi*rand(gas_Nt,num_modes)),[],1);
else % 'MPA'
    At_noise = repmat(fft(noise_amplitude.*randn(gas_Nt,1,num_modes).*exp(1i*2*pi*rand(gas_Nt,1,num_modes)),[],1),1,sim.MPA.M+1,1);
end

end
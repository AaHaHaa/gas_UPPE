function At_noise = shot_noise(gas_Nt,gas_dt,sim,num_modes)
%SHOT_NOISE It computes the shot noise included in the governing equation

h = 6.62607015e-34; % J*s

time_window = gas_Nt*gas_dt; % ps
f = ifftshift((-floor(gas_Nt/2):floor((gas_Nt-1)/2))'/time_window,1); % THz
real_f = (f+sim.f0)*1e12; % Hz

% I use analytical-signal representation for solving UPPE, so the field
% covers only the positive frequencies.
real_f(real_f<0) = 0; % no noise at negative frequencies

noise_amplitude = sqrt(h*real_f/(time_window*1e-12));

if isequal(sim.step_method,'RK4IP')
    At_noise = fft(noise_amplitude.*randn(gas_Nt,num_modes).*exp(1i*2*pi*rand(gas_Nt,num_modes)),[],1);
else % 'MPA'
    At_noise = repmat(fft(noise_amplitude.*randn(gas_Nt,1,num_modes).*exp(1i*2*pi*rand(gas_Nt,1,num_modes)),[],1),1,sim.MPA.M+1,1);
end

end
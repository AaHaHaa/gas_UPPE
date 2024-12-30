function sim_betas = calc_sim_betas(fiber,dt,omegas,fields)
%CALC_SIM_BETAS It computes the sim.betas used in UPPE, where sim.betas(1)
%is a free parameter to facilitate simulations and sim.betas(2) is the
%inverse velocity of the moving frame

% Obtain the betas of the input pulse
fftshift_omegas = fftshift(omegas,1);
spectrum = sum(abs(fftshift(ifft(fields,[],1),1)).^2,2);
omega0 = sum(fftshift_omegas.*spectrum)/sum(spectrum); % 2*pi*THz; the pulse center frequency (under shifted omega)
omega_range = 1/dt; % 2*pi*THz
omegas_idx_near_pulse = fftshift_omegas>omega0-omega_range/5 & fftshift_omegas<omega0+omega_range/5;% pick only the data near the pulse center frequency to find its beta0 and beta1
clear spectrum omega0 omega_range;

fit_order = max(2,min(7,sum(omegas_idx_near_pulse)-1)); % 2~7
[betas_Taylor_coeff,~,mu] = polyfit(fftshift_omegas(omegas_idx_near_pulse),real(fiber.betas(omegas_idx_near_pulse,1)),fit_order);
sim_betas = [betas_Taylor_coeff(end);betas_Taylor_coeff(end-1)];
sim_betas = [sim_betas(1)-sim_betas(2)*mu(1)/mu(2);...
             sim_betas(2)/mu(2)];

end


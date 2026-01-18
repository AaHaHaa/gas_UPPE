function beta = recompute_beta_gradient_pressure(eta,gas)
%RECOMPUTE_BETA_GRADIENT_PRESSURE It calculates the propagation constant
%with gas density eta (with the unit "amagat").
%
%   gas.permittivity_r is the gas permittivity with 1 amagat
%
%   gas.coeff_beta_with_eta = gas.k0.^2./gas.beta_mean_gas is a pre-calculated coefficient in calculating the propagation constant.

%% Refractive index
n_gas = sqrt((gas.permittivity_r - 1)*eta + 1); % refractive index of the gas

% pressure-induced absorption
n_gas = n_gas + 1i*gas.Raman_absorption*eta^2;

% beta = beta_mean_gas + gas.k0.^2./gas.beta_mean_gas.*(n_gas-1)
beta = gas.beta_mean_gas + gas.k0.*(n_gas-gas.n_gas_mean); % beta is propagation constant at the desired gas pressure

end
function beta = recompute_beta_gradient_pressure(wavelength,eta,gas)
%RECOMPUTE_BETA_GRADIENT_PRESSURE It calculates the propagation constant
%with gas density eta (with the unit "amagat").
%
%   gas.permittivity_r is the gas permittivity with 1 amagat
%   gas.permittiivty_r_113, gas.permittiivty_r_120 and gas.permittivity_r_400 are special values to be used to remove singularity in the Sellmeier equation.
%
%   gas.coeff_beta_with_eta = gas.k0.^2./gas.beta_no_gas is a pre-calculated coefficient in calculating the propagation constant.

%% Refractive index
switch gas.material
    case 'H2'
        n_gas = sqrt((gas.permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        n_gas_120 = sqrt((gas.permittivity_r_120 - 1)*eta + 1); %  % Sellmeier is valid only above 164nm
        n_gas(wavelength<120e-9) = n_gas_120;
        
        % pressure-induced absorption
        n_gas = n_gas + 1i*gas.Raman_absorption*eta^2;
    case 'O2'
        n_gas = sqrt((gas.permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        n_gas_400 = sqrt((gas.permittivity_r_400 - 1)*eta + 1); %  % Sellmeier is valid only above 400nm
        n_gas(wavelength<400e-9) = n_gas_400;
        
        % pressure-induced absorption
        n_gas = n_gas + 1i*gas.Raman_absorption*eta^2;
    case {'air','N2'}
        n_gas = sqrt((gas.permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        % pressure-induced absorption
        n_gas = n_gas + 1i*gas.Raman_absorption*eta^2;
    case {'Ar','Ne','He'}
        n_gas = sqrt((gas.permittivity_r - 1)*eta + 1); % refractive index of the gas
    case {'Xe','Kr'}
        n_gas = sqrt((gas.permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        n_gas_113 = sqrt((gas.permittivity_r_113 - 1)*eta + 1); %  % Sellmeier is valid only above 113nm
        n_gas(wavelength<113e-9) = n_gas_113;
    case 'CH4'
        n_gas = sqrt((gas.permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        % Avoid the singularity at resonances
        idx_resonance = n_gas < 1;
        n_gas(idx_resonance) = 1;
end

% beta = beta_no_gas + gas.k0.^2./gas.beta_no_gas.*(n_gas-1)
beta = gas.beta_no_gas + gas.coeff_beta_with_eta.*(n_gas-1); % beta is propagation constant at the desired gas pressure

end
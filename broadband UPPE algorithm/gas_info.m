function [fiber,sim,gas] = gas_info(fiber,sim,gas,wavelength)
%GAS_INFO It loads parameters related to gases
% It will load either different functions based on constant or gradient
% pressure considered.

if ~sim.include_Raman && ~sim.include_photoionization
    sim.include_heating = false;
end
if sim.include_heating % consider Raman/photoionization heating during propagation
    if isfield(gas,'pressure')
        if sum(gas.pressure) ~= 0
            [fiber,sim,gas] = gas_info_constant_pressure_heating(fiber,sim,gas,wavelength);
        else % no gas = vacuum, no heating
             % Just run with the regular constant-pressure model
            gas.include_heating = false; % turn thermal model off
            [fiber,sim,gas] = gas_info_constant_pressure(fiber,sim,gas,wavelength);
        end
    else
        error('gas_info:includeHeatingError',...
              'Current thermal model supports only constant pressure; otherwise, gas flow that takes away heat complicates the problem.');
    end
else
    if isfield(gas,'pressure_in') && isfield(gas,'pressure_out')
        [fiber,sim,gas] = gas_info_gradient_pressure(fiber,sim,gas,wavelength);
    elseif isfield(gas,'pressure') % constant gas pressure
        [fiber,sim,gas] = gas_info_constant_pressure(fiber,sim,gas,wavelength);
    else
        error('gas_info:gasPressureError',...
              '"gas" needs to have "pressure_in" and "pressure_out" for a gradient pressure or "pressure" only for a constrant pressure.');
    end
end

% There will be some situations where we want to modify, for example, the 
% Raman response functions, so gas_info() will be called before pulse
% propagation. In this case, we don't want the propagation function to call
% gas_info() again and rewrite the desired parameters set beforehand.
gas.info_called = true; % gas_info() has been called

end
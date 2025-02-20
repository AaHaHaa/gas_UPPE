function nonlinear_phase = accumulated_nonlinear_phase( fiber,sim,gas_species,prop_output )
%ACCUMULATED_NONLINEAR_PHASE Calculates the total accumulated nonlinear
%phase from the propagation.
%   
% fiber:
%   SR: 1/Aeff, where Aeff: a scalar (m^2)
% sim:
%   f0: a scalar; central frequency (THz)
%   mode_profiles.norms: norm of the mode profiles
%   X3
% gas_species:
%   R.preR: Raman strength of rotational SRS
%   R.omega: Raman transition frequency of rotational SRS
%   V.preR: Raman strength of vibrational SRS
%   V.omega: Raman transition frequency of vibrational SRS
% prop_output:
%   dt: a scalar; the difference between each time point (ps)
%   field: (Nt,1,save_num); fields during propagations (sqrt(W))
%   z: (1,save_num); the z position of each saved field (m)

prop_output.fields = prop_output.fields(:,1,:);

Nt = size(prop_output.fields,1);
time_window = Nt*prop_output.dt;
Aeff = 1/fiber.SR(1);

%% Nonlinear constant at the pulse's center frequency
Omega = 2*pi*(linspace(-floor(Nt/2), floor((Nt-1)/2), Nt))'/time_window; % in 1/ps, in the order that the ifft gives
omega = (Omega + 2*pi*sim.f0)*1e12; % Hz
permittivity0 = 8.85418782e-12;
prefactor = {omega/4./sim.mode_profiles.norms(:,1).^4,... % I use the 1st mode for norm^4 computation.
             ...                                                % This can create certain amount of deviation for higher-order modes, but it should be acceptable; otherwise, the computation is too heavy.
             3*permittivity0*fiber.X3/4}; % for Kerr term

% Find the nonlinear constant at the pulse center frequency
spectrum = sum(abs(fftshift(ifft(prop_output.fields),1)).^2,2);
omega0 = sum(omega.*spectrum,1)./sum(spectrum,1); % 2*pi*THz; the pulse center frequency (under shifted omega)
idx0 = zeros(1,1,length(prop_output.z));
for zi = 1:length(prop_output.z)
    idx0(zi) = find(omega>omega0(zi),1);
end
prefactor{1} = permute(prefactor{1}(idx0),[2,3,1]); % (general) nonlinear constant
prefactor{2} = permute(prefactor{2}(idx0),[2,3,1]); % nonlinear prefactor for the Kerr term

%% Peak intensity
Ip = max(abs(prop_output.fields).^2/Aeff,[],1);

%% Calculate accumulated nonlinear phase
nonlinear_phase_Kerr = prefactor{1}.*prefactor{2}.*Ip;
nonlinear_phase_Raman = prefactor{1}.*Ip.*(sum(gas_species.R.preR./(gas_species.R.omega*1e12),2) + sum(gas_species.V.preR./(gas_species.V.omega*1e12),2));

nonlinear_phase = trapz(prop_output.z',squeeze(nonlinear_phase_Kerr + nonlinear_phase_Raman));

end
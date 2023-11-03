function delta_permittivity = calc_Raman_delta_permittivity(sim,gas,gas_eqn,A_w,Nt)
%CALC_RAMAN_DELTA_PERMITTIVITY Summary of this function goes here
%   Detailed explanation goes here

A_t_upsampling = fft(cat(1,A_w(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,A_w(gas_eqn.n+1:end,:)));

R_delta_permittivity = permute(gas_eqn.R_delta_permittivity,[1,3,2]); % make it [Nt,num_modes,Raman_type=2]; The Raman_type dimension are R and V
switch gas.model
    case 0
        delta_permittivity = fft(R_delta_permittivity.*ifft(abs(A_t_upsampling).^2, 2*gas_eqn.Nt-1,1));
        delta_permittivity = delta_permittivity(gas_eqn.R_downsampling,:,:);
    case 1
        delta_permittivity = fft(R_delta_permittivity.*ifft(abs(A_t_upsampling).^2));
end
delta_permittivity = ifft(delta_permittivity.*(permute(max(max(real(sim.mode_profiles.mode_profiles),[],1),[],2),[1,3,2])./mean(sim.mode_profiles.norms)).^2); % find the max delta_permittivity of each mode
% Only the real part corresponds to the actual permittiviy contribution of each Raman response.
% The imaginary part is retained so that it's easier to visualize the "intensity" of the phonon strength by taking abs().
delta_permittivity = fft(delta_permittivity([1:gas_eqn.n,gas_eqn.Nt-(Nt-gas_eqn.n-1):gas_eqn.Nt],:,:)); % transform back to time domain

end


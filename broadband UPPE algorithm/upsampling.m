function [D_op_upsampling,prefactor,mode_profiles_norms,...
          upsampling_zeros] = upsampling(sim,...
                                         Nt,gas_Nt,...
                                         num_modes,...
                                         D_op,prefactor)
%UPSAMPLING Upsampling to avoid frequency aliasing
%
% Aliasing mostly comes from Raman shift. However, when the spectrum
% becomes broad, aliasing can come from Kerr effect as well due to
% four-wave mixing.
% Instead of applying upsampling to the Raman computation only, it's
% important to apply it to Kerr as well, especially when running
% supercontinuum generation or when your frequency window isn't large
% enough.
if sim.gpu_yes
    upsampling_zeros = complex(zeros(gas_Nt-Nt, num_modes, 'gpuArray'));
else
    upsampling_zeros = complex(zeros(gas_Nt-Nt, num_modes));
end
n = ceil(Nt/2);
% Zero-paddding in frequency domain for upsampling temporally
% Zeros are added at the low- and high-frequency side which is the center
% of the array after discrete Fourier transform.
if ~isempty(D_op) % not computing with gradient pressure
    D_op_upsampling = cat(1,D_op(1:n,:),upsampling_zeros,D_op(n+1:end,:));
else
    D_op_upsampling = [];
end
mode_profiles_norms = cat(1,sim.mode_profiles.norms(1:n),upsampling_zeros(:,1),sim.mode_profiles.norms(n+1:end)); % currently only for computing the nonlinear index change to save
mode_profiles_norms(mode_profiles_norms==0) = Inf; % As in computing index change, this term is in the denominator, we don't this 1/0 to blow up but give zero due to the upsampling operation.
prefactor{1} = cat(1,prefactor{1}(1:n),upsampling_zeros(:,1),prefactor{1}(n+1:end));
prefactor{2} = cat(1,prefactor{2}(1:n,:,:),real(upsampling_zeros),prefactor{2}(n+1:end,:,:));

if sim.include_photoionization
    prefactor{3} = cat(1,prefactor{3}(1:n),upsampling_zeros(:,1),prefactor{3}(n+1:end));
    prefactor{4} = cat(1,prefactor{4}(1:n,:),real(upsampling_zeros),prefactor{4}(n+1:end,:));
end

end


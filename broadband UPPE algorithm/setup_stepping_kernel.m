function [gpuDevice_Device,...
          cuda_SRSK,num_operations_SRSK,...
          cuda_sponRS,num_operations_sponRS,...
          cuda_MPA_psi_update] = setup_stepping_kernel(sim,Nt,num_modes)
%SETUP_STEPPING_KERNEL It sets cuda for computing sums of SR and SK terms, 
%and spontaneous Raman terms.

%% Use the specified GPU
% This needs to run at the beginning; otherwise, the already-stored values
% in GPU will be unavailable in a new GPU if the GPU device is switched.
try
    gpuDevice_Device = gpuDevice(sim.gpuDevice.Index); % use the specified GPU device
catch
    error('Please set the GPU you''re going to use by setting "sim.gpuDevice.Index".');
end

%% SR, SK
% Polarization modes
if sim.scalar
    polar_str = '';
else % polarized fields
    polar_str = '_with_polarization';
end

% Nonlinear term
SRSK_filename = ['UPPE_nonlinear_sum' polar_str];

switch SRSK_filename
    case 'UPPE_nonlinear_sum'
        num_operations_SRSK = 2;
    case 'UPPE_nonlinear_sum_with_polarization'
        num_operations_SRSK = 3;
end

% The number of blocks is set based on the total number of threads
if isequal(sim.step_method,'MPA')
    M = sim.MPA.M+1;
else % RK4IP
    M = 1;
end

% Kernal for computing the nonlinear term of UPPE
cuda_SRSK = setup_kernel_SRSK(SRSK_filename,sim.cuda_dir_path,Nt,M,num_operations_SRSK,num_modes^2);

%% MPA
% Kernal for updating psi in MPA stepping algorithm
if isequal(sim.step_method,'MPA')
    cuda_MPA_psi_update = setup_kernel_MPA_coeff('MPA_psi_update',sim.cuda_dir_path,Nt,M,num_modes);
else
    cuda_MPA_psi_update = [];
end

%% Spontaneous Raman scattering
% Nonlinear term
sponRS_filename = ['UPPE_sponRS_sum' polar_str];

switch sponRS_filename
    case 'UPPE_sponRS_sum'
        num_operations_sponRS = 1;
    case 'UPPE_sponRS_sum_with_polarization'
        num_operations_sponRS = 2;
end

% Kernal for computing the spontaneous Raman term of UPPE
cuda_sponRS = setup_kernel_SRSK(sponRS_filename,sim.cuda_dir_path,Nt,M,num_operations_sponRS,num_modes^2);

end
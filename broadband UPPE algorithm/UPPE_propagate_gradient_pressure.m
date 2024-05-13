function foutput = UPPE_propagate_gradient_pressure(fiber, initial_condition, sim, gas)
%UPPE_PROPAGATE_GRADIENT_PRESSURE Propagate an initial multimode pulse through 
%an arbitrary distance of a hollow-core fiber with an adaptive-step method.
%The HCF is filled with a gradient gas pressure.
%
% -------------------------------------------------------------------------
%
%   "fiber" is a structure with the fields:
%
%       Basic properties -->
%
%           betas - a (?,nm) matrix; "nm" = num_modes if under scalar fields;
%                                    otherwise, "nm" can be both num_modes or 2*num_modes depending on whether there's birefringence.
%                   betas(i, :) = (i-1)th order dispersion coefficient for each mode, in ps^n/m
%
%           L0 - length of fiber, in m
%
% -------------------------------------------------------------------------
%
%   "initial_condition" is a structure with the fields:
%
%       dt - time step
%       fields - initial field, in W^1/2, (N-by-num_modes).
%                If the size is (N-by-num_modes-by-S), then it will take the last S.
%
%                num_modes = num_modes
%
% -------------------------------------------------------------------------
%
%   "sim" is a structure with the fields:
%
%       Basic settings -->
%
%           betas - the betas for the slowly varying approximation and the moving frame, 
%                   that is to say, fiber.betas([1 2],:) = fiber.betas([1 2],:) - sim.betas;
%                   (2,1) column vector;
%                   if not set, no "sim.betas", the simulation will be run relative to the first mode
%           f0 - center frequency, in THz
%           save_period - spatial period between saves, in m. 0 = only save input and output
%
%       Adaptive method -->
%
%           adaptive_deltaZ.threshold - a scalar;
%                                       the accuracy used to determined whether to increase or decrease the step size.
%           adaptive_deltaZ.max_deltaZ - a scalar; the maximum adaptive step size
%
%       Algorithms to use -->
%
%           gpu_yes - 1(true) = GPU, 0(false) = CPU
%                     Whether or not to use the GPU. Using the GPU is HIGHLY recommended, as a speedup of 50-100x should be possible.
%           Raman_model - 0 = ignore Raman effect
%                         1 = include gas Raman
%
%       Others -->
%
%           pulse_centering - 1(true) = center the pulse according to the time window, 0(false) = do not
%                             The time delay will be stored in time_delay after running UPPE_propagate().
%           num_photon_noise_per_bin - a scalar; include photon noise (typically one photon per spectral discretization bin)
%           parallel_yes - show the simulation progress under the parallel "parfor".
%                          The progress bar will be ignored within parfor.
%           parallel_idx - the index of the session
%           gpuDevice.Index - a scalar; the GPU to use
%           gpuDevice.Device - the output of MATLAB "gpuDevice(gpu_index)"
%           cuda_dir_path - path to the cuda directory into which ptx files will be compiled and stored
%           progress_bar - 1(true) = show progress bar, 0(false) = do not
%                          It'll slow down the code slightly. Turn it off for performance.
%           progress_bar_name - the name of the UPPE propagation shown on the progress bar.
%                               If not set (no "sim.progress_bar_name"), it uses a default empty string, ''.
%
% =========================================================================
%
% Output:
% foutput.fields - (N, num_modes, num_save_points) matrix with the multimode field at each save point
% foutput.dt - time grid point spacing, to fully identify the field

%% Check the precision of the input fields and make it with "double"-precision
% Consider only the last fields of the initial condition
initial_condition.fields = double(initial_condition.fields(:,:,end));
if isreal(initial_condition.fields)
    initial_condition.fields = complex(initial_condition.fields);
end

%% Pick the stepping algorithm to use depending on the number of modes
if length(sim.midx) == 1 % single mode
    sim.step_method = 'RK4IP';
else % multimode
    sim.step_method = 'MPA';
    sim.MPA.M = ceil(abs(sim.MPA.M)/2)*2; % force M to be even; if odd, make it a larger even number
    sim.MPA.coeff = permute(MPA_AM_coeff(sim.MPA.M),[1,3,2]);
end

%% Preset some parameters
dt = initial_condition.dt;

% Get the numerical parameters from the initial condition.
size_field = size(initial_condition.fields);
Nt = size_field(1);
num_modes = size_field(2);

num_spatial_modes = size(fiber.SR,1);

% During the nonlinear computation, the frequency window is extended to
% three times to prevent aliasing from Raman or FWM. It's shrinked back to
% the original size after the operation.
gas_Nt = Nt*3;
gas_dt = dt/3;

time_window = Nt*dt;

permittivity0 = 8.8541878176e-12; % F/m

%% Check the validity of input parameters
if sim.save_period == 0
    sim.save_period = fiber.L0;
end

num_saves_total = fiber.L0/sim.save_period;
if rem(num_saves_total,1) && rem(num_saves_total+eps(num_saves_total),1) && rem(num_saves_total-eps(num_saves_total),1)
    error('UPPE_propagate:SizeIncommemsurateError',...
          'The save period is %f m and the fiber length is %f m, which are not commensurate', sim.save_period, fiber.L0)
else
    num_saves_total = round(num_saves_total);
end

% Error check on the dimensions (num_modes) of matrices
check_nummodes(sim, fiber, initial_condition.fields);

if sim.Raman_model ~= 0
    Fmax = 1/2/dt;
    % Below, only the available Raman (whose preR~=0) are considered.
    % max_Omega and max_T2 are used to determine which Raman model to use (see below).
    switch gas.gas_material
        case 'CH4'
            max_Omega = max(gas.(gas.gas_material).V.omega.*(gas.(gas.gas_material).V.preR~=0));
            max_T2 = max(gas.(gas.gas_material).V.T2.*(gas.(gas.gas_material).V.preR~=0));
        case {'H2','N2','O2'}
            max_Omega = max([gas.(gas.gas_material).R.omega.*(gas.(gas.gas_material).R.preR~=0),...
                             gas.(gas.gas_material).V.omega.*(gas.(gas.gas_material).V.preR~=0)]);
            max_T2 = max([gas.(gas.gas_material).R.T2.*(gas.(gas.gas_material).R.preR~=0),...
                          gas.(gas.gas_material).V.T2.*(gas.(gas.gas_material).V.preR~=0)]);
        case 'air'
            max_Omega = max([gas.N2.R.omega.*(gas.N2.R.preR~=0),...
                             gas.N2.V.omega.*(gas.N2.V.preR~=0),...
                             gas.O2.R.omega.*(gas.O2.R.preR~=0),...
                             gas.O2.V.omega.*(gas.O2.V.preR~=0)]);
            max_T2 = max([gas.N2.R.T2.*(gas.N2.R.preR~=0),...
                          gas.N2.V.T2.*(gas.N2.V.preR~=0),...
                          gas.O2.R.T2.*(gas.O2.R.preR~=0),...
                          gas.O2.V.T2.*(gas.O2.V.preR~=0)]);
        otherwise
            max_Omega = 0;
            max_T2 = 0;
    end
    
    % Check the frequency window
    % It needs to cover 2*f_R to include at least the first-order Stokes 
    % and anti-Stokes Raman scattering.
    if max_Omega/2/pi > Fmax
        error('UPPE_propagate:FrequencyWindowError',...
              ['The frequency window is too small.\n',...
               'It needs to cover 2*f_R to include at least the first-order Stokes and anti-Stokes Raman scattering.']);
    end
    
    % Select the computational model to use
    % 0: Do acyclic convolution for the Raman term; this requires expansion of the time window to prevent aliasing
    % 1: Do cyclic convolution for the Raman term
    %
    % The acyclic-convolution model is implemented because Raman has very 
    % long dephasing time (nanoseconds) where most physics under 
    % investigation happens around fs or ps time scale. Using a large time 
    % window in simulations is a waste of time by 1000 times.
    if time_window < max_T2*10
        gas.model = 0;
    else
        gas.model = 1;
    end
end

%% Polarized propagation constant (betas) and X3
% For polarized fields, the dimension of the input betas is allowed to be
% "num_spatial_modes", so it needs to be expanded into
% "2*num_spatial_modes" in the rest of the computation.
[sim,fiber] = expansion_including_polarization_modes(sim,fiber,num_modes);

%% Set up the GPU details
if sim.gpu_yes
    [sim.gpuDevice.Device,...
     sim.cuda_SRSK,  sim.cuda_num_operations_SRSK,...
     sim.cuda_sponRS,sim.cuda_num_operations_sponRS,...
     sim.cuda_MPA_psi_update] = setup_stepping_kernel(sim,gas_Nt,num_modes);
end

%% Pre-calculate sim.betas for the dispersion term
% The "omegas" here is actually (omega - omega0), omega: true angular frequency
%                                                 omega0: central angular frequency (=2*pi*f0)
if sim.gpu_yes
    dt = gpuArray(dt);
    time_window = gpuArray(time_window);
end
omegas = 2*pi*ifftshift(linspace(-floor(Nt/2), floor((Nt-1)/2), Nt))'/time_window; % in 1/ps, in the order that the fft gives
real_omegas = (omegas + 2*pi*sim.f0)*1e12; % Hz
c = 299792458; % m/s
wavelength = 2*pi*c./real_omegas; % m

if ~isfield(sim,'betas')
    sim.betas = calc_sim_betas(fiber,dt,omegas,initial_condition.fields);
end

%% Pre-calculate the factor used in UPPE (the nonlinear constant)
prefactor = {1i*real_omegas/4./ifftshift(sim.mode_profiles.norms(:,1),1).^4,... % I use the 1st mode for norm^4 computation.
             ...                                                                % This can create certain amount of deviation for higher-order modes, but it should be acceptable; otherwise, the computation is too heavy.
             3*permittivity0*ifftshift(sim.X3_prefactor,1)/4}; % for Kerr term

% Photoionization prefactor
if sim.photoionization_model
    me = 9.1093837e-31; % kg; electron mass
    e = 1.60217663e-19; % Coulomb; electron charge
    c = 299792458;
    prefactor = [prefactor,...
                 {-1i/4./ifftshift(sim.mode_profiles.norms(:,1),1).^2*e^2/me./real_omegas,... % ionization-related loss
                  -gas.ionization_energy/8/fiber.SR(1)*permittivity0*c./ifftshift(sim.mode_profiles.norms(:,1),1).^2}]; % ionization-related phase contribution
end

if sim.gpu_yes
    for i = 1:length(prefactor)
        prefactor{i} = gpuArray(prefactor{i});
    end
end

%% Zero-padding for upsampling computation of Kerr and Raman effect
% Upsampling to avoid frequency aliasing
%
% Aliasing mostly comes from Raman shift. However, when the spectrum
% becomes broad, aliasing can come from Kerr effect as well due to
% four-wave mixing.
% Instead of applying upsampling to the Raman computation only, it's
% important to apply it to Kerr as well, especially when running
% supercontinuum generation or when your frequency window isn't large
% enough.
[~,prefactor,...
 upsampling_zeros] = upsampling(sim,Nt,gas_Nt,num_modes,[],prefactor);

%% Set up some parameters for the gas Raman generation equations
% gas_eqn is a container of useful values that will be continuously used
% during propagation.
gas_func = Raman_model_for_UPPE_gradient_pressure();

gas_eqn = gas_func.precalc_gas_params(sim,gas,Nt,...
                                      gas_Nt,gas_dt,upsampling_zeros);

%% Include spontaneous Raman scattering
% Current model of spontaneous Raman scattering assumes only isotropic
% Raman scattering or in the case of conv(R,|A|^2) for the Raman index
% modulation. This form can be achieved when there is only isotropic Raman
% scattering or scalar fields. In other situations, to compute Raman
% correctly, include one photon per frequency band in the electric field of
% each mode instead.
sim.include_sponRS = sim.Raman_model ~= 0;
if sim.include_sponRS
    sponRS_prefactor = spontaneous_Raman(Nt,dt,sim,gas,gas_eqn);
else
    sponRS_prefactor = []; % dummay variable due to no Raman
end

%% Finalize a few parameters for MPA computation
% Its size needs to be modified according to the stepping algorithm
if isequal(sim.step_method,'MPA')
    prefactor{2} = permute(prefactor{2},[1,3,2]); % size: (Nt, M+1, num_modes)
    if sim.Raman_model ~= 0
        sponRS_prefactor{1} = permute(sponRS_prefactor{1},[1,3,2]); % size: (Nt, M+1, num_modes)
    end
end

%% Work out the overlap tensor details
[SK_info, SRa_info, SRb_info] = calc_SRSK(fiber,sim,num_spatial_modes);

%% Include the shot noise: "sim.num_photon_noise_per_bin" photons per mode
initial_condition.fields = include_shot_noise(sim,real_omegas,Nt,time_window,initial_condition.fields);

%% Setup the exact save points
% We will always save the initial condition as well
save_points = int64(num_saves_total + 1);
save_z = double(0:save_points-1)*sim.save_period;

%% Photoionization - erfi() lookup table
% Because calculating erfi is slow, it's faster if I create a lookup table
% and use interp1. The range of input variable for erfi is 0~sqrt(2) only.
if sim.photoionization_model ~= 0
    n_Am = 10; % the number of summation of Am term in photoionization
    gas_eqn.erfi_x = linspace(0,sqrt(2*(n_Am+1)),1000)';
    gas_eqn.erfi_y = erfi(gas_eqn.erfi_x);
end

%% Run the step function over each step
run_start = tic;
% -------------------------------------------------------------------------
[A_out,...
 save_z,save_deltaZ,...
 T_delay_out,...
 delta_permittivity,relative_Ne] = SteppingCaller_gradient_pressure(fiber, sim, gas, gas_eqn,...
                                                                    save_z, save_points,...
                                                                    initial_condition,...
                                                                    SK_info, SRa_info, SRb_info,...
                                                                    prefactor, sponRS_prefactor,...
                                                                    time_window,...
                                                                    omegas,wavelength,...
                                                                    gas_func);
% -------------------------------------------------------------------------
% Just to get an accurate timing, wait before recording the time
if sim.gpu_yes
    sim.betas = gather(sim.betas);
    wait(sim.gpuDevice.Device);
end
fulltime = toc(run_start);

%% Save the results in a struct
foutput = struct('z',save_z,...
                 'deltaZ',save_deltaZ,...
                 'fields',A_out,...
                 'dt',initial_condition.dt,...
                 'betas',sim.betas,...
                 'seconds',fulltime,...
                 't_delay',T_delay_out);
if sim.Raman_model ~= 0 && sim.scalar
    foutput.delta_permittivity = delta_permittivity;
end
if sim.photoionization_model ~= 0
    foutput.relative_Ne = relative_Ne;
end

end
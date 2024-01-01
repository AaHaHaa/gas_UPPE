function [fiber,sim] = load_default_UPPE_propagate( input_fiber,input_sim )
%LOAD_DEFAULT_UPPE_PROPAGATE It loads the default settings for "fiber"
%and "sim" for different types of modes used.
%
%   If a user has specified some of the parameters of "fiber" and "sim",
%   user-defined one will be chosen instead of default ones.
%
%
% Example Use:
%
%    % User-defined parameters
%    fiber.L0 = 3;
%    
%    % Incorporate default settings
%    [fiber,sim] = load_default_UPPE_propagate(fiber,[]);
%
%    % If there are "sim" settings
%    sim.Raman_model = 0;
%    [fiber,sim] =  load_default_UPPE_propagate(fiber,sim);
%
%    % Use only user-defined "sim", not "fiber"
%    [fiber,sim] = load_default_UPPE_propagate([],sim);
%
% -------------------------------------------------------------------------
%
%	Additional parameters:
%
%       input_sim.lambda0 - central wavelength, in m
%       input_sim.midx - an array of the mode index of modes in calculations
%
%       -------------------------------------------------------------------
%       --     Explanation of "midx":
%       --         If I want only mode 2 and 4 in simulations, they should be set as
%       -- 
%       --             sim.midx = [2 4];
%       -- 
%       -------------------------------------------------------------------
%
% -------------------------------------------------------------------------
%
%   If both "lambda0" and "f0" are set by users, the final value will depend on "f0".
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
%           SR - SR tensor, in m^-2
%           L0 - length of fiber, in m
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
%           midx - mode index; an integer array
%           f0 - center frequency, in THz
%           save_period - spatial period between saves, in m. 0 = only save input and output
%           X3 - electronic nonlinearity (Nt,num_modes)
%
%       Mode info -->
%
%           mode_profiles.mode_profiles - 2D mode profiles of each mode
%           mode_profiles.norms - the norms of the mode profiles
%           mode_profiles.r: (m)
%           mode_profiles.dr: (m)
%           mode_profiles.dtheta: (rad)
%
%       MPA -->
%
%           MPA.M - parallel extent for MPA;
%                   1 is no parallelization,
%                   5-20 is recommended; there are strongly diminishing returns after 5-10.
%           MPA.n_tot_max - maximum number of iterations for MPA
%                           This doesn't really matter because if the step size is too large, the algorithm will diverge after a few iterations.
%           MPA.n_tot_min - minimum number of iterations for MPA
%           MPA.tol - tolerance for convergence for MPA
%                     Value of the average NRMSE between consecutive itertaions in MPA at which the step is considered converged.
%
%       Polarization included -->
%
%           scalar - 0(false) = consider polarization mode coupling
%                    1(true) = don't consider polarization mode coupling
%
%           *If under scalar field, the input field takes only the scalar fields, e.g., [mode1, mode2, mode3......].
%           *Otherwise, the input field of each polarization needs to be specified in the order of [mode1_+ mode1_- mode2_+ mode2_-......], where (+,-) can be (x,y) or any orthogonal modes.
%           *SRSK is always loaded in a dimension of num_spatial_modes^4. It's automatically calculated to its polarized version in the code.
%
%           ellipticity - the ellipticity of the polarization modes; Please refer to "Nonlinear Fiber Optics, eq (6.1.18) Agrawal" for the equations.
%                         0: linear polarization   -> (+,-)=(x,y)
%                         1: circular polarization -> (+,-)=(right,left)
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
%           photoionization_model - 0 = ignore the photoionization effect
%                                   1 = include the photoionization effect
%                                   (Photoionization model is implemented currently in single-mode scenarios.)
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

%% Current path (or the folder where this "load_default_UPPE_propagate.m" is)
if ispc
    sep = '\';
else % unix
    sep = '/';
end
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep);
upper_folder = current_path(1:sep_pos(end-1));

%% Default settings below:

% Supperss warnings generated by the function, "catstruct", due to there
% are definitely duplicate elements between default and input.
warning('off','catstruct:DuplicatesFound');

if ~exist('input_fiber','var')
    input_fiber = [];
end
if ~exist('input_sim','var')
    input_sim = [];
end

% -------------------------------------------------------------------------
% Set some default parameters early here because these parameters will be
% used for loading files if multimode.
% If single-mode, only "fiber.lambda0" or "fiber.f0" is important.
% -------------------------------------------------------------------------
c = 2.99792458e-4; % speed of ligth, m/ps

% Get lambda0 from input f0 or lambda0
if isfield(input_sim,'f0')
    default_sim.lambda0 = c/input_sim.f0;
else
    if isfield(input_sim,'lambda0')
        default_sim.lambda0 = input_sim.lambda0;
    else
        default_sim.lambda0 = 1030e-9;
    end
end

% -------------------------------------------------------------------------
% capillary
% -------------------------------------------------------------------------
% Basic properties
if ~isfield(input_sim,'midx') || isempty(input_sim.midx)
    default_sim.midx = 1;
else
    default_sim.midx = input_sim.midx;
end

if isfield(input_fiber,'L0')
    default_fiber.L0 = input_fiber.L0;
else
    default_fiber.L0 = 1; % m
end

% -------------------------------------------------------------------------
% sim
% -------------------------------------------------------------------------
% Basic settings
default_sim.f0 = c/default_sim.lambda0;
default_sim.save_period = 0;
default_sim.ellipticity = 0; % linear polarization

% MPA
default_sim.MPA.M = 10;
default_sim.MPA.n_tot_max = 20;
default_sim.MPA.n_tot_min = 2;
default_sim.MPA.tol = 1e-6;

% Polarization modes
default_sim.scalar = true;

% Adaptive-step algorithm
if length(default_sim.midx) == 1 % single mode
    % Threshold error for adaptive RK4IP
    default_sim.adaptive_deltaZ.threshold = 1e-8; % the threshold of the adaptive method
                                                  % Recommended value is less than 1e-5.
                                                  % Values larger than 1e-3 are too large.
else % multimode
    % Threshold error for adaptive Adams-Moulton method in MPA
    default_sim.adaptive_deltaZ.threshold = 1e-6;
end

% Algorithms to use
default_sim.gpu_yes = true;
default_sim.Raman_model = 1; % consider Raman
default_sim.parallel_yes = false;
default_sim.parallel_idx = 1;
default_sim.photoionization_model = 0;

% Others
default_sim.pulse_centering = true; % center the pulse according to the time window
default_sim.num_photon_noise_per_bin = 1;
default_sim.gpuDevice.Index = 1; % the gpuDevice to use
default_sim.progress_bar = true;
default_sim.progress_bar_name = '';
default_sim.cuda_dir_path = fullfile(upper_folder,'cuda');

%%
% =========================================================================
% Merge settings with the input, which have higher priorities than the
% default ones.
% =========================================================================
if isempty(input_fiber)
    fiber = default_fiber;
elseif isstruct(input_fiber)
    fiber = catstruct(default_fiber, input_fiber);
else
    error('LoadDefaultUPPEPropagate:InputFiberError',...
            '"input_fiber" should be a "structure".');
end
if isempty(input_sim)
    sim = default_sim;
elseif isstruct(input_sim)
    sim = catstruct(default_sim, input_sim);
else
    error('LoadDefaultUPPEPropagate:InputSimError',...
            '"input_sim" should be a "structure".');
end

end
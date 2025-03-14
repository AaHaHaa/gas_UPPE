function [A_out,...
          save_z,save_dz,...
          T_delay_out,...
          delta_permittivity,relative_Ne] = SteppingCaller_gradient_pressure(fiber, sim, gas, gas_eqn,...
                                                                             save_z, save_points,...
                                                                             initial_condition,...
                                                                             SK_info, SRa_info, SRb_info,...
                                                                             prefactor,...
                                                                             At_noise,...
                                                                             time_window,...
                                                                             Omega, wavelength,...
                                                                             gas_func)
%STEPPINGCALLER_GRADIENT_PRESSURE It starts the pulse propagation.

Nt = size(initial_condition.fields,1);
num_modes = size(initial_condition.fields,2);
dt = initial_condition.dt;

save_dz = zeros(save_points,1);
T_delay_out = zeros(save_points,1);
if sim.include_Raman && sim.scalar
    delta_permittivity = zeros(Nt,num_modes,gas_eqn.num_Raman,save_points);
else
    delta_permittivity = []; % dummay variable for output
end
if sim.photoionization_model ~= 0
    relative_Ne = zeros(Nt,1,save_points); % excited electrons due to photoionization
else
    relative_Ne = []; % dummy variable for output
end

% Pulse centering based on the moment of its intensity
if sim.pulse_centering
    % Center the pulse
    temporal_profile = abs(initial_condition.fields).^2;
    temporal_profile(temporal_profile<max(temporal_profile,[],1)/10) = 0;
    TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*temporal_profile,[1,2])/sum(temporal_profile,[1,2]));
    % Because circshift is slow on GPU, I discard it.
    %last_result = ifft(circshift(initial_condition.fields,-tCenter),[],1);
    if TCenter ~= 0
        if TCenter > 0
            initial_condition.fields = [initial_condition.fields(1+TCenter:end,:);initial_condition.fields(1:TCenter,:)];
        elseif TCenter < 0
            initial_condition.fields = [initial_condition.fields(end+1+TCenter:end,:);initial_condition.fields(1:end+TCenter,:)];
        end

        if sim.gpu_yes
            TCenter = gather(TCenter);
        end
        T_delay = TCenter*initial_condition.dt;
    else
        T_delay = 0;
    end
else
    T_delay = 0;
end
T_delay_out(1) = T_delay;

A_out = zeros(Nt, num_modes, save_points);

% Start by saving the initial condition
A_out(:,:,1) = initial_condition.fields;

last_A = ifft(initial_condition.fields,[],1);
if sim.gpu_yes
    last_A = gpuArray(last_A);
end

% Create a progress bar first
if sim.progress_bar
    if ~isfield(sim,'progress_bar_name')
        sim.progress_bar_name = '';
    elseif ~ischar(sim.progress_bar_name)
        error('UPPE_propagate:ProgressBarNameError',...
            '"sim.progress_bar_name" should be a string.');
    end
    h_progress_bar = waitbar(0,sprintf('%s   0.00%%',sim.progress_bar_name),...
        'Name',sprintf('Running gas MM-UPPE: %s...',sim.progress_bar_name),...
        'CreateCancelBtn',...
        'setappdata(gcbf,''canceling'',1)');
    setappdata(h_progress_bar,'canceling',0);

    % Create the cleanup object
    cleanupObj = onCleanup(@()cleanMeUp(h_progress_bar));

    % Use this to control the number of updated time for the progress bar below 1000 times.
    num_progress_updates = 1000;
    progress_bar_z = (1:num_progress_updates)*save_z(end)/num_progress_updates;
    progress_bar_i = 1;
end

% Use this to control the number of updated time for the gas pressure below 100 times.
num_gas_updates = 100;
gas_pressure_steps = (1:num_gas_updates)*(gas.pressure_out-gas.pressure_in)/num_gas_updates + gas.pressure_in;
gas_i = 0;
if gas.pressure_out > gas.pressure_in % increasing gas pressure
    gas_stepping_condition = @(P,i) P - gas_pressure_steps(i);
elseif gas.pressure_out < gas.pressure_in % decreasing gas pressure
    gas_stepping_condition = @(P,i) gas_pressure_steps(i) - P;
else % same input and output pressures
    gas_stepping_condition = @(P,i) false; % There's no need to update propagating parameters due to the unchanged pressure
end

z = 0;
save_i = 2; % the 1st one is the initial field
a5 = [];
if ~isfield(sim.adaptive_dz,'max_dz')
    sim.adaptive_dz.max_dz = sim.save_period/10;
end
sim.dz = min(1e-6,sim.adaptive_dz.max_dz); % m; start with a small value to avoid initial blowup
save_dz(1) = sim.dz;
sim.last_dz = 1; % randomly put a number, 1, for initialization

time_to_update_D = 100; D_idx = time_to_update_D;
pressure0 = 1.01325e5; % Pa
temperature0 = 273.15; % 0 degree Celsius

UPPE_func = str2func(['stepping_' sim.step_method '_gradient_pressure']);

% Then start the propagation
while z+eps(z) < save_z(end) % eps(z) here is necessary due to the numerical error
    % Check for Cancel button press
    if sim.progress_bar && getappdata(h_progress_bar,'canceling')
        error('SteppingCaller_gradient_pressure:ProgressBarBreak',...
              'The "cancel" button of the progress bar has been clicked.');
    end

    % Update Raman parameters and dispersion based on a different gradient pressure
    % Calculate the required parameters based on the gas pressure at the position z
    gas_pressure = sqrt( gas.pressure_in^2 + z/fiber.L0*(gas.pressure_out^2-gas.pressure_in^2) );
    % eta is calculated with the unit, amagat
    % Ideal gas law is used here.
    eta = gas_pressure/pressure0*temperature0/gas.temperature;
    if z == 0 || gas_stepping_condition(gas_pressure,gas_i) > 0
        [gas_i,...
         gas,gas_eqn,...
         Raw,Rbw,...
         D_op,D_op_upsampling] = gas_func.update_R_D(fiber,sim,gas,gas_eqn,...
                                                     gas_pressure_steps,...
                                                     gas_pressure,eta,...
                                                     time_window,...
                                                     Omega,wavelength);
    end

    % Start the steppping
    ever_fail = false;
    previous_A = last_A;
    previous_a5 = a5;

    success = false;
    while ~success
        if ever_fail
            last_A = previous_A;
            a5 = previous_a5;
        end

        [last_A, a5,...
         opt_dz, success] = UPPE_func(last_A, a5,...
                                      sim, gas, gas_eqn,...
                                      SK_info, SRa_info, SRb_info,...
                                      Raw, Rbw,...
                                      D_op_upsampling,...
                                      prefactor,...
                                      At_noise,...
                                      dt, fiber.SR(1),...
                                      eta);

        if ~success
            if opt_dz < 1e-10
                error('SteppingCaller_gradient_pressure:adaptiveRK4IPError',...
                      'Adaptive RK4IP continues to fail.\nCheck simulation parameters.');
            end

            ever_fail = true;

            sim.dz = opt_dz;
        end
    end
    sim.last_dz = sim.dz; % previous dz
    
    % Check for any NaN elements
    if any(any(isnan(last_A))) %any(isnan(last_A),'all')
        error('SteppingCaller_gradient_pressure:NaNError',...
              'NaN field encountered, aborting.\nReduce the step size or increase the temporal or frequency window.');
    end

    % Center the pulse in the time window
    % Important:
    % In the modified shot-noise approach, the noise cannot be changed, so it needs to be translated as well.
    % This took me more than 12 hours of debugging to realize it.
    % Otherwise, the output field, if it has a strong frequency shift and shifts a lot in time relative to the time window, 
    % the noise without translation cannot match with the translated field,
    % resulting in a different noise field overlapped with the coherent pulse.
    % This will artificially create a noisy output.
    if sim.pulse_centering
        last_A_in_time = fft(last_A,[],1);
        temporal_profile = abs(last_A_in_time).^2;
        temporal_profile(temporal_profile<max(temporal_profile,[],1)/10) = 0;
        TCenter = floor(sum((-floor(Nt/2):floor((Nt-1)/2))'.*temporal_profile,[1,2])/sum(temporal_profile,[1,2]));
        % Because circshift is slow on GPU, I discard it.
        %last_result = ifft(circshift(last_A_in_time,-tCenter),[],1);
        if TCenter ~= 0
            if ~isempty(a5) % RK4IP reuses a5 from the previous step
                a5 = fft(a5,[],1);
            end
            if TCenter > 0
                if ~isempty(a5) % RK4IP reuses a5 from the previous step
                    a5 = ifft([a5(1+TCenter:end,:,:,:);a5(1:TCenter,:,:,:)],[],1);
                end
                last_A = ifft([last_A_in_time(1+TCenter:end,:);last_A_in_time(1:TCenter,:)],[],1);
                At_noise = cat(1,At_noise(1+TCenter:end,:,:),At_noise(1:TCenter,:,:));
            elseif TCenter < 0
                if ~isempty(a5) % RK4IP reuses a5 from the previous step
                    a5 = ifft([a5(end+1+TCenter:end,:,:,:);a5(1:end+TCenter,:,:,:)],[],1);
                end
                last_A = ifft([last_A_in_time(end+1+TCenter:end,:);last_A_in_time(1:end+TCenter,:)],[],1);
                At_noise = cat(1,At_noise(end+1+TCenter:end,:,:),At_noise(1:end+TCenter,:,:));
            end
            if sim.gpu_yes
                TCenter = gather(TCenter);
            end
            T_delay = T_delay + TCenter*dt;
        end
    end

    % Update z
    z = z + sim.dz;
    % Because the adaptive-step algorithm determines the step size by 
    % checking the error of the spectral intensities from RK3 and RK4, it
    % might ignore changes at the weakest part of the spectrum. This
    % happens in cases of noise-seeded four-wave mixing and noise-seeded
    % Raman scattering far from the pulse central frequency.
    %
    % To account for this effect, I limit the step size to be 10x smaller 
    % than the effective maximum beat length which is
    % 2*pi/(max(eff_betas)-min(eff_betas)).
    % eff_betas is from betas, propagation constants, throughout the time 
    % window but scaled according to the spectral intensity to prevent 
    % taking into account betas where there's no light at all or where 
    % there's some light starting to grow.
    if D_idx == time_to_update_D
        eff_range_D = find_range_D(abs(last_A).^2,imag(D_op));
        min_beat_length = 2*pi/eff_range_D;
        if strcmp(sim.step_method,'MPA')
            dz_resolve_beat_length = min_beat_length/4*sim.MPA.M;
        else
            dz_resolve_beat_length = min_beat_length/4;
        end
        
        D_idx = 0;
    end
    D_idx = D_idx + 1;

    sim.dz = min([opt_dz,save_z(end)-z,sim.adaptive_dz.max_dz,dz_resolve_beat_length]);

    % If it's time to save, get the result from the GPU if necessary,
    % transform to the time domain, and save it
    if z == sim.last_dz
        if sim.include_Raman && sim.scalar
            delta_permittivity(:,:,:,1) = calc_permittivity(sim,gas,gas_eqn,last_A,Nt);
        end
        if sim.photoionization_model ~= 0
            A_out_ii = fft(last_A,[],1);
            relative_Ne(:,:,1) = calc_Ne(A_out_ii, dt, fiber.SR(1), gas, gas_eqn, sim, eta);
        end
    end
    if z >= save_z(save_i)-eps(z)
        A_out_ii = fft(last_A,[],1);
        if sim.gpu_yes
            save_dz(save_i) = gather(sim.last_dz);
            save_z(save_i) = gather(z);
            A_out(:, :, save_i) = gather(A_out_ii);
        else
            save_dz(save_i) = sim.last_dz;
            save_z(save_i) = z;
            A_out(:, :, save_i) = A_out_ii;
        end
        if sim.include_Raman && sim.scalar
            delta_permittivity(:,:,:,save_i) = calc_permittivity(sim,gas,gas_eqn,last_A,Nt);
        end
        if sim.photoionization_model ~= 0
            relative_Ne(:,:,save_i) = calc_Ne(A_out_ii, dt, fiber.SR(1), gas, gas_eqn, sim, eta);
        end

        T_delay_out(save_i) = T_delay;

        save_i = save_i + 1;
    end

    % Report current status in the progress bar's message field
    if sim.progress_bar
        if z >= progress_bar_z(progress_bar_i)
            if  sim.parallel_yes
                fprintf('Lab %u: %s%5.1f%%\n',sim.parallel_idx,sim.progress_bar_name,z/fiber.L0*100);
            else
                waitbar(gather(z/fiber.L0),h_progress_bar,sprintf('%s%5.1f%%',sim.progress_bar_name,z/fiber.L0*100));
            end
            progress_bar_i = find(z<progress_bar_z,1);
        end
    end
end

end

%% Helper functions
function eff_range_D = find_range_D(spectrum,D)
%FIND_RANGE_D
%
% For an adaptive-dz method, the maximum dz is also limited by the 
% range of the propagation constant, beta0.
% If the FWM, Raman, or anything else happens for multiple frequencies 
% where dz can't resolve their beta0 difference, the outcome can be 
% wrong.
% Here, I multiply the (intensity)^(1/5) of the normalized spectrum to the 
% beta0 to consider beta0 difference of the pulse and exclude those without
% the pulse but within the frequency window. (1/5) is to maximize the 
% contribution of the weak part of the spectrum.

nonzero_field = max(spectrum)~=0;
spectrum = spectrum./max(spectrum(:));

eff_D = D(:,nonzero_field).*(spectrum(:,nonzero_field)).^(1/5); % I use ^(1/5) to emphasize the weak part
eff_range_D = max(eff_D(:)) - min(eff_D(:));

end
% -------------------------------------------------------------------------
function delta_permittivity = calc_permittivity(sim,gas,gas_eqn,A_w,Nt)
%CALC_PERMITTIVITY finds the Raman-induced permittivity change
%
% Only the imaginary part corresponds to the actual permittiviy contribution of each Raman response.
% The real part is retained so that it's easier to visualize the "intensity" of the phonon strength by taking abs().

if sim.ellipticity == 0 % linear polarization
    gas_eqn.R_delta_permittivity(:,1:2:end) = gas_eqn.R_delta_permittivity(:,1:2:end)*4;
end

A_t_upsampling = fft(cat(1,A_w(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,A_w(gas_eqn.n+1:end,:)),[],1);

R_delta_permittivity = permute(gas_eqn.R_delta_permittivity,[1,3,2]); % make it [Nt,num_modes,Raman_type]; The Raman_type dimension are R and V
switch gas.model
    case 0
        delta_permittivity = fft(R_delta_permittivity.*ifft(abs(A_t_upsampling).^2, gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1),[],1);
        delta_permittivity = delta_permittivity(gas_eqn.R_downsampling,:,:);
    case 1
        delta_permittivity = fft(R_delta_permittivity.*ifft(abs(A_t_upsampling).^2,[],1),[],1);
end
% Below follows the computational order as the electric field; however,
% it creates strong aliasing during the second "fft()" operation due to the
% downsampling process for an temporally-elongated delta permittivity.
% To do this correctly, we need to apply the temporal and spectral
% downsampling in the opposite order:
%    1. spectral downsampling
%    2. temporal downsampling with a already-zero-padding delta permittivity
%
% Reverse order creates "spectral downsampling with a non-zero-padding
% delta permittivity" that induces aliasing because delta permittivity
% doesn't reduce to zero at the temporal edges.
% 
% If spectral downsampling is by a factor of an integer, we can simply
% pick data points every "integer" points temporally for spectral
% downsampling. Therefore, only with an integer spectral downsampling
% ratio, both downsampling can be done in the order of temporal-then-
% spectral, downsampling.
%delta_permittivity = ifft(delta_permittivity,[],1).*(permute(max(max(real(sim.mode_profiles.mode_profiles),[],1),[],2),[1,3,2])./mean(sim.mode_profiles.norms,1)).^2; % find the max delta_permittivity of each mode
%delta_permittivity = fft(delta_permittivity([1:gas_eqn.n,gas_eqn.Nt-(Nt-gas_eqn.n-1):gas_eqn.Nt],:,:),[],1); % transform back to time domain

delta_permittivity = delta_permittivity(1:floor(gas_eqn.Nt/Nt):end,:,:).*(permute(max(max(real(sim.mode_profiles.mode_profiles),[],1),[],2),[1,3,2])./mean(sim.mode_profiles.norms,1)).^2; % find the max delta_permittivity of each mode

if sim.gpu_yes
    delta_permittivity = gather(delta_permittivity);
end

end
% -------------------------------------------------------------------------
function relative_Ne = calc_Ne(A_t, dt, inverse_Aeff, gas, gas_eqn, sim, eta)

Ne = photoionization_PPT_model(A_t, inverse_Aeff, gas.ionization.energy, sim.f0, dt, gas.Ng*eta,...
                               gas.ionization.l, gas.ionization.Z,...
                               gas_eqn.erfi_x, gas_eqn.erfi_y,...
                               sim.ellipticity);
relative_Ne = Ne/(gas.Ng*eta);

if sim.gpu_yes
    relative_Ne = gather(relative_Ne);
end

end
% =========================================================================
function cleanMeUp(h_progress_bar)
%CLEANMEUP It deletes the progress bar.

% DELETE the progress bar; don't try to CLOSE it.
delete(h_progress_bar);
    
end
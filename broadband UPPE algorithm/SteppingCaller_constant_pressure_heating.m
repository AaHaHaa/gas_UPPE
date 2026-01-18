function [A_out,...
          save_z,save_dz,...
          T_delay_out,...
          sim_betas_out,...
          delta_permittivity_Raman,delta_permittivity_electronic,...
          relative_Ne,...
          temperature] = SteppingCaller_constant_pressure_heating(fiber, sim, gas, gas_eqn,...
                                                                  save_z, save_points,...
                                                                  initial_condition,...
                                                                  SK_info, SRa_info, SRb_info,...
                                                                  prefactor,...
                                                                  At_noise,...
                                                                  mode_profiles_norms,...
                                                                  time_window,...
                                                                  Omega,...
                                                                  gas_func)
%STEPPINGCALLER_CONSTANT_PRESSURE_HEATING It starts the pulse propagation.

room_temperature = gas.temperature; % the initial temperature

Nt = size(initial_condition.fields,1);
num_modes = size(initial_condition.fields,2);
dt = initial_condition.dt;

save_dz = zeros(save_points,1);
T_delay_out = zeros(save_points,1);
sim_betas_out = zeros(1,2,save_points);
delta_permittivity_Raman = []; % dummay variable for output
delta_permittivity_electronic = []; % dummay variable for output
if sim.scalar
    delta_permittivity_electronic = zeros(Nt,num_modes,save_points);
    if sim.include_Raman
        delta_permittivity_Raman = zeros(Nt,num_modes,gas_eqn.cumsum_num_Raman(end),save_points);
    end
end
if sim.include_photoionization
    relative_Ne = zeros(Nt,1,save_points); % excited electrons due to photoionization
else
    relative_Ne = []; % dummy variable for output
end
temperature = zeros(save_points,1); % temperaturea along the fiber

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

z = 0;
save_i = 2; % the 1st one is the initial field
a5 = [];
if ~isfield(sim.adaptive_dz,'max_dz')
    sim.adaptive_dz.max_dz = min(1e-3,sim.save_period/10); % debug
end
sim.dz = min(1e-6,sim.adaptive_dz.max_dz); % m; start with a small value to avoid initial blowup
save_dz(1) = sim.dz;
sim.last_dz = sim.dz; % randomly put a number, for initialization

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

    % As rising temperature can change the Raman behavior, so iteration
    % between optical propagation and finding temperature is needed. This
    % iteration is triggered typically "only" in the first z propagation.
    % In the latter propagation, as dz is small, temperature variation is
    % so small and the iteration isn't needed.
    previous_temperature = gas.temperature;
    need_temperature_iteration = true;
    while need_temperature_iteration
        % Update Raman parameters and dispersion based on a different gradient pressure
        % eta is calculated with the unit, amagat
        % Ideal gas law is used here.
        eta = gas.pressure/pressure0*temperature0/gas.temperature;
        
        [~,...
         gas,gas_eqn,...
         Raw,Rbw,...
         sim_betas,...
         D_op,D_op_upsampling] = gas_func.update_R_D(fiber,sim,gas,gas_eqn,...
                                                     [],...
                                                     [],eta,...
                                                     time_window,...
                                                     Omega,...
                                                     dt,last_A);
    
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
    
        % Find the temperature
        heat_source = sum(abs(last_A).^2 - abs(previous_A).^2)*(Nt*dt)*1e-12*sim.rep_rate; % W
        if heat_source > 0
            error('SteppingCaller_constant_pressure_heating:HeatingError',...
                  'Error in simulations, leading to increasing energy and thus cooling of the system.');
        end
        Q = heat_source*fiber.SR/sim.dz; % W/m^3
        % Conductivity
        gas.heating.conductivity.core = gas_conductivity(gas.pressure_ratio,gas.material,gas.temperature);
        if gas.temperature < 320 % K; the minimm measured temperature of the experimental data of our silica data used here; see gas_info_constant_pressure_heating() for details
            conductivity_silica = 1.15; % the value at the minimum temperature of the experimental data of our silica data used here; see gas_info_constant_pressure_heating() for details
        else
            conductivity_silica = polyval(gas.heating.silica_p,gas.temperature); % see gas_info_constant_pressure_heating() for details
        end
        gas.heating.conductivity.clad = conductivity_silica*ones(1,length(gas.heating.radius.clad)); % W/m/K
        gas.heating.conductivity.clad(1:length(gas.heating.conductivity.clad)-1) = gas.heating.conductivity.core;
        % Equation is from Eq.(20) in
        % Chen et al., "Theoretical Analysis of Heat Distribution in Raman Fiber Lasers and Amplifiers Employing Pure Passive Fiber," IEEE Photonics J. 12(6), 1-13 (2020)
        % Its equation misses a negative for Q.
        gas.temperature = room_temperature - Q*gas.heating.radius.core^2*(1/2/gas.heating.convection/gas.heating.radius.clad(end) +...
                                                                          1/4/gas.heating.conductivity.core*(1-0.32^2) + ...
                                                                          sum(log(gas.heating.radius.clad/gas.heating.radius.core)/2./gas.heating.conductivity.clad)); % I pick r=0.32*core_radius

        % Re-run the propagation if the temperature fluctuates a lot
        if abs(gas.temperature - previous_temperature)/previous_temperature < 1e-3
            need_temperature_iteration = false;
        else % continue with iterations to obtain the correct temperature
            previous_temperature = gas.temperature; % remember the previous one for convergence check

            % To recover the parameters to re-run the propagation
            sim.dz = sim.last_dz;
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
        last_A_in_time = last_A_in_time.*sim.damped_window.t;
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
        % Use the first z point as the data for z=0. This should be at
        % z=1e-6 very close to the fiber input.
        if sim.gpu_yes
            temperature(1) = gather(gas.temperature);
        else
            temperature(1) = gas.temperature;
        end
        if sim.scalar
            delta_permittivity_electronic(:,:,1) = calc_electronic_permittivity(fiber,gas_eqn,mode_profiles_norms,last_A,Nt,prefactor{2},eta);
            if sim.include_Raman
                delta_permittivity_Raman(:,:,:,1) = calc_Raman_permittivity(fiber,sim,gas,gas_eqn,mode_profiles_norms,last_A,Nt,eta);
            end
        end
        if sim.include_photoionization
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
            temperature(save_i) = gather(gas.temperature);
        else
            save_dz(save_i) = sim.last_dz;
            save_z(save_i) = z;
            A_out(:, :, save_i) = A_out_ii;
            temperature(save_i) = gas.temperature;
        end
        if sim.scalar
            delta_permittivity_electronic(:,:,save_i) = calc_electronic_permittivity(fiber,gas_eqn,mode_profiles_norms,last_A,Nt,prefactor{2},eta); %#ok
            if sim.include_Raman
                delta_permittivity_Raman(:,:,:,save_i) = calc_Raman_permittivity(fiber,sim,gas,gas_eqn,mode_profiles_norms,last_A,Nt,eta); %#ok
            end
        end
        if sim.include_photoionization
            relative_Ne(:,:,save_i) = calc_Ne(A_out_ii, dt, fiber.SR(1), gas, gas_eqn, sim, eta);
        end

        T_delay_out(save_i) = T_delay;
        if sim.gpu_yes
            sim_betas_out(:,:,save_i) = gather(sim_betas);
        else
            sim_betas_out(:,:,save_i) = sim_betas;
        end

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
function delta_permittivity_Raman = calc_Raman_permittivity(fiber,sim,gas,gas_eqn,mode_profiles_norms,A_w,Nt,eta)
%CALC_RAMAN_PERMITTIVITY finds the Raman-induced permittivity change
%
% Only the imaginary part corresponds to the actual permittiviy contribution of each Raman response.
% The real part is retained so that it's easier to visualize the "intensity" of the phonon strength by taking abs().

num_gas = length(gas.material);
if sim.ellipticity == 0 % linear polarization
    for gas_i = 1:num_gas
        if ismember(gas.material{gas_i},{'H2','D2','N2','O2','air','N2O','CO2'}) % including rotational Raman
            gas_eqn.R_delta_permittivity(:,gas_eqn.cumsum_num_Raman(gas_i)+1) = gas_eqn.R_delta_permittivity(:,gas_eqn.cumsum_num_Raman(gas_i)+1)*4;
        end
    end
end

A_w_upsampling = cat(1,A_w(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,A_w(gas_eqn.n+1:end,:));

SR = zeros(1,size(A_w,2));
for num_idx = 1:size(A_w,2)
    SR(:,num_idx) = fiber.SR(num_idx,num_idx,num_idx,num_idx);
end
R_delta_permittivity = permute(gas_eqn.R_delta_permittivity,[1,3,2]); % make it [Nt,num_modes,Raman_type]; The Raman_type dimension are R and V
switch gas.model
    case 0
        delta_permittivity_Raman = fft(R_delta_permittivity.*ifft(abs(fft(sqrt(SR)./mode_profiles_norms.*A_w_upsampling,[],1)).^2, gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1),[],1);
        delta_permittivity_Raman = delta_permittivity_Raman(gas_eqn.R_downsampling,:,:);
    case 1
        delta_permittivity_Raman = fft(R_delta_permittivity.*ifft(abs(fft(sqrt(SR)./mode_profiles_norms.*A_w_upsampling,[],1)).^2,[],1),[],1);
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
delta_permittivity_Raman = delta_permittivity_Raman(1:floor(gas_eqn.Nt/Nt):end,:,:);

if sim.gpu_yes
    delta_permittivity_Raman = gather(delta_permittivity_Raman);
end

delta_permittivity_Raman = delta_permittivity_Raman*eta;

end
% -------------------------------------------------------------------------
function delta_permittivity_electronic = calc_electronic_permittivity(fiber,gas_eqn,mode_profiles_norms,A_w,Nt,prefactor,eta)
%CALC_ELECTRONIC_PERMITTIVITY finds the electronically-induced permittivity change

A_w_upsampling = cat(1,A_w(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,A_w(gas_eqn.n+1:end,:));

SR = zeros(1,size(A_w,2));
for num_idx = 1:size(A_w,2)
    SR(:,num_idx) = fiber.SR(num_idx,num_idx,num_idx,num_idx);
end

% delta_permittivity(t) = 3/4*epsilon0*X3*|A(t)|^2 in narrowband case where X3 is a constant
% In general, X3 is frequency-dependent, so we need to modify this into
%
% delta_permittivity = 3/4*epsilon0* invF[ sqrt(X3(omega)) * A(omega)  ]
%
% to incorporate frequency-dependent electronic effect with the field.
delta_permittivity_electronic = abs(fft(sqrt(SR)./mode_profiles_norms.*sqrt(eta*prefactor).*A_w_upsampling)).^2;
delta_permittivity_electronic = delta_permittivity_electronic(1:floor(gas_eqn.Nt/Nt):end,:);

end
% -------------------------------------------------------------------------
function relative_Ne = calc_Ne(A_t, dt, inverse_Aeff, gas, gas_eqn, sim, eta)

num_gas = length(gas.material);
Ne = 0; % initialization
for gas_i = 1:num_gas
    Ne_i = photoionization_PPT_model(A_t, inverse_Aeff, gas.(gas.material{gas_i}).ionization.energy, sim.f0, dt, gas.Ng(gas_i)*eta,...
                                     gas.(gas.material{gas_i}).ionization.l, gas.(gas.material{gas_i}).ionization.Z,...
                                     gas_eqn.erfi_x, gas_eqn.erfi_y,...
                                     sim.ellipticity);
    Ne = Ne + Ne_i;
end
relative_Ne = Ne/sum(gas.Ng*eta);

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
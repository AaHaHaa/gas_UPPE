function func = Raman_model_for_UPPE_gradient_pressure()
%RAMAN_MODEL_FOR_UPPE_GRADIENT_PRESSURE It is a function container for
%computing Raman under gradient pressure

func = struct('precalc_gas_params',@precalc_gas_params,...
              'update_R_D',@update_R_D);

end

%%
function gas_eqn = precalc_gas_params(sim,gas,Nt,...
                                      gas_Nt,gas_dt,upsampling_zeros)
%PRECALC_GAS_PARAMS It computes the parameters required for the Raman
%computation in gases.
%
%   sim: a structure containing
%       sim.include_Raman
%       sim.gpu_yes
%   gas: a structure containing
%       gas.material
%       gas.model
%       gas.(H2, D2, N2, O2, air, CH4, N2O, CO2)
%   time_window: the size of the time window (ps)
%   Nt: the number of time points in the simulation
%   gas_Nt: the number of time points in the simulation
%   gas_dt: the time interval

acyclic_conv_stretch = @(x) 2*x-1;

n = ceil(Nt/2);
if sim.include_Raman
    switch gas.model
        case 0
            gas_eqn = struct('Nt', gas_Nt,'dt',gas_dt,...
                             'acyclic_conv_stretch',acyclic_conv_stretch,...
                             'R_downsampling', [true(gas_Nt,1);false(acyclic_conv_stretch(gas_Nt)-gas_Nt,1)],...
                             'upsampling_zeros', upsampling_zeros,...
                             'n',n,'m',ceil(acyclic_conv_stretch(Nt)/2));
            if sim.gpu_yes
                gas_eqn.R_downsampling = gpuArray(gas_eqn.R_downsampling);
            end
        case 1
            gas_eqn = struct('Nt', gas_Nt,'dt',gas_dt,...
                             'upsampling_zeros', upsampling_zeros,...
                             'n',n);
    end
    
    % Consider only the important part during upsampling in time (zero-padding in frequencies)
    % Because of zero-padding spectrally, convolution operation in
    % frequencies includes a lot of multiplications of these zeros. They're
    % removed from the Raman response, gas_eqn.R.
    % Because the time window can be much smaller than the Raman-response
    % decay time, acyclic convolution is implemented here. This results in
    % rearrangement of gas_eqn.R and the corresponding phase shift for a
    % correct result, gas_eqn.phase_shift and gas_eqn.phase_shift_acyclic.
    % Rearrangement of gas_eqn.R is due to the use of MATLAB fft(X,n,2)
    % which zero-pads to n points at the trailing edge; otherwise, we need 
    % to zero-pad ourselves and might be slower than MATLAB internal
    % function, fft. However, an extra phase shift needs to be considered.
    % This might sound confusing. Think carefully about this.
    switch gas.model
        case 0
            m2 = acyclic_conv_stretch(Nt)-gas_eqn.m;
            gas_eqn.phase_shift_acyclic = exp(1i*2*pi*m2/acyclic_conv_stretch(gas_Nt)*(0:gas_Nt-1)');
            if sim.gpu_yes
                gas_eqn.phase_shift_acyclic = gpuArray(gas_eqn.phase_shift_acyclic);
            end
        case 1
            m2 = Nt-n;
            gas_eqn.phase_shift = exp(1i*2*pi*m2/gas_Nt*(0:gas_Nt-1)');
            if sim.gpu_yes
                gas_eqn.phase_shift = gpuArray(gas_eqn.phase_shift);
            end
    end
    gas_eqn.m2 = m2; % record the index for rearrangement which will be useful later
    
    num_gas = length(gas.material);
    num_Raman = zeros(1,length(gas.material));
    for gas_i = 1:num_gas
        switch gas.material{gas_i}
            case {'H2','D2','N2','O2','air'}
                num_Raman_i = 2;
            case {'N2O','CO2','CH4'}
                num_Raman_i = 1;
            otherwise
                num_Raman_i = 0;
        end
        num_Raman(gas_i) = num_Raman_i;
    end
    cumsum_num_Raman = [0,cumsum(num_Raman)];
    gas_eqn.cumsum_num_Raman = cumsum_num_Raman; % for outputing the delta_permittivity in scalar situations
else % no Raman
    gas_eqn = struct('Nt',gas_Nt,'dt',gas_dt,'upsampling_zeros', upsampling_zeros,'n',n);
end

end

%% UPDATE_R_D
function [gas_i,...
          gas,gas_eqn,...
          Raw,Rbw,...
          sim_betas,...
          D_op,D_op_upsampling] = update_R_D(fiber,sim,gas,gas_eqn,...
                                             gas_pressure_steps,...
                                             gas_pressure,eta,...
                                             time_window,...
                                             Omega,...
                                             dt,fields)
%UPDATE_R_D It updates Raman and dispersion at a specified gas pressure
%(eta)

T = (0:gas_eqn.Nt-1)'*gas_eqn.dt; % ps
if sim.gpu_yes
    T = gpuArray(T);
end

if isfield(gas,'pressure_in') % gradient pressure
    if gas.pressure_out > gas.pressure_in % increasing gas pressure
        gas_i = find(gas_pressure - gas_pressure_steps < 0,1);
    else % decreasing gas pressure
        gas_i = find(gas_pressure_steps - gas_pressure < 0,1);
    end
else % scenarios under constant pressure considering heating
    gas_i = [];
end
% ---------------------------------------------------------------------
num_gas = length(gas.material);
if sim.include_Raman
    % Set up some parameters for the gas Raman generation equations
    if sim.include_heating
        % Recompute Raman strength, preR, under a different temperature
        if sim.include_photoionization % "gas" will be regenerated below, so "ionization" needs to be pre-saved and put back to "gas" later
            gas0 = gas;
        end
        gas = Raman_model(gas,eta);
        if sim.include_photoionization
            for gas_ii = 1:num_gas
                gas.(gas.material{gas_ii}).ionization = gas0.(gas.material{gas_ii}).ionization;
            end
        end
    end
    % Because the dephasing time depends on the gas pressure, Raman
    % response needs to be re-calculated.
    gas = Raman_T2( gas,eta*gas.pressure_ratio );
    switch gas.model
        case 0
            R = []; % initialization
            for gas_ii = 1:num_gas
                switch gas.material{gas_ii}
                    case {'H2','D2','N2','O2'} % rotational + vibrational Raman
                        R_i = [sum(gas.(gas.material{gas_ii}).R.preR.*exp(-T./gas.(gas.material{gas_ii}).R.T2).*exp(1i*gas.(gas.material{gas_ii}).R.omega.*T),2),...
                               sum(gas.(gas.material{gas_ii}).V.preR.*exp(-T./gas.(gas.material{gas_ii}).V.T2).*exp(1i*gas.(gas.material{gas_ii}).V.omega.*T),2)]*gas_eqn.acyclic_conv_stretch(gas_eqn.Nt)*(gas_eqn.dt*1e-12); % acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12) is due to the DFT-version convolution theorem
                    case 'air'
                        R_i = [sum(gas.N2.R.preR.*exp(-T./gas.N2.R.T2).*exp(1i*gas.N2.R.omega.*T),2) + sum(gas.O2.R.preR.*exp(-T./gas.O2.R.T2).*exp(1i*gas.O2.R.omega.*T),2),...
                               sum(gas.N2.V.preR.*exp(-T./gas.N2.V.T2).*exp(1i*gas.N2.V.omega.*T),2) + sum(gas.O2.V.preR.*exp(-T./gas.O2.V.T2).*exp(1i*gas.O2.V.omega.*T),2)]*gas_eqn.acyclic_conv_stretch(gas_eqn.Nt)*(gas_eqn.dt*1e-12); % acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12) is due to the DFT-version convolution theorem
                    case 'CH4' % only vibrational Raman
                        R_i = sum(gas.(gas.material{gas_ii}).V.preR.*exp(-T./gas.(gas.material{gas_ii}).V.T2).*exp(1i*gas.(gas.material{gas_ii}).V.omega.*T),2)*gas_eqn.acyclic_conv_stretch(gas_eqn.Nt)*(gas_eqn.dt*1e-12); % acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12) is due to the DFT-version convolution theorem
                    case {'N2O','CO2'} % only rotational Raman
                        R_i = sum(gas.(gas.material{gas_ii}).R.preR.*exp(-T./gas.(gas.material{gas_ii}).R.T2).*exp(1i*gas.(gas.material{gas_ii}).R.omega.*T),2)*gas_eqn.acyclic_conv_stretch(gas_eqn.Nt)*(gas_eqn.dt*1e-12); % acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12) is due to the DFT-version convolution theorem
                    otherwise
                        R_i = [];
                end

                R = [R,R_i];
            end
            R(isnan(R)) = 0; % in case that some T2=0 such that -T/T2 has a 0/0 term (this happens when the gas pressure is zero)
            % zero-padding in time for the acyclic convolution theorem to avoid time-domain aliasing
            % X = ifft(Y,n,dim) returns the n-point inverse Fourier transform of Y by padding Y with trailing zeros along the dimension "dim" to length n.
            gas_eqn.R_delta_permittivity = ifft(R,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % Raman-induced permittivity change
            Rw = ifft(imag(R),gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % Raman response; only "sin()" part matters in Raman computations
        case 1
            R = []; % initialization
            for gas_ii = 1:num_gas
                switch gas.material{gas_ii}
                    case {'H2','D2','N2','O2'} % rotational + vibrational Raman
                        R_i = [sum(gas.(gas.material{gas_ii}).R.preR.*exp(-T./gas.(gas.material{gas_ii}).R.T2).*exp(1i*gas.(gas.material{gas_ii}).R.omega.*T),2),...
                               sum(gas.(gas.material{gas_ii}).V.preR.*exp(-T./gas.(gas.material{gas_ii}).V.T2).*exp(1i*gas.(gas.material{gas_ii}).V.omega.*T),2)]*(time_window*1e-12); % (time_window*1e-12) is due to the DFT-version convolution theorem
                    case 'air'
                        R_i = [sum(gas.N2.R.preR.*exp(-T./gas.N2.R.T2).*exp(1i*gas.N2.R.omega.*T),2) + sum(gas.O2.R.preR.*exp(-T/gas.O2.R.T2).*exp(1i*gas.O2.R.omega.*T),2),...
                               sum(gas.N2.V.preR.*exp(-T./gas.N2.V.T2).*exp(1i*gas.N2.V.omega.*T),2) + sum(gas.O2.V.preR.*exp(-T/gas.O2.V.T2).*exp(1i*gas.O2.V.omega.*T),2)]*(time_window*1e-12); % (time_window*1e-12) is due to the DFT-version convolution theorem
                    case 'CH4'
                        R_i = sum(gas.(gas.material{gas_ii}).V.preR.*exp(-T./gas.(gas.material{gas_ii}).V.T2).*exp(1i*gas.(gas.material{gas_ii}).V.omega.*T),2)*(time_window*1e-12); % (time_window*1e-12) is due to the DFT-version convolution theorem
                    case {'N2O','CO2'} % only rotational Raman
                        R_i = sum(gas.(gas.material{gas_ii}).R.preR.*exp(-T./gas.(gas.material{gas_ii}).R.T2).*exp(1i*gas.(gas.material{gas_ii}).R.omega.*T),2)*(time_window*1e-12); % (time_window*1e-12) is due to the DFT-version convolution theorem
                end

                R = [R,R_i];
            end
            R(isnan(R)) = 0; % in case that some T2=0 such that -T/T2 has a 0/0 term (this happens when the gas pressure is zero)
            gas_eqn.R_delta_permittivity = ifft(R,[],1); % Raman-induced permittivity change
            Rw = ifft(imag(R),[],1); % Raman response; only "sin()" part matters in Raman computations
    end
    
    % Raman response (under frequency domain for the convolution operation later)
    Rw_rot = 0; % initialization
    Rw_vib = 0; % initialization
    for gas_ii = 1:num_gas
        switch gas.material{gas_ii}
            case {'H2','D2','N2','O2','air'} % rotational + vibrational Raman
                Rw_rot_i = Rw(:,gas_eqn.cumsum_num_Raman(gas_ii)+1);
                Rw_vib_i = Rw(:,gas_eqn.cumsum_num_Raman(gas_ii)+2);
            case 'CH4' % only vibrational Raman
                Rw_rot_i = 0;
                Rw_vib_i = Rw(:,gas_eqn.cumsum_num_Raman(gas_ii)+1);
            case {'N2O','CO2'} % only rotational Raman
                Rw_rot_i = Rw(:,gas_eqn.cumsum_num_Raman(gas_ii)+1);
                Rw_vib_i = 0;
            otherwise
                Rw_rot_i = 0;
                Rw_vib_i = 0;
        end

        Rw_rot = Rw_rot + Rw_rot_i;
        Rw_vib = Rw_vib + Rw_vib_i;
    end
    clear Rw;
    Raw = Rw_vib - 2*Rw_rot; % isotropic Raman response
    Rbw = 6*Rw_rot; % anisotropic Raman response
    
    % For scalar fields, anisotropic Raman, Rbw, is incorporated into Raw.
    if sim.scalar
        if sim.ellipticity == 0 % linear polarization
            Raw = Raw + Rbw;
        else % circular polarization: its SRb=SRa/2, so the factor 1/2 is included here
            Raw = Raw + Rbw/2;
        end
        Rbw = []; % unused dummy variable now
    end

    % Consider only the important part during upsampling in time (zero-padding in frequencies)
    % Because of zero-padding spectrally, convolution operation in
    % frequencies includes a lot of multiplications of these zeros. They're
    % removed from the Raman response, gas_eqn.R.
    % Because the time window can be much smaller than the Raman-response
    % decay time, acyclic convolution is implemented here. This results in
    % rearrangement of gas_eqn.R and the corresponding phase shift for a
    % correct result, gas_eqn.phase_shift and gas_eqn.phase_shift_acyclic.
    % Rearrangement of gas_eqn.R is due to the use of MATLAB fft(X,n,2)
    % which zero-pads to n points at the trailing edge; otherwise, we need 
    % to zero-pad ourselves and might be slower than MATLAB internal
    % function, fft. However, an extra phase shift needs to be considered.
    % This might sound confusing. Think carefully about this.
    switch gas.model
        case 0
            Raw = [Raw(end-(gas_eqn.m2-1):end);Raw(1:gas_eqn.m)];
            if ~isempty(Rbw)
                Rbw = [Rbw(end-(gas_eqn.m2-1):end);Rbw(1:gas_eqn.m)];
            end
        case 1
            n = ceil(length(Omega)/2);
            Raw = [Raw(end-(gas_eqn.m2-1):end);Raw(1:n)];
            if ~isempty(Rbw)
                Rbw = [Rbw(end-(gas_eqn.m2-1):end);Rbw(1:n)];
            end
    end
else
    Raw = []; Rbw = [];
end
% ---------------------------------------------------------------------
% Calculate the propagation constant based on the current gas pressure
fiber.betas = recompute_beta_gradient_pressure(eta,gas);
sim_betas = calc_sim_betas(fiber,dt,Omega,fields);
D_op = 1i*(ifftshift(fiber.betas,1)-(sim_betas(1)+sim_betas(2)*Omega));
% Zero-padding for upsampling computation in RK4IP
D_op_upsampling = cat(1,D_op(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,D_op(gas_eqn.n+1:end,:));

end
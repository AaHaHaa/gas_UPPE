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
%       gas.gas_material
%       gas.model
%       gas.(H2, N2, O2, CH4, or air)
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
    
    switch gas.gas_material
        case {'H2','N2','O2'}
            gas_eqn.num_Raman = 2;
        case 'air'
            gas_eqn.num_Raman = 4;
        case 'CH4'
            gas_eqn.num_Raman = 1;
    end
else % no Raman
    gas_eqn = struct('Nt',gas_Nt,'dt',gas_dt,'upsampling_zeros', upsampling_zeros,'n',n);
end

% Create a damped frequency window to kill the peaks around the edges of the window
gas_eqn.damped_freq_window = create_damped_freq_window(Nt);
if sim.gpu_yes
    gas_eqn.damped_freq_window = gpuArray(gas_eqn.damped_freq_window);
end

end

%% Helper function
function damped_freq_window = create_damped_freq_window(Nt)
% This function creates a sharp damped frequency window to remove anything
% around the edges, especially the high-frequency edge of the window.
% I use a super-Gaussian function here.

f = fftshift((1:Nt)',1);
fc = floor(Nt/2)+1;
ffwhm = Nt*0.8;
f0 = ffwhm/(2*sqrt(log(2))); % 2*sqrt(log(2))=1.665
gexpo = 2*10; % 10 is to make it a sharp window

% A much sharper damped window is used to for the low-frequency side;
% otherwise, it'll easily remove the long-wavelength signal we want to see.
damped_freq_window = exp(-(f-fc).^gexpo/(2*f0^gexpo)).^20; % 20 is to make it a sharp window
%damped_freq_window_low = exp(-(f-fc).^(gexpo*2)/(2*f0^(gexpo*2))).^4; % (gexpo*2) and 4 are to make it a much sharper window
%damped_freq_window(fc:end) = damped_freq_window_low(fc:end);
damped_freq_window(fc:end) = 1;

end

%% UPDATE_R_D
function [gas_i,...
          gas,gas_eqn,...
          Raw,Rbw,...
          D_op,D_op_upsampling] = update_R_D(fiber,sim,gas,gas_eqn,...
                                             gas_pressure_steps,...
                                             gas_pressure,eta,...
                                             time_window,...
                                             omegas,wavelength)
%UPDATE_R_D It updates Raman and dispersion at a specified gas pressure
%(eta)

T = (0:gas_eqn.Nt-1)'*gas_eqn.dt; % ps
if sim.gpu_yes
    T = gpuArray(T);
end

if gas.pressure_out > gas.pressure_in % increasing gas pressure
    gas_i = find(gas_pressure - gas_pressure_steps < 0,1);
else % decreasing gas pressure
    gas_i = find(gas_pressure_steps - gas_pressure < 0,1);
end
% ---------------------------------------------------------------------
if sim.include_Raman
    % Set up some parameters for the gas Raman generation equations
    % Because the dephasing time depends on the gas pressure, Raman
    % response needs to be re-calculated.
    gas = Raman_T2( gas,eta );
    switch gas.model
        case 0
            switch gas.gas_material
                case 'H2'
                    R = [sum(gas.H2.R.preR.*exp(-T/gas.H2.R.T2).*exp(1i*gas.H2.R.omega.*T),2),...
                         sum(gas.H2.V.preR.*exp(-T/gas.H2.V.T2).*exp(1i*gas.H2.V.omega.*T),2)]*gas_eqn.acyclic_conv_stretch(gas_eqn.Nt)*(gas_eqn.dt*1e-12); % acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12) is due to the DFT-version convolution theorem
                case 'N2'
                    R = [sum(gas.N2.R.preR.*exp(-T/gas.N2.R.T2).*exp(1i*gas.N2.R.omega.*T),2),...
                         sum(gas.N2.V.preR.*exp(-T/gas.N2.V.T2).*exp(1i*gas.N2.V.omega.*T),2)]*gas_eqn.acyclic_conv_stretch(gas_eqn.Nt)*(gas_eqn.dt*1e-12); % acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12) is due to the DFT-version convolution theorem
                case 'O2'
                    R = [sum(gas.O2.R.preR.*exp(-T/gas.O2.R.T2).*exp(1i*gas.O2.R.omega.*T),2),...
                         sum(gas.O2.V.preR.*exp(-T/gas.O2.V.T2).*exp(1i*gas.O2.V.omega.*T),2)]*gas_eqn.acyclic_conv_stretch(gas_eqn.Nt)*(gas_eqn.dt*1e-12); % acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12) is due to the DFT-version convolution theorem
                case 'air'
                    R = [sum(gas.N2.R.preR.*exp(-T/gas.N2.R.T2).*exp(1i*gas.N2.R.omega.*T),2),...
                         sum(gas.N2.V.preR.*exp(-T/gas.N2.V.T2).*exp(1i*gas.N2.V.omega.*T),2),...
                         sum(gas.O2.R.preR.*exp(-T/gas.O2.R.T2).*exp(1i*gas.O2.R.omega.*T),2),...
                         sum(gas.O2.V.preR.*exp(-T/gas.O2.V.T2).*exp(1i*gas.O2.V.omega.*T),2)]*gas_eqn.acyclic_conv_stretch(gas_eqn.Nt)*(gas_eqn.dt*1e-12); % acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12) is due to the DFT-version convolution theorem
                case 'CH4'
                    R = gas.CH4.V.preR.*exp(-T/gas.CH4.V.T2).*exp(1i*gas.CH4.V.omega.*T)*gas_eqn.acyclic_conv_stretch(gas_eqn.Nt)*(gas_eqn.dt*1e-12); % acyclic_conv_stretch(gas_Nt)*(gas_dt*1e-12) is due to the DFT-version convolution theorem
            end
            R(isnan(R)) = 0; % in case that some T2=0 such that -T/T2 has a 0/0 term (this happens when the gas pressure is zero)
            % zero-padding in time for the acyclic convolution theorem to avoid time-domain aliasing
            % X = ifft(Y,n,dim) returns the n-point inverse Fourier transform of Y by padding Y with trailing zeros along the dimension "dim" to length n.
            gas_eqn.R_delta_permittivity = ifft(R,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % Raman-induced permittivity change
            Rw = ifft(imag(R),gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % Raman response; only "sin()" part matters in Raman computations
        case 1
            switch gas.gas_material
                case 'H2'
                    R = ifft([sum(gas.H2.R.preR.*exp(-T/gas.H2.R.T2).*exp(1i*gas.H2.R.omega.*T),2),...
                              sum(gas.H2.V.preR.*exp(-T/gas.H2.V.T2).*exp(1i*gas.H2.V.omega.*T),2)])*(time_window*1e-12); % (time_window*1e-12) is due to the DFT-version convolution theorem
                case 'N2'
                    R = ifft([sum(gas.N2.R.preR.*exp(-T/gas.N2.R.T2).*exp(1i*gas.N2.R.omega.*T),2),...
                              sum(gas.N2.V.preR.*exp(-T/gas.N2.V.T2).*exp(1i*gas.N2.V.omega.*T),2)])*(time_window*1e-12); % (time_window*1e-12) is due to the DFT-version convolution theorem
                case 'O2'
                    R = ifft([sum(gas.O2.R.preR.*exp(-T/gas.O2.R.T2).*exp(1i*gas.O2.R.omega.*T),2),...
                              sum(gas.O2.V.preR.*exp(-T/gas.O2.V.T2).*exp(1i*gas.O2.V.omega.*T),2)])*(time_window*1e-12); % (time_window*1e-12) is due to the DFT-version convolution theorem
                case 'air'
                    R = ifft([sum(gas.N2.R.preR.*exp(-T/gas.N2.R.T2).*exp(1i*gas.N2.R.omega.*T),2),...
                              sum(gas.N2.V.preR.*exp(-T/gas.N2.V.T2).*exp(1i*gas.N2.V.omega.*T),2),...
                              sum(gas.O2.R.preR.*exp(-T/gas.O2.R.T2).*exp(1i*gas.O2.R.omega.*T),2),...
                              sum(gas.O2.V.preR.*exp(-T/gas.O2.V.T2).*exp(1i*gas.O2.V.omega.*T),2)])*(time_window*1e-12); % (time_window*1e-12) is due to the DFT-version convolution theorem
                case 'CH4'
                    R = ifft(gas.CH4.V.preR.*exp(-T/gas.CH4.V.T2).*exp(1i*gas.CH4.V.omega.*T))*(time_window*1e-12); % (time_window*1e-12) is due to the DFT-version convolution theorem
            end
            R(isnan(R)) = 0; % in case that some T2=0 such that -T/T2 has a 0/0 term (this happens when the gas pressure is zero)
            gas_eqn.R_delta_permittivity = ifft(R); % Raman-induced permittivity change
            Rw = ifft(imag(R)); % Raman response; only "sin()" part matters in Raman computations
    end
    
    % Raman response (under frequency domain for the convolution operation later)
    switch gas.gas_material
        case {'H2','N2','O2','air'}
            Rw_vib = sum(Rw(:,2:2:end),2);
            Rw_rot = sum(Rw(:,1:2:end),2);
        case 'CH4' % only vib Raman
            Rw_vib = Rw;
            Rw_rot = 0; % no rotational Raman for CH4 due to its molecular symmetry
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
fiber.betas = recompute_beta_gradient_pressure(wavelength,eta,gas);
D_op = 1i*(ifftshift(fiber.betas,1)-(sim.betas(1)+sim.betas(2)*omegas));
% Zero-padding for upsampling computation in RK4IP
D_op_upsampling = cat(1,D_op(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,D_op(gas_eqn.n+1:end,:));

end
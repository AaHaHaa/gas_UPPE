function [A1w,a5,...
          opt_dz,success] = stepping_RK4IP_gradient_pressure(A0w, a5_1,...
                                                             sim, gas, gas_eqn,...
                                                             SK_info, SRa_info, SRb_info,...
                                                             Raw, Rbw,...
                                                             D_op,...
                                                             prefactor, sponRS_prefactor,...
                                                             dt, inverse_Aeff,...
                                                             eta)
%STEPPING_RK4IP_GRADIENT_PRESSURE Take one step according to the UPPE, using the 
%Runge-kutta under the interaction picture
% A0w - initial field, (Nt, num_modes) matrix, in the frequency domain in W^1/2
% dt - time grid point spacing, in ps
%
% sim.f0 - center frequency, in THz
% sim.dz - step size, in m
% sim.singe_yes - 1 = single, 0 = double
% sim.gpu_yes - 1 = GPU, 0 = CPU
% sim.include_Raman - whether Raman is included
%
% nonlin_const - n2*w0/c, in W^-1 m
% SRa_info.SRa - SRa tensor, in m^-2
% SRa_info.nonzero_midx1234s - required SR indices in total
% SRa_info.nonzero_midx34s - required (SR) indices for partial Raman term (only for CPU computation)
%
% omegas - angular frequencies in 1/ps, in the fft ordering
% D - dispersion term for all modes in m^-1, with size (Nt, num_modes)
%
% Output:
% A1w - (Nt, num_modes) matrix with the field after the step, for each mode, in the frequency domain

[Nt,num_modes] = size(A0w);

% Setup the matrices
if sim.gpu_yes
    Kerr = complex(zeros(gas_eqn.Nt, num_modes, num_modes, 'gpuArray'));
    Ra = complex(zeros(gas_eqn.Nt, num_modes, num_modes, 'gpuArray'));
    Rb = complex(zeros(gas_eqn.Nt, num_modes, num_modes, 'gpuArray'));

    Ra_sponRS = complex(zeros(gas_eqn.Nt, num_modes, num_modes, 'gpuArray'));
    Rb_sponRS = complex(zeros(gas_eqn.Nt, num_modes, num_modes, 'gpuArray'));
    
    At_noise = fft(sponRS_prefactor{1}.*sqrt(abs(randn(gas_eqn.Nt,num_modes,'gpuArray'))).*exp(1i*2*pi*rand(gas_eqn.Nt,num_modes,'gpuArray')));
else
    Kerr = complex(zeros(gas_eqn.Nt, num_modes, num_modes));
    Ra = complex(zeros(gas_eqn.Nt, num_modes, num_modes));
    Rb = complex(zeros(gas_eqn.Nt, num_modes, num_modes));

    Ra_sponRS = complex(zeros(gas_eqn.Nt, num_modes, num_modes));
    Rb_sponRS = complex(zeros(gas_eqn.Nt, num_modes, num_modes));
    
    At_noise = fft(sponRS_prefactor{1}.*sqrt(abs(randn(gas_eqn.Nt,num_modes))).*exp(1i*2*pi*rand(gas_eqn.Nt,num_modes)));
end

% Upsampling to avoid frequency aliasing
%
% Aliasing mostly comes from Raman shift. However, when the spectrum
% becomes broad, aliasing can come from Kerr effect as well due to
% four-wave mixing.
% Instead of applying upsampling to the Raman computation only, it's
% important to apply it to Kerr as well, especially when running
% supercontinuum generation or when your frequency window isn't large
% enough.
A0w_upsampling = cat(1,A0w(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,A0w(gas_eqn.n+1:end,:));
if ~isempty(a5_1)
    a5_1_upsampling = cat(1,a5_1(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,a5_1(gas_eqn.n+1:end,:));
end

D = D_op*sim.dz/2;
expD = exp(D);

% 1) Represented under the interaction picture
A_IP = expD.*A0w_upsampling;

% 2) Propagate through the nonlinearity
if isempty(a5_1)
    a5_1_upsampling = N_op(       A0w_upsampling,...
                           sim, gas, gas_eqn,...
                           SK_info, SRa_info, SRb_info,...
                           Kerr, Ra, Rb,...
                           Raw, Rbw,...
                           At_noise, sponRS_prefactor,...
                           Ra_sponRS, Rb_sponRS,...
                           prefactor,...
                           num_modes,...
                           inverse_Aeff, dt/3,...
                           eta);
end
a1 = expD.*a5_1_upsampling;
a2 =                  N_op(       A_IP+a1*(sim.dz/2),...
                           sim, gas, gas_eqn,...
                           SK_info, SRa_info, SRb_info,...
                           Kerr, Ra, Rb,...
                           Raw, Rbw,...
                           At_noise, sponRS_prefactor,...
                           Ra_sponRS, Rb_sponRS,...
                           prefactor,...
                           num_modes,...
                           inverse_Aeff, dt/3,...
                           eta);
a3 =                  N_op(       A_IP+a2*(sim.dz/2),...
                           sim, gas, gas_eqn,...
                           SK_info, SRa_info, SRb_info,...
                           Kerr, Ra, Rb,...
                           Raw, Rbw,...
                           At_noise, sponRS_prefactor,...
                           Ra_sponRS, Rb_sponRS,...
                           prefactor,...
                           num_modes,...
                           inverse_Aeff, dt/3,...
                           eta);
a4 =                  N_op(expD.*(A_IP+a3*(sim.dz)),...
                           sim, gas, gas_eqn,...
                           SK_info, SRa_info, SRb_info,...
                           Kerr, Ra, Rb,...
                           Raw, Rbw,...
                           At_noise, sponRS_prefactor,...
                           Ra_sponRS, Rb_sponRS,...
                           prefactor,...
                           num_modes,...
                           inverse_Aeff, dt/3,...
                           eta);

A1w = expD.*(A_IP + (a1+2*a2+2*a3)*(sim.dz/6)) + a4*(sim.dz/6);

% 3) Local error estimate
a5 =                  N_op(       A1w,...
                           sim, gas, gas_eqn,...
                           SK_info, SRa_info, SRb_info,...
                           Kerr, Ra, Rb,...
                           Raw, Rbw,...
                           At_noise, sponRS_prefactor,...
                           Ra_sponRS, Rb_sponRS,...
                           prefactor,...
                           num_modes,...
                           inverse_Aeff, dt/3,...
                           eta);
err = sum(abs((a4-a5)*(sim.dz/10)).^2,1);

% 4) Stepsize control
normA = sum(abs(A1w).^2,1);
err = sqrt(err./normA);
err = max(err(normA~=0));
if isnan(err)
    opt_dz = 0.5*sim.dz;
    success = false;
else
    opt_dz = max(0.5,min(2,0.8*(sim.adaptive_dz.threshold/err)^(1/4)))*sim.dz;

    success = err < sim.adaptive_dz.threshold;
end

% Downsample them back
A1w = cat(1,A1w(1:gas_eqn.n,:),A1w(end-(Nt-gas_eqn.n-1):end,:));
a5 = cat(1,a5(1:gas_eqn.n,:),a5(end-(Nt-gas_eqn.n-1):end,:));

end

function dAdz = N_op(Aw,...
                     sim, gas, gas_eqn,...
                     SK_info, SRa_info, SRb_info,...
                     Kerr, Ra, Rb,...
                     Raw, Rbw,...
                     At_noise, sponRS_prefactor,...
                     Ra_sponRS, Rb_sponRS,...
                     prefactor,...
                     num_modes,...
                     inverse_Aeff, dt,...
                     eta)
%N_op Calculate dAdz

At = fft(Aw);
At_wNoise = At + At_noise;

% Calculate the large num_modes^4 sum term
if sim.gpu_yes
    % If using the GPU, do the computation with fast CUDA code
    if sim.scalar % scalar fields
        [Kerr,...
         Ra,...
         Ra_sponRS] = feval(sim.cuda_SRSK,...
                            Kerr, Ra, Ra_sponRS,...
                            complex(At), At_noise,...
                            SK_info.SK, SRa_info.SRa,...
                            SRa_info.nonzero_midx1234s,...
                            SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                            sim.include_Raman,...
                            int32(gas_eqn.Nt), 1,...
                            num_modes,...
                            sim.cuda_num_operations_SRSK);
    else % polarized fields
        [Kerr,...
         Ra, Rb,...
         Ra_sponRS, Rb_sponRS] = feval(sim.cuda_SRSK,...
                                       Kerr, Ra, Rb, Ra_sponRS, Rb_sponRS,...
                                       complex(At), At_noise,...
                                       SK_info.SK,   SK_info.nonzero_midx1234s,  SK_info.beginning_nonzero,  SK_info.ending_nonzero,...
                                       SRa_info.SRa, SRa_info.nonzero_midx1234s, SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                                       SRb_info.SRb, SRb_info.nonzero_midx1234s, SRb_info.beginning_nonzero, SRb_info.ending_nonzero,...
                                       sim.include_Raman, ~sim.scalar,...
                                       int32(gas_eqn.Nt), 1,...
                                       num_modes,...
                                       sim.cuda_num_operations_SRSK);
    end
    Kerr = sum(Kerr,3);
else
    % If using the CPU, first precompute SR_mn.
    if sim.include_Raman
        midx34s_sub2ind = @(x)...
            cellfun(@(xx)...
                feval(@(sub) sub2ind(num_modes*ones(1,2),sub{:}), num2cell(xx)),... % this extra "feval" is to get "xx", which is of the size 2x1, into the input arguments of "sub2ind", so transforming "xx" into a 2x1 cell, each containing an integer, and using {:} expansion is necessary
            mat2cell(x,2,ones(1,size(x,2)))); %#ok transform (2,num_nonzero34) midx34s into linear indices of a num_modes-by-num_modes matrix
            % What "midx34s_sub2ind" does (e.g.):
            %
            %   x = [1 3;
            %        5 4]
            %
            %   After "mat2cell": {[1;  {[3;  (2x1 cells, each having 2x1 array)
            %                       5]}   4]}
            %
            %   First,
            %
            %   xx = {[1;  , then after "num2cell": {{1}; (1 cell with 2x1 cell)
            %          5]}                           {5}}
            %
            %   The purpose of separating 1 and 5 into cells is to use
            %   index expansion, {:}, to put them into the input
            %   arguments of "sub2ind" function.
            %
            %   For 6 modes and thus for 6x6 matrix, sub2ind([6 6],1,5) = 25
            %
            %   Do the same for xx = {[3;  and get sub2ind([6 6],3,4) = 21
            %                          4]}
            %   Finally, midx34s_sub2ind = [25 21] (1x2 array)

        SRa_nonzero_midx34s = midx34s_sub2ind(SRa_info.nonzero_midx34s); % the corresponding linear indices of the 3rd-dimensional "num_nonzero34" above
        Ra_mn = At(:, SRa_info.nonzero_midx34s(1,:)).*conj(At(:, SRa_info.nonzero_midx34s(2,:))); % (Nt,num_nonzero34)
        if ~sim.scalar
            SRb_nonzero_midx34s = midx34s_sub2ind(SRb_info.nonzero_midx34s); % the corresponding linear indices of the 3rd-dimensional "num_nonzero34" above
            Rb_mn = At(:, SRb_info.nonzero_midx34s(1,:)).*conj(At(:, SRb_info.nonzero_midx34s(2,:))); % (Nt,num_nonzero34)
        end
        
        % spontaneous Raman scattering
        Ra_mn_sponRS =      At(:, SRa_info.nonzero_midx34s(1,:)).*conj(At_noise(:, SRa_info.nonzero_midx34s(2,:))) +... % (Nt,num_nonzero34)
                       At_noise(:, SRa_info.nonzero_midx34s(1,:)).*conj(     At(:, SRa_info.nonzero_midx34s(2,:))) +...
                       At_noise(:, SRa_info.nonzero_midx34s(1,:)).*conj(At_noise(:, SRa_info.nonzero_midx34s(2,:)));
        if ~sim.scalar
            Rb_mn_sponRS =      At(:, SRb_info.nonzero_midx34s(1,:)).*conj(At_noise(:, SRb_info.nonzero_midx34s(2,:))) +... % (Nt,num_nonzero34)
                           At_noise(:, SRb_info.nonzero_midx34s(1,:)).*conj(     At(:, SRb_info.nonzero_midx34s(2,:))) +...
                           At_noise(:, SRb_info.nonzero_midx34s(1,:)).*conj(At_noise(:, SRb_info.nonzero_midx34s(2,:)));
        end
    end

    % Then calculate Kerr, Ra, and Rb.
    for midx1 = 1:num_modes
        % Kerr
        nz_midx1 = find( SK_info.nonzero_midx1234s(1,:)==midx1 );
        midx2 = SK_info.nonzero_midx1234s(2,nz_midx1);
        midx3 = SK_info.nonzero_midx1234s(3,nz_midx1);
        midx4 = SK_info.nonzero_midx1234s(4,nz_midx1);
        Kerr(:,midx1) = sum(permute(SK_info.SK(nz_midx1),[2 1]).*At_wNoise(:, midx2).*At_wNoise(:, midx3).*conj(At_wNoise(:, midx4)),2);
        if sim.include_Raman
            % Ra
            for midx2 = 1:num_modes
                nz_midx1 = find( SRa_info.nonzero_midx1234s(1,:)==midx1 );
                nz_midx = nz_midx1( SRa_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                midx3 = SRa_info.nonzero_midx1234s(3,nz_midx);
                midx4 = SRa_info.nonzero_midx1234s(4,nz_midx);
                idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                idx = arrayfun(@(i) find(SRa_nonzero_midx34s==i,1), idx); % the indices connecting to the 2nd-dimensional "num_nonzero34" of Ra_mn
                Ra(:, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[2 1]).*Ra_mn(:, idx),2);
                Ra_sponRS(:, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[2 1]).*Ra_mn_sponRS(:, idx),2);
            end
            % Rb
            if ~sim.scalar
                for midx2 = 1:num_modes
                    nz_midx1 = find( SRb_info.nonzero_midx1234s(1,:)==midx1 );
                    nz_midx = nz_midx1( SRb_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                    midx3 = SRb_info.nonzero_midx1234s(3,nz_midx);
                    midx4 = SRb_info.nonzero_midx1234s(4,nz_midx);
                    idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                    idx = arrayfun(@(i) find(SRb_nonzero_midx34s==i,1), idx); % the indices connecting to the 3rd-dimensional "num_nonzero34" of Rb_mn
                    Rb(:, midx1, midx2) = sum(permute(SRb_info.SRb(nz_midx),[2 1]).*Rb_mn(:, idx),2);
                    Rb_sponRS(:, midx1, midx2) = sum(permute(SRb_info.SRb(nz_midx),[2 1]).*Rb_mn_sponRS(:, idx),2);
                end
            end
        end
    end
    if sim.include_Raman
        clear Ra_mn Ra_mn_sponRS;
        if ~sim.scalar
            clear Rb_mn Rb_mn_sponRS;
        end
    end
end

% Calculate h*Ra as F-1(h F(Ra))
% The convolution using Fourier Transform is faster if both arrays are
% large. If one of the array is small, "conv" can be faster.
% Please refer to
% "https://blogs.mathworks.com/steve/2009/11/03/the-conv-function-and-implementation-tradeoffs/"
% for more information.
if sim.include_Raman
    switch gas.model
        case 0
            Ra_upsampling = ifft(Ra,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % zero-padding for the acyclic convolution theorem to avoid time-domain aliasing
            Ra_upsampling = cat(1,Ra_upsampling(end-(gas_eqn.m2-1):end,:,:),Ra_upsampling(1:gas_eqn.m,:,:)); % only the acyclic_conv_stretch(Nt) points contain useful information; others are from the zero-padding result
            
            RaAA = fft(Raw.*Ra_upsampling,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1);
            RaAA = RaAA(gas_eqn.R_downsampling,:,:).*gas_eqn.phase_shift_acyclic; % select only the valid part
            
            % spontaneous Raman scattering
            Ra_sponRS = ifft(Ra_sponRS,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % zero-padding for the acyclic convolution theorem to avoid time-domain aliasing
            Ra_sponRS = cat(1,Ra_sponRS(end-(gas_eqn.m2-1):end,:,:),Ra_sponRS(1:gas_eqn.m,:,:)); % only the acyclic_conv_stretch(Nt) points contain useful information; others are from the zero-padding result

            RaAA_sponRS = fft(Raw.*Ra_sponRS.*sponRS_prefactor{2},gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1);
            RaAA_sponRS = RaAA_sponRS(gas_eqn.R_downsampling,:,:).*gas_eqn.phase_shift_acyclic; % select only the valid part
            
            if ~sim.scalar % polarized fields with an anisotropic Raman contribution from rotational Raman
                Rb_upsampling = ifft(Rb,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % zero-padding for the acyclic convolution theorem to avoid time-domain aliasing
                Rb_upsampling = cat(1,Rb_upsampling(end-(gas_eqn.m2-1):end,:,:),Rb_upsampling(1:gas_eqn.m,:,:)); % only the acyclic_conv_stretch(Nt) points contain useful information; others are from the zero-padding result

                RbAA = fft(Rbw.*Rb_upsampling,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1);
                RbAA = RbAA(gas_eqn.R_downsampling,:,:).*gas_eqn.phase_shift_acyclic; % select only the valid part
                
                % spontaneous Raman scattering
                Rb_sponRS = ifft(Rb_sponRS,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % zero-padding for the acyclic convolution theorem to avoid time-domain aliasing
                Rb_sponRS = cat(1,Rb_sponRS(end-(gas_eqn.m2-1):end,:,:),Rb_sponRS(1:gas_eqn.m,:,:)); % only the acyclic_conv_stretch(Nt) points contain useful information; others are from the zero-padding result

                RbAA_sponRS = fft(Rbw.*Rb_sponRS.*sponRS_prefactor{2},gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1);
                RbAA_sponRS = RbAA_sponRS(gas_eqn.R_downsampling,:,:).*gas_eqn.phase_shift_acyclic; % select only the valid part
            end
        case 1
            Ra = ifft(Ra); % transform into the frequency domain for upsampling
            Ra_upsampling = cat(1,Ra(end-(gas_eqn.m2-1):end,:,:),Ra(1:gas_eqn.n,:,:));
            
            RaAA = fft(Raw.*Ra_upsampling,gas_eqn.Nt,1).*gas_eqn.phase_shift;
            
            % spontaneous Raman scattering
            Ra_sponRS = ifft(Ra_sponRS); % transform into the frequency domain for upsampling
            Ra_sponRS = cat(1,Ra_sponRS(end-(gas_eqn.m2-1):end,:,:),Ra_sponRS(1:gas_eqn.n,:,:));
            
            RaAA_sponRS = fft(Raw.*Ra_sponRS.*sponRS_prefactor{2},gas_eqn.Nt,1).*gas_eqn.phase_shift;
            
            if ~sim.scalar % polarized fields with an anisotropic Raman contribution from rotational Raman
                Rb = ifft(Rb); % transform into the frequency domain for upsampling
                Rb_upsampling = cat(1,Rb(end-(gas_eqn.m2-1):end,:,:),Rb(1:gas_eqn.n,:,:));

                RbAA = fft(Rbw.*Rb_upsampling,gas_eqn.Nt,1).*gas_eqn.phase_shift;
                
                % spontaneous Raman scattering
                Rb_sponRS = ifft(Rb_sponRS); % transform into the frequency domain for upsampling
                Rb_sponRS = cat(1,Rb_sponRS(end-(gas_eqn.m2-1):end,:,:),Rb_sponRS(1:gas_eqn.n,:,:));

                RbAA_sponRS = fft(Rbw.*Rb_sponRS.*sponRS_prefactor{2},gas_eqn.Nt,1).*gas_eqn.phase_shift;
            end
    end
    
    if sim.scalar
        nonlinear = prefactor{2}.*ifft(Kerr) + ifft(sum((RaAA+RaAA_sponRS).*permute(At,[1 3 2]),3));
    else % polarized fields with an anisotropic Raman contribution from rotational Raman
        nonlinear = prefactor{2}.*ifft(Kerr) + ifft(sum((RaAA+RbAA+RaAA_sponRS+RbAA_sponRS).*permute(At,[1 3 2]),3));
    end
else
    nonlinear = prefactor{2}.*ifft(Kerr);
end

if sim.photoionization_model
    [Ne,DNeDt] = photoionization_PPT_model(At, inverse_Aeff, gas.ionization.energy, sim.f0, dt, gas.Ng*eta,...
                                           gas.ionization.l, gas.ionization.Z,...
                                           gas_eqn.erfi_x, gas_eqn.erfi_y,...
                                           sim.ellipticity);
    
    inverse_A2 = abs(At).^2;
    inverse_A2(inverse_A2<max(inverse_A2)/1e5) = max(inverse_A2)/1e5;
    nonlinear_photoionization = prefactor{3}.*ifft(Ne.*At) + prefactor{4}.*ifft(DNeDt./inverse_A2.*At);
else
    nonlinear_photoionization = 0;
end

% Finalized
dAdz = eta*prefactor{1}.*nonlinear + nonlinear_photoionization;

end
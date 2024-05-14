function [A1,a5,...
          opt_deltaZ,success] = stepping_RK4IP_constant_pressure(A0, a5_1,...
                                                                 sim, gas, gas_eqn,...
                                                                 SK_info, SRa_info, SRb_info,...
                                                                 Raw, Rbw,...
                                                                 Raw_sponRS, Rbw_sponRS,...
                                                                 D_op,...
                                                                 prefactor, sponRS_prefactor,...
                                                                 dt, inverse_Aeff)
%STEPPING_RK4IP_CONSTANT_PRESSURE Take one step according to the UPPE, using the 
%Runge-kutta under the interaction picture
%
% Input:
%    A0 - initial field, (N, num_modes) matrix, in the frequency domain in W^1/2
%    a5_1 - the previous RK4 term that can be reused
%
%    sim.scalar - scalar or polarized fields
%    sim.gpu_yes - true = GPU, false = CPU
%
%    sim.cuda_SRSK - the cuda for computing SR and SK values
%
%    sim.Raman_model - which Raman model is used
%
%    sim.f0 - center frequency, in THz
%    sim.deltaZ - step size, in m
%
%    sim.adaptive_deltaZ.threshold - the accuracy used to determined whether to increase or decrease the step size
%
%    gas
%    gas_eqn
%
%    SRa_info.SRa - SRa tensor; m^-2
%    SRa_info.nonzero_midx1234s - required SRa indices in total
%    SRa_info.nonzero_midx34s - required (SRa) indices for partial Raman term (only for CPU computation)
%    SRb_info.SRb - SRb tensor; m^-2
%    SRb_info.nonzero_midx1234s - required SRb indices in total
%    SRb_info.nonzero_midx34s - required (SRb) indices for partial Raman term (only for CPU computation)
%    SK_info.SK - SK tensor; m^2 (unempty if considering polarizaton modes)
%    SK_info.nonzero_midx1234s - required SK indices in total (unempty if considering polarizaton modes)
%
%    Raw - isotropic Raman response in the frequency domain
%    Rbw - anisotropic Raman response in the frequency domain
%
%    D_op - dispersion term Dz/2 (N, num_modes)
%
%    prefactor
%    sponRS_prefactor - prefactor for the spontaneous Raman scattering
%
%    dt - time grid point spacing, in ps
%    inverse_Aeff - 1/Aeff for calculating photoionization (1/m^2)
%
% Output:
%    A1 - the field (in the frequency domain) after one step size (N, num_modes)
%    a5 - the RK4 term that can be reused in the next step
%    opt_deltaZ - recommended step size
%    success - whether the current step size is sufficiently small for the required tolerance

[Nt,num_modes] = size(A0);

% Setup the matrices
if sim.gpu_yes
    Kerr = complex(zeros(gas_eqn.Nt, num_modes, num_modes, 'gpuArray'));
    Ra = complex(zeros(gas_eqn.Nt, num_modes, num_modes, 'gpuArray'));
    Rb = complex(zeros(gas_eqn.Nt, num_modes, num_modes, 'gpuArray'));
    if sim.include_sponRS % spontaneous Raman scattering
        Ra_sponRS = complex(zeros(gas_eqn.Nt, num_modes, num_modes, 'gpuArray'));
        Rb_sponRS = complex(zeros(gas_eqn.Nt, num_modes, num_modes, 'gpuArray'));
        
        A_sponRS = fft(sponRS_prefactor{1}.*sqrt(abs(randn(gas_eqn.Nt,num_modes,'gpuArray'))).*exp(1i*2*pi*rand(gas_eqn.Nt,num_modes,'gpuArray')));
    else
        Ra_sponRS = [];
        Rb_sponRS = [];
        A_sponRS = [];
    end
else
    Kerr = complex(zeros(gas_eqn.Nt, num_modes));
    Ra = complex(zeros(gas_eqn.Nt, num_modes, num_modes));
    Rb = complex(zeros(gas_eqn.Nt, num_modes, num_modes));
    if sim.include_sponRS % spontaneous Raman scattering
        Ra_sponRS = complex(zeros(gas_eqn.Nt, num_modes, num_modes));
        Rb_sponRS = complex(zeros(gas_eqn.Nt, num_modes, num_modes));
        
        A_sponRS = fft(sponRS_prefactor{1}.*sqrt(abs(randn(gas_eqn.Nt,num_modes))).*exp(1i*2*pi*rand(gas_eqn.Nt,num_modes)));
    else
        Ra_sponRS = [];
        Rb_sponRS = [];
        A_sponRS = [];
    end
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
A0_upsampling = cat(1,A0(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,A0(gas_eqn.n+1:end,:));
if ~isempty(a5_1)
    a5_1_upsampling = cat(1,a5_1(1:gas_eqn.n,:),gas_eqn.upsampling_zeros,a5_1(gas_eqn.n+1:end,:));
end

D = D_op*sim.deltaZ/2;
expD = exp(D);

% 1) Represented under the interaction picture
A_IP = expD.*A0_upsampling;

% 2) Propagate through the nonlinearity
if isempty(a5_1)
    a5_1_upsampling = N_op(       A0_upsampling,...
                           sim, gas, gas_eqn,...
                           SK_info, SRa_info, SRb_info,...
                           Kerr, Ra, Rb,...
                           Raw, Rbw,...
                           Raw_sponRS, Rbw_sponRS,...
                           A_sponRS, sponRS_prefactor,...
                           Ra_sponRS, Rb_sponRS,...
                           prefactor,...
                           num_modes,...
                           inverse_Aeff, dt/3);
end
a1 = expD.*a5_1_upsampling;
a2 =                  N_op(       A_IP+a1*(sim.deltaZ/2),...
                           sim, gas, gas_eqn,...
                           SK_info, SRa_info, SRb_info,...
                           Kerr, Ra, Rb,...
                           Raw, Rbw,...
                           Raw_sponRS, Rbw_sponRS,...
                           A_sponRS, sponRS_prefactor,...
                           Ra_sponRS, Rb_sponRS,...
                           prefactor,...
                           num_modes,...
                           inverse_Aeff, dt/3);
a3 =                  N_op(       A_IP+a2*(sim.deltaZ/2),...
                           sim, gas, gas_eqn,...
                           SK_info, SRa_info, SRb_info,...
                           Kerr, Ra, Rb,...
                           Raw, Rbw,...
                           Raw_sponRS, Rbw_sponRS,...
                           A_sponRS, sponRS_prefactor,...
                           Ra_sponRS, Rb_sponRS,...
                           prefactor,...
                           num_modes,...
                           inverse_Aeff, dt/3);
a4 =                  N_op(expD.*(A_IP+a3*(sim.deltaZ)),...
                           sim, gas, gas_eqn,...
                           SK_info, SRa_info, SRb_info,...
                           Kerr, Ra, Rb,...
                           Raw, Rbw,...
                           Raw_sponRS, Rbw_sponRS,...
                           A_sponRS, sponRS_prefactor,...
                           Ra_sponRS, Rb_sponRS,...
                           prefactor,...
                           num_modes,...
                           inverse_Aeff, dt/3);

A1 = expD.*(A_IP + (a1+2*a2+2*a3)*(sim.deltaZ/6)) + a4*(sim.deltaZ/6);

% 3) Local error estimate
a5 =                  N_op(       A1,...
                           sim, gas, gas_eqn,...
                           SK_info, SRa_info, SRb_info,...
                           Kerr, Ra, Rb,...
                           Raw, Rbw,...
                           Raw_sponRS, Rbw_sponRS,...
                           A_sponRS, sponRS_prefactor,...
                           Ra_sponRS, Rb_sponRS,...
                           prefactor,...
                           num_modes,...
                           inverse_Aeff, dt/3);
err = sum(abs((a4-a5)*(sim.deltaZ/10)).^2,1);

% 4) Stepsize control
normA = sum(abs(A1).^2,1);
err = sqrt(err./normA);
err = max(err(normA~=0));
if isnan(err)
    opt_deltaZ = 0.5*sim.deltaZ;
    success = false;
else
    opt_deltaZ = max(0.5,min(2,0.8*(sim.adaptive_deltaZ.threshold/err)^(1/4)))*sim.deltaZ;

    success = err < sim.adaptive_deltaZ.threshold;
end

% Downsample them back
A1 = cat(1,A1(1:gas_eqn.n,:),A1(end-(Nt-gas_eqn.n-1):end,:));
a5 = cat(1,a5(1:gas_eqn.n,:),a5(end-(Nt-gas_eqn.n-1):end,:));

end

function dAdz = N_op(A_w,...
                     sim, gas, gas_eqn,...
                     SK_info, SRa_info, SRb_info,...
                     Kerr, Ra, Rb,...
                     Raw, Rbw,...
                     Raw_sponRS, Rbw_sponRS,...
                     A_sponRS, sponRS_prefactor,...
                     Ra_sponRS, Rb_sponRS,...
                     prefactor,...
                     num_modes,...
                     inverse_Aeff, dt)
%N_op Calculate dAdz

A_t = fft(A_w);

% Calculate large num_modes^4 Kerr, Ra, and Rb terms.
% If not using the GPU, we will precompute Ra_mn and Rb_mn before the num_modes^4 sum
if sim.gpu_yes
    % If using the GPU, do the computation with fast CUDA code
    if sim.scalar % scalar fields
        [Kerr,...
         Ra] = feval(sim.cuda_SRSK,...
                     Kerr, Ra,...
                     complex(A_t),...
                     SK_info.SK, SRa_info.SRa,...
                     SRa_info.nonzero_midx1234s,...
                     SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                     sim.Raman_model~=0,...
                     int32(gas_eqn.Nt), 1,...
                     num_modes,...
                     sim.cuda_num_operations_SRSK);
        if sim.include_sponRS % spontaneous Raman scattering
            Ra_sponRS = feval(sim.cuda_sponRS,...
                              Ra_sponRS,...
                              complex(A_t), A_sponRS,...
                              SRa_info.SRa,...
                              SRa_info.nonzero_midx1234s,...
                              SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                              int32(gas_eqn.Nt), 1,...
                              num_modes);
        end
    else % polarized fields
        [Kerr,...
         Ra, Rb] = feval(sim.cuda_SRSK,...
                         Kerr, Ra, Rb,...
                         complex(A_t),...
                         SK_info.SK,   SK_info.nonzero_midx1234s,  SK_info.beginning_nonzero,  SK_info.ending_nonzero,...
                         SRa_info.SRa, SRa_info.nonzero_midx1234s, SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                         SRb_info.SRb, SRb_info.nonzero_midx1234s, SRb_info.beginning_nonzero, SRb_info.ending_nonzero,...
                         sim.Raman_model~=0, ~sim.scalar,...
                         int32(gas_eqn.Nt), 1,...
                         num_modes,...
                         sim.cuda_num_operations_SRSK);
        if sim.include_sponRS % spontaneous Raman scattering
            [Ra_sponRS, Rb_sponRS] = feval(sim.cuda_sponRS,...
                                           Ra_sponRS, Rb_sponRS,...
                                           complex(A_t), A_sponRS,...
                                           SRa_info.SRa, SRa_info.nonzero_midx1234s, SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                                           SRb_info.SRb, SRb_info.nonzero_midx1234s, SRb_info.beginning_nonzero, SRb_info.ending_nonzero,...
                                           ~sim.scalar,...
                                           int32(gas_eqn.Nt), 1,...
                                           num_modes,...
                                           sim.cuda_num_operations_sponRS);
        end
    end
    Kerr = sum(Kerr,3);
else
    % If using the CPU, first precompute Ra_mn and Rb_mn.
    if sim.Raman_model ~= 0
        midx34s_sub2ind = @(x)...
            cellfun(@(xx)...
                feval(@(sub) sub2ind(num_modes*ones(1,2),sub{:}), num2cell(xx)),... % this extra "feval" is to get "xx", which is of the size 2x1, into the input arguments of "sub2ind", so transforming "xx" into a 2x1 cell, each containing an integer, and using {:} expansion is necessary
            mat2cell(x,2,ones(1,size(x,2)))); % transform (2,num_nonzero34) midx34s into linear indices of a num_modes-by-num_modes matrix
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
        Ra_mn = A_t(:, SRa_info.nonzero_midx34s(1,:)).*conj(A_t(:, SRa_info.nonzero_midx34s(2,:))); % (N,num_nonzero34)
        if ~sim.scalar
            SRb_nonzero_midx34s = midx34s_sub2ind(SRb_info.nonzero_midx34s); % the corresponding linear indices of the 3rd-dimensional "num_nonzero34" above
            Rb_mn = A_t(:, SRb_info.nonzero_midx34s(1,:)).*conj(A_t(:, SRb_info.nonzero_midx34s(2,:))); % (N,num_nonzero34)
        end
        if sim.include_sponRS % spontaneous Raman scattering
            Ra_mn_sponRS =      A_t(:, SRa_info.nonzero_midx34s(1,:)).*conj(A_sponRS(:, SRa_info.nonzero_midx34s(2,:))) +... % (N,num_nonzero34)
                           A_sponRS(:, SRa_info.nonzero_midx34s(1,:)).*conj(     A_t(:, SRa_info.nonzero_midx34s(2,:))) +...
                           A_sponRS(:, SRa_info.nonzero_midx34s(1,:)).*conj(A_sponRS(:, SRa_info.nonzero_midx34s(2,:)));
            if ~sim.scalar
                Rb_mn_sponRS =      A_t(:, SRb_info.nonzero_midx34s(1,:)).*conj(A_sponRS(:, SRb_info.nonzero_midx34s(2,:))) +... % (N,num_nonzero34)
                               A_sponRS(:, SRb_info.nonzero_midx34s(1,:)).*conj(     A_t(:, SRb_info.nonzero_midx34s(2,:))) +...
                               A_sponRS(:, SRb_info.nonzero_midx34s(1,:)).*conj(A_sponRS(:, SRb_info.nonzero_midx34s(2,:)));
            end
        end
    end

    % Then calculate Kerr, Ra, and Rb.
    for midx1 = 1:num_modes
        % Kerr
        nz_midx1 = find( SK_info.nonzero_midx1234s(1,:)==midx1 );
        midx2 = SK_info.nonzero_midx1234s(2,nz_midx1);
        midx3 = SK_info.nonzero_midx1234s(3,nz_midx1);
        midx4 = SK_info.nonzero_midx1234s(4,nz_midx1);
        Kerr(:,midx1) = sum(permute(SK_info.SK(nz_midx1),[2 1]).*A_t(:, midx2).*A_t(:, midx3).*conj(A_t(:, midx4)),2);
        if sim.Raman_model ~= 0
            % Ra
            for midx2 = 1:num_modes
                nz_midx1 = find( SRa_info.nonzero_midx1234s(1,:)==midx1 );
                nz_midx = nz_midx1( SRa_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                midx3 = SRa_info.nonzero_midx1234s(3,nz_midx);
                midx4 = SRa_info.nonzero_midx1234s(4,nz_midx);
                idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                idx = arrayfun(@(i) find(SRa_nonzero_midx34s==i,1), idx); % the indices connecting to the 2nd-dimensional "num_nonzero34" of Ra_mn
                Ra(:, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[2 1]).*Ra_mn(:, idx),2);
                if sim.include_sponRS % spontaneous Raman scattering
                    Ra_sponRS(:, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[2 1]).*Ra_mn_sponRS(:, idx),2);
                end
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
                    if sim.include_sponRS % spontaneous Raman scattering
                        Rb_sponRS(:, midx1, midx2) = sum(permute(SRb_info.SRb(nz_midx),[2 1]).*Rb_mn_sponRS(:, idx),2);
                    end
                end
            end
        end
    end
    if sim.Raman_model ~= 0
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
if sim.Raman_model ~= 0
    switch gas.model
        case 0
            Ra_upsampling = ifft(Ra,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % zero-padding for the acyclic convolution theorem to avoid time-domain aliasing
            Ra_upsampling = cat(1,Ra_upsampling(end-(gas_eqn.m2-1):end,:,:),Ra_upsampling(1:gas_eqn.m,:,:)); % only the acyclic_conv_stretch(Nt) points contain useful information; others are from the zero-padding result
            
            RaAA = fft(Raw.*Ra_upsampling,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1);
            RaAA = RaAA(gas_eqn.R_downsampling,:,:).*gas_eqn.phase_shift_acyclic; % select only the valid part
            
            if sim.include_sponRS % spontaneous Raman scattering
                Ra_sponRS = ifft(Ra_sponRS,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % zero-padding for the acyclic convolution theorem to avoid time-domain aliasing
                Ra_sponRS = cat(1,Ra_sponRS(end-(gas_eqn.m2-1):end,:,:),Ra_sponRS(1:gas_eqn.m,:,:)); % only the acyclic_conv_stretch(Nt) points contain useful information; others are from the zero-padding result

                RaAA_sponRS = fft(Raw_sponRS.*Ra_sponRS.*sponRS_prefactor{2},gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1);
                RaAA_sponRS = RaAA_sponRS(gas_eqn.R_downsampling,:,:).*gas_eqn.phase_shift_acyclic; % select only the valid part
            end
            
            if ~sim.scalar % polarized fields with an anisotropic Raman contribution from rotational Raman
                Rb_upsampling = ifft(Rb,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % zero-padding for the acyclic convolution theorem to avoid time-domain aliasing
                Rb_upsampling = cat(1,Rb_upsampling(end-(gas_eqn.m2-1):end,:,:),Rb_upsampling(1:gas_eqn.m,:,:)); % only the acyclic_conv_stretch(Nt) points contain useful information; others are from the zero-padding result

                RbAA = fft(Rbw.*Rb_upsampling,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1);
                RbAA = RbAA(gas_eqn.R_downsampling,:,:).*gas_eqn.phase_shift_acyclic; % select only the valid part
                
                if sim.include_sponRS % spontaneous Raman scattering
                    Rb_sponRS = ifft(Rb_sponRS,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % zero-padding for the acyclic convolution theorem to avoid time-domain aliasing
                    Rb_sponRS = cat(1,Rb_sponRS(end-(gas_eqn.m2-1):end,:,:),Rb_sponRS(1:gas_eqn.m,:,:)); % only the acyclic_conv_stretch(Nt) points contain useful information; others are from the zero-padding result

                    RbAA_sponRS = fft(Rbw_sponRS.*Rb_sponRS.*sponRS_prefactor{2},gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1);
                    RbAA_sponRS = RbAA_sponRS(gas_eqn.R_downsampling,:,:).*gas_eqn.phase_shift_acyclic; % select only the valid part
                end
            end
        case 1
            Ra = ifft(Ra); % transform into the frequency domain for upsampling
            Ra_upsampling = cat(1,Ra(end-(gas_eqn.m2-1):end,:,:),Ra(1:gas_eqn.n,:,:));
            
            RaAA = fft(Raw.*Ra_upsampling,gas_eqn.Nt,1).*gas_eqn.phase_shift;
            
            if sim.include_sponRS % spontaneous Raman scattering
                Ra_sponRS = ifft(Ra_sponRS); % transform into the frequency domain for upsampling
                Ra_sponRS = cat(1,Ra_sponRS(end-(gas_eqn.m2-1):end,:,:),Ra_sponRS(1:gas_eqn.n,:,:));
                
                RaAA_sponRS = fft(Raw_sponRS.*Ra_sponRS.*sponRS_prefactor{2},gas_eqn.Nt,1).*gas_eqn.phase_shift;
            end
            
            if ~sim.scalar % polarized fields with an anisotropic Raman contribution from rotational Raman
                Rb = ifft(Rb); % transform into the frequency domain for upsampling
                Rb_upsampling = cat(1,Rb(end-(gas_eqn.m2-1):end,:,:),Rb(1:gas_eqn.n,:,:));

                RbAA = fft(Rbw.*Rb_upsampling,gas_eqn.Nt,1).*gas_eqn.phase_shift;
                
                if sim.include_sponRS % spontaneous Raman scattering
                    Rb_sponRS = ifft(Rb_sponRS); % transform into the frequency domain for upsampling
                    Rb_sponRS = cat(1,Rb_sponRS(end-(gas_eqn.m2-1):end,:,:),Rb_sponRS(1:gas_eqn.n,:,:));

                    RbAA_sponRS = fft(Rbw_sponRS.*Rb_sponRS.*sponRS_prefactor{2},gas_eqn.Nt,1).*gas_eqn.phase_shift;
                end
            end
    end
    
    if sim.scalar
        if sim.include_sponRS % spontaneous Raman scattering
            nonlinear = prefactor{2}.*ifft(Kerr) + ifft(sum((RaAA+RaAA_sponRS).*permute(A_t,[1 3 2]),3));
        else
            nonlinear = prefactor{2}.*ifft(Kerr) + ifft(sum(RaAA.*permute(A_t,[1 3 2]),3));
        end
    else % polarized fields with an anisotropic Raman contribution from rotational Raman
        if sim.include_sponRS % spontaneous Raman scattering
            nonlinear = prefactor{2}.*ifft(Kerr) + ifft(sum((RaAA+RbAA+RaAA_sponRS+RbAA_sponRS).*permute(A_t,[1 3 2]),3));
        else
            nonlinear = prefactor{2}.*ifft(Kerr) + ifft(sum((RaAA+RbAA).*permute(A_t,[1 3 2]),3));
        end
    end
else
    nonlinear = prefactor{2}.*ifft(Kerr);
end

if sim.photoionization_model
    [Ne,DNeDt] = photoionization_PPT_model(A_t, inverse_Aeff, gas.ionization.energy, sim.f0, dt, gas.Ng,...
                                           gas.ionization.l, gas.ionization.Z,...
                                           gas_eqn.erfi_x, gas_eqn.erfi_y,...
                                           sim.ellipticity);
    
    inverse_A2 = abs(A_t).^2;
    inverse_A2(inverse_A2<max(inverse_A2)/1e5) = max(inverse_A2)/1e5;
    nonlinear_photoionization = prefactor{3}.*ifft(Ne.*A_t) + prefactor{4}.*ifft(DNeDt./inverse_A2.*A_t);
else
    nonlinear_photoionization = 0;
end

% Finalized
dAdz = prefactor{1}.*nonlinear + nonlinear_photoionization;

end
function [A1w, dummy_a5_1,...
          opt_dz,success]= stepping_MPA_gradient_pressure(A0w, dummy_a5_1,...
                                                          sim, gas, gas_eqn,...
                                                          SK_info, SRa_info, SRb_info,...
                                                          Raw, Rbw,...
                                                          D_op,...
                                                          prefactor,...
                                                          At_noise,...
                                                          dummy_dt, dummy_inverse_Aeff,...
                                                          eta)
%STEPPING_MPA_GRADIENT _PRESSURE Take one step with MPA
%
% Input:
%    A0w - initial field in the frequency domain (Nt, num_modes); sqrt(W)
%    dummy_dt - time grid point spacing; ps; but this isn't used in this function
%
%    sim.dz - step size; m
%
%    sim.MPA.M - parallel extent, 1 is no parallelization
%    sim.MPA.n_tot_max - maximum number of iterations
%    sim.MPA.n_tot_min - minimum number of iterations
%    sim.MPA.tol - tolerance for convergence at each iteration
%
%    sim.scalar - scalar or polarized fields
%    sim.gpu_yes - true = GPU, false = CPU
%
%    sim.cuda_SRSK - the cuda for computing SR and SK values
%
%    sim.include_Raman - whether Raman is included
%
%    sim.adaptive_dz.threshold - the accuracy used to determined whether to increase or decrease the step size
%
%    prefactor - 1i*n2*omega/c; m/W
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
%    D.pos - dispersion term exp(Dz) (Nt, num_modes, M+1)
%    D.neg - dispersion term exp(-Dz) (Nt, num_modes, M+1)
%    haw - isotropic Raman response in the frequency domain
%    hbw - anisotropic Raman response in the frequency domain
%
%    dummy_a5_1 - unused variable
%
%    dummy_G - gain term (Nt, 1); m^-1; but this isn't used in this function
%    dummy_saturation - for gain saturation; but this isn't used in this function
%    dummy_fr - Raman fraction; but this isn't used in this function
%
% Output:
%    A1w - the field (in the frequency domain) after one step size (Nt, num_modes)
%    dummy_a5 - unused variable
%    opt_dz - recommended step size
%    success - whether the current step size is sufficiently small for the required tolerance

small_dz = sim.dz/sim.MPA.M; % the step size between each parallelization plane in MPA

[Nt,num_modes] = size(A0w);

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

% We can pre-compute exp(D_op*z) and exp(-D_op*z) for all z.
z = permute((0:sim.MPA.M)',[2,3,1]);% to make "Dz" into (Nt,num_modes,M+1)
Dz = D_op*small_dz.*z;
D = struct('pos', exp(Dz), 'neg', exp(-Dz));

% Set initial values for psi
psi = repmat(A0w_upsampling, 1, 1, sim.MPA.M+1); % M copies of psi(w,z) = A(w,z), in the frequency domain
current_psi_half_parallelization = psi(:,:,sim.MPA.M+1);

% Iterate to solve for psi
for n_it = 1:sim.MPA.n_tot_max
    % Calculate A(w,z) at all z:
    %    A_w = D.pos*psi
    Aw = permute(D.pos.*psi, [1 3 2]); % (Nt, M+1, num_modes)
    % Calculate A(t,z) at all z
    At = fft(Aw,[],1);
    At_wNoise = At + At_noise;

    % Calculate Kerr, SRa, and SRb terms
    % If not using the GPU, we will precompute Ra_mn and Rb_mn before the num_modes^4 sum

    % Set up matrices for the following Kerr, Ra, and Rb computations
    if sim.gpu_yes
        Kerr = complex(zeros(gas_eqn.Nt, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));
        Ra = complex(zeros(gas_eqn.Nt, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));
        Rb = complex(zeros(gas_eqn.Nt, sim.MPA.M+1, num_modes, num_modes, 'gpuArray'));
    else
        Kerr = complex(zeros(gas_eqn.Nt, sim.MPA.M+1, num_modes));
        Ra = complex(zeros(gas_eqn.Nt, sim.MPA.M+1, num_modes, num_modes));
        Rb = complex(zeros(gas_eqn.Nt, sim.MPA.M+1, num_modes, num_modes));
    end

    % Calculate large num_modes^4 Kerr, Ra, and Rb terms.
    % If not using the GPU, we will precompute Ra_mn and Rb_mn before the num_modes^4 sum
    if sim.gpu_yes
        % If using the GPU, do the computation with fast CUDA code
        if sim.scalar % scalar fields
            [Kerr,...
             Ra] = feval(sim.cuda_SRSK,...
                         Kerr, Ra,...
                         complex(At_wNoise),...
                         SK_info.SK, SRa_info.SRa,...
                         SRa_info.nonzero_midx1234s,...
                         SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                         sim.include_Raman,...
                         int32(gas_eqn.Nt), sim.MPA.M+1,...
                         num_modes,...
                         sim.cuda_num_operations_SRSK);
        else % polarized fields
            [Kerr,...
             Ra, Rb] = feval(sim.cuda_SRSK,...
                             Kerr, Ra, Rb,...
                             complex(At_wNoise),...
                             SK_info.SK,   SK_info.nonzero_midx1234s,  SK_info.beginning_nonzero,  SK_info.ending_nonzero,...
                             SRa_info.SRa, SRa_info.nonzero_midx1234s, SRa_info.beginning_nonzero, SRa_info.ending_nonzero,...
                             SRb_info.SRb, SRb_info.nonzero_midx1234s, SRb_info.beginning_nonzero, SRb_info.ending_nonzero,...
                             sim.include_Raman, ~sim.scalar,...
                             int32(gas_eqn.Nt), sim.MPA.M+1,...
                             num_modes,...
                             sim.cuda_num_operations_SRSK);
        end
        Kerr = sum(Kerr,4);
    else
        % If using the CPU, first precompute Ra_mn and Rb_mn.
        if sim.include_Raman
            if n_it == 1
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
                if ~sim.scalar
                    SRb_nonzero_midx34s = midx34s_sub2ind(SRb_info.nonzero_midx34s); % the corresponding linear indices of the 3rd-dimensional "num_nonzero34" above
                end
            end
            Ra_mn = At_wNoise(:, :, SRa_info.nonzero_midx34s(1,:)).*conj(At_wNoise(:, :, SRa_info.nonzero_midx34s(2,:))); % (Nt,M+1,num_nonzero34)
            if ~sim.scalar
                Rb_mn = At_wNoise(:, :, SRb_info.nonzero_midx34s(1,:)).*conj(At_wNoise(:, :, SRb_info.nonzero_midx34s(2,:))); % (Nt,M+1,num_nonzero34)
            end
        end

        % Then calculate Kerr, Ra, and Rb.
        for midx1 = 1:num_modes
            % Kerr
            nz_midx1 = find( SK_info.nonzero_midx1234s(1,:)==midx1 );
            midx2 = SK_info.nonzero_midx1234s(2,nz_midx1);
            midx3 = SK_info.nonzero_midx1234s(3,nz_midx1);
            midx4 = SK_info.nonzero_midx1234s(4,nz_midx1);
            Kerr(:,:,midx1) = sum(permute(SK_info.SK(nz_midx1),[3 2 1]).*At_wNoise(:, :, midx2).*At_wNoise(:, :, midx3).*conj(At_wNoise(:, :, midx4)),3);
            if sim.include_Raman
                % Ra
                for midx2 = 1:num_modes
                    nz_midx1 = find( SRa_info.nonzero_midx1234s(1,:)==midx1 );
                    nz_midx = nz_midx1( SRa_info.nonzero_midx1234s(2,nz_midx1)==midx2 ); % all the [midx1;midx2;?;?]
                    midx3 = SRa_info.nonzero_midx1234s(3,nz_midx);
                    midx4 = SRa_info.nonzero_midx1234s(4,nz_midx);
                    idx = midx34s_sub2ind([midx3;midx4]); % the linear indices
                    idx = arrayfun(@(i) find(SRa_nonzero_midx34s==i,1), idx); % the indices connecting to the 3rd-dimensional "num_nonzero34" of Ra_mn
                    Ra(:, :, midx1, midx2) = sum(permute(SRa_info.SRa(nz_midx),[3 2 1]).*Ra_mn(:, :, idx),3);
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
                        Rb(:, :, midx1, midx2) = sum(permute(SRb_info.SRb(nz_midx),[3 2 1]).*Rb_mn(:, :, idx),3);
                    end
                end
            end
        end
        if sim.include_Raman
            clear Ra_mn;
            if ~sim.scalar
                clear Rb_mn;
            end
        end
    end

    % Apply the convolution for each part of the Raman sum:
    % The convolution using Fourier transform is faster if both arrays are
    % large. If one of the array is small, "conv" can be faster.
    % Please refer to
    % "https://blogs.mathworks.com/steve/2009/11/03/the-conv-function-and-implementation-tradeoffs/"
    % for more information.
    if sim.include_Raman
        switch gas.model
            case 0
                Ra_upsampling = ifft(Ra,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % zero-padding for the acyclic convolution theorem to avoid time-domain aliasing
                Ra_upsampling = cat(1,Ra_upsampling(end-(gas_eqn.m2-1):end,:,:,:),Ra_upsampling(1:gas_eqn.m,:,:,:)); % only the (2*Nt-1) points contain useful information; others are from the zero-padding result

                RaAA = fft(Raw.*Ra_upsampling,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1);
                RaAA = RaAA(gas_eqn.R_downsampling,:,:,:).*gas_eqn.phase_shift_acyclic; % select only the valid part
                
                if ~sim.scalar % polarized fields with an anisotropic Raman contribution from rotational Raman
                    Rb_upsampling = ifft(Rb,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1); % zero-padding for the acyclic convolution theorem to avoid time-domain aliasing
                    Rb_upsampling = cat(1,Rb_upsampling(end-(gas_eqn.m2-1):end,:,:,:),Rb_upsampling(1:gas_eqn.m,:,:,:)); % only the (2*Nt-1) points contain useful information; others are from the zero-padding result

                    RbAA = fft(Rbw.*Rb_upsampling,gas_eqn.acyclic_conv_stretch(gas_eqn.Nt),1);
                    RbAA = RbAA(gas_eqn.R_downsampling,:,:,:).*gas_eqn.phase_shift_acyclic; % select only the valid part
                end
            case 1
                Ra = ifft(Ra,[],1); % transform into the frequency domain for upsampling
                Ra_upsampling = cat(1,Ra(end-(gas_eqn.m2-1):end,:,:,:),Ra(1:gas_eqn.n,:,:,:));

                RaAA = fft(Raw.*Ra_upsampling,gas_eqn.Nt,1).*gas_eqn.phase_shift;
                
                if ~sim.scalar % polarized fields with an anisotropic Raman contribution from rotational Raman
                    Rb = ifft(Rb,[],1); % transform into the frequency domain for upsampling
                    Rb_upsampling = cat(1,Rb(end-(gas_eqn.m2-1):end,:,:,:),Rb(1:gas_eqn.n,:,:,:));

                    RbAA = fft(Rbw.*Rb_upsampling,gas_eqn.Nt,1).*gas_eqn.phase_shift;
                end
        end

        if sim.scalar
            nonlinear = prefactor{2}.*ifft(Kerr,[],1) + ifft(sum(RaAA.*permute(At_wNoise,[1 2 4 3]),4),[],1);
        else % polarized fields with an anisotropic Raman contribution from rotational Raman
            nonlinear = prefactor{2}.*ifft(Kerr,[],1) + ifft(sum((RaAA+RbAA).*permute(At_wNoise,[1 2 4 3]),4),[],1);
        end
    else
        nonlinear = prefactor{2}.*ifft(Kerr,[],1);
    end
    
    % Finalize by multiplying with the prefactor
    nonlinear = eta*prefactor{1}.*nonlinear;

    % Incorporate dz and D.neg term for the integration.
    nonlinear = permute(nonlinear,[1,3,2]); % (Nt,num_modes,M+1)
    nonlinear = small_dz*D.neg.*nonlinear; % (Nt, num_modes, M+1)

    % Save the previous psi at the right end, then compute the new psi's
    last_psi_half_parallelization = current_psi_half_parallelization;
    % ---------------------------------------------------------------------
    % "psi" is calculated from trapezoidal integrals, which results in
    % (1/2,1,1,......,1,1/2) coefficients.
    % Note that for each z plane, it has a (1/2, 1, 1...,1/2) factor in "nonlinear".
    %nonlinear(:,:,1) = nonlinear(:,:,1)/2;
    %psi(:,:,2:end) = psi(:,:,1) + cumsum(nonlinear(:,:,1:sim.MPA.M),3) + nonlinear(:,:,2:end)/2;
    % ---------------------------------------------------------------------
    % "psi" is calculated with a multistep Adams-Moulton method for each parallelization plane.
    % A M-step Adams-Moulton method (with the summation of M's derivatives) has a (M+1)th-order accuracy.
    % With it, the latter psi will not have a bigger error.
    % The error builds up for latter parallization planes in the typical trapezoidal quadrature.
    % The highest-order accuracy is 2nd order O(h^2) from
    %   y(n+1)=y(n)+1/2*h*( f(t(n+1),y(n+1)) + f(t(n),y(n)) ),
    % so the error estimate of adaptive-step-size control uses cubic root: error = O(h^3).
    if sim.gpu_yes
        psi = feval(sim.cuda_MPA_psi_update,psi,...
                                            nonlinear,...
                                            sim.MPA.coeff,...
                                            gas_eqn.Nt,sim.MPA.M+1,num_modes);
    else
        for Midx = 1:sim.MPA.M
            psi(:,:,Midx+1) = psi(:,:,1) + sum(nonlinear(:,:,1:Midx+1).*sim.MPA.coeff(Midx,:,1:Midx+1),3);
        end
    end

    % Calculate the average NRMSE = take the RMSE between the previous psi
    % and the current psi at the right edge, normalize by the absolute max,
    % and average over all modes
    current_psi_half_parallelization = psi(:,:,1) + 2*sum(nonlinear(:,:,1:2:end).*sim.MPA.coeff(sim.MPA.M/2,:,1:sim.MPA.M/2+1),3);
    energy_current_psi_half_parallelization = sum(abs(current_psi_half_parallelization).^2);
    weight = energy_current_psi_half_parallelization/sum(energy_current_psi_half_parallelization);
    NRMSE_p = sqrt(sum(abs(current_psi_half_parallelization-last_psi_half_parallelization).^2)./energy_current_psi_half_parallelization).*weight;
    NRMSE_p(isnan(NRMSE_p)) = 0; % exclude modes with all zero fields
    avg_NRMSE = sum(NRMSE_p);

    % If it has converged to within tol, then quit
    if avg_NRMSE < sim.MPA.tol && n_it >= sim.MPA.n_tot_min
        break
    end

    if n_it == sim.MPA.n_tot_max
        error('Error in stepping_MPA_gradient_pressure: The step did not converge after %d iterations, aborting.', sim.MPA.n_tot_max);
    end
end

% Get back to A from psi
psi_full_parallelization = psi(:,:,sim.MPA.M+1);
A1w = D.pos(:,:,sim.MPA.M+1).*psi_full_parallelization;

% Stepsize control
err = sqrt(sum(abs(current_psi_half_parallelization-psi_full_parallelization).^2)./energy_current_psi_half_parallelization).*weight;
err(isnan(err)) = 0; % exclude modes with all zero fields
err = sum(err);
if isnan(err) % the computation is just so wrong, so we reduce the step size and do it again
    opt_dz = 0.5*sim.dz;
    success = false;
else
    opt_dz = max(0.5,min(2,0.8*(sim.adaptive_dz.threshold/err)^(1/3)))*sim.dz; % limit the optimal step size to 0.5~2 times as the current step size

    success = err < sim.adaptive_dz.threshold;
end

% Downsample them back
A1w = cat(1,A1w(1:gas_eqn.n,:),A1w(end-(Nt-gas_eqn.n-1):end,:));

end
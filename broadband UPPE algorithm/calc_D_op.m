function [D_op,sim] = calc_D_op(fiber,sim,dt,Omega,fields)
%CALC_D_OP It computes the dispersion operator used in UPPE

% If the code applies narrowband transformation to the coherent fields,
% fields are downsampled in UPPE_propagate() already. If
% fiber.betas is a function of frequency, it will then contain too many
% points. Its downsampling operation is applied below.
if sim.cs.cs > 1
    fiber.betas = fiber.betas(1:sim.cs.cs:end,:);
end

if ~isfield(sim,'betas')
    if ~any(fields) % fields is all-zero
        sim.betas = [0;0];
    else
        % Obtain the betas of the input pulse
        sim.betas = calc_sim_betas(fiber,dt,Omega,ifft(fields,[],1));
    end
end

D_op = 1i*(ifftshift(fiber.betas,1)-(sim.betas(1)+sim.betas(2)*Omega));

% Scaled dispersion due to the narrowband transformation (scaled Fourier transform)
if sim.cs.cs > 1
    D_op = D_op/sim.cs.cs;
end

end
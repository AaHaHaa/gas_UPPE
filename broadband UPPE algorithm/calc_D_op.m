function [D_op,sim] = calc_D_op(fiber,sim,dt,omegas,fields)
%CALC_D_OP It computes the dispersion operator used in UPPE

if ~isfield(sim,'betas')
    if ~any(fields) % fields is all-zero
        sim.betas = [0;0];
    else
        % Obtain the betas of the input pulse
        sim.betas = calc_sim_betas(fiber,dt,omegas,fields);
    end
end

D_op = 1i*(ifftshift(fiber.betas,1)-(sim.betas(1)+sim.betas(2)*omegas));

end


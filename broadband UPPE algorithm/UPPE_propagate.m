function foutput = UPPE_propagate(fiber, initial_condition, sim, gas)
%UPPE_PROPAGATE Propagate an initial multimode pulse through an arbitrary distance of a hollow-core fiber
%   This is a caller function for UPPE_propagate_constant_pressure() and UPPE_propagate_gradient_pressure().

%%
if ispc
    sep_char = '\';
else % unix
    sep_char = '/';
end

% Load the folder
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep_char);
current_folder = current_path(1:sep_pos(end));

sim.cuda_dir_path = [current_folder '../cuda'];

%% Apply the narrowband transformation (due to the scaled Fourier transform)
scaledFT_func = narrowband_scaledFT();

if sim.cs.cs > 1
    Nt = size(initial_condition.fields,1);
    num_modes = size(initial_condition.fields,2);
    transformed_At = zeros(round(Nt/sim.cs.cs),num_modes); % I use round() here to prevent error. The error check for "cs" will be done in scaledFT_func.convert later, so we just let it pass here.
    for midx = 1:num_modes
        transformed_At(:,midx) = scaledFT_func.convert(initial_condition.fields(:,midx,end),sim.cs.cs);
    end

    initial_condition.dt = initial_condition.dt*sim.cs.cs;
    initial_condition.fields = transformed_At;
end

%% Run the pulse propagation
if sim.include_heating % consider Raman/photoionization heating during propagation
    % Some pre-checks should be done when calling gas_info() beforehand, so
    % I don't repeat them here.
    foutput = UPPE_propagate_constant_pressure_heating(fiber, initial_condition, sim, gas);
else
    if isfield(gas,'pressure_in') && isfield(gas,'pressure_out')
        foutput = UPPE_propagate_gradient_pressure(fiber, initial_condition, sim, gas);
    elseif isfield(gas,'pressure') % constant gas pressure
        foutput = UPPE_propagate_constant_pressure(fiber, initial_condition, sim, gas);
    else
        error('UPPE_propagate:gasPressureError',...
              '"gas" needs to have "pressure_in" and "pressure_out" for a gradient pressure or "pressure" only for a constrant pressure.');
    end
end

%% Recover from the narrowband transformation
if sim.cs.cs > 1
    foutput.dt = foutput.dt/sim.cs.cs;
    foutput.fields = scaledFT_recover_field(scaledFT_func,foutput.fields,sim.cs.cs);
end

end

%% Helper functions for recovering the field from the narrowband transformation due to the scaled Fourier transform
function rA = scaledFT_recover_field(scaledFT_func,A,cs)

Nt = size(A,1);
num_modes = size(A,2);
Nz = size(A,3);

rA = zeros(Nt*cs,num_modes,Nz);
for zi = 1:Nz
    for midx = 1:num_modes
        rA(:,midx,zi) = scaledFT_func.recover(A(:,midx,zi),cs);
    end
end

end
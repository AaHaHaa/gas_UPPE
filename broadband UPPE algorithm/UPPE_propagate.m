function foutput = UPPE_propagate(fiber, initial_condition, sim, gas)
%UPPE_PROPAGATE Propagate an initial multimode pulse through an arbitrary distance of a hollow-core fiber
%   This is a caller function for UPPE_propagate_constant_pressure().

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

%% Run the pulse propagation
if isfield(gas,'pressure_in') && isfield(gas,'pressure_out')
    foutput = UPPE_propagate_gradient_pressure(fiber, initial_condition, sim, gas);
elseif isfield(gas,'pressure') % constant gas pressure
    foutput = UPPE_propagate_constant_pressure(fiber, initial_condition, sim, gas);
else
    error('UPPE_propagate:gasPressureError',...
          '"gas" needs to have "pressure_in" and "pressure_out" for a gradient pressure or "pressure" only for a constrant pressure.');
end

end
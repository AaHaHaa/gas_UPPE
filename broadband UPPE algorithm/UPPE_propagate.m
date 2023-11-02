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
foutput = UPPE_propagate_constant_pressure(fiber, initial_condition, sim, gas);

end
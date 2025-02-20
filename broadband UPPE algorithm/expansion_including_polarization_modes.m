function [sim,fiber] = expansion_including_polarization_modes(sim,fiber,num_modes)
%EXPANSION_INCLUDING_POLARIZATION_MODES It extends betas into 2*num_spatial_modes if necessary.

num_modes_betas = size(fiber.betas,2);
if isfield(fiber,'X3') % constant pressure
    num_modes_X3 = size(fiber.X3,2);
else % gradient pressure
    num_modes_X3 = size(fiber.X3_prefactor,2);
end

if ~sim.scalar
    % betas
    betas = complex(zeros(size(fiber.betas,1),num_modes));
    if num_modes_betas == num_modes/2 % num_modes = 2*num_spatial_modes
        betas(:,2:2:num_modes) = fiber.betas;
        betas(:,1:2:num_modes-1) = fiber.betas;
        fiber.betas = betas;
    end

    % X3
    X3 = complex(zeros(size(fiber.X3,1),num_modes));
    if num_modes_X3 == num_modes/2 % num_modes = 2*num_spatial_modes
        X3(:,2:2:num_modes) = fiber.X3;
        X3(:,1:2:num_modes-1) = fiber.X3;
        fiber.X3 = X3;
    end
    
    % mode profiles
    mode_profiles = complex(zeros([size(sim.mode_profiles.mode_profiles,[1,2]),num_modes]));
    if num_modes_X3 == num_modes/2 % num_modes = 2*num_spatial_modes
        mode_profiles(:,:,2:2:num_modes) = sim.mode_profiles.mode_profiles;
        mode_profiles(:,:,1:2:num_modes-1) = sim.mode_profiles.mode_profiles;
        sim.mode_profiles.mode_profiles = mode_profiles;
    end
    
    % mode profile norms
    mode_profile_norms = complex(zeros(size(sim.mode_profiles.norms,1),num_modes));
    if num_modes_X3 == num_modes/2 % num_modes = 2*num_spatial_modes
        mode_profile_norms(:,2:2:num_modes) = sim.mode_profiles.norms;
        mode_profile_norms(:,1:2:num_modes-1) = sim.mode_profiles.norms;
        sim.mode_profiles.norms = mode_profile_norms;
    end
end

end
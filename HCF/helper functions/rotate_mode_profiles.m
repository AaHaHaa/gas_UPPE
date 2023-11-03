function mode_profiles = rotate_mode_profiles(mode_profiles, rot, dtheta)
%ROTATE_MODE_PROFILES Summary of this function goes here
%   Detailed explanation goes here

theta_shift = floor(rot/dtheta);

% rotate the coordinates
mode_profiles = circshift(mode_profiles,theta_shift,4);

% rotate the vector
tmp = mode_profiles(:,:,:,:,:,:);
mode_profiles(:,:,:,:,1,:) = tmp(:,:,:,:,1,:)*cos(rot) - tmp(:,:,:,:,2,:)*sin(rot);
mode_profiles(:,:,:,:,2,:) = tmp(:,:,:,:,1,:)*sin(rot) + tmp(:,:,:,:,2,:)*cos(rot);

end
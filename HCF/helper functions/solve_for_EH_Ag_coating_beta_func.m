function [beta,SR,mode_profile] = solve_for_EH_Ag_coating_beta_func(wavelength,eta,sim,gas)
%   
%   wavelength: wavelength points in simulations (m)
%   eta: gas density (amagats)
%   gas: core_radius: the core radius of the hollow-core fiber (m)
%        wavelength: wavelengths to calculate the propagation constant and the loss
%                    It contains two fields, range and num
%        order: make sure that the propagation constant have smooth 
%               higher-order derivatives up this order
%        gas_material: 'H2','air','N2','O2','Ar','Xe','Kr','Ne','He'

% The order of modes, based on propagation constant (beta), is found by
% running "find_order_of_EH_modes" first.

if ispc
    sep_char = '\';
else % unix
    sep_char = '/';
end

% Load the folder
current_path = mfilename('fullpath');
sep_pos = strfind(current_path,sep_char);
current_folder = current_path(1:sep_pos(end));
upper_folder = current_path(1:sep_pos(end-2));
addpath(current_folder);
addpath(fullfile(upper_folder,'gas absorption spectra'));

Nf = uint32(length(wavelength));

c = 299792.458; % nm/ps
gas_wavelength = wavelength; % m
gas_f = c./gas_wavelength*1e-9; % THz

%% Refractive index

% Reference:
% 1. Walter G., et el, "On the Dependence of the Refractive Index of Gases on Temperature" (1903)
% 2. Arthur L. Ruoff and Kouros Ghandehari, "THE REFRACTIVE INDEX OF HYDROGEN AS A FUNCTION OF PRESSURE" (1993)

[a,b] = Sellmeier_coefficients(gas.gas_material); % Sellmeier coefficients
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
switch gas.gas_material
    case 'H2'
        %n_gas = calc_n_H2(gas_wavelength*1e9,sim.cuda_dir_path,gas.wavelength_order);
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        permittivity_r = n_from_Sellmeier(gas_wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        permittivity_r = n_from_Sellmeier(0.120).^2;
        n_gas_120 = sqrt((permittivity_r - 1)*eta + 1); %  % Sellmeier is valid only above 164nm
        n_gas(gas_wavelength<120e-9) = n_gas_120;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas.gas_material,gas_wavelength,eta);
        n_gas = n_gas + 1i*Raman_absorption./(2*pi./gas_wavelength);
    case 'O2'
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        permittivity_r = n_from_Sellmeier(gas_wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        permittivity_r = n_from_Sellmeier(0.4).^2;
        n_gas_400 = sqrt((permittivity_r - 1)*eta + 1); %  % Sellmeier is valid only above 400nm
        n_gas(gas_wavelength<120e-9) = n_gas_400;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas.gas_material,gas_wavelength,eta);
        n_gas = n_gas + 1i*Raman_absorption./(2*pi./gas_wavelength);
    case {'air','N2'}
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        permittivity_r = n_from_Sellmeier(gas_wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas.gas_material,gas_wavelength,eta);
        n_gas = n_gas + 1i*Raman_absorption./(2*pi./gas_wavelength);
    case {'Ar','Xe','Kr','Ne','He'}
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        permittivity_r = n_from_Sellmeier(gas_wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
    case 'CH4'
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        permittivity_r = n_from_Sellmeier(gas_wavelength*1e6).^2;
        n_gas = sqrt((permittivity_r - 1)*eta + 1); % refractive index of the gas
        
        % Avoid the singularity at resonances
        idx_resonance = n_gas < 1;
        n_gas(idx_resonance) = 1;
end

%%

% (n,m) modes to solve for
user_midx = [1,4,9,17,28,40]; % a maximum of 6 circular symmetric modes included here
user_midx = user_midx(sim.midx);
num_modes = length(user_midx);

n_out = calc_n_Ag(gas_wavelength*1e9,sim.gpu_yes,sim.cuda_dir_path,gas.wavelength_order);

%%
load(fullfile(current_folder,'nm_order.mat'),'sorted_nm');
nm = sorted_nm(:,user_midx);
% For EH modes,
% EH_(-|n|,m) is degenerate (same progation constant) with EH_(|n|+2,m)
nm = [nm,[2-nm(1,:);nm(2,:)]]; % [2-nm(1,:);nm(2,:)] are the degenerate modes

k0 = 2*pi./gas_wavelength;

ZdYd = zeros(length(wavelength),num_modes*2);
unm = zeros(1,num_modes*2);
mode = cell(1,num_modes*2);
for midx = 1:num_modes*2
    % Zd, Yd
    if nm(1,midx) == 0 % only TE is considered in this code
        mode{midx} = 'TE';
    else
        mode{midx} = 'EH';
    end
    
    ZdYd(:,midx) = calc_ZdYd(n_gas,n_out,mode{midx});
    
    % unm: zeros of the Bessel function of the first kind
    u = besselzero(nm(1,midx)-1,nm(2,midx),1);
    unm(midx) = u(end);
end

if sim.gpu_yes
    ZdYd = gpuArray(ZdYd);
    k0 = gpuArray(k0);
    unm = gpuArray(unm);
    n_gas = gpuArray(n_gas);
end

%ki = complex(unm/gas.core_radius);
%ki = unm./(1+1i*ZdYd./(k0.*n_gas)/gas.core_radius)/gas.core_radius;
ki = (1i*ZdYd.*nm(1,:)./k0/gas.core_radius-1).*unm./(1i*ZdYd./k0/gas.core_radius.*(nm(1,:)-1)-1)/gas.core_radius;
gamma = sqrt((k0.*n_gas).^2 - ki.^2);

%% Propagation constant
chosen_midx = [];
next_midx = 1;
for midx = 1:num_modes
    if nm(1,midx) == 1
        acum_midx = next_midx + 0;
    else
        acum_midx = next_midx + [0,1];
    end
    next_midx = next_midx + 2;
    chosen_midx = [chosen_midx,acum_midx]; %#ok
end

if sim.gpu_yes
    beta = complex(zeros(Nf,num_modes*2,'gpuArray'));
else
    beta = complex(zeros(Nf,num_modes*2));
end
beta(:,1:2:end) = gamma(:,1:num_modes);
beta(:,2:2:end) = gamma(:,1:num_modes);
beta = beta(:,chosen_midx);
if sim.gpu_yes
    beta = gather(beta); % gather the output back from GPU
end

%% Spatial profile
r = permute(linspace(gas.core_radius/gas.xy_sampling*1e-3,gas.core_radius,gas.xy_sampling),[1,3,2]);
theta = permute(linspace(0,2*pi,gas.xy_sampling+1),[1,3,4,2]);  theta = theta(1:end-1); % 0 and 2*pi are the same
dr = gas.core_radius/gas.xy_sampling;
dtheta = 2*pi/gas.xy_sampling;

num_polarized = 2; % two polarizations (r,theta)
mode_profiles = complex(zeros(1,num_modes*2,gas.xy_sampling,gas.xy_sampling,num_polarized));

%target_wavelength = wavelength(floor(length(wavelength)/2)+1); % m
target_wavelength = gas.mode_profile_wavelength; % m

abs_n_gas = interp1(gas_f,abs(n_gas),c*1e-9/target_wavelength);
ang_n_gas = interp1(gas_f,unwrap(angle(n_gas),[],1),c*1e-9/target_wavelength);
target_n_gas = abs_n_gas.*exp(1i*ang_n_gas);

abs_ki = interp1(gas_f,abs(ki),c*1e-9/target_wavelength);
ang_ki = interp1(gas_f,unwrap(angle(ki),[],1),c*1e-9/target_wavelength);
target_ki = abs_ki.*exp(1i*ang_ki);

target_n_out = interp1(gas_f,n_out,c*1e-9/target_wavelength); % for real n_out
target_k0 = interp1(gas_f,k0,c*1e-9/target_wavelength);

% "GPU besselj" doesn't allow complex double" input, so I need to gather 
% the data back from GPU.
if sim.gpu_yes
    mode_profiles = gpuArray(mode_profiles);
    target_ki = gather(target_ki);
    target_k0 = gather(target_k0);
    target_n_gas = gather(target_n_gas);
    target_n_out = gather(target_n_out);
    unm = gather(unm);
    theta = gather(theta);
end
% (r,theta) basis from Marcatili's paper
for midx = 1:num_modes*2
    switch mode{midx}
        case 'TE'
            mode_profiles(:,midx,:,:,2) = repmat(besselj(1,target_ki(:,midx).*r),[1,1,1,gas.xy_sampling,1]);
        case 'EH'
            theta0 = 0;
            Dbesselj = @(n,z) -besselj(n+1,z) + n./z.*besselj(n,z);
            mode_profiles(:,midx,:,:,1) = (besselj(nm(1,midx)-1,target_ki(:,midx).*r)+1i*unm(midx)./(2*target_k0.*target_n_gas.*r).*sqrt(target_n_out.^2-target_n_gas.^2).*besselj(nm(1,midx),target_ki(:,midx).*r)).*sin(nm(1,midx)*(theta+theta0));
            mode_profiles(:,midx,:,:,2) = (besselj(nm(1,midx)-1,target_ki(:,midx).*r)+1i*unm(midx)^2./(2*nm(1,midx)*target_k0.*target_n_gas*gas.core_radius).*sqrt(target_n_out.^2-target_n_gas.^2).*Dbesselj(nm(1,midx),target_ki(:,midx).*r)).*cos(nm(1,midx)*(theta+theta0));
    end
end
norm = sqrt(sum(sum(sum(abs(mode_profiles).^2,5).*r*dr*dtheta,3),4));
nonzero_idx = norm~=0;
mode_profiles(nonzero_idx) = mode_profiles(nonzero_idx)./norm(nonzero_idx);

% Transform into (x,y) basis
if sim.gpu_yes
    mode_profiles_xy = zeros(size(mode_profiles),'gpuArray');
else
    mode_profiles_xy = zeros(size(mode_profiles));
end
mode_profiles_xy(:,:,:,:,1) = mode_profiles(:,:,:,:,1).*cos(theta) - mode_profiles(:,:,:,:,2).*sin(theta);
mode_profiles_xy(:,:,:,:,2) = mode_profiles(:,:,:,:,1).*sin(theta) + mode_profiles(:,:,:,:,2).*cos(theta);

%% Form a basis of linear polarizations (scalar modes)
% Arrange degenerate modes into the 6th dimension
linear_mode_profiles = 1/sqrt(2)*(mode_profiles_xy(:,1:num_modes,:,:,:) + mode_profiles_xy(:,num_modes+1:num_modes*2,:,:,:));
linear_mode_profiles = cat(6,linear_mode_profiles,...
                       1/sqrt(2)*(mode_profiles_xy(:,1:num_modes,:,:,:) - mode_profiles_xy(:,num_modes+1:num_modes*2,:,:,:))...
                       ); % degenerate mode

% Because it may not align with x or y axis, I rotate the profiles.
xy_contribution = sum(sum(abs(linear_mode_profiles(1,:,:,:,:,:)).^2.*r*dr*dtheta,3),4);
xy_ratio = xy_contribution(:,:,:,:,1,:)./xy_contribution(:,:,:,:,2,:);
degenerate_mode_ratio = xy_ratio(:,:,:,:,:,1)./xy_ratio(:,:,:,:,:,2);
tol = 0.1;
rot = pi/4;
for midx = 1:num_modes
    if abs(degenerate_mode_ratio(:,midx,:,:) - 1) < tol
        linear_mode_profiles(:,midx,:,:,:,:) = rotate_mode_profiles(linear_mode_profiles(:,midx,:,:,:,:), rot, dtheta);
    end
end

% Now rotate the polarization only to make it x-polarized
xy_contribution = sum(sum(abs(linear_mode_profiles(1,:,:,:,:,:)).^2.*r*dr*dtheta,3),4);
xy_ratio = xy_contribution(:,:,:,:,1,:)./xy_contribution(:,:,:,:,2,:);
for di = 1:2
    for midx = 1:num_modes
        if xy_ratio(:,midx,:,:,:,di) < 1
            linear_mode_profiles(:,midx,:,:,1,di) = -conj(linear_mode_profiles(:,midx,:,:,2,di));
            linear_mode_profiles(:,midx,:,:,2,di) =  conj(linear_mode_profiles(:,midx,:,:,1,di));
        end
    end
end

% Discard the y-polarized components and eliminate the polarization dimension
linear_mode_profiles = -permute(linear_mode_profiles(:,:,:,:,1,:),[1,2,3,4,6,5]);

%% Refine the information
mode_profiles = complex(zeros(1,num_modes*2,gas.xy_sampling,gas.xy_sampling));
if sim.gpu_yes
    mode_profiles = gpuArray(mode_profiles);
end
mode_profiles(:,1:2:end,:,:) = linear_mode_profiles(:,:,:,:,1); % the 1st degenerate mode
mode_profiles(:,2:2:end,:,:) = linear_mode_profiles(:,:,:,:,2); % the 2nd degenerate mode

mode_profiles = mode_profiles(:,chosen_midx,:,:);

% Orthogonality
% Make the modes orthogonal to each other
norm = sqrt(sum(sum(abs(mode_profiles(:,1,:,:)).^2.*r*dr*dtheta,3),4));
mode_profiles(:,1,:,:) = mode_profiles(:,1,:,:)./norm;
for midx1 = 2:length(chosen_midx)
    for midx2 = 1:(midx1-1)
        mode_profiles(:,midx1,:,:) = mode_profiles(:,midx1,:,:) - sum(sum(mode_profiles(:,midx1,:,:).*conj(mode_profiles(:,midx2,:,:)).*r*dr*dtheta,3),4).*mode_profiles(:,midx2,:,:);
    end
    norm = sqrt(sum(sum(abs(mode_profiles(:,midx1,:,:)).^2.*r*dr*dtheta,3),4));
    mode_profiles(:,midx1,:,:) = mode_profiles(:,midx1,:,:)./norm;
end

%% Mode-field computation: QK, QRa
permittivity0 = 8.8541878176e-12; % F/m
norm_spatial_modes = sqrt(permittivity0*real(beta./(2*pi./wavelength))*(c*1e3)/2);

SR = calc_SR_tensors(mode_profiles,r,dr,dtheta,sim,current_folder);

if sim.gpu_yes
    mode_profiles = gather(mode_profiles);
    norm_spatial_modes = gather(norm_spatial_modes);
end
mode_profile = struct('mode_profiles',permute(mode_profiles,[3,4,2,1]),'norms',norm_spatial_modes,...
                      'r',r,'dr',dr,'dtheta',dtheta);

end
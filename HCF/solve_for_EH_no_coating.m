% This code solves for the propagation constant of a silica capillary.
% The loss of the silica is considered.
%
% The order of modes, based on propagation constant (beta), is found by
% running "find_order_of_EH_modes" first.
%
% To calculate the spatial profile, I follow Marcitili-Schmeltzer model
% (check their paper) which gives EH modes. I did some calculation to
% transform them into LP modes.

clearvars; close all;

addpath('helper functions','../gas absorption spectra/');

use_gpu = false;%true; % GPU

gas_material = 'N2';

pressure = 0; % atm
temperature = 273.15 + 25; % 15 degree Celsius
core_radius = 150e-6; % core radius; m

% Don't change "wavelength"!
num_wavelength = 1e5;
wavelength_range = [100,30e3]; % nm
c = 299792.458; % nm/ps
f = linspace(c/wavelength_range(1),c/wavelength_range(2),num_wavelength)'; % THz
wavelength = c./f; % nm
target_wavelength = 1030; % nm

% the number of sampling points of the spatial fields
r_sampling = 101;
theta_sampling = 101;

% (n,m) modes to solve for
% [1,4,9,17,28,40] are the first six radial EH0m modes
user_midx = 1;%[1,4,9,17,28,40];
num_modes = length(user_midx);

% refractive index
% Reference:
% 1. Walter G., et el, "On the Dependence of the Refractive Index of Gases on Temperature" (1903)
% 2. Arthur L. Ruoff and Kouros Ghandehari, "THE REFRACTIVE INDEX OF HYDROGEN AS A FUNCTION OF PRESSURE" (1993)
pressure = pressure*1.01325e5; % Pa
pressure0 = 1.01325e5; % Pa
temperature0 = 273.15; % 0 degree Celsius
[a,b] = Sellmeier_coefficients(gas_material); % Sellmeier coefficients
Sellmeier_terms = @(lambda,a,b) a.*lambda.^2./(lambda.^2 - b.^2);
cuda_dir_path = '../cuda';
diff_order = 6;
switch gas_material
    case 'H2'
        %n_gas = calc_n_H2(wavelength,cuda_dir_path,diff_order);
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        permittivity_r = n_from_Sellmeier(wavelength*1e-3).^2;
        n_gas = sqrt((permittivity_r - 1)*(pressure/temperature)/(pressure0/temperature0) + 1); % refractive index of the gas
        
        permittivity_r = n_from_Sellmeier(0.12).^2;
        n_gas_120 = sqrt((permittivity_r - 1)*(pressure/temperature)/(pressure0/temperature0) + 1); % Sellmeier is valid only above 164nm
        n_gas(wavelength<120) = n_gas_120;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas_material,wavelength*1e-9,(pressure/temperature)/(pressure0/temperature0));
        n_gas = n_gas + 1i*Raman_absorption./(2*pi./(wavelength*1e-9));
    case 'O2'
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        permittivity_r = n_from_Sellmeier(wavelength*1e-3).^2;
        n_gas = sqrt((permittivity_r - 1)*(pressure/temperature)/(pressure0/temperature0) + 1); % refractive index of the gas
        
        permittivity_r = n_from_Sellmeier(0.4).^2;
        n_gas_400 = sqrt((permittivity_r - 1)*(pressure/temperature)/(pressure0/temperature0) + 1); % Sellmeier is valid only above 400nm
        n_gas(wavelength<400) = n_gas_400;
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas_material,wavelength*1e-9,(pressure/temperature)/(pressure0/temperature0));
        n_gas = n_gas + 1i*Raman_absorption./(2*pi./(wavelength*1e-9));
    case {'air','N2'}
        n_from_Sellmeier = @(lambda) sum(Sellmeier_terms(lambda,a,b),2) + 1;
        
        permittivity_r = n_from_Sellmeier(wavelength*1e-3).^2;
        n_gas = sqrt((permittivity_r - 1)*(pressure/temperature)/(pressure0/temperature0) + 1); % refractive index of the gas
        
        % pressure-induced absorption
        Raman_absorption = read_absorption(gas_material,wavelength*1e-9,(pressure/temperature)/(pressure0/temperature0));
        n_gas = n_gas + 1i*Raman_absorption./(2*pi./(wavelength*1e-9));
    case {'Ar','Ne','He'}
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        permittivity_r = n_from_Sellmeier(wavelength*1e-3).^2;
        n_gas = sqrt((permittivity_r - 1)*(pressure/temperature)/(pressure0/temperature0) + 1); % refractive index of the gas
    case {'Xe','Kr'}
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        permittivity_r = n_from_Sellmeier(wavelength*1e-3).^2;
        n_gas = sqrt((permittivity_r - 1)*(pressure/temperature)/(pressure0/temperature0) + 1); % refractive index of the gas
        
        permittivity_r = n_from_Sellmeier(0.15).^2;
        n_gas_150 = sqrt((permittivity_r - 1)*(pressure/temperature)/(pressure0/temperature0) + 1); % Sellmeier is valid only above ~150nm
        n_gas(wavelength<150) = n_gas_150;
    case 'CH4'
        n_from_Sellmeier = @(lambda) sqrt(1+sum(Sellmeier_terms(lambda,a,b),2));
        
        permittivity_r = n_from_Sellmeier(wavelength*1e-3).^2;
        n_gas = sqrt((permittivity_r - 1)*(pressure/temperature)/(pressure0/temperature0) + 1); % refractive index of the gas
        
        % Avoid the singularity at resonances
        idx_resonance = n_gas < 1;
        n_gas(idx_resonance) = 1;
end

n_out = calc_n_silica(wavelength,use_gpu,cuda_dir_path,diff_order);

if abs(round(pressure/1.01325e5)-pressure/1.01325e5) < eps(1) % pressure is an integer number of atm
    saved_filename = sprintf('info_no_coating_%s_%dum_%datm.mat',gas_material,core_radius*2*1e6,pressure/1.01325e5);
else
    saved_filename = sprintf('info_no_coating_%s_%dum_%.1fatm.mat',gas_material,core_radius*2*1e6,pressure/1.01325e5);
end

%%
% the mode index
load('nm_order.mat');
nm = sorted_nm(:,user_midx);
% For EH modes,
% EH_(-|n|,m) is degenerate (same progation constant) with EH_(|n|+2,m)
nm = [nm,[2-nm(1,:);nm(2,:)]]; % [2-nm(1,:);nm(2,:)] are the degenerate modes

wavelength = wavelength*1e-9; % m
target_wavelength = target_wavelength*1e-9; % m
if size(wavelength,1) == 1 % make it column vector
    wavelength = wavelength.';
end
wavelength_sampling = length(wavelength);
k0 = 2*pi./wavelength;

r = permute(linspace(core_radius/r_sampling*1e-3,core_radius,r_sampling),[1,3,2]);
theta = permute(linspace(0,2*pi,theta_sampling+1),[1,3,4,2]);  theta = theta(1:end-1); % 0 and 2*pi are the same
dr = core_radius/r_sampling;
dtheta = 2*pi/theta_sampling;

ZdYd = zeros(wavelength_sampling,num_modes*2);
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

if use_gpu
    ZdYd = gpuArray(ZdYd);
    k0 = gpuArray(k0);
    unm = gpuArray(unm);
    n_gas = gpuArray(n_gas);
    wavelength = gpuArray(wavelength);
end

target_wavelength_sampling = length(target_wavelength);
%ki = complex(unm/core_radius).*ones(size(k0));
%ki = unm./(1+1i*ZdYd./(k0.*n_gas)/core_radius)/core_radius;
ki = (1i*ZdYd.*nm(1,:)./k0/core_radius-1).*unm./(1i*ZdYd./k0/core_radius.*(nm(1,:)-1)-1)/core_radius;
gamma = sqrt((k0.*n_gas).^2 - ki.^2);

% Check validity of this Marcatili's equation
figure;
h = plot(wavelength*1e6,abs(n_out./n_gas).*unm./(k0.*core_radius));
set(h,'linewidth',2); set(gca,'fontsize',18);
xlabel('Wavelength (\mum)');
ylabel('|n_{out}^{relative}|u_{nm}/(k_0a)');
title('Validaity check 1 (<<1)');

figure;
h = plot(wavelength*1e6,abs(gamma./k0-1));
set(h,'linewidth',2); set(gca,'fontsize',18);
xlabel('Wavelength (\mum)');
ylabel('|\gamma/k_0-1|');
title('Validaity check 2 (<<1)');

%% Spatial profiles below
% Interpolate to the target wavelength
if use_gpu
    cuda_interp1_D2 = setup_kernel('interp1_D2',cuda_dir_path,target_wavelength_sampling*num_modes*2);
    target_ki = complex(zeros(target_wavelength_sampling,num_modes*2,'gpuArray'));
    target_ki = feval(cuda_interp1_D2,...
                                     wavelength,complex(ki),false,uint32(wavelength_sampling),...
                                     target_wavelength,target_ki,false,uint32(target_wavelength_sampling),uint32(length(target_ki)),...
                                     0.5);
    target_k0 = interp1(wavelength,k0,target_wavelength);
else
    abs_ki = interp1(wavelength,abs(ki),target_wavelength);
    ang_ki = interp1(wavelength,unwrap(angle(ki),[],1),target_wavelength);
    target_ki = abs_ki.*exp(1i*ang_ki);
    
    target_k0 = interp1(wavelength,k0,target_wavelength);
end
target_n_in = interp1(wavelength,n_gas,target_wavelength);
target_n_out = interp1(wavelength,n_out,target_wavelength);

num_polarized = 2; % two polarizations (r,theta)
mode_profiles = complex(zeros(target_wavelength_sampling,num_modes*2,r_sampling,theta_sampling,num_polarized));
% "GPU besselj" doesn't allow complex double" input, so I need to gather 
% the data back from GPU.
if use_gpu
    mode_profiles = gpuArray(mode_profiles);
    target_ki = gather(target_ki);
    target_k0 = gather(target_k0);
    target_n_in = gather(target_n_in);
    unm = gather(unm);
    theta = gather(theta);
end
% (r,theta) basis from Marcatili's paper
for midx = 1:num_modes*2
    switch mode{midx}
        case 'TE'
            mode_profiles(:,midx,:,:,2) = repmat(besselj(1,target_ki(:,midx).*r),[1,1,1,theta_sampling,1]);
        case 'EH'
            theta0 = 0;
            Dbesselj = @(n,z) -besselj(n+1,z) + n./z.*besselj(n,z);
            mode_profiles(:,midx,:,:,1) = (besselj(nm(1,midx)-1,target_ki(:,midx).*r)+1i*unm(midx)./(2*target_k0.*target_n_in.*r).*sqrt(target_n_out.^2-target_n_in.^2).*besselj(nm(1,midx),target_ki(:,midx).*r)).*sin(nm(1,midx)*(theta+theta0));
            mode_profiles(:,midx,:,:,2) = (besselj(nm(1,midx)-1,target_ki(:,midx).*r)+1i*unm(midx)^2./(2*nm(1,midx)*target_k0.*target_n_in*core_radius).*sqrt(target_n_out.^2-target_n_in.^2).*Dbesselj(nm(1,midx),target_ki(:,midx).*r)).*cos(nm(1,midx)*(theta+theta0));
    end
end
norm = sqrt(sum(sum(sum(abs(mode_profiles).^2,5).*r*dr*dtheta,3),4));
nonzero_idx = norm~=0;
mode_profiles(nonzero_idx) = mode_profiles(nonzero_idx)./norm(nonzero_idx);

% Transform into (x,y) basis
if use_gpu
    mode_profiles_xy = zeros(size(mode_profiles),'gpuArray');
else
    mode_profiles_xy = zeros(size(mode_profiles));
end
mode_profiles_xy(:,:,:,:,1) = mode_profiles(:,:,:,:,1).*cos(theta) - mode_profiles(:,:,:,:,2).*sin(theta);
mode_profiles_xy(:,:,:,:,2) = mode_profiles(:,:,:,:,1).*sin(theta) + mode_profiles(:,:,:,:,2).*cos(theta);

%% Form a basis of linear polarizations
% Arrange degenerate modes into the 6th dimension
linear_mode_profiles = 1/sqrt(2)*(mode_profiles_xy(:,1:num_modes,:,:,:) + mode_profiles_xy(:,num_modes+1:num_modes*2,:,:,:));
linear_mode_profiles = cat(6,linear_mode_profiles,...
                       1/sqrt(2)*(mode_profiles_xy(:,1:num_modes,:,:,:) - mode_profiles_xy(:,num_modes+1:num_modes*2,:,:,:))...
                       ); % degenerate mode

% Because it may not align with x or y axis, I rotate the profiles.
center_wavelength_idx = ceil(target_wavelength_sampling/2);
xy_contribution = sum(sum(abs(linear_mode_profiles(center_wavelength_idx,:,:,:,:,:)).^2.*r*dr*dtheta,3),4);
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
xy_contribution = sum(sum(abs(linear_mode_profiles(center_wavelength_idx,:,:,:,:,:)).^2.*r*dr*dtheta,3),4);
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
mode_profiles = complex(zeros(target_wavelength_sampling,num_modes*2,r_sampling,theta_sampling));
if use_gpu
    mode_profiles = gpuArray(mode_profiles);
end
mode_profiles(:,1:2:end,:,:) = linear_mode_profiles(:,:,:,:,1); % the 1st degenerate mode
mode_profiles(:,2:2:end,:,:) = linear_mode_profiles(:,:,:,:,2); % the 2nd degenerate mode

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
mode_profiles = mode_profiles(:,chosen_midx,:,:);

% Orthogonality
% Make the modes are orthogonal to each other
norm = sqrt(sum(sum(abs(mode_profiles(:,1,:,:)).^2.*r*dr*dtheta,3),4));
mode_profiles(:,1,:,:) = mode_profiles(:,1,:,:)./norm;
for midx1 = 2:length(chosen_midx)
    for midx2 = 1:(midx1-1)
        mode_profiles(:,midx1,:,:) = mode_profiles(:,midx1,:,:) - sum(sum(mode_profiles(:,midx1,:,:).*conj(mode_profiles(:,midx2,:,:)).*r*dr*dtheta,3),4).*mode_profiles(:,midx2,:,:);
    end
    norm = sqrt(sum(sum(abs(mode_profiles(:,midx1,:,:)).^2.*r*dr*dtheta,3),4));
    mode_profiles(:,midx1,:,:) = mode_profiles(:,midx1,:,:)./norm;
end

%% Mode-field diameter
SR = calc_SR_tensors(mode_profiles,r,dr,dtheta,struct('gpu_yes',use_gpu),'./');
fprintf('MFD=%4.2f um\n',sqrt(1/SR/pi)*2*1e6);

%% Save data
r = squeeze(r);
theta = squeeze(theta);

% change the notation into "beta"
if use_gpu
    beta = complex(zeros(wavelength_sampling,num_modes*2,'gpuArray'));
else
    beta = complex(zeros(wavelength_sampling,num_modes*2));
end
beta(:,1:2:end) = gamma(:,1:num_modes);
beta(:,2:2:end) = gamma(:,1:num_modes);
beta = beta(:,chosen_midx);

if use_gpu
    wavelength = gather(wavelength);
    beta = gather(beta);
    mode_profiles = gather(mode_profiles);
end

save(saved_filename,'wavelength','beta','mode_profiles','r','theta','pressure','temperature','core_radius','n_out');

%% Plot it!
plot_theta = [theta;theta(1)];
plot_mode_profiles = cat(4,mode_profiles,mode_profiles(:,:,:,1));
for midx = 1:length(chosen_midx)
    figure; polarPcolor(r',plot_theta'*180/pi,squeeze(real(plot_mode_profiles(ceil(target_wavelength_sampling/2),midx,:,:)))); colormap(jet);
end

% Plot the loss of each mode
figure;
h = plot(wavelength*1e6,imag(beta)*20*log10(exp(1)));
xlabel('Wavelength (\mum)');
ylabel('Loss (dB/m)');
set(h,'linewidth',2);
set(gca,'fontsize',20);
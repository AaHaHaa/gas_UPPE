% Calculate the Taylor series coefficients of the propagation constant,
% which represents the phase velocity, group velocity, and group
% dispersion, etc.

%clearvars; close all;

addpath('../user_helpers');

num_disp_orders = 3;
use_gpu = false;%true;

load('info_AR_HC_PCF_N2_30um_100atm_300nm.mat','beta','wavelength');

c = 2.99792458e-4; % speed of ligth; m/ps

Nf = size(beta,1);
num_modes = size(beta,2);

wavelength_min = 0.1; % um
wavelength_max = 10; % um

% Show the Stokes,pump,anti-Stokes phase matching condition 
wavelength_pump = 1.03; % um
Stokes_shift_V = 125;
Stokes_shift_R = 17.6;

%% Calculate the propagation constants
wavelength = wavelength*1e6; % um
f_calc = c./wavelength*1e6; % THz; loaded from the file

f = linspace(f_calc(end),f_calc(1),Nf)'; % THz; resample for applying Taylor series expansion
if use_gpu
    beta = myPchip( flipud(f_calc),flipud(beta),f,6,'../cuda/' );
else
    beta = interp1(f_calc,beta,f,'pchip','extrap');
end

omega = 2*pi*f; % angular frequencies in 1/ps
df = f(2)-f(1);
domega = 2*pi*df;
beta = real(beta); % beta in 1/m

%% Calculate the derivatives

% We need to use cell arrays because the higher orders are calculated from
% finite differences. This means that each order has one less data point
% than the previous one.
omega_full = cell(num_disp_orders+1, 1); % 2*pi*THz
wavelength_full = cell(num_disp_orders+1, 1); % um

central_diff_coeff = {[   0,    0, -1/2,   0, 1/2,   0,   0],...
                      [   0,    0,    1,  -2,   1,   0,   0],...
                      [   0, -1/2,    1,   0,  -1, 1/2,   0],...
                      [   0,    1,   -4,   6,  -4,   1,   0],...
                      [-1/2,    2, -5/2,   0, 5/2,  -2, 1/2],...
                      [   1,   -6,   15, -20,  15,  -6,   1]};

omega_full{1} = omega;
wavelength_full{1} = 2*pi*c./omega_full{1}*1e6;
for disp_order = 1:num_disp_orders
    switch disp_order
        case {1,2}
            omega_full{disp_order+1} = omega(2:end-1); % 2*pi*THz
        case {3,4}
            omega_full{disp_order+1} = omega(3:end-2);
        case {5,6}
            omega_full{disp_order+1} = omega(4:end-3);
    end
        
    wavelength_full{disp_order+1} = 2*pi*c./omega_full{disp_order+1}*1e6; % in um
end

% beta_full will have all of the orders, for each mode, as a function of wavlength
beta_full = cell(num_disp_orders+1, 1);
beta_full{1} = beta;
for disp_order = 1:num_disp_orders
    switch disp_order
        case {1,2}
            beta_full{disp_order+1} = central_diff_coeff{disp_order}(3)*beta(1:end-2,:) + ...
                                      central_diff_coeff{disp_order}(4)*beta(2:end-1,:) + ...
                                      central_diff_coeff{disp_order}(5)*beta(3:end,:);
        case {3,4}
            beta_full{disp_order+1} = central_diff_coeff{disp_order}(2)*beta(1:end-4,:) + ...
                                      central_diff_coeff{disp_order}(3)*beta(2:end-3,:) + ...
                                      central_diff_coeff{disp_order}(4)*beta(3:end-2,:) + ...
                                      central_diff_coeff{disp_order}(5)*beta(4:end-1,:) + ...
                                      central_diff_coeff{disp_order}(6)*beta(5:end,:);
        case {5,6}
            beta_full{disp_order+1} = central_diff_coeff{disp_order}(1)*beta(1:end-6,:) + ...
                                      central_diff_coeff{disp_order}(2)*beta(2:end-5,:) + ...
                                      central_diff_coeff{disp_order}(3)*beta(3:end-4,:) + ...
                                      central_diff_coeff{disp_order}(4)*beta(4:end-3,:) + ...
                                      central_diff_coeff{disp_order}(5)*beta(5:end-2,:) + ...
                                      central_diff_coeff{disp_order}(6)*beta(6:end-1,:) + ...
                                      central_diff_coeff{disp_order}(7)*beta(7:end,:);
    end
    beta_full{disp_order+1} = beta_full{disp_order+1}/domega^disp_order;
end

% ps^n/m to fs^n/mm
for disp_order = 1:num_disp_orders+1
    beta_full{disp_order} = beta_full{disp_order}*10^(3*disp_order-6);
end

%% Display the results
coo = distinguishable_colors(num_modes);

ylabels = cell(num_disp_orders+1, 1);
ylabels{1} = '1/mm';
ylabels{2} = 'fs/mm';
for disp_order = 2:num_disp_orders
    ylabels{disp_order+1} = ['fs^' num2str(disp_order) '/mm'];
end

% Plot the results
figure;
lowest_order = 1;
for disp_order = lowest_order+1:num_disp_orders+1
    subplot(1,num_disp_orders+1-lowest_order,disp_order-lowest_order)
    hold on
    for midx = 1:num_modes
        h = plot(wavelength_full{disp_order}, beta_full{disp_order}(:, midx), 'Color', coo(midx,:));
        set(h,'linewidth',2);
        xlim([wavelength_min,wavelength_max]);
    end
    hold off
    set(gca,'fontsize',20);
    ylabel(ylabels{disp_order});
    xlabel('\mum')
    title(['\beta_' num2str(disp_order-1)]);
end

%% Plot beta2 just for presentation
figure;
for midx = 1:num_modes
    h = plot(wavelength_full{3}, beta_full{3}(:,midx)); hold on;
    set(h,'linewidth',2);
end
hold off;
xlim([wavelength_min,wavelength_max]);
set(gca,'fontsize',20);
xlabel('Wavelength (\mum)');
ylabel('\beta_2 (fs^2/mm)');
%{
figure;
for midx = 1:num_modes
    h = plot(wavelength_full{3}, beta_full{3}(:,midx).*(-2*pi*3e8./wavelength_full{3}.^2/1e9)); hold on;
    set(h,'linewidth',2);
end
hold off;
xlim([wavelength_min,wavelength_max]);
set(gca,'fontsize',20);
xlabel('Wavelength (\mum)');
ylabel('D (ps/nm/km)');
%}
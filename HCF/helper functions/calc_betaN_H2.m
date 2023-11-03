clearvars; %close all;

addpath('SLMtools');

num_disp_orders = 4;

c = 2.99792458e-4; % speed of ligth m/ps

load('n_H2.mat','slm');
wavelength_i = 100e-9; % m
wavelength_f = 15e-6;
w = linspace(c/wavelength_f,c/wavelength_i,10000)'; % THz
wavelength = c./w; % m
n_H2 = slmeval(wavelength*1e6,slm,0); % refractive index of H2
beta = n_H2*2*pi./wavelength;

Nf = size(beta,1);
num_modes = 1;

%% Calculate the propagation constants
l = wavelength*1e6; % um
fo = c./l*1e6; % THz

f = linspace(fo(end),fo(1),Nf)';
abs_beta = interp1(fo,abs(beta),f,'pchip');
ang_beta = interp1(fo,unwrap(angle(beta),[],1),f,'pchip');
beta = abs_beta.*exp(1i*ang_beta);

w = 2*pi*f; % angular frequencies in 1/ps
df = f(2)-f(1);
dw = 2*pi*df;
beta_calc = real(beta); % beta in 1/m

%% Display the results

% We need to use cell arrays because the higher orders are calculated from
% finite differences. This means that each order has one less data point
% than the previous one.
w_vectors = cell(num_disp_orders+1, 1); % omegas, in 1/ps
l_vectors = cell(num_disp_orders+1, 1); % lambdas, in um
w_vectors{1} = w;
l_vectors{1} = 2*pi*c./w_vectors{1}*1e6;
for disp_order = 1:num_disp_orders
    w_prev = w_vectors{disp_order};
    w_vectors{disp_order+1} = dw/2 + w_prev(1:length(w_prev)-1); % in 1/ps
    l_vectors{disp_order+1} = 2*pi*c./w_vectors{disp_order+1}*1e6; % in um
end

% beta_full will have all of the orders, for each mode, as a function of
% wavlength
beta_full = cell(num_disp_orders+1, 1);
beta_full{1} = beta_calc/1000;
for disp_order = 1:num_disp_orders
    beta_full{disp_order+1} = zeros(Nf-disp_order, num_modes);
end

% Take the differences to calculate the higher orders
for disp_order = 1:num_disp_orders
    beta_full{disp_order+1} = diff(beta_full{disp_order})/dw*1000;
end

coo=hsv(num_modes);

ylabels = cell(num_disp_orders+1, 1);
ylabels{1} = '1/mm';
ylabels{2} = 'fs/mm';
for disp_order = 2:num_disp_orders
    ylabels{disp_order+1} = ['fs^' num2str(disp_order) '/mm'];
end

% Plot the results
figure;
for disp_order = 1:num_disp_orders+1
    subplot(1,num_disp_orders+1,disp_order)
    hold on
    for midx = 1:num_modes
        plot(l_vectors{disp_order}, beta_full{disp_order}(:, midx), 'Color', coo(midx,:))
    end
    hold off
    ylabel(ylabels{disp_order})
    xlabel('\mum')
    axis tight
end
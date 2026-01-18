% This script plots both the Raman-gain spectrum [imag(ifft(R))] and the
% Raman coefficient.

clearvars; close all;

addpath('../broadband UPPE algorithm/')

%% Configurations
gas.material = {'D2'};
pressure0 = 1.01325e5; % Pa
temperature0 = 273.15; % 0 degree Celsius
gas.pressure = pressure0;
gas.temperature = temperature0;
% eta is calculated with the unit, amagat
% Ideal gas law is used here.
eta = gas.pressure/pressure0*temperature0/gas.temperature;

% Number density of the gas
k = 1.38064852e-23; % Boltzmann constant (MKS unit)
gas.Ng = gas.pressure/k/gas.temperature; % m^(-3)

time_window = 5e3; % ps
Nt = 2^22;
dt = time_window/Nt; % ps

T = (0:Nt-1)'*dt; % ps

% For linear polarization, the rotational Raman contribution is 4 times the
% Raman strength R_rot here.
%
% See Eq.(39) in
% Chen and Wise, "Unified and vector theory of Raman scattering in 
% gas-filled hollow-core fiber across temporal regimes," APL Photonics
% 9(3), 030902 (2024).
%
% For circular polarization, the factor is 1.
polarization_factor = 4; % 4 for linear polarization and 1 for circular polarization

gas = Raman_model(gas,eta);
if isfield(gas.(gas.material{1}),'V')
    max_f_to_plot = max(gas.(gas.material{1}).V.omega/2/pi*1.1);
else
    max_f_to_plot = max(gas.(gas.material{1}).R.omega/2/pi*1.1);
end

%% Raman-gain spectrum
Raman_gain_R = [];
Raman_gain_V = [];
if isfield(gas.(gas.material{1}),'R')
    Raman_gain_R = sum(gas.(gas.material{1}).R.preR.*exp(-T./gas.(gas.material{1}).R.T2).*sin(gas.(gas.material{1}).R.omega.*T),2)*polarization_factor;
end
if isfield(gas.(gas.material{1}),'V')
    Raman_gain_V = sum(gas.(gas.material{1}).V.preR.*exp(-T./gas.(gas.material{1}).V.T2).*sin(gas.(gas.material{1}).V.omega.*T),2);
end
Raman_gain = [Raman_gain_R,Raman_gain_V];

%% Raman coefficient
Raman_f = zeros(100,1);
Raman_preR = zeros(size(Raman_f));

if isfield(gas.(gas.material{1}),'R')
    Raman_f(1:length(gas.(gas.material{1}).R.omega),1) = gas.(gas.material{1}).R.omega/2/pi;
    Raman_preR(1:length(gas.(gas.material{1}).R.preR),1) = gas.(gas.material{1}).R.preR*polarization_factor;
    if isfield(gas.(gas.material{1}),'V')
        Raman_f(length(gas.(gas.material{1}).R.omega)+1:length(gas.(gas.material{1}).R.omega)+length(gas.(gas.material{1}).V.omega),1) = gas.(gas.material{1}).V.omega/2/pi;
        Raman_preR(length(gas.(gas.material{1}).R.omega)+1:length(gas.(gas.material{1}).R.omega)+length(gas.(gas.material{1}).V.omega),1) = gas.(gas.material{1}).V.preR;
    end
else
    Raman_f(1:length(gas.(gas.material{1}).V.omega),1) = gas.(gas.material{1}).V.omega/2/pi;
    Raman_preR(1:length(gas.(gas.material{1}).V.preR),1) = gas.(gas.material{1}).V.preR;
end

final_Raman_f = linspace(0,max_f_to_plot,30000)';
final_Raman_preR = zeros(30000,1);
for j = 1:100
    if Raman_f(j) > 0
        idx = find(final_Raman_f(:)>Raman_f(j),1);
        final_Raman_f(idx) = Raman_f(j);
        final_Raman_preR(idx) = Raman_preR(j);
    end
end

%%
f = (-Nt/2:Nt/2-1)/time_window;

% Plot Raman-gain spectrum
figure;
plot(f,imag(fftshift(ifft(Raman_gain),1)),'linewidth',2,'Color','b');
xlabel('Raman frequency (THz)');
ylabel('Raman gain');
set(gca,'fontsize',25);

% Plot Raman coefficient
figure;
fp = get(gcf,'position');
screen_size = get(0,'ScreenSize');
original_top = screen_size(4)-fp(2)-fp(4);
set(gcf,'position',[fp(1) screen_size(4)-original_top-fp(4) fp(3)*6/4 fp(4)]);
plot(final_Raman_f,final_Raman_preR,'linewidth',2,'Color','b');
xlabel('Raman frequency (THz)');
ylabel('R^{coeff}');
%xlim([-5,20]);
%ylim([0,1e-24]);
set(gca,'fontsize',25);
%print(gcf,'Raman coeff (rot).pdf','-dpdf');
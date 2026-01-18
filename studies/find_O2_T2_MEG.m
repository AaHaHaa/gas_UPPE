% This code computes the O2's pressure-dependent T2 from the MEG model and
% saves it to a MATLAB .mat file for the main script to use.
%
% See find_Raman_linewidth() for detail.

addpath('../broadband UPPE algorithm/');

clearvars; close all;

all_eta = 0.1:0.1:500;
length_J = 33; % for O2; check Raman_model()
T2 = zeros(length(all_eta),length_J);

for i = 1:length(all_eta)
    eta = all_eta(i); % gas density in amagat
    pressure0 = 1.01325e5; % Pa
    temperature0 = 273.15; % 0 degree Celsius
    gas.temperature = temperature0 + 25; % K
    gas.pressure = eta*gas.temperature/temperature0*pressure0; % Pa
    c = 299792458; % m/s
    h = 6.62607015e-34; % J*s
    hbar = h/(2*pi); % J*s
    k = 1.38064852e-23; % Boltzmann constant (SI unit)
    au_polarizability = 1.64878e-41; % F*m^2; atomic unit
    permittivity0 = 8.85e-12;
    gas.Ng = gas.pressure/k/gas.temperature; % m^(-3)
    
    %% The calculated Raman gain of my derived equation
    gas.material = {'O2'};
    if i == 1
        gas = Raman_model( gas,eta );
    else
        gas = Raman_model( gas,eta,T2(i-1,:) );
    end

    T2(i,:) = gas.(gas.material{1}).V.T2;
end

eta = all_eta';
save(sprintf('%s_T2_vs_P.mat',gas.material{1}),'eta','T2');
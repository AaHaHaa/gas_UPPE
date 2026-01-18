% This code plots the vibrational Raman gain spectrum of O2.
% Collisional narrowing is observed.
%
% See find_Raman_linewidth() for detail.

addpath('../broadband UPPE algorithm/');

clearvars; close all;

all_eta = 20;%0.1:0.1:30;
length_J = 33; % for O2; check Raman_model()
T2 = zeros(length(all_eta),length_J);

save_point = length(all_eta);
Frame(save_point) = struct('cdata',[],'colormap',[]);
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
    gas = Raman_model( gas,eta );

    max_J = length(gas.(gas.material{1}).V.omega)-1; % minus one is due to Ji beginning from zero

    Ji = 0:max_J; % rotational quantum number
    
    ignore_Ji = (gas.(gas.material{1}).R.gJ(Ji) == 0);
    Ji(ignore_Ji) = []; % ignore transitions that has zero strength due to no population for it to happen
    
    vib_omega = gas.(gas.material{1}).V.omega;
    vib_omega(ignore_Ji) = []; % ignore transitions that has zero strength due to no population for it to happen
    
    % Energies for ground and first-excited vibrational levels
    EJ_low = @(J) (gas.(gas.material{1}).R.B0*J.*(J+1) - gas.(gas.material{1}).R.D0*J.^2.*(J+1).^2)*1e2*h*c;
    EJ_up = @(J) EJ_low(J) + hbar*gas.(gas.material{1}).V.omega(J+1)*1e12; % in the unit "J"
    
    % Population
    Z = sum([gas.(gas.material{1}).R.gJ(Ji).*(2*Ji+1).*exp(-EJ_low(Ji)/k/gas.temperature),...
             gas.(gas.material{1}).R.gJ(Ji).*(2*Ji+1).*exp( -EJ_up(Ji)/k/gas.temperature)]); % partition function considering only the ground vibrational state
    rho_low = @(J) gas.(gas.material{1}).R.gJ(J).*(2*J+1).*exp(-EJ_low(J)/k/gas.temperature)/Z; % population
    rho_up = @(J) gas.(gas.material{1}).R.gJ(J).*(2*J+1).*exp(-EJ_up(J)/k/gas.temperature)/Z; % population
    
    % Raman strength is determined by the population in the initial state
    strength_ratio = (2*Ji+1).*rho_low(Ji);

    domega = max(vib_omega) - min(vib_omega); % 2*pi*THz
    Nf = 2^10;
    omega = linspace(min(vib_omega)-domega,max(vib_omega)+domega,Nf)'; % THz

    Raman_gain = calc_Raman_gain(strength_ratio,omega,gas.(gas.material{1}).V.T2(~ignore_Ji),vib_omega);

    f = (omega/2/pi*1e12)/(100*c); % cm^(-1)

    fig = figure;
    plot(f,Raman_gain/max(Raman_gain),'LineWidth',2,'Color','b');
    xlabel('Wavenumber (cm^{-1})');
    ylabel('Raman gain (norm.)');
    set(gca,'fontsize',25);
    xlim([min(f),max(f)]);
    ylim([0,1]);
    text(min(f)+(max(f)-min(f))*0.1,0.9,sprintf('%4.1f bar',gas.pressure/pressure0),'FontSize',25);

    %{
    gas.(gas.material{1}).V.MEG = struct('alpha', 0.01528*1e2*c*1e-12*2*pi,... % 2*pi*THz/atm
                     'beta', 1.362,...
                     'delta', 1.786,...
                     'a', 0.50);
    T2 = find_Raman_linewidth(gas,gas.material{1},eta(1),gas.(gas.material{1}).R,gas.(gas.material{1}).V);
    %}
    
    set(fig,'Color',[1,1,1]);

    Frame(i) = getframe(fig);
    close(fig);
end
% Movie
implay(Frame,30);

exportVideo = VideoWriter(sprintf('%s_T2',gas.material{1}));
exportVideo.FrameRate = 30;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);

%% Helper functions
function output = F_damped_sine(omega,T2,omega_R)

output = 1./(omega_R.^2 + (1./T2 + 1i*omega).^2);

end

function Raman_gain = calc_Raman_gain(strength_ratio,omega,T2,omega_R)

Raman_gain = -imag(sum(strength_ratio.*F_damped_sine(omega,T2,omega_R),2));

end

function loss = fit_MEG(x,...
                        strength_ratio,...
                        omega,omega_R,...
                        MEG_gain)

% Make it row vectors (and omega is a column vector)
if size(x,1) ~= 1
    x = x.';
end
T2 = abs(x);

Raman_gain = calc_Raman_gain(strength_ratio,omega,T2,omega_R);
Raman_gain = Raman_gain/max(Raman_gain);
loss = sum(abs(MEG_gain - Raman_gain).^2);

end
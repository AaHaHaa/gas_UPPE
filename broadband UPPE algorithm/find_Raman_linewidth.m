function T2 = find_Raman_linewidth(gas,gas_material,eta,rot,vib)
%FIND_RAMAN_LINEWIDTH It finds the Raman linewidth/T2 with the modified
%exponential gap (MEG) model. Due to the collisional narrowing, increasing
%pressure leads to a reduce linewidth for the vibrational Raman spectral
%line.
%
% Because of collisional narrowing, Raman transitions are coupled. Here, I
% assume that they reasonably maintain a sinusoidal response such that my
% derived Raman equation is still valid. I assume that the interaction of
% line mixing results only in varied dephasing of transitions while
% maintaining the (temporal) Raman strength, preR in this code. For
% example, increasing pressure leads to a drop of dephasing of some
% transitions such that they dominate the process (e.g. Raman gain). The
% empirical MEG model approximates the collisional-narrowing effect by
% finding the cross-dephasing terms, which validates my assumption.
%
% The equation here follows
% Lempert et al., "Stimulated Raman scattering and coherent anti-Stokes
% Raman spectroscopy in high-pressure oxygen," J. Opt. Soc. Am. B 7(5),
% 715-721 (1990).
%
% Also check
% Rahn and Palmer, "Studies of nitrogen self-broadening at high temperature
% with inverse Raman spectroscopy," J. Opt. Soc. Am. B 3(9), 1164-1169
% (1986).

c = 299792458; % m/s
h = 6.62607015e-34; % J*s
hbar = h/(2*pi); % J*s
k = 1.38064852e-23; % Boltzmann constant (SI unit)

%% The calculated Raman gain of my derived equation
max_J = length(vib.omega)-1; % minus one is due to Ji beginning from zero

Ji = 0:max_J; % rotational quantum number

ignore_Ji = (rot.gJ(Ji) == 0);
Ji(ignore_Ji) = []; % ignore transitions that has zero strength due to no population for it to happen

vib_omega = vib.omega;
vib_omega(ignore_Ji) = []; % ignore transitions that has zero strength due to no population for it to happen

if isscalar(vib.T2)
    T2 = vib.T2.*ones(1,length(Ji));
else % if there is an array of initial guess for T2 optimization process
    T2 = vib.T2;
    T2(ignore_Ji) = [];
    T2 = T2 + (rand(1,length(T2))-1/2)*std(T2)/5;
end

% Energies for ground and first-excited vibrational levels
EJ_low = @(J) (rot.B0*J.*(J+1) - rot.D0*J.^2.*(J+1).^2)*1e2*h*c; % in the unit "J"
% vib.omega is the Q transition with J fixed, so the energy of (v=1,J) is
% the summation of (v=0,J) from EJ_low(J) and Q(J) transition.
EJ_up = @(J) EJ_low(J) + hbar*vib.omega(J+1)*1e12; % in the unit "J"

% Population
Z = sum([rot.gJ(Ji).*(2*Ji+1).*exp(-EJ_low(Ji)/k/gas.temperature),...
         rot.gJ(Ji).*(2*Ji+1).*exp( -EJ_up(Ji)/k/gas.temperature)]); % partition function considering only the ground vibrational state
rho_low = @(J) rot.gJ(J).*(2*J+1).*exp(-EJ_low(J)/k/gas.temperature)/Z; % population
rho_up = @(J) rot.gJ(J).*(2*J+1).*exp(-EJ_up(J)/k/gas.temperature)/Z; % population

% Raman strength is determined by the population in the initial state
strength_ratio = rho_low(Ji);

%% Plotting parameters to show the Raman gain spectra
domega = max(vib_omega) - min(vib_omega); % 2*pi*THz
Nf = 2^10;
omega = permute(linspace(min(vib_omega)-domega,max(vib_omega)+domega,Nf)',[2,3,1]); % THz

%% MEG model
pressure = eta*gas.temperature/273.15; % eta is in amagat, so we need to consider the temperature, and the reference temperature 273.15 K used in amagat, to transform into pressure
Tref = 295; % K; reference temperature for the fitting data
% The empirical MEG formula doesn't have "n" for O2 to correct for its
% temperature dependence but only for N2. Since the computations of
% ultrafast gas nonlinear optics is mostly applied in room temperature, I
% keep N2's n-including MEG formula even for O2 as it shouldn't affect its
% result too much. Although I haven't significantly studied its origin, I
% speculate that this "n" might be applicable to O2 as well.
n = 1.346; % From (N2 data) Rahn and Palmer, "Studies of nitrogen self-broadening at high temperature with inverse Raman spectroscopy," J. Opt. Soc. Am. B 3(9), 1164-1169 (1986)
% Cross-dephasing term gamma_ij: line mixing effect
% With its off-diagonal terms being zero, rovibrational Raman transitions
% are treated as independent and there is no line-mixing effect.
%
% From Eq.(1) in Rahn and Palmer's paper, also from Eq.(5) in Lempert's paper
% Eq.(1) is the more-general version.
gamma_up = vib.MEG.alpha*pressure*(1-exp(-vib.MEG.m))/(1-exp(-vib.MEG.m*gas.temperature/Tref))*(Tref/gas.temperature)^n*((1+vib.MEG.a*EJ_low(Ji)/k/gas.temperature/vib.MEG.delta)./(1+vib.MEG.a*EJ_low(Ji)/k/gas.temperature)).^2.*exp(-vib.MEG.beta*abs(EJ_low(Ji)'-EJ_low(Ji))/k/gas.temperature);
gamma_down = (gamma_up.*rho_low(Ji)./rho_low(Ji)').'; % from the detailed balance; Eq.(2) in Rahn and Palmer's

% The final gamma_ij matrix
idx = 1:length(Ji);
gamma = gamma_up;
idx2 = idx' < idx;
gamma(idx2) = gamma_down(idx2);
idx3 = idx' == idx;

% The diagonal terms are made zero because later we use it to compute the
% G-matrix whose diagonal term doesn't come from gamma_ij but from the
% resonance feature of each transition, (1i*(omega(ii)-vib_omega) +
% Gamma2).
gamma(idx3) = 0;

% Its diagonal terms correspond to the dephasing of each transition
Gamma2 = -sum(gamma,1); % Gamma/2 = 1/2/T2; Eq.(3) in Rahn and Palmer's

% Compute the G-matrix [Eq.(4) in Lempert's paper]
G = repmat(gamma,1,1,Nf); % extend to include the frequency dimension
for ii = 1:Nf
    tmp = eye(length(Ji)).*(1i*(omega(ii)-vib_omega) + Gamma2);
    G(:,:,ii) = G(:,:,ii) + tmp;
end
invG = pageinv(G);
% Then compute the MEG gain
X3_factor = squeeze(sum(invG.*(rho_low(Ji)-rho_up(Ji)),[1,2])); % Eq.(3) in Lempert's paper
MEG_gain = -real(X3_factor);
% Normalize the MEG's gain for optimization to find the corresponding T2
MEG_gain = MEG_gain/max(MEG_gain);

%% Fitting
% Find dominant transitions that determine the Raman gain shape
min_dominant_omega_idx = find(MEG_gain > max(MEG_gain)/4,1);
max_dominant_omega_idx = find(MEG_gain > max(MEG_gain)/4,1,'last');
dominant_idx = min_dominant_omega_idx:max_dominant_omega_idx;
dominant_omega = omega(dominant_idx); % omega in the dominant regime
max_dominant_Jidx = find(vib_omega < min(dominant_omega),1); % vib_omega is a monotonically-decreasing array
if isempty(max_dominant_Jidx)
    max_dominant_Jidx = length(vib_omega);
end
min_dominant_Jidx = find(vib_omega > max(dominant_omega),1,'last'); % vib_omega is a monotonically-decreasing array
if isempty(min_dominant_Jidx)
    min_dominant_Jidx = 1;
end
dominant_Jidx = min_dominant_Jidx:max_dominant_Jidx; % dominant J's


omega = squeeze(omega);

% Rough optimization with dominant transitions
find_opt = @(x) fit_MEG(x,strength_ratio(dominant_Jidx),squeeze(dominant_omega),vib_omega(dominant_Jidx),MEG_gain(dominant_idx));
option = optimset('TolFun',0.01,'MaxFunEvals',300*length(dominant_Jidx),'MaxIter',300*length(dominant_Jidx),'Display','off');
%option = optimset('TolFun',0.1,'PlotFcns',@optimplotfval,'MaxFunEvals',300*length(dominant_Jidx),'MaxIter',300*length(dominant_Jidx)); % plot the process of optimization
y = fminsearch(find_opt,T2(dominant_Jidx),option);
y = smoothing(y);

% Fine optimization with all transitions
T2(dominant_Jidx) = y;
find_opt = @(x) fit_MEG(x,strength_ratio,omega,vib_omega,MEG_gain);
option = optimset('TolX',0.01,'TolFun',0.01,'MaxFunEvals',500*length(vib_omega),'MaxIter',500*length(vib_omega),'Display','off');
%option = optimset('TolX',0.1,'TolFun',0.1,'PlotFcns',@optimplotfval,'MaxFunEvals',500*length(vib_omega),'MaxIter',500*length(vib_omega)); % plot the process of optimization
y = fminsearch(find_opt,T2,option);
y = smoothing(y);

T2 = vib.T2.*ones(1,length(vib.omega)); % initialization with a non-zero value; otherwise, T2 at ignore_Ji will be zero and will make the computation wrong in later pulse propagation due to the included 1/T2 computation
T2(~ignore_Ji) = y;

%% Comparison during debugging
debug_plot_yes = false;
if debug_plot_yes
    f = (omega/2/pi*1e12)/(100*c); % cm^(-1)

    if isscalar(vib.T2)
        Raman_gain0 = calc_Raman_gain(strength_ratio,omega,vib.T2,vib_omega);
    else
        Raman_gain0 = calc_Raman_gain(strength_ratio,omega,vib.T2(~ignore_Ji),vib_omega);
    end
    Raman_gain = calc_Raman_gain(strength_ratio,omega,y,vib_omega);

    figure;
    plot(f,-real(X3_factor)/max(-real(X3_factor)),'linewidth',2,'Color','b');
    hold on;
    plot(f,Raman_gain/max(Raman_gain),'linewidth',2,'Color','r');
    %plot(f,Raman_gain0/max(Raman_gain0),'linewidth',2,'Color','c');
    hold off;
    xlabel('Wavenumber (cm^{-1})');
    ylabel('Raman gain (a.u.)');
    set(gca,'fontsize',25);
    %legend('MEG','Fitted','Isolated');
    legend('MEG','Fitted');
    print(gcf,sprintf('%s_vib_Raman.jpg',gas_material),'-djpeg');
end

end

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
T2 = smoothing(x);

Raman_gain = calc_Raman_gain(strength_ratio,omega,T2,omega_R);
Raman_gain = Raman_gain/max(Raman_gain);
loss = sum(abs(MEG_gain - Raman_gain).^2);

end

function x = smoothing(x)
% Smoothing is required to remove significant variations of T2

x = abs(x); % optimization might give negative values, but T2 cannot be negative
%x = ( x + circshift(x,1) + circshift(x,-1) )/3; % smoothing
x(x > mean(x) + std(x)) = mean(x); % remove outliers
x = smooth(x)'; % smoothing

end
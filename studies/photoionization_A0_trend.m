% This code plots the A0 w.r.t. Keldysh parameter for H2.
%
% Please see the supplement of the following paper for details about the
% photoionization computation and what A0 is.
%
% Chen et al., "Femtosecond long-wave-infrared generation in
% hydrogen-filled hollow-core fiber," J. Opt. Soc. Am. B 40(4), 796-806
% (2023)

clearvars; close all;

lambda = 800e-9; % m
c = 299792458; % m/s
E_b = 15.43/(2*13.6); % ionization energy under a.u.
omega = 2*pi*c/lambda/4.13414e16; % omega under a.u.

Z = 1;

r = linspace(0,10,1000)'; r = r(2:end); % Keldysh parameter
n = Z*floor(E_b/omega+1);
v = E_b/omega*(1+1/2./r.^2);
beta = 2*r./sqrt(1+r.^2);

%% Photoionization - erfi() lookup table
% Because calculating erfi is slow, it's faster if I create a lookup table
% and use interp1. The range of input variable for erfi is 0~sqrt(2) only.
n_Am = 1000; % the number of summation of Am term in photoionization
erfi_x = linspace(0,sqrt(2*(n_Am+1)),1000)';
erfi_y = erfi(erfi_x);

%%
A0 = zeros(size(r));
for i = 0:200
    n_v = ceil(v) - v+i;
    erfix = interp1(erfi_x,erfi_y,sqrt(beta.*n_v));
    
    A0 = A0 + 2/sqrt(3)*r.^2./(1+r.^2).*exp(-2*asinh(r).*n_v).*erfix;
    
    if i == 0
        A = A0;
    end
    if i == 10
        A2 = A0;
    end
end

figure;
h = plot(r,[A0,A2,A],'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r'); set(h(3),'Color','k');
set(gca,'fontsize',20);
xlabel('Keldysh parameter \gamma');
ylabel('A_0');
legend('S=0~200','S=0~10','S=0');
%print(gcf,'A0_trend_summation_number.pdf','-dpdf');
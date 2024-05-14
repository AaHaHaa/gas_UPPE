% This code plots the A1 w.r.t. Keldysh parameter for Ar.
%
% Please see the supplement of the following paper for details about the
% photoionization computation and what A1 is.
%
% Chen et al., "Femtosecond long-wave-infrared generation in
% hydrogen-filled hollow-core fiber," J. Opt. Soc. Am. B 40(4), 796-806
% (2023)

clearvars; close all;

lambda = 800e-9; % m
c = 299792458; % m/s
E_b = 15.759/(2*13.6); % ionization energy under a.u.
omega = 2*pi*c/lambda/4.13414e16; % omega under a.u.

Z = 1;

r = linspace(0,10,1000)'; % Keldysh parameter
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
A1 = zeros(size(r));
for i = 0:50
    n_v = ceil(v) - v+i;
    erfix = interp1(erfi_x,erfi_y,sqrt(beta.*n_v));
    
    A1 = A1 + 4/sqrt(3*pi)*r.^2./(1+r.^2).*exp(-2*n_v.*asinh(r)).*...
              ( sqrt(pi)/2*beta.*n_v.*erfix + ...
                sqrt(pi)/4*erfix - ...
                sqrt(beta.*n_v)/2.*exp(beta.*n_v) );
    
    if i == 0
        A = A1;
    end
    if i == 10
        A2 = A1;
    end
end

figure;
h = plot(r,[A1,A2,A],'linewidth',2);
set(h(1),'Color','b'); set(h(2),'Color','r'); set(h(3),'Color','k');
set(gca,'fontsize',20);
xlabel('Keldysh parameter \gamma');
ylabel('A_0');
legend('S=0~200','S=0~10','S=0');
%print(gcf,'A1_trend_summation_number.pdf','-dpdf');
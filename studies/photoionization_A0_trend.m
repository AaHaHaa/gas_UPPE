clearvars; close all;

lambda = 800e-9; % m
c = 299792458; % m/s
E_b = 15.43/(2*13.6); % ionization energy under a.u.
omega = 2*pi*c/lambda/4.13414e16; % omega under a.u.

r = linspace(0,10,1000)'; % Keldysh parameter
n = floor(E_b/omega+1);
v = E_b/omega*(1+1/2./r.^2);

A0 = zeros(size(r));
for i = 0:200
    n_v = ceil(v) - v+i;
    A0 = A0 + 2/sqrt(3)*r.^2./(1+r.^2).*exp(-2*asinh(r).*n_v).*erfi(sqrt(2*r./sqrt(1+r.^2).*n_v));
    
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
print(gcf,'A0_trend_summation_number.pdf','-dpdf');
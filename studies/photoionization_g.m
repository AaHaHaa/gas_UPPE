clearvars; close all;

r = linspace(0,100,1000)';
g = 3/2./r.*((1+1/2./r.^2).*asinh(r) - sqrt(1+r.^2)/2./r);
g_at_large_r = 3*log(2*r)/2./r;


figure;
plot(r,[g,g_at_large_r]);
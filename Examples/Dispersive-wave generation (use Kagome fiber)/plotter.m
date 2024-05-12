clearvars; close all;

addpath('../../../user_helpers/');

load('DW.mat');

spectrum = sum(abs(fftshift(ifft(prop_output.fields),1)).^2.*3e2./(3e2./f).^2,2);
spectrum = squeeze(spectrum).';
log_spectrum = 10*log10(spectrum); log_spectrum = log_spectrum - max(log_spectrum(:));

figure;
pcolor(f/1e3,prop_output.z*100,log_spectrum); shading interp;
cmap = whitejet_lower; colormap(cmap); colormap(jet);
caxis([-60,0]); xlim([0.1,1.5]);
c = colorbar; ylabel(c,'Intensity (dB)');
set(gca,'fontsize',20);
xlabel('Frequency (THz)');
ylabel('Propagation distance (cm)');
title('Spectrum');

% Draw white dashed lines
hold on;
x1 = [0.9,0.9,1.03,1.03,0.9]';
y1 = [10,20,20,10,10]';
plot(x1,y1,'--','Color','w','Linewidth',2);
text(0.9,8,'HE_{11}','Color','w','fontsize',12);
x2 = [1.095,1.095,1.224,1.224,1.095]';
y2 = [10,20,20,10,10]';
plot(x2,y2,'--','Color','w','Linewidth',2);
text(1.095,8,'HE_{12}','Color','w','fontsize',12);
x3 = [1.3,1.3,1.374,1.374,1.3]';
y3 = [10,20,20,10,10]';
plot(x3,y3,'--','Color','w','Linewidth',2);
text(1.3,8,'HE_{13}','Color','w','fontsize',12);
hold off;
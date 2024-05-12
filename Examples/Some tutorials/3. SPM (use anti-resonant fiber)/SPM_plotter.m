close all; clearvars;

load('SPM.mat');

%% Plot
% Time
figure;
h = plot(t,abs(prop_output.fields(:,:,end)).^2);
xlabel('t');
ylabel('Power');
title('Field');
set(h,'linewidth',2);
set(gca,'fontsize',14);
xlim([-1.2,1.2]);

% Spectrum
figure;
h = plot(f,abs(fftshift(ifft(prop_output.fields(:,:,end)),1)).^2);
xlabel('Frequency (THz)');
ylabel('PSD');
title('Spectrum');
set(h,'linewidth',2);
set(gca,'fontsize',14);
xlim([200,360]);

% Comparison of time
figure;
[x,y] = meshgrid(t,prop_output.z);
pcolor(x,y,permute(abs(prop_output.fields(:,1,:)).^2,[3 1 2]));
shading interp; colormap(jet);
xlabel('t');
ylabel('z');
title('Field during propagation');
set(gca,'fontsize',14);
xlim([-1.2,1.2]);

% Comparison of spectra
figure;
[x,y] = meshgrid(f,prop_output.z(2:end));
pcolor(x,y,permute(abs(fftshift(ifft(prop_output.fields(:,1,2:end)),1)).^2,[3 1 2]));
shading interp; colormap(jet);
xlabel('Frequency (THz)');
ylabel('z');
title('Spectrum during propagation');
set(gca,'fontsize',14);
xlim([200,360]);
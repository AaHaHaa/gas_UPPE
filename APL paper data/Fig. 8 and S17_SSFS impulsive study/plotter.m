%clearvars;
close all;

addpath('../../../MMTools/gas_UPPE_vector/user_helpers');

load('SSFS_H2.mat');

% pulse duration
pulse_FWHM = zeros(1,length(prop_output.z));
for zi = 1:length(prop_output.z)
    threshold = max(abs(prop_output.fields(:,:,zi)).^2)/1.01;
    [~,~,tmp_pulse_width,~] = findpeaks(abs(prop_output.fields(:,:,zi)).^2,t*1e3,'MinPeakHeight',threshold,'WidthReference','halfheight');
    pulse_FWHM(zi) = tmp_pulse_width(1);
end
diff_FWHM = diff(pulse_FWHM);
ratio_FWHM = diff_FWHM./pulse_FWHM(1:end-1);
idx = find(ratio_FWHM<0.01,1);
transition_z = prop_output.z(idx);
figure;
plot(prop_output.z,pulse_FWHM,'linewidth',2,'Color','b');
this_ylim = get(gca,'YLim');
set(gca,'fontsize',25);
xlabel('Propagation (m)');
ylabel('Duration (fs)');
hold on;
plot(transition_z*ones(2,1),[0;max(pulse_FWHM)*2],'linewidth',2,'Color','k','LineStyle','--');
hold off;
ylim([0,this_ylim(2)]);
print(gcf,'duration (H2).pdf','-dpdf');
disp(['duration (H2)=' num2str(pulse_FWHM(idx))]);

de = imag(prop_output.delta_permittivity(:,1,1,idx));
de = de/max(de);
P = abs(prop_output.fields(:,:,idx)).^2;
P = P/max(P);
figure;
yyaxis right;
plot(t,de,'linewidth',10);
ylabel('\Delta\epsilon');
ylim([-0.5,1.15]);
set(gca,'YTick',[]);
yyaxis left;
plot(t,P,'linewidth',10,'Color','b');
ylabel('Power');
ylim([0,1.3]);
set(gca,'fontsize',40);
xlim([4.25,4.65]);
set(gca,'Color','None','YTick',[],'XTick',[4.3,4.5],'YColor','b');
xlabel('Time (ps)');
print(gcf,'SSSFS (H2)_paper_closeview.pdf','-dpdf');

%%
load('SSFS_N2.mat');

% pulse duration
pulse_FWHM = zeros(1,length(prop_output.z));
for zi = 1:length(prop_output.z)
    threshold = max(abs(prop_output.fields(:,:,zi)).^2)/1.01;
    [~,~,tmp_pulse_width,~] = findpeaks(abs(prop_output.fields(:,:,zi)).^2,t*1e3,'MinPeakHeight',threshold,'WidthReference','halfheight');
    pulse_FWHM(zi) = tmp_pulse_width(1);
end
diff_FWHM = diff(pulse_FWHM);
ratio_FWHM = diff_FWHM./pulse_FWHM(1:end-1);
idx = find(ratio_FWHM<0.01,1);
transition_z = prop_output.z(idx);
figure;
plot(prop_output.z,pulse_FWHM,'linewidth',2,'Color','b');
this_ylim = get(gca,'YLim');
set(gca,'fontsize',25);
xlabel('Propagation (m)');
ylabel('Duration (fs)');
hold on;
plot(transition_z*ones(2,1),[0;max(pulse_FWHM)*2],'linewidth',2,'Color','k','LineStyle','--');
hold off;
ylim(this_ylim);
print(gcf,'duration (N2).pdf','-dpdf');
disp(['duration (N2)=' num2str(pulse_FWHM(idx))]);

de = imag(prop_output.delta_permittivity(:,1,1,idx));
de = de/max(de);
P = abs(prop_output.fields(:,:,idx)).^2;
P = P/max(P);
figure;
yyaxis right;
plot(t,de,'linewidth',10);
ylabel('\Delta\epsilon');
ylim([-0.5,1.15]);
set(gca,'YTick',[]);
yyaxis left;
plot(t,P,'linewidth',10,'Color','b');
ylabel('Power');
ylim([0,1.3]);
set(gca,'fontsize',40);
xlim([1.2,2.8]);
set(gca,'Color','None','YTick',[],'XTick',[1.5,2,2.5],'YColor','b');
xlabel('Time (ps)');
print(gcf,'SSSFS (N2)_paper_closeview.pdf','-dpdf');
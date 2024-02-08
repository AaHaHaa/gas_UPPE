clearvars; close all;

%% H2
filename = 'H2.mat';

Raman_f = zeros(100,1);
Raman_preR = zeros(size(Raman_f));

load(filename);

Raman_f(1:length(gas.H2.R.omega),1) = gas.H2.R.omega/2/pi;
Raman_preR(1:length(gas.H2.R.preR),1) = gas.H2.R.preR*4;
Raman_f(length(gas.H2.R.omega)+1:length(gas.H2.R.omega)+length(gas.H2.V.omega),1) = gas.H2.V.omega/2/pi;
Raman_preR(length(gas.H2.R.omega)+1:length(gas.H2.R.omega)+length(gas.H2.V.omega),1) = gas.H2.V.preR;

final_Raman_f = linspace(0,150,30000)';
final_Raman_preR_H2 = zeros(30000,1);
for j = 1:100
    if Raman_f(j) > 0
        idx = find(final_Raman_f(:)>Raman_f(j),1);
        final_Raman_f(idx) = Raman_f(j);
        final_Raman_preR_H2(idx) = Raman_preR(j);
    end
end

%% N2
filename = 'N2.mat';

Raman_f = zeros(100,1);
Raman_preR = zeros(size(Raman_f));

load(filename);

Raman_f(1:length(gas.N2.R.omega),1) = gas.N2.R.omega/2/pi;
Raman_preR(1:length(gas.N2.R.preR),1) = gas.N2.R.preR*4;
Raman_f(length(gas.N2.R.omega)+1:length(gas.N2.R.omega)+length(gas.N2.V.omega),1) = gas.N2.V.omega/2/pi;
Raman_preR(length(gas.N2.R.omega)+1:length(gas.N2.R.omega)+length(gas.N2.V.omega),1) = gas.N2.V.preR;

%final_Raman_f = linspace(0,150,30000)';
final_Raman_preR_N2 = zeros(30000,1);
for j = 1:100
    if Raman_f(j) > 0
        idx = find(final_Raman_f(:)>Raman_f(j),1);
        final_Raman_f(idx) = Raman_f(j);
        final_Raman_preR_N2(idx) = Raman_preR(j);
    end
end

%% O2
filename = 'O2.mat';

Raman_f = zeros(100,1);
Raman_preR = zeros(size(Raman_f));

load(filename);

Raman_f(1:length(gas.O2.R.omega),1) = gas.O2.R.omega/2/pi;
Raman_preR(1:length(gas.O2.R.preR),1) = gas.O2.R.preR*4;
Raman_f(length(gas.O2.R.omega)+1:length(gas.O2.R.omega)+length(gas.O2.V.omega),1) = gas.O2.V.omega/2/pi;
Raman_preR(length(gas.O2.R.omega)+1:length(gas.O2.R.omega)+length(gas.O2.V.omega),1) = gas.O2.V.preR;

%final_Raman_f = linspace(0,150,30000)';
final_Raman_preR_O2 = zeros(30000,1);
for j = 1:100
    if Raman_f(j) > 0
        idx = find(final_Raman_f(:)>Raman_f(j),1);
        final_Raman_f(idx) = Raman_f(j);
        final_Raman_preR_O2(idx) = Raman_preR(j);
    end
end

%%
final_Raman_preR = [final_Raman_preR_H2,final_Raman_preR_O2,final_Raman_preR_N2];

%%
figure;
fp = get(gcf,'position');
screen_size = get(0,'ScreenSize');
original_top = screen_size(4)-fp(2)-fp(4);
set(gcf,'position',[fp(1) screen_size(4)-original_top-fp(4) fp(3)*1 fp(4)]);
cc = distinguishable_colors(2);
h = plot(final_Raman_f,final_Raman_preR,'linewidth',2);
set(h(1),'Color',cc(1,:)); set(h(2),'Color',cc(2,:)); set(h(3),'Color','k');
xlabel('Raman frequency (THz)');
ylabel('R^{coeff}');
xlim([0,40]);
ylim([0,8.5e-24]);
set(gca,'fontsize',25);
print(gcf,'Raman coeff (rot).pdf','-dpdf');

pos = get(gca,'Position');
pos2 = pos(2);
pos4 = pos(4);

figure;
plot(final_Raman_f,final_Raman_preR_H2,'linewidth',2,'Color',cc(1,:));
xlim([122,125]);
ylim([0,8.5e-24]);
set(gca,'fontsize',25,'YColor','None','YTick',[]);
pos = get(gca,'Position');
set(gca,'Position',[pos(1),pos2,pos(3)/3,pos4]);
print(gcf,'Raman coeff (H2 vib).pdf','-dpdf');

figure;
plot(final_Raman_f,final_Raman_preR_O2,'linewidth',2,'Color',cc(2,:));
xlim([46,47]);
ylim([0,8.5e-24]);
set(gca,'fontsize',25,'YColor','None','YTick',[]);
pos = get(gca,'Position');
set(gca,'Position',[pos(1),pos2,pos(3)/3,pos4]);
print(gcf,'Raman coeff (O2 vib).pdf','-dpdf');

figure;
plot(final_Raman_f,final_Raman_preR_N2,'linewidth',2,'Color','k');
xlim([69.5,70]);
ylim([0,8.5e-24]);
set(gca,'fontsize',25,'YColor','None','YTick',[]);
pos = get(gca,'Position');
set(gca,'Position',[pos(1),pos2,pos(3)/3,pos4]);
print(gcf,'Raman coeff (N2 vib).pdf','-dpdf');
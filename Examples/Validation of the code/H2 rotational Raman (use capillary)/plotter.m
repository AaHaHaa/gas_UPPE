clearvars; close all;

addpath('../../../user_helpers');

filename = {'S1_Raman_gain_250.mat',...
            'S1_Raman_gain_500.mat',...
            'S1_Raman_gain_750.mat',...
            'S1_Raman_gain_1000.mat',...
            'S1_Raman_gain_1250.mat',...
            'S1_Raman_gain_1500.mat',...
            'S1_Raman_gain_1750.mat',...
            'S1_Raman_gain_2000.mat',...
            'S1_Raman_gain_2250.mat',...
            'S1_Raman_gain_2500.mat'};
        
coupled_pulse_energy = 250:250:2500;
rep_rate = 2; % MHz

residual_energy_all = zeros(length(filename),1);
Raman_energy_all = zeros(length(filename),1);
for i = 1:length(filename)
    load(filename{i},'residual_energy','Raman_energy');
    
    residual_energy_all(i) = residual_energy;
    Raman_energy_all(i) = Raman_energy;
end

Li_pump_power = [0.3,0.6,1.45,2.2,3,3.7,4.6,5.2]';
Li_residual_power = [0.2,0.55,1.2,1.3,1.1,1,0.9,0.95]';
Li_Raman_power = [0.05,0.07,0.1,0.5,1.3,1.9,2.5,2.9]';

figure;
name = {'Color'}; value = {[0,0,1];[1,0,0]};
h = plot(coupled_pulse_energy*rep_rate/1e3,[residual_energy_all,Raman_energy_all]*rep_rate/1e3,'linewidth',2);
set(gca,'fontsize',20);
set(h,name,value);
xlim([min(Li_pump_power),max(Li_pump_power)]);
xlabel('Coupled power (W)');
ylabel('Output power (W)');
hold on;
h2 = plot(Li_pump_power,[Li_residual_power,Li_Raman_power],'linewidth',2);
name = {'Color','LineStyle'}; value = {[0,0,1],'--';[1,0,0],'--'};
set(h2,name,value);
legend('Residual','Stokes','location','northwest');
hold off;
print('H2 rot comparison.eps','-depsc');
% As the file name, this code draws the anti-Stokes Raman gain function and
% its corresponding amplitudes w.r.t. wave-vector mismatch.
% Please see Fig. S7 for details.

close all; clearvars;

addpath('../../../../user_helpers');

Im_gamma = -[1,0.001,0.1,3,5];
gamma = 1+1i*Im_gamma;

xLimits = 10;
delta = linspace(-xLimits*10,xLimits,100000)';

%% + sign in gain
y = abs(1i*(delta/2-gamma) + (1/2*sqrt((4*gamma-delta).*delta)))./abs(gamma);

ccc = distinguishable_colors(length(Im_gamma));
ccc = ccc([1:3,5:end],:);
figure;
h = plot(delta,y,'linewidth',2);
set(h(1),'Color','k','linewidth',5);
for i = 1:size(ccc,1)
    set(h(i+1),'Color',ccc(i,:));
end
xlabel('\Delta\beta (norm.)');
ylabel('|C_+|');
xlim([-1,1]*xLimits);
ylim([0,1]);
set(gca,'fontsize',25);
set(gca,'XTick',[-20,0,20],'YTick',[0,1],'YTickLabel',{'0','1'},'XTick',[-1,0,1]*xLimits,'XTickLabel',{'-1','0','1'});
print(gcf,'Raman suppression function (anti-Stokes_+).pdf','-dpdf');

% Compare with Stokes Raman gain to see how they align with each other
y_Stokes = abs(real(1/2*sqrt((4*gamma-delta).*delta)));

ccc = distinguishable_colors(length(Im_gamma));
ccc = ccc([1:3,5:end],:);
figure;
h = plot(delta,y_Stokes,'linewidth',2);
set(h(1),'Color','k','linewidth',5);
for i = 1:size(ccc,1)
    set(h(i+1),'Color',ccc(i,:));
end
xlabel('\Delta\beta (norm.)');
ylabel('g (norm.)');
xlim([-1,1]*xLimits);
ylim([0,6]);
set(gca,'fontsize',25);
set(gca,'XTick',[-20,0,20],'YTick',[0,6],'YTickLabel',{'0','1'},'XTick',[-1,0,1]*xLimits,'XTickLabel',{'-1','0','1'});
print(gcf,'Raman suppression function (Stokes for anti-Stokes comparison).pdf','-dpdf');

%% - sign in gain
y = abs(1i*(delta/2-gamma) - (1/2*sqrt((4*gamma-delta).*delta)))./abs(gamma);

[~,ymax_idx] = max(y,[],1);

ccc = distinguishable_colors(length(Im_gamma));
ccc = ccc([1:3,5:end],:);
figure;
h = plot(delta,y,'linewidth',2);
set(h(1),'Color','k','linewidth',5);
for i = 1:size(ccc,1)
    set(h(i+1),'Color',ccc(i,:));
end
xlabel('\Delta\beta (norm.)');
ylabel('|C_-|');
xlim([-1,1]*xLimits);
ylim([0,6]);
set(gca,'fontsize',25);
set(gca,'XTick',[-20,0,20],'YTick',[0,1,6],'YTickLabel',{'0','1','6'},'XTick',[-1,0,1]*xLimits,'XTickLabel',{'-1','0','1'});
print(gcf,'Raman suppression function (anti-Stokes_-).pdf','-dpdf');
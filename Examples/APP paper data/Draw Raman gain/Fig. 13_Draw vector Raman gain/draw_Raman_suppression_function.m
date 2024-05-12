% This code draws the cross-linear Stokes Raman gain function w.r.t. 
% wave-vector mismatch.
% Please see Fig. 13 for details.

close all; clearvars;

addpath('../../../../user_helpers');

Im_gamma = -[3/4,0.01,0.5,3,5];
gamma = -3/4+1i*Im_gamma;

xLimits = 8;
delta = linspace(-xLimits*10,xLimits,10000)';
y = abs(real(1/2*sqrt((gamma-delta).*(5/3+delta))));

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
ylabel('|g^{cross-linear}| (norm.)');
xlim([-1,1]*xLimits);
ylim([0,1.3]);
set(gca,'fontsize',25);
set(gca,'XTick',[-5,0,5],'YTick',[0,1.3],'YTickLabel',{'0','1'},'XTick',[-1,0,1]*xLimits,'XTickLabel',{'-1','0','1'});
print(gcf,'Raman suppression function.pdf','-dpdf');
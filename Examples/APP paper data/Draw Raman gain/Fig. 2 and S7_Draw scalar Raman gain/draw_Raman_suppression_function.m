% This code draws the Stokes Raman gain function w.r.t. wave-vector
% mismatch.
% Please see Fig. 2 for details.

close all; clearvars;

addpath('../../../../user_helpers');

Im_gamma = -[1,0.01,0.5,3,5];
gamma = 1+1i*Im_gamma;

xLimits = 30;
delta = linspace(-xLimits*10,xLimits,10000)';
y = abs(real(1/2*sqrt((4*gamma-delta).*delta)));

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
ylabel('g (norm.)');
xlim([-1,1]*xLimits);
ylim([0,6]);
set(gca,'fontsize',25);
set(gca,'XTick',[-20,0,20],'YTick',[0,6],'YTickLabel',{'0','1'},'XTick',[-1,0,1]*xLimits,'XTickLabel',{'-1','0','1'});
print(gcf,'Raman suppression function.pdf','-dpdf');

%% Draw separate lines
figure;
h = plot(delta,y(:,2),'Color',ccc(1,:),'linewidth',10);
xlim([-1,1]*xLimits);
ylim([0,6]);
set(gca,'Color','None','XTick',[],'YTick',[],'XColor','None','YColor','None');
print(gcf,'Raman suppression function (small Im).pdf','-dpdf');

figure;
h = plot(delta,y(:,1),'Color','k','linewidth',10);
xlim([-1,1]*xLimits);
ylim([0,6]);
set(gca,'Color','None','XTick',[],'YTick',[],'XColor','None','YColor','None');
print(gcf,'Raman suppression function (medium Im).pdf','-dpdf');

figure;
h = plot(delta,y(:,end),'Color',ccc(end,:),'linewidth',10);
xlim([-1,1]*xLimits);
ylim([0,6]);
set(gca,'Color','None','XTick',[],'YTick',[],'XColor','None','YColor','None');
print(gcf,'Raman suppression function (large Im).pdf','-dpdf');
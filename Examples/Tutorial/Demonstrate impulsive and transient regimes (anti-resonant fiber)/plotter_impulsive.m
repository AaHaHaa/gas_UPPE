clearvars; close all;

addpath('../../user_helpers');

load('impulsive.mat');

permittivity0 = 8.85e-12;

epsilon = squeeze(real(prop_output.delta_permittivity(:,1,2,:))/permittivity0);

fig_size = [570,413];
save_point = size(prop_output.fields,3);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    figs = figure;
    yyaxis left;
    hr = plot(t*1e3,abs(prop_output.fields(:,:,i)).^2/1e6,'linewidth',2);
    ylabel('Power (MW)');
    set(gca,'fontsize',20);
    ylim([-1,1]*300);
    yyaxis right;
    hl = plot(t*1e3,epsilon(:,i),'linewidth',2);
    set(gca,'fontsize',20);
    ylabel('\epsilon_r');
    xlabel('Time (fs)');
    ylim([-1,1]*2e-5);
    xlim([-1,1]*50);
    
    set(gca, 'SortMethod', 'depth')
    hl.ZData = zeros(size(hl.XData));
    hr.ZData = ones(size(hr.XData));
    
    drawnow;
    set(figs,'Color',[1,1,1], 'Units','pixels');
    drawnow;
    fig_pos = get(figs,'Position');
    set(gcf,'Position', [fig_pos(1:2),fig_size]);

    Frame(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame,20);

exportVideo = VideoWriter('impulsive');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);
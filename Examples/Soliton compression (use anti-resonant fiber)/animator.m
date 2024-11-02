% Please run any script first.

log_yes = false; % Use "log_yes = true" to see, under log scale, how spectral interference generates temporal fringes.
save_point = size(prop_output.fields,3);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    [~,~,~,figs,ax] = calc_spectrogram(t,f,prop_output.fields(:,1,i),[-0.2,0.2],[800,1400],400,400,true,true,log_yes);
    set(ax(2),'YLim',[0,8e8]);
    set(figs,'Color',[1,1,1]);
    
    Frame(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame,20);

exportVideo = VideoWriter('soliton_compression');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);

%%
permittivity0 = 8.85e-12;
epsilonR = squeeze(real(prop_output.delta_permittivity(:,1,1,:))/permittivity0)/2e-5;
epsilonV = squeeze(real(prop_output.delta_permittivity(:,1,2,:))/permittivity0)/4e-6;

save_point = size(prop_output.fields,3);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    figs = figure;
    yyaxis left;
    hr = plot(t,abs(prop_output.fields(:,:,i)).^2/1e6,'linewidth',2,'Color','k');
    ylabel('Power (MW)');
    set(gca,'fontsize',20,'YColor','k');
    ylim([0,800]);
    yyaxis right;
    hl = plot(t,[epsilonR(:,i),epsilonV(:,i)],'linewidth',2);
    set(hl(1),'Color','b'); set(hl(2),'Color','r');
    set(gca,'fontsize',20);
    ylabel('\epsilon_r (norm.)');
    xlabel('Time (ps)');
    ylim([-1,1]);
    xlim([-1,1]);
    legend(hl,'rot','vib');
    
    set(gca, 'SortMethod', 'depth');
    hl(1).ZData = zeros(size(hl(1).XData));
    hl(2).ZData = zeros(size(hl(2).XData));
    hr.ZData = ones(size(hr.XData));
    
    drawnow;
    set(figs,'Color',[1,1,1], 'Units','pixels');
    drawnow;

    Frame(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame,20);

exportVideo = VideoWriter('soliton_compression_index');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);
% Please run any script first.

spectrum = abs(fftshift(ifft(prop_output.fields),1)).^2./lambda.^2;
spectrum = spectrum/max(max(max(spectrum)));

log_yes = false; % Use "log_yes = true" to see, under log scale, how spectral interference generates temporal fringes.
save_point = size(prop_output.fields,3);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    figs = figure;
    subplot(1,2,1);
    h = plot(t,abs(prop_output.fields(:,:,i)).^2/1e6,'linewidth',2);
    set(h(1),'Color','b'); set(h(2),'Color','r');
    set(gca,'fontsize',20);
    legend('x','y');
    xlabel('Time (ps)');
    ylabel('Power (MW)');
    ylim([0,30]);
    xlim([-3,3]);

    subplot(1,2,2);
    h = plot(lambda,spectrum(:,:,i),'linewidth',2);
    set(h(1),'Color','b'); set(h(2),'Color','r');
    set(gca,'fontsize',20);
    legend('x','y');
    xlabel('Wavelength (nm)');
    ylabel('PSD (norms)');
    ylim([0,1]);
    xlim([700,1800]);

    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);

    set(figs,'Color',[1,1,1]);
    
    Frame(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame,20);

exportVideo = VideoWriter('SSFS');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);
% Please run any script first.

center_lambda = 1180;
bandwidth_lambda = 300;
gaussexpo = 4;
Stokes1 = gaussian_spectral_filter(prop_output{1}, sim.f0, center_lambda, bandwidth_lambda, gaussexpo);
pump1 = gaussian_spectral_filter(prop_output{1}, sim.f0, pump_wavelength*1e9, bandwidth_lambda, gaussexpo);
permittivity0 = 8.85e-12;
epsilonV1 = squeeze(real(prop_output{1}.delta_permittivity(:,1,2,:))/permittivity0)/2e-4;

Stokes2 = gaussian_spectral_filter(prop_output{2}, sim.f0, center_lambda, bandwidth_lambda, gaussexpo);
pump2 = gaussian_spectral_filter(prop_output{2}, sim.f0, pump_wavelength*1e9, bandwidth_lambda, gaussexpo);
epsilonV2 = squeeze(real(prop_output{2}.delta_permittivity(:,1,2,:))/permittivity0)/2e-4;

save_point = size(prop_output{1}.fields,3);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    figs = figure;
    subplot(1,2,1);
    yyaxis left;
    h = plot(t,abs([pump1.fields(:,:,i),Stokes1.fields(:,:,i)]).^2/1e6,'linewidth',2);
    set(h(1),'Color','k','linestyle','-'); set(h(2),'Color','r','linestyle','-');
    set(gca,'fontsize',20);
    xlabel('Time (ps)');
    ylabel('Power (W)');
    ylim([0,35]);
    xlim([-5,5]);
    yyaxis right;
    plot(t,epsilonV1(:,i),'linewidth',2);
    ylabel('\epsilon_r (norm.)');
    ylim([-1,1]);
    legend(h,'pump','Stokes');

    subplot(1,2,2);
    yyaxis left;
    h = plot(t,abs([pump2.fields(:,:,i),Stokes2.fields(:,:,i)]).^2/1e6,'linewidth',2);
    set(h(1),'Color','k','linestyle','-'); set(h(2),'Color','r','linestyle','-');
    set(gca,'fontsize',20);
    legend('pump','Stokes');
    xlabel('Time (ps)');
    ylabel('Power (W)');
    ylim([0,35]);
    xlim([-15,15]);
    yyaxis right;
    plot(t,epsilonV2(:,i),'linewidth',2);
    ylabel('\epsilon_r (norm.)');
    ylim([-1,1]);
    legend(h,'pump','Stokes');

    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1:2),pos(3)*2,pos(4)]);

    set(figs,'Color',[1,1,1]);
    
    Frame(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame,20);

exportVideo = VideoWriter('Stokes_generation');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);
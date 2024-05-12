clearvars; close all;

addpath('../../../user_helpers');

load('photoionization_blueshift_5.2uJ.mat');

% High-resolution spectrogram
log_yes = false; % Use "log_yes = true" to see, under log scale, how spectral interference generates temporal fringes.
save_point = size(prop_output.fields,3);
Frame(save_point) = struct('cdata',[],'colormap',[]);
for i = 1:save_point
    [~,~,~,figs,ax] = calc_spectrogram(t,f,prop_output.fields(:,1,i),[-0.05,0.05],[300,1100],50,50,true,true,log_yes);
    set(figs,'Color',[1,1,1]);

    Frame(i) = getframe(figs);
    close(figs);
end
% Movie
implay(Frame,20);

exportVideo = VideoWriter('pulse_evolution');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);
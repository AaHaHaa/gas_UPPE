save('Stokes_pulse_compression.mat');

% Stokes
analyze_field( t,f,Stokes.fields(:,:,end),'Treacy-t',pi/6,1e-3/600,true );

% Movie
implay(Frame_S,10);

exportVideo = VideoWriter('Stokes');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);

% -------------------------------------------------------------------------
% pump
% Movie
implay(Frame_p,10);

exportVideo = VideoWriter('pump');
exportVideo.FrameRate = 20;
open(exportVideo);
writeVideo(exportVideo, Frame);
close(exportVideo);
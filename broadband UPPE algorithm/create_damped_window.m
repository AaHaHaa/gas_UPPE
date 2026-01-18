function func = create_damped_window(Nt)

func.f = create_damped_freq_window(Nt);
func.t = create_damped_time_window(Nt);

end

function damped_time_window = create_damped_time_window(Nt)
% This function creates a sharp damped time window to remove anything
% around the edges.
% I use a super-Gaussian function here.

t = (1:Nt)';
tc = floor(Nt/2)+1;
tfwhm = Nt*0.85;
t0 = tfwhm/(2*sqrt(log(2))); % 2*sqrt(log(2))=1.665
gexpo = 2*20; % 20 is to make it a sharp window

damped_time_window = exp(-(t-tc).^gexpo/(2*t0^gexpo)).^20; % 20 is to make it a sharp window
damped_time_window(damped_time_window>0.9) = 1;

end

function damped_freq_window = create_damped_freq_window(Nt)
% This function creates a sharp damped frequency window to remove anything
% around the edges, especially the high-frequency edge of the window.
% I use a super-Gaussian function here.

f = ifftshift((1:Nt)',1);
fc = floor(Nt/2)+1;
ffwhm = Nt*0.85;
f0 = ffwhm/(2*sqrt(log(2))); % 2*sqrt(log(2))=1.665
gexpo = 2*20; % 20 is to make it a sharp window

% A much sharper damped window is used only for the low-frequency side;
% otherwise, it'll easily remove the long-wavelength signal we want to see.
damped_freq_window = exp(-(f-fc).^gexpo/(2*f0^gexpo)).^20; % 20 is to make it a sharp window
damped_freq_window(fc:end) = 1;

end
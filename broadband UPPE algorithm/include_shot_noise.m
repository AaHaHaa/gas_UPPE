function fields = include_shot_noise(sim,real_omegas,Nt,time_window,fields)
%INCLUDE_SHOT_NOISE It adds shot noise to the fields
%
% Input arguments:
%   sim.num_photon_noise_per_bin: the number of photon noise per spectral discretization bin
%   sim.gpu_yes: whether to use GPU or not
%   real_omegas: omegas for computing the photon energy (Hz)
%   Nt: the number of numerical sampling points
%   time_window: (ps)
%   fields: the electric field (sqrt(W))
%
% -------------------------------------------------------------------------
% Unit explanation:
%   intensity = abs(field).^2;
%   energy = trapz(t,intensity) = trapz(intensity)*dt;       % pJ
%   
%   spectrum_unknown_unit = abs(fftshift(ifft(field),1)).^2;
%
%   Parseval's theorem: sum(intensity) = sum(spectrum_unknown_unit)*N;
%                       * Note that spectrum_unknown_unit is from "ifft".
%   therefore sum(intensity)*dt = sum(spectrum_unknown_unit)*N*dt
%                               = sum(spectrum_unknown_unit)*(N*dt)^2/(N*dt)
%                               = sum(spectrum_unknown_unit)*(N*dt)^2*df
%                               = sum(spectrum_f)*df
%
%   spectrum_f = spectrum_unknown_unit*(N*dt)^2;
%   energy = trapz(f,spectrum_f) = trapz(spectrum_f)*df      % pJ
%                                = trapz(spectrum_f)/(N*dt);
%                                = trapz(spectrum_unknown_unit)*(N*dt)
%
% -------------------------------------------------------------------------
%   c = 299792.458;     % nm/ps
%   wavelength = c./f;  % nm
%   spectrum_wavelength = spectrum_f.*(c./wavelength.^2);
%   energy = -trapz(wavelength,spectrum_wavelength);         % pJ
% -------------------------------------------------------------------------
%   The noise is added as "one photon per frequency bin," so each
%   frequency bin adds one photon energy, hbar*omega, to the total energy.
%   This corresponds to its "spectrum_unknown_unit" counterpart as
%   "hbar*omega/(N*dt)," whose field amplitude is
%   "sqrt(hbar*omega/(N*dt))."
%
%   If the frequency window is small, assume that all omegas=omegas0,
%   adding photon noise adds N*hbar*omegas0 to the total energy,
%   proportional to the number of number of numerical sampling points.
% -------------------------------------------------------------------------
%   Note that the implementation of having photon noise being
%   hbar*omegas/(N*dt) relies on having "ifft" as Fourier Transform. The
%   relation might be different when Fourier Transform becomes "fft"
%   because they have different constants to divide in the integral
%   relations.
% -------------------------------------------------------------------------

if sim.num_photon_noise_per_bin ~= 0
    hbar = 6.62607015e-34/2/pi; % J*s
    photon_noise_intensity = hbar*real_omegas/(time_window*1e-12)*sim.num_photon_noise_per_bin;
    % I use analytical-signal representation for solving UPPE, so the field covers only the positive frequencies.
    photon_noise_intensity(photon_noise_intensity<0) = 0; % no noise at negative frequencies
    photon_noise_amplitude = sqrt(photon_noise_intensity);
    if sim.gpu_yes
        rand_phase = exp(1i*2*pi*rand(Nt,size(fields,2),'gpuArray'));
    else
        rand_phase = exp(1i*2*pi*rand(Nt,size(fields,2)));
    end
    fields = fft(ifft(fields) + photon_noise_amplitude.*rand_phase);

    clear hbar photon_noise_intensity photon_noise_amplitude rand_phase;
end

end


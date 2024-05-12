function filter_freq = find_filter_Raman_solitoin(freq,spectrum,varargin)
%FIND_FILTER_RAMAN_SOLITOIN Summary of this function goes here
%   Detailed explanation goes here

%%
% Accept only 2 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 2
    error('find_filter_Raman_solitoin:TooManyInputs', ...
          'It takes only at most 2 optional inputs');
end

% Set defaults for optional inputs
thresholding = 0.1;
smoothing = ceil(length(freq)/500);
optargs = {thresholding,smoothing};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[thresholding,smoothing] = optargs{:};

%%
if smoothing > 1
    smoothing_kernel = ones(smoothing,1)/100;
    spectrum = conv(spectrum,smoothing_kernel,'same');
end
spectrum(spectrum < max(spectrum)*thresholding) = 0;

first_nonzero = find(spectrum ~= 0,1);
first_zero_after_peak = find(spectrum(first_nonzero:end) == 0,1) + first_nonzero - 1;
first_nonzero_after_peak = find(spectrum(first_zero_after_peak:end) ~= 0,1) + first_zero_after_peak - 1;
if isempty(first_nonzero_after_peak) % there's only one peak in the spectrum
    first_nonzero_after_peak = length(spectrum);
end

filter_freq = freq(floor((first_zero_after_peak + first_nonzero_after_peak)/2));

end


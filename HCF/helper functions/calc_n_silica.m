function target_n = calc_n_silica(wavelength,use_gpu,varargin)
%
%   wavelength: (nm)

%% Default optional input arguments
% Accept only 4 optional inputs at most
numvarargs = length(varargin);
if numvarargs > 4
    error('calc_n_silica:TooManyInputs', ...
          'It takes only at most 4 optional inputs');
end

% Set defaults for optional inputs
cuda_dir_path = '';
diff_order = 0;
optargs = {cuda_dir_path,diff_order};

% Now put these defaults into the valuesToUse cell array, 
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[cuda_dir_path,diff_order] = optargs{:};

%% Correct "wavelength" dimension
flip_yes = false;
if wavelength(end) < wavelength(1)
    wavelength = flipud(wavelength);
    flip_yes = true;
end

%% Experimental refractive index
data = dlmread('n_silica (Palik Handbook).csv',',',1,0);
%data = dlmread('Yang.csv',',',1,0);
[~,ui] = unique(data(:,1));
data_wl = data(ui,1); % um
data_n = data(ui,2);
data_k = data(ui,3);

%% Interpolation
if use_gpu
    wl_all = cell(1,diff_order+1);
    n_all = cell(1,diff_order+1);
    k_all = cell(1,diff_order+1);
    wl_all{1} = data_wl;
    n_all{1} = data_n;
    k_all{1} = data_k;
    for i = 1:diff_order
        wl_all{i+1} = (wl_all{i}(1:end-1) + wl_all{i}(2:end))/2;
        n_all{i+1} = diff(n_all{i})./diff(wl_all{i});
        k_all{i+1} = diff(k_all{i})./diff(wl_all{i});

        [n_all{i+1},k_all{i+1}] = smooth_nk(false,'',wl_all{i+1},n_all{i+1},k_all{i+1},5,3);
    end

    data_wl = wl_all{end};
    for i = 1:diff_order+1
        n_all{i} = interp1(wl_all{i},n_all{i},data_wl,'pchip');
        k_all{i} = interp1(wl_all{i},k_all{i},data_wl,'pchip');
    end

    target_n = complex(zeros(length(wavelength),1,'gpuArray'));
    n_all = cell2mat(n_all);
    k_all = cell2mat(k_all);
    cuda_mySpline = setup_kernel('mySpline',cuda_dir_path,length(wavelength));
    target_n = feval(cuda_mySpline,...
                                  data_wl,n_all+1i*k_all,uint32(length(data_wl)),...
                                  wavelength*1e-3,target_n,uint32(length(wavelength)),...
                                  uint32(diff_order));
else
    % If no GPU, I use MATLAB internal interpolation "pchip".
    % I'm too lazy to write another CCU version of the higher-order spline interpolation
    target_n = interp1(data_wl,data_n+1i*data_k,wavelength*1e-3,'pchip');
end

%% Smooth the final data
% This might be ignorable.
% Just some random smoothing of the refractive index.

if use_gpu
    %target_n = gather(target_n);
    wavelength = gpuArray(wavelength);
end
for i = 1:3
    [Dn,Dk] = smooth_nk(use_gpu,cuda_dir_path,wavelength(2:end)*1e-3,diff(real(target_n))./diff(wavelength),diff(imag(target_n))./diff(wavelength),max(ceil(length(wavelength)/2000),3),1);
    target_n = cumtrapz(wavelength,[0;Dn+1i*Dk]) + target_n(1);
end
for i = 1:5
    real_n = smooth(real(target_n),ceil(length(wavelength)/1000));
    imag_n = smooth(imag(target_n),ceil(length(wavelength)/1000));
    target_n = real_n + 1i*imag_n;
end

smooth_wl = wavelength*1e-3 > 2 & wavelength*1e-3 < 6;
for i = 1:5
    real_n(smooth_wl) = smooth(real(target_n(smooth_wl)),ceil(length(wavelength)/1000));
    imag_n(smooth_wl) = smooth(imag(target_n(smooth_wl)),ceil(length(wavelength)/1000));
    target_n = real_n + 1i*imag_n;
end

idx = imag(target_n)<0;
target_n(idx) = real(target_n(idx)); % loss can't be negative

if flip_yes
    target_n = flipud(target_n);
end

end

function [n,k] = smooth_nk(use_gpu,cuda_dir_path,wl,n,k,smooth_range,smooth_rep)

smooth_wl = wl > 0.5;

if use_gpu
    cuda_mySmooth = setup_kernel('mySmooth',cuda_dir_path,length(wl));
end

if smooth_range ~= 0
    for i = 1:smooth_rep
        if use_gpu
            smoothing_nk = complex(zeros(sum(smooth_wl),1,'gpuArray'));
            smoothing_nk = feval(cuda_mySmooth,...
                                              smoothing_nk,...
                                              wl(smooth_wl),n(smooth_wl)+1i*k(smooth_wl),...
                                              abs(smooth_range/2/length(wl)*(max(wl)-min(wl))),uint32(smooth_range),uint32(sum(smooth_wl)));
            n(smooth_wl) = real(smoothing_nk);
            k(smooth_wl) = imag(smoothing_nk);
        else
            n(smooth_wl) = smooth(wl(smooth_wl),n(smooth_wl),smooth_range,'sgolay',2);
            k(smooth_wl) = smooth(wl(smooth_wl),k(smooth_wl),smooth_range,'sgolay',2);
        end
    end
end

end
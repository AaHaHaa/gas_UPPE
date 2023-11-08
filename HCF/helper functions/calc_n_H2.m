function target_n = calc_n_H2(wavelength,cuda_dir_path,diff_order)
%
%   wavelength: (nm)

flip_yes = false;
if wavelength(end) < wavelength(1)
    wavelength = flipud(wavelength);
    flip_yes = true;
end

%% Experimental refractive index
data = dlmread('n_H2 (Peck).csv',',',1,0);
[~,ui] = unique(data(:,1));
data_wl = data(ui,1); % um
data_n = data(ui,2);

%% Sellmeier formula
wavenumber_R = 1./data_wl;
refractivity = @(wavenumber) 14895.6./(180.7-wavenumber.^2) + 4903.7./(92-wavenumber.^2); % 10^6(n-1)
%n_Sellmeier = refractivity(wavenumber_R)/1e6 + 1;

%% Extended range of wavelengths
wavelength_L = linspace(0.01,0.167,1000)';
wavelength_R = [linspace(1.7,5,10),linspace(5.1,60,20)]'; % um
%wavenumber_L = 1./wavelength_L;
wavenumber_R = 1./wavelength_R;

%n_Sellmeier_calc_L = refractivity(wavenumber_L)/1e6 + 1;
n_Sellmeier_calc_R = refractivity(wavenumber_R)/1e6 + 1;

slm = slmengine([data_wl;wavelength_R],[data_n;n_Sellmeier_calc_R],'knots',[data_wl;wavelength_R],... 
   'decreasing','on','plot','off','extrapolation','cubic','jerk','negative');

n_data_calc_R = slmeval(wavelength_R,slm,0); % refractive index of H2
n_data_calc_L = slmeval(wavelength_L,slm,0); % refractive index of H2

% use the new range
data_wl = [wavelength_L;data_wl;wavelength_R];
data_n = [n_data_calc_L;data_n;n_data_calc_R];

%% Interpolation
wl_all = cell(1,diff_order+1);
n_all = cell(1,diff_order+1);
wl_all{1} = data_wl;
n_all{1} = data_n;
for i = 1:diff_order
    wl_all{i+1} = (wl_all{i}(1:end-1) + wl_all{i}(2:end))/2;
    n_all{i+1} = diff(n_all{i})./diff(wl_all{i});
    
    n_all{i+1} = smooth_n(false,'',wl_all{i+1},n_all{i+1},5,3);
end

data_wl = wl_all{end};
for i = 1:diff_order+1
    n_all{i} = interp1(wl_all{i},n_all{i},data_wl,'pchip');
end

target_n = complex(zeros(length(wavelength),1,'gpuArray'));
n_all = cell2mat(n_all);
cuda_mySpline = setup_kernel('mySpline',cuda_dir_path,length(wavelength));
target_n = feval(cuda_mySpline,...
                              data_wl,n_all,uint32(length(data_wl)),...
                              wavelength*1e-3,target_n,uint32(length(wavelength)),...
                              uint32(diff_order));

idx = imag(target_n)<0;
target_n(idx) = real(target_n(idx)); % loss can't be negative

%% Smooth the final data
for i = 1:3
    Dn = smooth_n(use_gpu,cuda_dir_path,wavelength(2:end)*1e-3,diff(real(target_n))./diff(wavelength),max(ceil(length(wavelength)/2000),3),1);
    target_n = cumtrapz(wavelength,[0;Dn]) + target_n(1);
end
for i = 1:5
    target_n = smooth(real(target_n),ceil(length(wavelength)/1000));
end

if flip_yes
    target_n = flipud(target_n);
end

end

function n = smooth_n(use_gpu,cuda_dir_path,wl,n,smooth_range,smooth_rep)

smooth_wl = wl > 0.5;

if use_gpu
    cuda_mySmooth = setup_kernel('mySmooth',cuda_dir_path,length(wl));
end

if smooth_range ~= 0
    for i = 1:smooth_rep
        if use_gpu
            smoothing_n = complex(zeros(sum(smooth_wl),1,'gpuArray'));
            smoothing_n = feval(cuda_mySmooth,...
                                              smoothing_n,...
                                              wl(smooth_wl),complex(n(smooth_wl)),...
                                              abs(smooth_range/2/length(wl)*(max(wl)-min(wl))),uint32(smooth_range),uint32(sum(smooth_wl)));
            n(smooth_wl) = real(smoothing_n);
        else
            n(smooth_wl) = smooth(wl(smooth_wl),n(smooth_wl),smooth_range,'sgolay',2);
        end
    end
end

end
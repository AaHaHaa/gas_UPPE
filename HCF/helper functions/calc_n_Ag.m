function target_n = calc_n_Ag(wavelength,use_gpu,cuda_dir_path,diff_order)
%
%   wavelength: (nm)

flip_yes = false;
if wavelength(end) < wavelength(1)
    wavelength = flipud(wavelength);
    flip_yes = true;
end

%% Experimental refractive index
data = dlmread('n_Ag (Yang).csv',',',1,0);
%data = dlmread('Yang.csv',',',1,0);
[~,ui] = unique(data(:,1));
data_wl = data(ui,1); % um
data_n = data(ui,2);
data_k = data(ui,3);

%% Smooth
[Ddata_n,Ddata_k] = smooth_nk(data_wl(2:end),diff(data_n)./diff(data_wl),diff(data_k)./diff(data_wl),150,3);
data_n = cumtrapz(data_wl,[0;Ddata_n]) + data_n(1);
data_k = cumtrapz(data_wl,[0;Ddata_k]) + data_k(1);

added_wl = 50;
data_n = interp1(data_wl,data_n,[data_wl;added_wl],'pchip','extrap');
data_k = interp1(data_wl,data_k,[data_wl;added_wl],'linear','extrap');
data_wl = [data_wl;50];

%% Interpolation
if use_gpu
    target_n = myPchip( data_wl,data_n+1i*data_k,...
                        wavelength*1e-3,...
                        diff_order,cuda_dir_path );
    target_n = gather(target_n);
else
    % If no GPU, I use MATLAB internal interpolation "pchip".
    % I'm too lazy to write another CPU version of the higher-order spline interpolation
    target_n = interp1(data_wl,data_n+1i*data_k,wavelength*1e-3,'pchip');
end

idx = imag(target_n)<0;
target_n(idx) = real(target_n(idx)); % loss can't be negative

if flip_yes
    target_n = flipud(target_n);
end

end

function [n,k] = smooth_nk(wl,n,k,smooth_range,smooth_rep)

smooth_wl = wl > 0.5;

for i = 1:smooth_rep
    n(smooth_wl) = smooth(wl(smooth_wl),n(smooth_wl),smooth_range,'sgolay',2);
    k(smooth_wl) = smooth(wl(smooth_wl),k(smooth_wl),smooth_range,'sgolay',2);
end

end
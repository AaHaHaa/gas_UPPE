clearvars; close all;

% (n,m) modes to solve for
k = 1;
for i = 1:100
    for j = 1:40
        nm(k,:) = [i,j];
        k = k+1;
    end
end

wavelength = 1000; % nm

a = 12.5e-6; % m

% refractive index
load('n_H2.mat');

%%
if size(wavelength,1) == 1 % make it column vector
    wavelength = wavelength.';
end
cuda_dir_path = '../cuda';
diff_order = 6;
n_in = calc_n_silica([wavelength-10:wavelength+10]',cuda_dir_path,diff_order); % refractive index of H2
n_in = n_in(11);
k0 = 2*pi./wavelength*1e9;

if size(nm,2) == 2
    nm = nm.';
end

num_modes = size(nm,2);

vn = zeros(1,num_modes);
unm = zeros(1,num_modes);
mode = cell(1,num_modes);
for i = 1:num_modes
    % vn
    if nm(1,i) == 0 % only TE is considered in this code
        mode{i} = 'TE';
    else
        mode{i} = 'EH';
    end
    
    % unm: zeros of the Bessel function of the first kind
    u = besselzero(nm(1,i)-1,nm(2,i),1);
    unm(i) = u(end);
end

ki = complex(unm/a);
gamma = sqrt((k0*n_in).^2 - ki.^2);

beta = real(gamma);
%loss = imag(gamma);

[~,idx] = sort(beta,'descend');
sorted_nm = nm(:,idx);

save('nm_order.mat','sorted_nm');
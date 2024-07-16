function SR = calc_SR_tensors(mode_profiles,r,dr,dtheta,sim,current_folder)

num_modes = size(mode_profiles,2);
Nx = size(mode_profiles,3);
fields = real(permute(mode_profiles,[3,4,2,1])); % (x,y,midx); the computation below required it to be real.
r = permute(r,[3,1,2]);

SR = zeros(num_modes^4,1);

if sim.gpu_yes
    SR = gpuArray(SR);
else % CPU
    if ispc
        [~,systemview] = memory;
        mem_planned_for_SR = systemview.PhysicalMemory.Available/3*2^20;
    elseif isunix
        [~,w] = unix('free -b | grep Mem'); % Display the memory info in Bytes
        stats = str2double(regexp(w, '[0-9]*', 'match'));
        %memsize = stats(1)/1e6;
        freemem = stats(end); % availabel memory
        mem_planned_for_SR = freemem/3;
    end
    
    fields = gather(fields);
end

%% Calculate the overlap integrals

% SR will hold the tensor. We only need to calculate SR

% If using the GPU, we need to do some things differently
if sim.gpu_yes
    fields = permute(fields, [3 1 2]); % The order needs to go (num_modes, Nx, Nx)

    
    specific_filename = 'calculate_tensors_double';
    
    kernel = setup_kernel(specific_filename,[current_folder '../cuda'],num_modes^4);

    % Run the CUDA code
    SR = feval(kernel, SR, fields, int32(num_modes), int32(Nx), r, dr, dtheta);
    %wait(gd);
    SR = gather(SR);
else
    % If we're not using the GPU, then do all the calculations directly in MATLAB
    [midx1,midx2,midx3,midx4] = ind2sub(num_modes*ones(1,4),1:num_modes^4);

    precision = 8;
    
    Nf = 1;
    num_segments = max(1,ceil(Nx^2*num_modes^4*precision*Nf/mem_planned_for_SR));
    num_each_segment = ceil(num_modes^4/num_segments);
    segments = [num_each_segment*ones(1,num_segments-1) num_modes^4-num_each_segment*(num_segments-1)];
    cs = [0 cumsum(segments)];

    for segi = 1:num_segments
        s = cs(segi)+1;
        e = cs(segi+1);
        SR(s:e) = permute(sum(sum(r.*fields(:, :, midx1(s:e)).*fields(:, :, midx2(s:e)).*fields(:, :, midx3(s:e)).*fields(:, :, midx4(s:e)),1),2),[3 1 2])*dr*dtheta;
    end
end

%% Eliminate the zero elements

thresholdzero = max(SR)/1e5; % This is fairly arbitrary

zero_SR = find(abs(SR)<thresholdzero);
SR(zero_SR) = 0;
cnt = num_modes^4 - length(zero_SR);

SR = reshape(SR, [num_modes, num_modes, num_modes, num_modes]);

end
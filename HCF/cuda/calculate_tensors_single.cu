__global__ void calculate_tensors(float* QR, const float* normalized_fields, const int num_modes, const int Nx, const float* r, const float dr, const float dtheta) {
    unsigned int full_thread_idx = threadIdx.x + blockIdx.x*blockDim.x;

    // Calculate the index
    unsigned int nmp4 = num_modes*num_modes*num_modes*num_modes;
    unsigned int Nxnm = Nx*num_modes;

    if (full_thread_idx >= nmp4) {
        return;
    }

    // Turn linear index into components
    unsigned int midx1 = full_thread_idx % num_modes;
    unsigned int midx2 = (full_thread_idx/num_modes) % num_modes;
    unsigned int midx3 = (full_thread_idx/num_modes/num_modes) % num_modes;
    unsigned int midx4 = (full_thread_idx/num_modes/num_modes/num_modes);

    // Compute the sum
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Nx; j++) {
            QR[full_thread_idx] += r[i]*normalized_fields[midx1+i*num_modes+j*Nxnm]*normalized_fields[midx2+i*num_modes+j*Nxnm]*normalized_fields[midx3+i*num_modes+j*Nxnm]*normalized_fields[midx4+i*num_modes+j*Nxnm];
        }
    }

    // Normalize
    QR[full_thread_idx] = QR[full_thread_idx]*dr*dtheta;
}

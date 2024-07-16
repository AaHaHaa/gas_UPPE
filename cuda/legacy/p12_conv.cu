__global__ void p12_conv(double2* p12,
                         const unsigned int Nt,
                         const double hbar, const double B, const double dt,
                         const double2* exp12, const double* A2) {

    if (blockIdx.x >= Nt*2) return;

    const unsigned int MaxThreadsPerBlock = 1024;
    __shared__ double exp12A2[MaxThreadsPerBlock];
    exp12A2[threadIdx.x] = 0; // initialize to 0
    __syncthreads();

    const unsigned int Ni = blockIdx.x/2;
    const bool iseven = ((blockIdx.x & 1)==0);
    const unsigned int num_each_thread = ceil((double)Nt/(double)blockDim.x);

    unsigned int ni;
    if (iseven) {
        for (int i = 0; i<num_each_thread; i++) {
            ni = i + threadIdx.x*num_each_thread;
            if (ni < Ni)
                exp12A2[threadIdx.x] = exp12A2[threadIdx.x] + exp12[ni].y*A2[Ni-ni];
        }
    } else {
        for (int i = 0; i<num_each_thread; i++) {
            ni = i + threadIdx.x*num_each_thread;
            if (ni < Ni)
                exp12A2[threadIdx.x] = exp12A2[threadIdx.x] - exp12[ni].x*A2[Ni-ni];
        }
    }
    __syncthreads();
    
    if (threadIdx.x == 0) {
        if (iseven) {
            for (int i = 0; i<blockDim.x; i++)
                p12[Ni].x = p12[Ni].x + exp12A2[i];

            p12[Ni].x = B/2/hbar*dt*p12[Ni].x;
        } else {
            for (int i = 0; i<blockDim.x; i++)
                p12[Ni].y = p12[Ni].y + exp12A2[i];

            p12[Ni].y = B/2/hbar*dt*p12[Ni].y;
        }
    }
}
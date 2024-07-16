__global__ void myConv(float2* output,
                       const float* v, const float2* A,
                       const unsigned int N, const unsigned int total_size) {
    const unsigned int thread_idx = threadIdx.x + blockIdx.x*blockDim.x;

    if (thread_idx >= total_size) return; // total_size is the number of elements in A

    const unsigned int higherDim_idx = thread_idx / N;
    const unsigned int Ni = thread_idx - higherDim_idx*N; // thread_idx % Nx

    float2 conv_sum; conv_sum.x = 0;
                      conv_sum.y = 0;
    for(unsigned int i = 0; i<=Ni; i++) {
        conv_sum.x += A[thread_idx-i].x*v[i];
        conv_sum.y += A[thread_idx-i].y*v[i];
    }

    output[thread_idx].x = conv_sum.x;
    output[thread_idx].y = conv_sum.y;
}
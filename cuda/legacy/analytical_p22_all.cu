__global__ void analytical_p22_all(double* p22_integrand,
                                   const double2* p12,
                                   const double* E2, const unsigned int NE,
                                   const double* t,
                                   const double T1_R1, const double T1_R2, const double T1_V,
                                   const double* a_R1, const double* a_R2, const double* a_V,
                                   const double hbar, const double unit_dt) {
    const unsigned int thread_idx0 = threadIdx.x + blockIdx.x*blockDim.x;
    unsigned int thread_idx;

    if (thread_idx0 >= NE*3) return;

    const unsigned int Ni = thread_idx0%NE;

    if (thread_idx0 < NE) { // R1
        p22_integrand[thread_idx0] = -a_R1[2]/hbar*exp(1/T1_R1*t[Ni])*p12[thread_idx0].y*E2[thread_idx0]*unit_dt; // p22_integrand_R1
    } else if (thread_idx0 < NE*2) { // R2
        thread_idx = thread_idx0 - NE;
        p22_integrand[thread_idx0] = -a_R2[2]/hbar*exp(1/T1_R2*t[Ni])*p12[thread_idx+NE].y*E2[thread_idx]*unit_dt; // p22_integrand_R2
    } else { // V
        thread_idx = thread_idx0 - NE*2;
        p22_integrand[thread_idx0] = -a_V[2]/hbar*exp(1/T1_V*t[Ni])*p12[thread_idx+NE*2].y*E2[thread_idx]*unit_dt; // p22_integrand_V
    }
}
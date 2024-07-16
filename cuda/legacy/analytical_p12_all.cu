#include "analytical_p12_helper.cu"

__global__ void analytical_p12_all(double2* p12_integrand,
                                   const double* E2, const unsigned int NE,
                                   const double T2_R1, const double T2_R2, const double T2_V,
                                   const double Omega_R1, const double Omega_R2, const double Omega_V,
                                   const double* a_R1, const double* a_R2, const double* a_V,
                                   const double hbar,
                                   const double unit_dt, const double* t) {
    const unsigned int thread_idx = threadIdx.x + blockIdx.x*blockDim.x;

    if (thread_idx >= NE*3) return;

    if (thread_idx < NE) { // R1
        analytical_p12(p12_integrand, // p12_integrand_R1
                          thread_idx,
                          E2,
                          T2_R1,Omega_R1,a_R1,hbar,unit_dt,
                          t);
    }
    else if (thread_idx < NE*2) { // R2
        analytical_p12(p12_integrand+(NE), // p12_integrand_R2
                          thread_idx - NE,
                          E2,
                          T2_R2,Omega_R2,a_R2,hbar,unit_dt,
                          t);

    } else  { // V
        analytical_p12(p12_integrand+(NE*2), // p12_integrand_V
                          thread_idx - NE*2,
                          E2,
                          T2_V,Omega_V,a_V,hbar,unit_dt,
                          t);
    }
}
#include "calc_trapz_p12.cu"

__device__ void analytical_p12(double2* p12_integrand,
                                  const unsigned int thread_idx,
                                  const double* E2,
                                  const double T2, const double Omega, const double* a, const double hbar, const double unit_dt,
                                  const double* t) {
    const unsigned int Ni = thread_idx;

/* ------------------------
   trapezoidal rule
 * ------------------------ */
    double2 p12_tmp;
    p12_tmp = calc_trapz_p12(thread_idx,
                             Ni,
                             T2,Omega,
                             unit_dt,
                             t,E2);

    p12_integrand[thread_idx].x =  p12_tmp.y*a[2]/2/hbar;
    p12_integrand[thread_idx].y = -p12_tmp.x*a[2]/2/hbar;
}
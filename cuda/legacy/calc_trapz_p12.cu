__device__ double2 calc_trapz_p12(const unsigned int thread_idx,
                                  const unsigned int Ni,
                                  const double T2, const double Omega,
                                  const double unit_dt,
                                  const double* t, const double* E2) {
    double2 p12_integrand;

    p12_integrand.x =  exp(1/T2*t[Ni])*cos(Omega*t[Ni])*E2[thread_idx]*unit_dt;
    p12_integrand.y = -exp(1/T2*t[Ni])*sin(Omega*t[Ni])*E2[thread_idx]*unit_dt;

    return p12_integrand;
}
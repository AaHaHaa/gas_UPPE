__global__ void get_GNLSE_nonlinearity_air(double2* PK, double2* P_p12,
                                            const double2* A_t,
                                            const double permittivity0,
                                            const double Ng_N2,
                                            const double B_N2, const double2* p12_R_N2,
                                            const double a_N2, const double2* p12_V_N2,
                                            const double Ng_O2,
                                            const double B_O2, const double2* p12_R_O2,
                                            const double a_O2, const double2* p12_V_O2,
                                            const unsigned int N,
                                            const double* DW) {
    const unsigned int thread_idx = threadIdx.x + blockIdx.x*blockDim.x;
    const unsigned int parallel_num = 10; // the number of parallelization
    const unsigned int Nf_in_block = 102;// floor(MaxThreadsPerBlock / parallel_num)
    const unsigned int Ni_in_block = threadIdx.x / parallel_num;
    const unsigned int parallel_idx = thread_idx % parallel_num;
    const unsigned int Ni = thread_idx / parallel_num;

    if (thread_idx >= parallel_num*N) return;

    __shared__ double I[Nf_in_block];
    __shared__ double2 A2[Nf_in_block];

    const double2 A = A_t[Ni];
    /* 
    *  I = abs(A)^2 + abs(A_p)^2
    *  A2 = A^2 + A_p^2
    */
    if (parallel_idx == 0) {
        I[Ni_in_block] = pow(A.x,2) + pow(A.y,2);
            A2[Ni_in_block].x = pow(A.x,2) - pow(A.y,2);
            A2[Ni_in_block].y = 2*A.x*A.y;
    }
    __syncthreads();

    // Kerr
    double2 conjA; conjA.x = A.x; conjA.y = -A.y;
    double2 K1;
    // Raman
    __shared__ double2 this_p12_R_N2[Nf_in_block];
    __shared__ double2 this_p12_V_N2[Nf_in_block];
    __shared__ double2 this_p12_R_O2[Nf_in_block];
    __shared__ double2 this_p12_V_O2[Nf_in_block];
    switch (parallel_idx) {
        // N2's R
        case 0: // p12_R_N2
            this_p12_R_N2[Ni_in_block] = p12_R_N2[Ni];
            break;

        // N2's V
        case 1: // p12_V_N2
            this_p12_V_N2[Ni_in_block] = p12_V_N2[Ni];
            break;

        // O2's R
        case 2: // p12_R_O2
            this_p12_R_O2[Ni_in_block] = p12_R_O2[Ni];
            break;

        // O2's V
        case 3: // p12_V_O2
            this_p12_V_O2[Ni_in_block] = p12_V_O2[Ni];
            break;
    }
    __syncthreads();

    __shared__ double2 PR_p12_N2[Nf_in_block];
    __shared__ double2 PV_p12_N2[Nf_in_block];
    __shared__ double2 PR_p12_O2[Nf_in_block];
    __shared__ double2 PV_p12_O2[Nf_in_block];

    switch (parallel_idx) {
        /* Kerr term */
        // X3 is taken out because it can be the function of frequency, so it's considered later.
        case 0: // x
            K1.x = conjA.x*A2[Ni_in_block].x - conjA.y*A2[Ni_in_block].y;
            PK[Ni].x = permittivity0/4*(K1.x + 2*A.x*I[Ni_in_block]);
            break;
        case 1: // y
            K1.y = conjA.x*A2[Ni_in_block].y + conjA.y*A2[Ni_in_block].x;
            PK[Ni].y = permittivity0/4*(K1.y + 2*A.y*I[Ni_in_block]);
            break;

        /* stimulated Raman term: N2's R */
        case 2:
            PR_p12_N2[Ni_in_block].x = Ng_N2*B_N2*(2*this_p12_R_N2[Ni_in_block].x)*A.x;
            break;
        case 3:
            PR_p12_N2[Ni_in_block].y = Ng_N2*B_N2*(2*this_p12_R_N2[Ni_in_block].x)*A.y;
            break;

        /* stimulated Raman term: N2's V */
        case 4:
            PV_p12_N2[Ni_in_block].x = Ng_N2*a_N2*(2*this_p12_V_N2[Ni_in_block].x)*A.x;
            break;
        case 5:
            PV_p12_N2[Ni_in_block].y = Ng_N2*a_N2*(2*this_p12_V_N2[Ni_in_block].x)*A.y;
            break;

        /* stimulated Raman term: O2's R */
        case 6:
            PR_p12_O2[Ni_in_block].x = Ng_O2*B_O2*(2*this_p12_R_O2[Ni_in_block].x)*A.x;
            break;
        case 7:
            PR_p12_O2[Ni_in_block].y = Ng_O2*B_O2*(2*this_p12_R_O2[Ni_in_block].x)*A.y;
            break;

        /* stimulated Raman term: O2's V */
        case 8:
            PV_p12_O2[Ni_in_block].x = Ng_O2*a_O2*(2*this_p12_V_O2[Ni_in_block].x)*A.x;
            break;
        case 9:
            PV_p12_O2[Ni_in_block].y = Ng_O2*a_O2*(2*this_p12_V_O2[Ni_in_block].x)*A.y;
            break;
    }
    __syncthreads();

    // Sum Raman terms up
    switch (parallel_idx) {
        /* p12 */
        case 0:
            P_p12[Ni].x = (PR_p12_N2[Ni_in_block].x + PV_p12_N2[Ni_in_block].x + \
                           PR_p12_O2[Ni_in_block].x + PV_p12_O2[Ni_in_block].x   )*DW[Ni];
            break;
        case 1:
            P_p12[Ni].y = (PR_p12_N2[Ni_in_block].y + PV_p12_N2[Ni_in_block].y + \
                           PR_p12_O2[Ni_in_block].y + PV_p12_O2[Ni_in_block].y   )*DW[Ni];
            break;
    }
}
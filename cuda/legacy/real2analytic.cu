// the following two lines are to use the M_PI constant=3.14...
#define _USE_MATH_DEFINES
#include <math.h>

__device__ void find_as(const unsigned int thread_idx,
                       double2* A, const double2* A0,
                       const double* f, const double* f0,
                       const unsigned int Nt, const unsigned int Nt0) {
    const unsigned int Ni = thread_idx % Nt;
    const double this_f = f[Ni];

    int base_idx0;
    if (this_f >= 0) {
    	base_idx0 = -(int)Nt0/2; // Nt0/2 is due to the frequency sequence (ifftshift).
                            // A of the positive frequency is at the first half of the matrix
                            // while f is sorted.
    } else {
        base_idx0 =  (int)Nt0/2;
    }

    unsigned int upper_idx, lower_idx;
    unsigned int idx;
    double target_f;
    unsigned int min_idx = 0;
    unsigned int max_idx = Nt0-1;
    if (this_f == f0[min_idx]) {
        A[thread_idx].x = A0[base_idx0+min_idx].x*2;
        A[thread_idx].y = A0[base_idx0+min_idx].y*2;
        return;
    } else if (this_f <f0[min_idx]) {
        A[thread_idx].x = 0;
        A[thread_idx].y = 0;
        return;
    }
    if (this_f == f0[max_idx]) {
        A[thread_idx].x = A0[base_idx0+max_idx].x*2;
        A[thread_idx].y = A0[base_idx0+max_idx].y*2;
        return;
    } else if (this_f > f0[max_idx]) {
        A[thread_idx].x = 0;
        A[thread_idx].y = 0;
        return;
    }
    while (true) {
        idx = floor(((double)max_idx+(double)min_idx)/2);
        target_f = f0[idx];
        if (target_f > this_f) {
            if (this_f > f0[idx-1]) {
                upper_idx = idx;
                lower_idx = idx-1;
                break;
            } else if (this_f == f0[idx-1]) {
                A[thread_idx].x = A0[base_idx0+idx-1].x*2;
                A[thread_idx].y = A0[base_idx0+idx-1].y*2;
                return;
            } else {
                max_idx = idx;
            }
        } else if (target_f < this_f) {
            if (this_f < f0[idx+1]) {
                upper_idx = idx+1;
                lower_idx = idx;
                break;
            } else if (this_f == f0[idx+1]) {
                A[thread_idx].x = A0[base_idx0+idx+1].x*2;
                A[thread_idx].y = A0[base_idx0+idx+1].y*2;
                return;
            } else {
                min_idx = idx;
            }
        } else { // target_f == this_f
            A[thread_idx].x = A0[base_idx0+idx].x*2;
            A[thread_idx].y = A0[base_idx0+idx].y*2;
            return;
        }
    }
    int base_idx0_upper = base_idx0;
    int base_idx0_lower = base_idx0;
    if (f0[upper_idx] >= 0 && f0[lower_idx] < 0) {
        base_idx0_upper = -(int)Nt0/2;
        base_idx0_lower =  (int)Nt0/2;
    }

    // coefficients for interpolation
    const double c1 = (this_f-f0[lower_idx])/(f0[upper_idx]-f0[lower_idx]);
    const double c2 = (f0[upper_idx]-this_f)/(f0[upper_idx]-f0[lower_idx]);

    const double r0_upper = sqrt(pow(A0[base_idx0_upper+upper_idx].x,2) + pow(A0[base_idx0_upper+upper_idx].y,2));
    const double r0_lower = sqrt(pow(A0[base_idx0_lower+lower_idx].x,2) + pow(A0[base_idx0_lower+lower_idx].y,2));
    double angle0_upper = atan2(A0[base_idx0_upper+upper_idx].y,A0[base_idx0_upper+upper_idx].x);
    double angle0_lower = atan2(A0[base_idx0_lower+lower_idx].y,A0[base_idx0_lower+lower_idx].x);
    const double tol = 0.5;
    if ((angle0_upper - angle0_lower) >= 2*M_PI*tol) {
        angle0_upper -= 2*M_PI;
    } else if ((angle0_lower - angle0_upper) >= 2*M_PI*tol) {
        angle0_upper += 2*M_PI;
    }

    const double target_r = r0_upper*c1 + r0_lower*c2;
    const double target_angle = angle0_upper*c1 + angle0_lower*c2;
    // Multiplication by two is due to analytical-signal representation
    A[thread_idx].x = target_r*cos(target_angle)*2;
    A[thread_idx].y = target_r*sin(target_angle)*2;
}

__global__ void real2analytic(double2* Pas, const double2* Pr,
                              const double* f, const double* f0,
                              const unsigned int Nt, const unsigned int Nt0) {
    const unsigned int thread_idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (thread_idx >= Nt) return;

	// p12
        find_as(thread_idx,
                Pas, Pr,
                f, f0,
                Nt, Nt0);
}
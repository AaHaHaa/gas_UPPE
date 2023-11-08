/* ------------------------------------------------------------------------
* ------------------------------- mySmooth --------------------------------
* -------------------------------------------------------------------------
* During calculation of propagation constants, material data is extracted from experiments.
* As a result, it is important to smooth the data; otherwise, the higher-order derivatives will oscillate strongly,
* making group velocity, dispersion, third-order dispersion, and so on, untrustable.
* 
* This function is the cuda implementation of the Gaussian smoothing.
* Smoothing is done by averaging the neighboring "smooth_n" points with
* their weight determined by the following Gaussian function:
* weight = exp(-(x-x0)^2/sigma^2)
*
* smooth_y[i] = ( sum of "j" from (-n) to n: weight*y[i+j] ) / total_weight of averaging neighbors
*
* Input related arguments:
*    smooth_y:  (N, 1); the smoothed output y
*    x: (N, 1); x coordinates
*    y: (N, 1); y values of each x coordinate
*    sigma: a scalar integer; the width of the Gaussian smoothing function
*    smooth_n: a scalar integer; the total number of averaging neighbors (=2*n+1)
*    num: a scalar integer; the length of array (=N)
* -----------------------------------------------------------------------*/
__global__ void mySmooth(double2* smooth_y,
                         const double* x, const double2* y,
                         const double sigma, const unsigned int smooth_n, const unsigned int num) {
    const unsigned int full_thread_idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (full_thread_idx >= num) return;

    double2 smooth_result;
    smooth_result.x = 0; smooth_result.y = 0;
    const unsigned int this_smooth_n = int(smooth_n/2)*2+1; // smooth_n needs to be odd
    const int n = (this_smooth_n-1)/2;
    const double sigma2 = sigma*sigma;
    double weight, total_weight; total_weight = 0;
    for(int i = -n; i <= n; i++) {
        if ( ((int(full_thread_idx)+i) >= 0) && ((int(full_thread_idx)+i) < num)) {
            weight = exp(-pow(x[full_thread_idx+i]-x[full_thread_idx],2.0)/sigma2);
            smooth_result.x += weight*y[full_thread_idx+i].x;
            smooth_result.y += weight*y[full_thread_idx+i].y;

            total_weight += weight;
        }
    }
    smooth_y[full_thread_idx].x = smooth_result.x/total_weight;
    smooth_y[full_thread_idx].y = smooth_result.y/total_weight;
}
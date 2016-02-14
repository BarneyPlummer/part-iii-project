// ======================================================================
// fourier.hpp
// ======================================================================

// include standard C/C++ headers
#include <cmath>
#include <iostream>

// include FFTW headers
#include <fftw3.h>

// ********************************************************************** //
// class FFTW
//
// class to call FFTW3 library
//
// Initialise the class with the number of points in the DFT and also
// the patience with which the FFTW3 planner will look for the optimal
// solution (FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT or
// FFTW_EXHAUSTIVE)
//
// ********************************************************************** //
class FFTW
{
  int npts_in;
  int npts_fft;
  
  // variables for forward fft (real time series -> complex FFT)
  double *data_fwd_r;
  fftw_complex *result_fwd_r;
  fftw_plan plan_fwd_r;

  // variables for inverse fft
  fftw_complex *data_inv_r;
  double *result_inv_r;
  fftw_plan plan_inv_r;

  // private functions
  double acorr_zero_lag_r(double *func, int npts);


public:
  FFTW();
  FFTW(int npts, int plan_type);
  ~FFTW();
  void init(int npts, int plan_type);
  fftw_complex * fft_real(double *ts, int npts, int *npts_cplx);
  double *ifft_real(fftw_complex *fft_in, int npts_cplx);
  double *cross_correlate(double *ts1, int npts1, double *ts2, int npts2);
  float *cross_correlate(float *ts1, int npts1, float *ts2, int npts2);

};

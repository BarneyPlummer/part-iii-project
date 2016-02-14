// ======================================================================
// fourier.cpp
//
// Author: Jamie Barron
// Date: February 2011
// ======================================================================

#include "fourier.hpp"

using namespace std;

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

// ---------------------------------------------------------------------- //
// constructor -- initialise memory and fft plans
// ---------------------------------------------------------------------- //
FFTW::FFTW()
{
  npts_in = 0;
  npts_fft = 0;
  
  data_fwd_r = NULL;
  result_fwd_r = NULL;
  plan_fwd_r = NULL;

  data_inv_r = NULL;
  result_inv_r = NULL;
  plan_inv_r = NULL;

}

FFTW::FFTW(int npts, int plan_type)
{
  init(npts,plan_type);
}

void
FFTW::init(int npts, int plan_type)
{
  
  // store npts for Fourier transform, and calculate output
  // number of points (the complex result if folded because there 
  // is redundancy with real input out[i] = *out[n-i], i.e. we only 
  // calculate the part of the FT giving out[i] up to (n/2+1) because 
  // we can calculate the rest from this Hermitian symmetry
  npts_in = npts;
  npts_fft = npts/2 + 1;

  // allocate forward arrays and plan
  data_fwd_r = (double *) fftw_malloc(sizeof(double)*npts_in);
  result_fwd_r = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*npts_fft);
  plan_fwd_r = fftw_plan_dft_r2c_1d(npts_in,data_fwd_r,result_fwd_r,plan_type);

  // allocate reverse arrays and plan
  data_inv_r = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*npts_fft);
  result_inv_r = (double *) fftw_malloc(sizeof(double)*npts_in);
  plan_inv_r = fftw_plan_dft_c2r_1d(npts_in,data_inv_r,result_inv_r,plan_type);


}

// ---------------------------------------------------------------------- //
// destructor -- clean up memory
// ---------------------------------------------------------------------- //
FFTW::~FFTW()
{

  fftw_free(data_fwd_r);
  fftw_free(result_fwd_r);
  fftw_destroy_plan(plan_fwd_r);

  fftw_free(data_inv_r);
  fftw_free(result_inv_r);
  fftw_destroy_plan(plan_inv_r);

}

// ---------------------------------------------------------------------- //
// forward transform from real time series
// (returns positive frequency Fourier transform -- negative frequencies
//  are the complex conjugate of the positive frequencies, i.e. 
//  FT(f) = FT(-f)* ).
// ---------------------------------------------------------------------- //
fftw_complex*
FFTW::fft_real(double *ts, int npts, int *npts_cplx)
{
  
  // check number of points does not exceed that set up for the FT
  if ( npts > npts_in )
    {
      cerr << "Error: number of points in input is less than that set up "
	"in FFT plan" << endl;
      return NULL;
    }

  // pad with zeroes to length of Fourier transform, otherwise
  for ( int i=0; i < npts_in; i++ )
    {
      if ( i < npts ) data_fwd_r[i] = ts[i];
      else data_fwd_r[i] = 0.0;
    }
  *npts_cplx = npts_fft;

  // do FFT
  fftw_execute(plan_fwd_r);

  // copy result to new array
  fftw_complex *out;
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*npts_fft);
  for ( int i=0; i < npts_fft; i++ )
    {
      out[i][0] = result_fwd_r[i][0];
      out[i][1] = result_fwd_r[i][1];
    }
  
  return out;

}


// ---------------------------------------------------------------------- //
// inverse transform from positive frequency half of FT which has 
// symmetry such that the negative frequency values are the complex
// conjugate of the corresponding positive frequency values, i.e
// FT(f) = FT(-f)*
// ---------------------------------------------------------------------- //
double *
FFTW::ifft_real(fftw_complex *fft_in, int npts_cplx)
{

  // check number of pts expected in output time series, argument npts
  // is equivalent to that expected by the Fourier transform routines
  if ( npts_cplx != npts_fft )
    {
      cerr << "Incorrect array length for inverse FFT" << endl;
      return NULL;
    }

  // copy fft data into data array
  for ( int i=0; i < npts_fft; i++ )
    {
      data_inv_r[i][0] = fft_in[i][0];
      data_inv_r[i][1] = fft_in[i][1];
    }
  
  // do FFT
  fftw_execute(plan_inv_r);

  // copy result to new array and return
  double *out;
  out = (double *) fftw_malloc(sizeof(double)*npts_in);
  for ( int i = 0; i < npts_in; i++)
    {
      out[i] = result_inv_r[i];
    }
  
  return out;

}

// ---------------------------------------------------------------------- //
// function to calculate the autocorrelation at zero lag of a function
// (real) defined at npts points.  Just integrates the function squared 
// using the trapezoidal rule
// ---------------------------------------------------------------------- //
double 
FFTW::acorr_zero_lag_r(double *func, int npts)
{

  double result;

  result = pow(func[0],2) + pow(func[npts-1],2);
  for ( int i=1; i < (npts-1); i++ )
    {
      result += 2*pow(func[i],2);
    }
  result = (double(npts)-1.0)/(2.0*double(npts))*result;

  return result;
}

// ---------------------------------------------------------------------- //
// cross-correlation
// ---------------------------------------------------------------------- //
double *
FFTW::cross_correlate(double *ts1, int npts1, double *ts2, int npts2)
{

  // calculate the autocorrelation at zero lag, i.e. the integral
  // of the time-series squared
  double cnorm1 = sqrt(acorr_zero_lag_r(ts1,npts1));
  double cnorm2 = sqrt(acorr_zero_lag_r(ts2,npts2));
  
  // take the ffts of the time series
  fftw_complex *fft1;
  fftw_complex *fft2;
  int npts_cplx; 
  fft1 = fft_real(ts1,npts1,&npts_cplx);
  fft2 = fft_real(ts2,npts2,&npts_cplx);

  // normalize each Fourier series by its autocorrelation
  for ( int i=0; i < npts_cplx; i++ )
    {
      fft1[i][0] = fft1[i][0] / cnorm1;
      fft1[i][1] = fft1[i][1] / cnorm1;

      fft2[i][0] = fft2[i][0] / cnorm2;
      fft2[i][1] = fft2[i][1] / cnorm2;
    }

  // do the cross correlation
  // in the frequency domain corr(g,h)_k = (G_k)(H_k)*
  fftw_complex *ft_ccorr;
  ft_ccorr = (fftw_complex*) fftw_malloc( sizeof(fftw_complex)*npts_cplx );
  for ( int i=0; i < npts_cplx; i++ )
    {
      ft_ccorr[i][0] = fft1[i][0]*fft2[i][0] + fft1[i][1]*fft2[i][1];
      ft_ccorr[i][1] = fft1[i][1]*fft2[i][0] - fft1[i][0]*fft2[i][1];

      // normalize for inverse FFT
      ft_ccorr[i][0] = ft_ccorr[i][0] / double(npts_in);
      ft_ccorr[i][1] = ft_ccorr[i][1] / double(npts_in);
    }

  // do the inverse transform
  double *ccorr;
  ccorr = ifft_real(ft_ccorr,npts_cplx);
  
  // clean up memory
  fftw_free(fft1);
  fftw_free(fft2);
  fftw_free(ft_ccorr);
  
  return ccorr;

}

float *
FFTW::cross_correlate(float *ts1, int npts1, float *ts2, int npts2)
{

  // convert floats to double to call the double precision cross-correlate
  // function
  double *ts1d;
  double *ts2d;
  ts1d = new double [npts1];
  ts2d = new double [npts2];
 
  for ( int i=0; i < npts1; i++ )
    ts1d[i] = double(ts1[i]);
  for ( int i=0; i < npts2; i++ )
    ts2d[i] = double(ts2[i]);

  // call cross correlation
  double *ccorrd;
  ccorrd = cross_correlate(ts1d,npts1,ts2d,npts2);
  
  // convert double array back to float
  float *ccorr;
  ccorr = new float[npts_in];
  for ( int i=0; i < npts_in; i++ )
    ccorr[i] = float(ccorrd[i]);

  // free memory
  delete [] ts1d;
  delete [] ts2d;
  fftw_free(ccorrd);

  return ccorr;

}

// **********************************************************************
// whiten.cpp
//
// Routine to spectrally whiten SAC file data
// calls:
//    fortran routine grbeamfep_ in grbeam.f
//
// Jamie Barron
// January 2010
// **********************************************************************

#include "whiten.hpp"

// 
// function whitens the signal in struct Sac *sac
// whitens whole signal based on specified noise window
//
// other arguments:
//  nstart - start of noise window used for whitening
//  nstop  - end of noise window used for whitening
//    N.B. these arguments should be the actual sample number, not 
//    the C array index, as we're using a Fortran code.
// 
int whiten_signal(struct Sac *sac,int nstrt,int nstop)
{

  // *** set up variables for grbeamfep_ Fortran routine *** //
  
  int nc = 1;                   // input number of channels 
                                //    (grbeamfep_ is multiplex)
  int np = sac->header.npts;    // number of points in time series
  int nc1 = 1;                  // number of channels in output
  int ncn[nc1];                 // array to tell subroutine which 
  ncn[0] = 1;                   //    channels to use

  int ip = 1;                   // initial point to process in input array
  int np1 = np - ip + 1;        // output number of points

  float x2[np1];                // output array section - not whitened
  float s1[nc1];                // mean values of channels
  float s[nc1];                 // variances of channels
  float s2;                     // average variance of channels
  float x[np1];                 // beam time series (when stacking channels,
                                //    not used here)
  float bm;                     // mean value of beam
  float bs;                     // variance of beam
  int m = 3;                    // order of model for autoregression whitening
  float cr[m];                  // correlation function of beam
  float af[m*(m+1)/2];          // correlation matrix of beam
  float bar[m+1];               // vector of autoregressive params of beam
  float disp = 0;               // dispersion of autoregression residuals

  float *xaw;                   // whitened data array
  xaw = new float [np1];
//  cout << endl << "We have got to here in whiten.cpp with np1=" << np1 << endl;
//  cout << endl << "We have got to here in whiten.cpp with xaw=" << xaw << endl;
  float *sw;                // whitened beam power
  sw = new float [np1];
//  cout << endl << "We DONT get to here in whiten.cpp with nstop=" << nstop << endl;
  float sw1;                    // mean value of whitened beam
  int iswit = 1;                // iswit = 1 calc and apply whitening filter
                                // iswit = 0 apply prev calc whitening params
  float ya[m+1];                // auxilliary array
  int mr = 0;                   // moving window for phase detection (n/a here)

  int iopt = 2;                 // iopt = 2 -> do beam whitening
  int ier;



  // call grbeamfep_
  grbeamfep_(&nc, &np, &nc1, ncn, &ip, &np1,
	     &m, &mr,
	     sac->y, x2, s1, s, &s2,
	     x, &bm, &bs,
	     cr, bar, &disp, &ier,
	     xaw, &sw1, sw,
	     &iswit,
	     ya, af, &nstrt, &nstop, &iopt);
  
  
  cout << "Dispersion of Ar residuals: " << disp << endl;

  // assign whitened signal to SAC array
  delete [] sac->y;
  sac->y = xaw;
  
  return 0;

}


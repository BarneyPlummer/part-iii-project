// whiten.hpp

#include <cmath>
#include <iostream>
#include "SacFile.hpp"

using namespace std;

// prototype for grbeamfep Fortran function
extern "C" 
{
  void grbeamfep_ ( int *nc, int *np, int *nc1, int *ncn, int *ip, int *np1,
		    int *m, int *mr,
		    float *x1, float *x2, float *s1, float *s, float *s2,
		    float *x, float *bm, float *bs,
		    float *cr, float *bar, float *disp, int *ier,
		    float *xaw, float *sw1, float *sw,
		    int *iswit,
		    float *ya, float *af, int *n1, int *n2, int *iopt );
}

// prototype for whiten function
int whiten_signal(struct Sac *sac, int nstrt, int nstop);

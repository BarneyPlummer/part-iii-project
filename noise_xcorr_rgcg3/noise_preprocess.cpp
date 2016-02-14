// **********************************************************************
// noise_preprocess.cpp
//
// Program to do preprocessing prior to noise cross-correlation
// 
// Will do the following:
//  1. bandpass filter the data
//  2. remove mean and trend
//  3. do temporal normalisation (one-bit, clipping, and running mean)
//  4. whiten the signal
//
// Jamie Barron
// January 2010
// **********************************************************************

#include <cmath>
#include <fstream>
#include <iostream>

// include gsl least-squares fitting
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>

// include header for SAC file read/write routines
#include "SacFile.hpp"

// include header for whitening routine
#include "whiten.hpp"

// include boost library headers
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace std;

const int NORM_ONE_BIT = 1;
const int NORM_CLIPPING = 2;
const int NORM_RUN_MEAN = 3;
const int NORM_RUN_MEAN_FILT = 4;

// function prototypes
int clipping_norm(float *x, int npts);
int hanning_taper(float *y, int npts, int width);
int one_bit_norm(float *x, int npts);
int run_abs_mean_norm(float *x, int npts, int norm_win_len);
int run_abs_mean_norm(float *x, int npts, int norm_win_len,
		      float delta, int npoles, float f1, float f2);
int remove_mean(float *x, int npts);
int remove_trend(float *x, int npts);
float root_mean_square(float *x,int npts);

// prototype for Fortran iir function
extern "C" {
  void iir_(float *of, float *af, float *dt, int *npts, int *nl, float *fl, 
	    int *nh, float *fh);
}

// ----------------------------------------------------------------------
// function to do clipping normalisation of data
// -- normalize by rms of signal
// -- clips normalized ampl. to 1 where it exceeds the rms of the signal
// 
// arguments:
//   *x   - pointer to array containing data to normalise
//   npts - number of data points in array
//
// returns:
//   0 on successful completion
// ----------------------------------------------------------------------
int 
clipping_norm(float *x, int npts)
{
  
  // calculate rms
  float rms = root_mean_square(x,npts);

  // do the normalisation
  for (int i=0; i < npts; i++)
    {
      x[i] = x[i] / rms;
      if ( fabs(x[i]) > 1 )
	{
	  x[i] = x[i] / fabs(x[i]);
	}
    }

  return 0;

}


// ----------------------------------------------------------------------
// function to do one-bit normalisation of data
//
// arguments:
//   *x   - pointer to array containing data to normalise
//   npts - number of data points in array
//
// returns: 
//   0 on successful completion
// ----------------------------------------------------------------------
int
one_bit_norm(float *x, int npts)
{
  
  for (int i=0; i < npts; i++) 
    {
      x[i] = x[i] / fabs(x[i]);
    }

  return 0;

}


// ----------------------------------------------------------------------
// function to normalise the data using a running mean window
// if called in the first version below, no filtering is carried out
// if the second version is called, then the normalisation weights
// are derived from a seismogram bandpass filtered as specified
//
// arguments:
//   *x            - pointer to array containing data to normalise
//   npts          - number of data points in array
//   norm_win_len  - length of window to use for running mean
//
// additional arguments for filtered version:
//   delta         - time series sample interval
//   npoles        - number of poles to use for the filtering
//                   (uses same number for both high and low pass)
//   f1            - freq limit 1
//   f2            - freq limit 2 
// (n.b. frequencies will be swapped to appear in the right order
//  to call the filtering routine)
//
// returns:
//   0 on successful completion
// ----------------------------------------------------------------------
int 
run_abs_mean_norm(float *x, int npts, int norm_win_len)
{

  // calculate normalisation weights
  //
  // at either end of the data series we truncate the normalisation
  // window so that it's asymmetrical about the point for which
  // we're calculating the weight
  //
  float *w;
  w = new float [npts];
  
  for (int i=0; i < npts; i++ )
    {
      int norm_win_start = max(i-norm_win_len,0);
      int norm_win_end = min(i+norm_win_len,npts-1);

      w[i] = 0;
      for ( int j=norm_win_start; j <= norm_win_end; j++ ) 
	{
	  w[i] += fabs(x[j]);
	}
      w[i] = w[i] / (norm_win_end - norm_win_start + 1);
      
    }

  // normalise data using calculated weights
  for (int i=0; i < npts; i++ )
    {
      x[i] = x[i] / w[i];
    }

  // clean up memory
  delete [] w;

  return 0;

}


int run_abs_mean_norm(float *x, int npts, int norm_win_len,
		      float delta, int npoles, float f1, float f2)
{

  // filter the data at specified parameters
  float *xfilt;
  xfilt = new float [npts];
  
  if ( f1 > f2 )
    {
      iir_(x,xfilt,&delta,&npts,&npoles,&f1,&npoles,&f2);
    }
  else
    {
      iir_(x,xfilt,&delta,&npts,&npoles,&f2,&npoles,&f1);
    }

  // calculate normalisation weights based on filtered seismogram
  // 
  // at either end of the data series we truncate the normalisation
  // window so that it's asymmetrical about the point for which we're
  // calculating the weight
  float *w;
  w = new float [npts];
  
  for (int i=0; i < npts; i++ )
    {
      int norm_win_start = max(i-norm_win_len,0);
      int norm_win_end = min(i+norm_win_len,npts-1);
      
      w[i] = 0;
      for ( int j=norm_win_start; j <= norm_win_end; j++ ) 
	{
	  w[i] += fabs(xfilt[j]);
	}
      w[i] = w[i] / (norm_win_end - norm_win_start + 1);
      
    }

  // normalise data using calculated weights
  for (int i=0; i < npts; i++ )
    {
      x[i] = x[i] / w[i];
    }


  // clean up memory
  delete [] xfilt;
  delete [] w;

  return 0;

}


// ----------------------------------------------------------------------
// function to calculate root mean square of an array
//
// arguments:
//   *x   - pointer to array containing data
//   npts - number of data points in array
//
// returns: 
//   rms value of array
// ----------------------------------------------------------------------
float
root_mean_square(float *x, int npts)
{
  
  float rms = 0;
  
  for (int i=0; i < npts; i++ )
    {
      rms += x[i]*x[i];
    }

  rms = rms / npts;
  rms = sqrt(rms);
  return rms;

}


// ----------------------------------------------------------------------
// function to calculate and remove mean of an array
//
// arguments:
//   *x    - pointer to array containing data
//   npts  - number of data points in array
// 
// returns:
//   0 on successful completion
// ----------------------------------------------------------------------
int
remove_mean(float *x, int npts)
{
  float m = 0;

  cout << "Removing mean" << endl ;
  for (int i=0; i < npts; i++)
    m += x[i];
  m = m/npts;

  for (int i=0; i < npts; i++)
    x[i] -= m;

  return 0;
}

// ----------------------------------------------------------------------
// function to calculate and remove the trend of a seismogram
//
// arguments:
//   *y    - pointer to array containing data
//   npts  - number of data points in array
// 
// returns:
//   0 on successful completion
// ----------------------------------------------------------------------
int
remove_trend(float *y, int npts)
{

  double c0, c1, cov00, cov01, cov11, sumsq;
  
  // create x array for fit
  double *x;
  x = new double [npts];
  for (int i = 0; i < npts; i++)
    {
      x[i] = double(i);
    }

  // convert y array to double yy
  double *yy;
  yy = new double [npts];
  for (int i=0; i < npts; i++)
    {
      yy[i] = double(y[i]);
    }

  // do linear fit
  gsl_fit_linear(x,1,yy,1,npts,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
  delete [] x;
  delete [] yy;

  cout << "Removing trend - y = c0 + c1*x: " << endl;
  cout << "Intercept c0 = " << c0 << " std dev = " << sqrt(cov00) << endl;
  cout << "Gradient  c1 = " << c1 << " std dev = " << sqrt(cov11) << endl << endl;

  // remove trend from array y
  for ( int i=0; i < npts; i++ )
    {
      y[i] = y[i] - ( float(c0) + float(c1)*float(i) );
    }
  
  return 0;

}

// ----------------------------------------------------------------------
// function to apply Hanning taper to seismogram
//
// arguments:
//   *y    - pointer to array containing data
//   npts  - number of data points in array
//   width - width of Hanning taper at each end (in samples)
//   
// returns:
//   0 on successful completion
// ----------------------------------------------------------------------
int
hanning_taper(float *y, int npts, int width)
{
  // generate taper function
  float *taper;
  taper = new float [npts];
  
  for ( int i=0; i<npts; i++ )
    {
      if ( i < width )
	{
	  taper[i] = 0.5*(1-cos(M_PI*i/width));
	}
      else if ( i > (npts-width) )
	{
	  taper[i] = 0.5*(1-cos(M_PI*(npts-i-1)/width));
	}
      else
	{
	  taper[i] = 1;
	}
    }

  // multiply the data by the taper
  cout << "Applying Hanning taper of width " << width << " samples" << endl;
  
  for ( int i=0; i<npts; i++ )
    {
      y[i] = y[i]*taper[i];
    }

  delete [] taper;

  return 0;
}

// **********************************************************************
// MAIN PROGRAM
// **********************************************************************
int
main ( int argv, char* argc [] )
{

  // ** set up parsers for command line and configuration files
  string usage = "noise_preprocess [options] input_file.SAC";
  string cfgfile;
  string normtypestr;
  int taperlen;

  // command line options
  po::options_description cmdline("Command line options");
  cmdline.add_options()
    ("config-file,c",po::value<string>(&cfgfile)->default_value("noise.cfg"),
     "load specified configuration file")
    ("output-dir,o",po::value<string>(),
     "set output directory (default preproc_files)")
    ("help,h", "output help message")
    ("normalization,n",po::value<string>(&normtypestr)->default_value("runmean_filt"),
     "set time-domain normalisation type\n"
     "   clipping\n"
     "   onebit\n"
     "   runmean\n"
     "   runmean_filt (default)\n")
    ("no-band-pass","turn off initial band-pass filtering of record")
    ("no-remove-mean", "turn off removal of mean")
    ("no-remove-trend", "turn off removal of trend")
    ("no-taper","turn off tapering of signal prior to bandpass")
    ("no-whiten", "turn off whitening of signal")
    ("overwrite","overwrite existing files (instead of writing to separate "
     "directory")
    ("taper-window,t",po::value<int>(&taperlen)->default_value(100),
     "taper window length (number of samples each end of data)")
    ;

  // positional options (add to hidden options group so doesn't appear in 
  // options list)
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file",po::value< vector<string> >(),"input file")
    ;
  po:: positional_options_description p;
  p.add("input-file",-1);

  // cmdline_pos combines both command line and hidden positional options, to
  // use to parse the command line
  po::options_description cmdline_pos;
  cmdline_pos.add(cmdline).add(hidden);
  po::variables_map vm;

  // parse command line
  try {
    po::store(po::command_line_parser(argv,argc).options(cmdline_pos)
	      .positional(p).run(), vm);
    po::notify(vm);

    // validate input
    if ( vm.count("help") )
      {
	cerr << endl;
	cerr << usage << endl;
	cerr << endl;
	cerr << cmdline << endl;
	return 1;
      }

    if (! vm.count("input-file") )
      {
	cerr << "Must supply input SAC filename(s)" << endl;
	cerr << usage << endl;
	cerr << endl;
	cerr << cmdline << endl;
	return 1;
      }
    
    if ( (normtypestr !=  "onebit")  && 
	 (normtypestr != "clipping") && 
	 (normtypestr != "runmean")  && 
	 (normtypestr != "runmean_filt") )
      {
	cerr << "Unknown normalization type specified: " << normtypestr << endl;
	cerr << usage << endl;
	cerr << endl;
	cerr << cmdline << endl;
	return 1;
      }

    // ensure we haven't both specified an output directory and 
    // to overwrite existing files
    if ( vm.count("output-dir") && vm.count("overwrite") )
      {
	cerr << "Cannot specify output-dir and overwrite options "
	  "simultaneously" << endl;
	cerr << usage << endl;
	cerr << endl;
	cerr << cmdline << endl;
	return 1;
      }
    

  }
  catch(exception &e)
    {
      cerr << "Error: " << e.what() << endl;
      cerr << usage << endl;
      return 1;
    }

  // assign parsed values to variables (cfgfile already assigned in processing)
  vector<string> sacFilenames = vm["input-file"].as< vector<string> >();
  bool rmmean = vm.count("no-remove-mean") ? false : true;
  bool rmtrend = vm.count("no-remove-trend") ? false: true;
  bool bandpass = vm.count("no-band-pass") ? false : true;
  bool taper = vm.count("no-taper") ? false : true;
  bool whiten = vm.count("no-whiten") ? false : true;
  
  int normtype;
  if ( normtypestr == "onebit" ) 
    { normtype = NORM_ONE_BIT; }
  else if ( normtypestr == "clipping" ) 
    { normtype = NORM_CLIPPING;}
  else if ( normtypestr == "runmean" ) 
    { normtype = NORM_RUN_MEAN;}
  else if ( normtypestr == "runmean_filt" ) 
    { normtype = NORM_RUN_MEAN_FILT; }
  else {
    cerr << "Error: Unknown normalization type: " << normtypestr << endl;
    cerr << endl;
    return 1;
  }

  // setup output directory (if applicable) and overwrite flag
  string outdir;
  bool overwrite = false;
  if ( vm.count("output-dir") )
    {
      outdir = vm["output-dir"].as<string>();
    }
  else if ( vm.count("overwrite") )
    {
      overwrite = true;
    }
  else
    {
      outdir = "preproc_files";
    }
      
  // now set up config file parser and parse the configuration file
  po::options_description cfg_parser("Config file options");
  cfg_parser.add_options()
    ("prefilter.bp_low_corner",po::value<float>(),"low corner freq of bp filter")
    ("prefilter.bp_high_corner",po::value<float>(),"high corner freq of bp filter")
    ("prefilter.bp_npoles",po::value<int>(),"number of poles of bp filter")
    ("runmean.window_length",po::value<int>(),"window length for running mean norm.")
    ("runmean.bp_low_corner",po::value<float>(),"low corner freq of bp filt for norm.")
    ("runmean.bp_high_corner",po::value<float>(),"high corner freq of bp filt for norm.")
    ("runmean.bp_npoles",po::value<int>(),"number of poles of bp  filt for norm.")
    ;
  
  ifstream ifscfg;
  try
    {
      if ( bandpass || (normtype==NORM_RUN_MEAN) ||
	   normtype==NORM_RUN_MEAN_FILT)
	{
	  ifscfg.open(cfgfile.c_str());
	  if (! ifscfg.is_open())
	    {
	      cerr << "Unable to open configuration file: " << cfgfile << endl;
	      return 1;
	    }
	  po::store(po::parse_config_file(ifscfg,cfg_parser),vm);
	  po::notify(vm);
	  ifscfg.close();
	}

      // validate input
      if ( bandpass ) 
	{
	  if ( ! ( vm.count("prefilter.bp_low_corner") 
		   && vm.count("prefilter.bp_high_corner")
		   && vm.count("prefilter.bp_npoles") ))
	    {
	      cerr << "Config file " << cfgfile << " does not specify one or more of"
		   << endl;
	      cerr << "  bp_low_corner" << endl;
	      cerr << "  bp_high_corner" << endl;
	      cerr << "  bp_npoles" << endl;
	      cerr << "in the options group [prefilter]" << endl;
	      cerr << endl;
	      cerr << "These options are required unless --no-bandpass commondline"
		" option given" << endl;
	      return 1;
	    }
	}

      if ( (normtype==NORM_RUN_MEAN) || (normtype==NORM_RUN_MEAN_FILT) )
	{
	  if ( ! vm.count("runmean.window_length") )
	    {
	      cerr << "Config file " << cfgfile << " does not specify the window length "
		"for running mean normalisation: " << endl ;
	      cerr << "section [runmean], variable window_length" << endl;
	      cerr << endl;
	      return 1;
	    }
	}

      if ( normtype == NORM_RUN_MEAN_FILT )
	{
	  if ( ! ( vm.count("runmean.bp_low_corner")
		   && vm.count("runmean.bp_high_corner")
		   && vm.count("runmean.bp_npoles") ))
	    {
	      cerr << "Config file " << cfgfile << " does not specify one or more of"
		   << endl;
	      cerr << "  bp_low_corner" << endl;
	      cerr << "  bp_high_corner" << endl;
	      cerr << "  bp_npoles" << endl;
	      cerr << "in the options group [runmean]" << endl;
	      cerr << endl;
	      cerr << "These options are required for runmean_filt normalization"<< endl;
	      return 1;
	    }
	}

    }
  catch(exception &e)
    {
      cerr << "Error: " << e.what() << endl;
      cerr << "Problem parsing configuration file " << cfgfile << endl;
      cerr << endl;
      if (ifscfg.is_open()) ifscfg.close();
      return 1;
    }

  // assign variables parsed from config file
  // --pre-filter parameters
  int pre_filt_npoles;
  float pre_filt_f1;
  float pre_filt_f2;
  if (bandpass)
    {
      pre_filt_npoles = vm["prefilter.bp_npoles"].as<int>();
      pre_filt_f2 = vm["prefilter.bp_low_corner"].as<float>();
      pre_filt_f1 = vm["prefilter.bp_high_corner"].as<float>();
    }

  // --running mean normalization parameters
  int win_len;
  float norm_filt_npoles;
  float norm_filt_f1;
  float norm_filt_f2;
  
  if ( ( normtype == NORM_RUN_MEAN ) || ( normtype == NORM_RUN_MEAN_FILT ) )
    {
      win_len = vm["runmean.window_length"].as<int>();
    }

  if ( normtype == NORM_RUN_MEAN_FILT )
    {
      norm_filt_npoles = vm["runmean.bp_npoles"].as<int>();
      norm_filt_f1 = vm["runmean.bp_low_corner"].as<float>();
      norm_filt_f2 = vm["runmean.bp_high_corner"].as<float>();

    }
  
  // iterate over input sac files
  for ( vector<string>::iterator it = sacFilenames.begin(); 
	it != sacFilenames.end();
	++it)

    {
      // read in the SAC file
      string sacFilename = *it;
      struct Sac *sac_in;
      cout << endl << "****************************************"
	"**************************************" << endl;
      cout << "Reading file: " << sacFilename << endl;
      sac_in = read_sacfile(sacFilename);
      if ( sac_in == NULL )
	{
	  return -1;
	}

      // remove mean and trend prior to filtering
      if ( rmmean ) remove_mean(sac_in->y,sac_in->header.npts);
      if ( rmtrend ) remove_trend(sac_in->y, sac_in->header.npts);

      // do band-pass filtering of record
      if ( bandpass ) 
	{
	  if (taper) hanning_taper(sac_in->y,sac_in->header.npts,taperlen);

	  float *yfilt;
	  yfilt = new float [sac_in->header.npts];
  
	  // ensure filter arguments are in the right order for the iir_ subroutine
	  if ( pre_filt_f1 < pre_filt_f2 )
	    {
	      float temp = pre_filt_f1;
	      pre_filt_f1 = pre_filt_f2;
	      pre_filt_f2 = temp;
	    }


	  printf("Filtering data between %f and %f Hz\n",pre_filt_f2,pre_filt_f1);
	  iir_(sac_in->y,yfilt,
	       &sac_in->header.delta,&sac_in->header.npts,
	       &pre_filt_npoles, &pre_filt_f1, 
	       &pre_filt_npoles, &pre_filt_f2);
   
	  // assign filtered record to array
	  delete [] sac_in->y;
	  sac_in->y = yfilt;

	  // remove mean and trend of filtered signal
	  if ( rmmean ) remove_mean(sac_in->y,sac_in->header.npts);
	  if ( rmtrend ) remove_trend(sac_in->y, sac_in->header.npts);
	}


      // do normalisation
      switch (normtype)
	{

	case NORM_ONE_BIT:
	  cout << "Applying one-bit normalization to data" << endl;
	  one_bit_norm(sac_in->y,sac_in->header.npts);
	  break;
      
	case NORM_CLIPPING:
	  cout << "Applying normalization: clipping data to rms value" << endl;
	  clipping_norm(sac_in->y,sac_in->header.npts);
	  break;

	case NORM_RUN_MEAN:
	  cout << "Applying running mean normalization: window length "<<win_len<<endl;
	  run_abs_mean_norm(sac_in->y,sac_in->header.npts,win_len);
	  break;

	case NORM_RUN_MEAN_FILT:
	  cout << "Applying running mean normalization with filter" << endl;
	  cout << "   window length = " << win_len << endl;
	  cout << "   normalizing by record filtered between " << norm_filt_f1 << 
	    " and " << norm_filt_f2 << " (" << norm_filt_npoles << " poles)." << endl;
	  run_abs_mean_norm(sac_in->y,sac_in->header.npts,win_len,
			    sac_in->header.delta,norm_filt_npoles,
			    norm_filt_f1,norm_filt_f2);
	  break;
      
	default:
	  break;

	}

      // do spectral whitening
      if (whiten) 
	{
	  cout << endl << "Applying spectral whitening to signal" << endl;
	  cout << endl << "Npts is " << sac_in->header.npts << endl;
	  whiten_signal(sac_in,1,sac_in->header.npts);
	}


  
      // write out SAC file
      string outfilename;
      if ( overwrite )
	{
	  // set output filename to input filename, i.e. overwrite
	  outfilename = sacFilename;
	}
      else
	{
	  // create output directory if it doesn't exist
	  // (and check no file exists with the same name!)
	  fs::path outpath(outdir);
	  if ( ! fs::is_directory(outpath) ){
	    if ( fs::exists(outpath) ) 
	      {
		cerr << "File with name of output directory already exists" << endl;
		return 1;
	      }
	    printf("Making directory: %s\n",outdir.c_str());
	    fs::create_directory(outpath);
	  }
      
	  // make output filename
	  fs::path inpath(sacFilename);
	  fs::path outfile = operator /(outpath,inpath.filename());
	  outfilename = outfile.string();
	}

      cout << "Writing: " << outfilename << endl;

      int stat = write_sacfile(outfilename,sac_in);
      if ( stat )
	{
	  fprintf(stderr,"Error writing file: %s", outfilename.c_str());
	  return -1;
	}

      // clean up memory
      delete [] sac_in->y;
      delete sac_in;

    }

  cout << endl << "*****************************************************"
    "*************************" << endl;

  return 0;

}

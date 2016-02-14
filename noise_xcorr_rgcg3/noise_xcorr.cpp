// **********************************************************************
// noise_xcorr.cpp
//
// Program to do cross-correlation between signals
//
// Author: Jamie Barron
// Date: February 2011
//
// Modifications
// Author: Rob Green
// Modifications to run cross-correlation and autocorrelation (Lines 177-208)
//
// Modifications to output shorter Noise Correlation Functions with an 
// odd npts for functionality in AFTAN. (Lines 359-386)
//
// Necessary parsing modifications are made in lines 793-820
//
// **********************************************************************

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <vector>

// include boost headers
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

// include other headers
#include "SacFile.hpp"
#include "string_utils.hpp"
#include "time_utils.hpp"
#include "constants.hpp"
#include "fourier.hpp"

using namespace std;
namespace fs = boost::filesystem;
namespace po = boost::program_options;

// ---------------------------------------------------------------------- //
// StationPair structure holds matched paths for station pairs
// ---------------------------------------------------------------------- //
struct StationPair {
  string stn1;
  string stn2;
  vector<fs::path> stn1_files;
  vector<fs::path> stn2_files;
};


// ====================================================================== //
// function prototypes
// ====================================================================== //
bool comp_filenames(fs::path path1, fs::path path2);
bool comp_file_times_day(fs::path path1, fs::path path2);
bool comp_file_times_hour(fs::path path1, fs::path path2);
vector<StationPair> construct_station_pairs(string filename, int tslen, int cortyp, int jset);
float *pad_zeros_front(float *inarray, int npts, int npad);
int correlate_station_pair(StationPair sp, int tslen,float delta);
int output_pair_files(vector<StationPair> station_pairs);
int output_correlation_files(string filenamebase,float *ccorr, 
			     SacHeader ccorr_header);


// ====================================================================== //
// FUNCTIONS
// ====================================================================== //

// ---------------------------------------------------------------------- //
// function to compare filenames - used to sort vector of paths by 
// filename
// ---------------------------------------------------------------------- //
bool 
comp_filenames(fs::path path1, fs::path path2)
{
  string filename1 = path1.filename().string();
  string filename2 = path2.filename().string();

  return (filename1 < filename2);
}


// ---------------------------------------------------------------------- //
// function to compare file times based on year/day -- used to match
// up files recorded at the same time 
// ---------------------------------------------------------------------- //
bool
comp_file_times_day(fs::path path1, fs::path path2)
{
  string filename1 = path1.filename().string();
  string filename2 = path2.filename().string();
  
  // compare by day
  string time1 = filename1.substr(0,8);
  string time2 = filename2.substr(0,8);
  
  return (time1 < time2);
}

// ---------------------------------------------------------------------- //
// function to compare file times based on year/day/hr -- used to match
// up files recorded at the same time 
// ---------------------------------------------------------------------- //
bool
comp_file_times_hour(fs::path path1, fs::path path2)
{
  string filename1 = path1.filename().string();
  string filename2 = path2.filename().string();
  
  // compare by hour
  string time1 = filename1.substr(0,11);
  string time2 = filename2.substr(0,11);
  
  return (time1 < time2);
}

// ---------------------------------------------------------------------- //
// function to pair stations and match up files recorded at same time
// returns a vector of StationPair structures... this will be returned
// empty in the case of a problem within the function
// ---------------------------------------------------------------------- //
vector<StationPair> 
construct_station_pairs(string filename, int tslen, int cortyp, int jset)
{
  vector<StationPair> station_pairs;

  // read in details of files to match, store in map object indexed by station
  map<string,vector<fs::path> > stationfiles;
  ifstream infile;
  infile.open(filename.c_str());
  if ( ! infile.is_open() )
    {
      cerr << "Error opening input file list: " << filename << endl;
      return station_pairs;
    }

  string filestring;
  getline(infile,filestring);
  while (!infile.eof())
    {
      // check file exists
      fs::path filepath(filestring);
      if ( ! fs::exists(filepath) )
	{
	  cerr << "File: " << filestring << "does not exist" << endl;
	  return station_pairs;
	}
      // get station from filename
      string filename = filepath.filename().string();
      vector<string> split = tokenize(filename,".");

      // store in stationfiles and stationfilenames maps
      stationfiles[split[7]].push_back(filepath);
      getline(infile,filestring);
    }
  infile.close();

  // sort vectors in map
  for ( map<string,vector<fs::path> >::iterator m = stationfiles.begin();
	m != stationfiles.end(); ++m )
    {
      sort(m->second.begin(),m->second.end(),comp_filenames);
    }

  // get list of unique stations from map
  vector<string> stationlist;
  for ( map<string,vector<fs::path> >::iterator station = stationfiles.begin();
	station != stationfiles.end(); ++station)
    {
      stationlist.push_back(station->first);
    }

  // Set maximum i for station loop (RG)
  unsigned int maxi=stationlist.size();
  
  // loop through and create vector of StationPair structs for each station pair
  for ( unsigned int i=0; i < maxi; i++ )
    {

      // Set maximum j for station loop (RG)
      unsigned int maxj=stationlist.size();
      if ( jset == 1)
      {
      // For xcorr and bothcorr
        maxj=stationlist.size();
      }
      else if ( jset == 0 )
      {
      // For autocorr
        maxj=i+1;
      }
      else
      {
      cerr << "Error: unknown value for jset in construct_station_pairs" << endl;
      station_pairs.clear();
      return station_pairs;
      }

      // Autocorr/Xcorr sets the chosen station pairings (RG)
      // cortype sets where j loop starts from (RG)
      for ( unsigned int j=i+cortyp; j < maxj; j++ )
	{
	  StationPair station_pair;

	  station_pair.stn1 = stationlist[i];
	  station_pair.stn2 = stationlist[j];

	  // note -- sort above means these lists are sorted by filename
	  vector<fs::path> st1paths = stationfiles[stationlist[i]];
	  vector<fs::path> st2paths = stationfiles[stationlist[j]];

	  // iterate over one list of files and use binary search to 
	  // match against file times in second list
	  pair<vector<fs::path>::iterator,vector<fs::path>::iterator> bounds;
	  for ( vector<fs::path>:: iterator it = st1paths.begin();
		it != st1paths.end(); ++it )
	    {
	      if ( tslen == DAY_BLOCKS )
		{
		  bounds = equal_range( st2paths.begin(), st2paths.end(), 
					*it, comp_file_times_day);
		}
	      else if ( tslen == HOUR_BLOCKS )
		{
		  bounds = equal_range( st2paths.begin(), st2paths.end(), 
					*it, comp_file_times_hour);
		}
	      else
		{
		  cerr << "Error: unknown value for tslen in construct_station_pairs" << endl;
		  cerr << "Acceptable values are:" << endl;
		  cerr << "DAY_BLOCKS: " << DAY_BLOCKS;
		  cerr << "HOUR_BLOCKS: " << HOUR_BLOCKS;
		  station_pairs.clear();
		  return station_pairs;
		}

	      // calculate number of matches
	      // warn if more than one match, but take the first match
	      int nmatch = int(bounds.second - st2paths.begin()) - 
		int(bounds.first - st2paths.begin());
	      if ( nmatch > 1 )
		{
		  cout << endl;
		  cout << "Warning: more than one match for file at path: " << endl;
		  cout << *it << endl;
		  cout << "Using first match: " << endl;
		  cout << *bounds.first;
		  cout << "Other matches are: " << endl;
		  for ( vector<fs::path>::iterator it2 = bounds.first+1;
			it2 < bounds.second; ++it2 )
		    {
		      cout << *it2 << endl;
		    }
		  cout << endl;
		}
	      // otherwise, if we have a match, then add to pair.
	      else if ( nmatch == 1 )
		{
		  station_pair.stn1_files.push_back(*it);
		  station_pair.stn2_files.push_back(*bounds.first);
		}

	    }
	  // if pair is non-empty, then assign to station_pairs vector
	  if ( station_pair.stn1_files.size() != 0 ) 
	    station_pairs.push_back(station_pair);
	}
    }
  
  // return the station_pairs vector
  return station_pairs;
}


// ---------------------------------------------------------------------- //
// function to output files containing details of the paired up files 
// for each station pair
// ---------------------------------------------------------------------- //
int
output_pair_files(vector<StationPair> station_pairs)
{
  for ( vector<StationPair>::iterator sp = station_pairs.begin();
	sp != station_pairs.end(); ++sp)
    {
      // generate file name based on station names
      string spfile = "files_" + sp->stn1 + "_" + sp->stn2;
      ofstream outfile;
      outfile.open(spfile.c_str());
      if ( ! outfile.is_open() ) 
	{
	  cerr << "Unable to open file for output: " << spfile << endl;
	  return -1;
	}

      cout << "Writing pairs to file: " << spfile << endl;
      for ( unsigned int i = 0; i < sp->stn1_files.size(); i++ )
	{
	  string file1 = sp->stn1_files[i].string();
	  string file2 = sp->stn2_files[i].string();
	  outfile << file1 << " " << file2 << endl;
	}

      outfile.close();
    }
  
  return 0;

}


// ---------------------------------------------------------------------- //
// function to output SAC files with the correlations in
// outputs 2 files - one with +ve and -ve lags, and another with the 
// +ve and -ve lags averaged to get a one-sided correlation function 
// ---------------------------------------------------------------------- //
int
output_correlation_files(string filenamebase,float *ccorr, 
			 SacHeader ccorr_header, int long_out)
{

  // split the ccorr array up into +ve and -ve lag sections
  // (the -ve lag being flipped, so the arrivals closest to zero
  // are first in the array).   Also join the complete signal
  // (+ve and -ve lags) in the correct order.
  float *ccorr_pos;
  float *ccorr_neg_rev;
  float *ccorr_shifted;
  float *ccorr_short;
  float *ccorr_short_pos;

  int n_half = ccorr_header.npts / 2;
  int n_pos = n_half + 1;
  int n_neg;
  if ( (ccorr_header.npts % 2)  == 0 )
    n_neg = n_half - 1;
  else
    n_neg = n_half;

  ccorr_pos = new float [n_pos];
  ccorr_neg_rev = new float [n_neg + 1];
  //ccorr_neg_rev = new float [n_pos];
  ccorr_shifted = new float [ccorr_header.npts];

  for (int i=0; i < ccorr_header.npts; i++)
    {
      if ( i < n_pos )
	{
	  ccorr_pos[i] = ccorr[i];
	  ccorr_shifted[n_neg + i] = ccorr[i];
	}
      else
	{
	  ccorr_neg_rev[i-n_pos+1] = ccorr[ccorr_header.npts-i+n_pos-1];
	  ccorr_shifted[i-n_pos] = ccorr[i];
	}
    }
  ccorr_neg_rev[0] = ccorr[0];


  // Initiate SAC file header
  Sac sacf;
  sacf.header = ccorr_header;
  // write long SAC file for xcorr with +ve and -ve lags (depends on usage-RG)
  if ( long_out < 1 )
    {
    sacf.header.b = -(n_neg)*sacf.header.delta;
    sacf.header.e = (n_pos-1)*sacf.header.delta;
    sacf.y = ccorr_shifted;
    write_sacfile(filenamebase+"_full_long.SAC",&sacf);
    cout << "Written long xcorr with 131072 npts: " << endl;
    }

  // Make a shorter length xcorr to output (RG)
  int out_length = 14401; // MUST BE ODD! / Length of short fullNCF (14401)
  int out_half_length = out_length / 2; // Float rounded down to an int (7200)
  int out_neg = out_half_length;
  int out_length_avg = out_half_length + 1; // Length of short avgNCF (7201)
  ccorr_short = new float [out_length];
  for ( int i=0; i < out_length; i++ )
    {
      ccorr_short[i] = ccorr_shifted[i+n_neg-out_neg];
    }

  // write shorter SAC file with +ve and -ve lags (RG)
  sacf.header.npts = out_length;
  sacf.header.b = -(out_neg)*sacf.header.delta;
  sacf.header.e = out_neg*sacf.header.delta;
  sacf.y= ccorr_short;
  write_sacfile(filenamebase+"_full.SAC",&sacf);
  cout << "Written xcorr with npts: " << out_length << endl;


  // calculate average of +ve and -ve lags
  for ( int i=0; i < n_pos; i++ )
    {
      ccorr_pos[i] = (ccorr_pos[i] + ccorr_neg_rev[i])/2;
    }

  // Make a shorter length avg xcorr to output (RG)
  ccorr_short_pos = new float [out_length_avg];
  for ( int i=0; i < out_length_avg; i++ )
    {
      ccorr_short_pos[i] = ccorr_pos[i];
    }

  // write short SAC file for average of +ve and -ve lags (RG)
  sacf.header.b = 0;
  sacf.header.e = (out_length_avg-1)*sacf.header.delta;
  sacf.header.npts = (out_length_avg);
  sacf.y = ccorr_short_pos;
  write_sacfile(filenamebase+"_avg.SAC",&sacf);


  // write long SAC file for average of +ve and -ve lags (depends on usage-RG)
  if ( long_out < 1 )
    {
    sacf.header.b = 0;
    sacf.header.e = (n_pos-1)*sacf.header.delta;
    sacf.header.npts = n_pos;
    sacf.y = ccorr_pos;
    write_sacfile(filenamebase+"_avg_long.SAC",&sacf);
    }

  // clean up memory
  delete [] ccorr_pos;
  delete [] ccorr_neg_rev;
  delete [] ccorr_shifted;
  
  return 0;

}

float *
pad_zeros_front(float *inarray, int32_t *npts, int npad)
{
  
  int32_t newnpts;
  float *newarray;
  
  if ( npad == 0 )
    return inarray;

  newnpts = *npts + npad;
  newarray = new float [newnpts];
  for ( int i = 0; i < npad; i++ )
    newarray[i] = 0.;
  for ( int i = npad; i < newnpts; i++ )
    newarray[i] = inarray[i-npad];
  
  inarray = newarray;
  *npts = newnpts;
  
  return newarray;

}

// ---------------------------------------------------------------------- //
// function to do the cross correlation for a particular station pair
// (cross correlates and stacks all recordings for that pair, and 
// outputs the cross correlation files 
// ---------------------------------------------------------------------- //
int 
correlate_station_pair(StationPair sp, int tslen, float delta, int longoutfile)
{
  // calculate the length of Fourier transform to use
  int overlap_sec = 3600; 
  int fft_min_npts;
  if ( tslen == HOUR_BLOCKS ) 
    {
      fft_min_npts = 3600/delta + overlap_sec/delta + 1;
    }
  else if ( tslen == DAY_BLOCKS )
    {
      fft_min_npts = 86400/delta + overlap_sec/delta + 1;
    }
  else
    {
      cerr << "Unknown value for tslen in correlate_station_pair" << endl;
      return 1;
    }

  // make length of FT the next-largest power of 2
  // (FFTW doesn't require this but it seems to be faster like that)
  int fft_power_2 = ceil(log2((double)fft_min_npts));
//  cout << "fft_min_npts is: " << endl;
//  cout << fft_min_npts << endl;
//  cout << "fft_power_2 is: " << endl;
//  cout << fft_power_2 << endl;
  int fft_npts = pow(2,fft_power_2);
//  cout << endl << "fft_npts is: " << fft_npts << endl;

  // initialise FFTW class
  FFTW fft(fft_npts,FFTW_ESTIMATE);

  // ccorr_header (header for Xcorr output) built from first pair of SAC files
  SacHeader ccorr_header = sac_null;
  
  // set up stack
  float *ccorr_stack;
  ccorr_stack = new float[fft_npts];
  for ( int i = 0; i < fft_npts; i++ )
    ccorr_stack[i] = 0;
  
  // loop over pairs of files -- read in and do the cross-correlation
  // between the timeseries
  cout << endl << "Correlating station pair: " << sp.stn1 << " " << sp.stn2 << endl;
  for ( unsigned int i = 0; i < sp.stn1_files.size(); i++ )
    {
      string file1 = sp.stn1_files[i].string();
      string file2 = sp.stn2_files[i].string();

      struct Sac *sac1;
      struct Sac *sac2;

      sac1 = read_sacfile(file1);
      if ( sac1 == NULL )
	{
	  cerr << "Error reading SAC file: " << file1 << endl;
	  return 2;
	}
      sac2 = read_sacfile(file2);
      if ( sac2 == NULL )
	{
	  cerr << "Error reading SAC file: " << file2 << endl;
	  return 2;
	}

      // data must match specified sample rate (so that we have a 
      // consistent sample rate across the entire dataset).
      float eps = 1e-06;
      if ( ( abs(sac1->header.delta - delta) > eps ) ||
	   ( abs(sac2->header.delta - delta) > eps ) )
	{
	  cerr << "Data not at specified sample rate of " << delta 
	       << " samples per second" << endl;
	  return 3;
	}
	
      // set up cross-correlogram header on first correlation
      if ( i == 0 ) 
	{
	  ccorr_header.npts = fft_npts;
	  ccorr_header.delta = sac1->header.delta;
	  ccorr_header.iftype = ITIME;
	  ccorr_header.leven = TRUE;
	  ccorr_header.idep = IUNKN;
	  ccorr_header.lovrok = TRUE;
	  ccorr_header.lcalda = TRUE;
	  
	  ccorr_header.nzyear = 2000;
	  ccorr_header.nzjday = 1;
	  ccorr_header.nzhour = 0;
	  ccorr_header.nzmin = 0;
	  ccorr_header.nzsec = 0;
	  ccorr_header.nzmsec = 0;

	  // n.b. have to use loop rather than strcpy to copy over
	  // station names because strings in header are not null-terminated
	  for (int j=0; j<8; j++)
	    {
	      ccorr_header.kstnm[j] = sac1->header.kstnm[j];
	      ccorr_header.kevnm[j] = sac2->header.kstnm[j];
	    }
		
	  ccorr_header.stla = sac1->header.stla;
	  ccorr_header.stlo = sac1->header.stlo;
	  ccorr_header.stel = sac1->header.stel;
	  ccorr_header.stdp = sac1->header.stdp;
	  ccorr_header.evla = sac2->header.stla;
	  ccorr_header.evlo = sac2->header.stlo;
	  ccorr_header.evel = sac2->header.stel;
	  ccorr_header.evdp = sac2->header.stdp;

	}

      // do correlation
      cout << "Correlating: " << endl;
      cout << file1 << endl;
      cout << file2 << endl;
      // Write xcorr length to screen
      //cout << endl << "fft_npts is: " << fft_npts << endl;

      // zero pad time series that are not complete hour/day blocks.
      // (only need to worry about the start -- end is done automatically
      // when Fourier transforming).
      time_t sac1_start_timet, sac2_start_timet;
      time_t sac1_block_timet, sac2_block_timet;
      sac1_start_timet = get_time_t(sac1->header.nzyear, sac1->header.nzjday,
				    sac1->header.nzhour, sac1->header.nzmin,
				    sac1->header.nzsec);
      sac2_start_timet = get_time_t(sac2->header.nzyear, sac2->header.nzjday,
				    sac2->header.nzhour, sac2->header.nzmin,
				    sac2->header.nzsec);
				    
      if ( tslen == DAY_BLOCKS )
	{
	  sac1_block_timet = get_time_t(sac1->header.nzyear,
					sac1->header.nzjday,0,0,0);
	  sac2_block_timet = get_time_t(sac2->header.nzyear,
					sac2->header.nzjday,0,0,0);
	}
      else if ( tslen == HOUR_BLOCKS )
	{
	  sac1_block_timet = get_time_t(sac1->header.nzyear,sac1->header.nzjday,
					sac1->header.nzhour,0,0);
	  sac2_block_timet = get_time_t(sac2->header.nzyear,sac2->header.nzjday,
					sac2->header.nzhour,0,0);
	}
      else
	{
	  // n.b. shouldn't get here -- put this in in case an additional
	  // block length is added at the top of the program but not here
	  cerr << "Inconsistency in program noise_xcorr.cpp" << endl;
	  cerr << "Block length code " << tslen << " not completely supported"
	       << endl;
	  return 1;
	}
		
      // add the padding (note that this makes the SAC headersccin_1_npts inconsistent,
      // but given that we're finished with them after this it doesn't matter!
      int sac1_npad = difftime(sac1_start_timet,sac1_block_timet)/delta;
      int sac2_npad = difftime(sac2_start_timet,sac2_block_timet)/delta;
      
      // check for negative padding -- this is an error.
      if ( (sac1_npad < 0) || (sac2_npad < 0) )
	{
	  cerr << "Error: SAC files not split up correctly..." << endl;
	  cerr << "You should probably split up your files using split_data" << endl;
	  cerr << "and start again..." << endl;
	  cerr << "If you've used split_data then there's a bug somewhere..." << endl;
	  return 4;
	}

      float *ccin_1, *ccin_2;
      int ccin_1_npts, ccin_2_npts;
      ccin_1 = pad_zeros_front(sac1->y,&sac1->header.npts,sac1_npad);
      ccin_2 = pad_zeros_front(sac2->y,&sac2->header.npts,sac2_npad);
      ccin_1_npts = sac1->header.npts + sac1_npad;
      ccin_2_npts = sac2->header.npts + sac2_npad;

      // check the data -- make sure there are no NaNs etc.
      bool data_1_good = true;
      bool data_2_good = true;
      for ( int j = 0; j < ccin_1_npts; j++ )
	{
	  if ( isinf(ccin_1[j]) || isnan(ccin_1[j]) )
	    {
	      data_1_good = false;
	    }
	  if ( isinf(ccin_2[j]) || isnan(ccin_2[j]))
	    {
	      data_2_good = false;
	    }
	}

      float *ccorr;      
      if ( data_1_good && data_2_good )
	{
	  // actually do the cross-correlation!
	  ccorr = fft.cross_correlate(ccin_1,ccin_1_npts,ccin_2,ccin_2_npts);
	  
	  // add to stack
	  for ( int j = 0; j < fft_npts; j++ )
	    ccorr_stack[j] += ccorr[j];

	  delete [] ccorr;
	}
      else
	{
	  if ( ! data_1_good ) 
	    cerr << "Bad file: " << file1 << endl;
	  if ( ! data_2_good )
	    cerr << "Bad file: " << file2 << endl;
	}
 
      // free the ccin_* structures if we've used padding
      // (in this case they're separate from the sac*->y arrays
      //  but otherwise they're just pointers to existing arrays)
      if ( sac1_npad != 0 )
	{
	  delete [] ccin_1;
	}
      if ( sac2_npad != 0 )
	{
	  delete [] ccin_2;
	}

      // free memory
      delete [] sac1->y;
      delete sac1;
      delete [] sac2->y;
      delete sac2;


    }

  // output stack to file
  string outfilebase = "ccorr_" + sp.stn1 + "_" + sp.stn2;
  output_correlation_files(outfilebase,ccorr_stack,ccorr_header,longoutfile);
  
  delete [] ccorr_stack;

  return 0;
}

int
main(int argc, char *argv[])
{

  // ************************* parser ************************* //
  string usage = "noise_xcorr [options] input_file_list";

  po::options_description cmdline("Command line options");
  cmdline.add_options()
    ("help,h", "output this help message")
    ("dry-run","do the file pairing but don't do the correlation")
    ("day,d","correlate day blocks (default - use with day option in split_data)")
    ("hour,u","correlate hour blocks (use with hour option in split_data)")
    ("xcorr,x","run cross-correlations between station pairs (default)")
    ("autocorr,a","run auto-correlations on stations")
    ("bothcorr,b","run both auto and cross-correlations")
    ("sample-rate,s",po::value<float>(),"specify sample rate of input data "
     "in samples/sec (default 1)")
    ("longoutput,l","output longer records as well (default is only short)")
    ;

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file",po::value<string>(),"input file")
    ;

  po::positional_options_description p;
  p.add("input-file",1);

  po::options_description cmdline_pos;
  cmdline_pos.add(cmdline).add(hidden);
  po::variables_map vm;

  // run the parser
  try
    {
      po::store(po::command_line_parser(argc,argv).options(cmdline_pos)
		.positional(p).run(), vm);
      po::notify(vm);

      if ( vm.count("help") )
	{
	  cerr << endl;
	  cerr << usage << endl;
	  cerr << endl;
	  cerr << cmdline << endl;
	  return 1;
	}

      if ( ! vm.count("input-file") )
	{
	  cerr << "Must supply input file list" << endl;
	  cerr << usage << endl;
	  cerr << endl;
	  cerr << cmdline << endl;
	  return 1;
	}
	
      if ( vm.count("day") && vm.count("hour") )
	{
	  cerr << "Can only specify one block length for cross-correlation" << endl;
	  cerr << usage << endl;
	  cerr << endl;
	  cerr << cmdline << endl;
	  return 1;
	}

    }
  
  catch (exception& e)
    {
      cerr << "Error: " << e.what() << endl;
      cerr << usage << endl;
      cerr << endl;
      cerr << cmdline << endl;
      return 1;
    }

  // assign command line arguments to variables
  int ts_len = DAY_BLOCKS;
  if ( vm.count("day") )
    ts_len = DAY_BLOCKS;
  else if ( vm.count("hour") )
    ts_len = HOUR_BLOCKS;

  // assign variables for differing auto/xcorr pair selection (RG)
  int cor_typ = 1;
  int j_set = 1;
  if ( vm.count("xcorr") )
  {
    cor_typ = 1;
    j_set = 1;
  }
  else if ( vm.count("autocorr") )
  {
    cor_typ = 0;
    j_set = 0;
  }
  else if ( vm.count("bothcorr") )
  {
    cor_typ = 0;
    j_set = 1;
  }

  // assign variable for long file output (RG)
  // the variable changes from (1)longout > (2)longoutfile > (3)long_out
  // in the   (1)parser to (2)correlatepairfunction to (3)outputfunction
  int longout = 1;
  if ( vm.count("longoutput") )
  {
    longout = 0;
  }

  string filelist = vm["input-file"].as<string>();

  float delta = 1.0;
  if ( vm.count("sample-rate") )
    {
      float sample_rate = vm["sample-rate"].as<float>();
      if ( sample_rate != 0 )
	{
	  delta = 1 / sample_rate;
	}
      else
	{
	  cerr << "Sample rate cannot be zero!" << endl;
	  return -1;
	}
    }

  bool dryrun = vm.count("dry-run") ? true : false;

  // ********************* end of parser ********************* //


  // read in details of data from file, create station pairs
  vector<StationPair> station_pairs;
  station_pairs = construct_station_pairs(filelist,ts_len,cor_typ,j_set);

  // if vector is empty then there was a problem, so return.
  if ( station_pairs.size() == 0 )
    return 1;
  
  // output files required for each pairing (so as to be able to collect
  // these stations so cluster jobs can be set up with just one pair)
  output_pair_files(station_pairs);

  // if we're doing a dry-run, return here
  if ( dryrun )
    return 0;

  // for each pair of stations, cross-correlate pairs of time series
  // and stack (correlate_station_pair outputs SAC files containing
  // the correlations).
  for ( vector<StationPair>::iterator sp = station_pairs.begin();
	sp != station_pairs.end(); ++sp )
    {
      int stat = correlate_station_pair(*sp,ts_len,delta,longout);

      if ( stat != 0 )
	return 1;
    }

  return 0;

}

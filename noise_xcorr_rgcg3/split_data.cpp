// **********************************************************************
// split_data.cpp 
//
// Program to split up long records into day/hour blocks
// 
// Author: Jamie Barron
// Date: January 2011
// **********************************************************************

// include standard system headers
#include <ctime>
#include <iostream>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

// include boost library headers
// -- filesystem utilities & program option parser
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

// include other headers
#include "SacFile.hpp"
#include "string_utils.hpp"
#include "time_utils.hpp"
#include "constants.hpp"

namespace po = boost::program_options;
using namespace std;

struct Block {
  int start_sample;
  int end_sample;
  time_t start_time;
};  
  
// ---------------------------------------------------------------------- //
// function prototypes
// ---------------------------------------------------------------------- //
vector<Block> split_into_blocks(time_t begin_time_t,float delta,
				int nsamples, int mode);
int write_out_blocks(string infilename,string outdir,vector<Block> blocks, 
		     int msec, struct Sac *sac_in);

// ---------------------------------------------------------------------- //
// function to split the data into blocks 
// 
// arguments:
//   time_t begin_time_t -- begin time for the seismogram
//          (ignoring millisecs -- need to keep track of these separately)
//   float delta -- sample rate
//   int nsamples -- number of samples in complete record
//   int mode -- sets whether to split into hour or day blocks
//
// Note that if the ctime functions take account of leapseconds, etc.,
// this function should too...
// ---------------------------------------------------------------------- //
vector<Block> split_into_blocks(time_t begin_time_t,float delta,
				int nsamples, int mode)
{

  float eps = 1E-06;
  float frac = 0.8;

  printf("delta: %f\n", delta);
  printf("npts: %i\n", nsamples);

  struct Block block;
  vector<Block> blocks;

  // set the number of seconds per block (hour/day etc.)
  int sec_per_block;
  if ( mode == HOUR_BLOCKS) 
    { 
      sec_per_block = 3600; 
    }
  else if ( mode == DAY_BLOCKS ) 
    { 
      sec_per_block = 3600 * 24 ;
    }
  else 
    {
    printf("Unrecognized mode for function split_into_blocks: %d",mode);
    return blocks;
    }

  // check that we haven't got some strange sample rate that will
  // cause problems when breaking into consecutive hour/day blocks
  float fsamples_per_block = sec_per_block / delta ;
  int nsamples_per_block = sec_per_block / delta ;

    printf("Mode is %d\n",mode);
    printf("nsamples per block= %d\n",nsamples_per_block);

  if ( abs(fsamples_per_block - nsamples_per_block) > eps ) {
    printf("Error: can't split into blocks with sample rate %f\n",delta);
    return blocks;
  }
  
  // get end time for the potential first block, which is the first day/hour
  // after the start of the data.  This block will only be used if it
  // is of length > frac*sec_per_block
  struct tm end_time = *gmtime(&begin_time_t);
  end_time.tm_sec = 0;
  end_time.tm_min = 0;
  if ( mode == HOUR_BLOCKS )
    {
      end_time.tm_hour++;
    }
  else if ( mode == DAY_BLOCKS )
    {
      end_time.tm_hour = 0;
      end_time.tm_mday++;
    }
  end_time.tm_isdst = 0;
  time_t end_time_t = mktime(&end_time);
  int time_diff = difftime(end_time_t,begin_time_t);

  // break up into blocks and add details to vector<Block> blocks
  
  // first check that the data doesn't end before the end of the first
  // block -- if so, we return here with an empty vector blocks
  if ( time_diff > int(nsamples*delta) )
    {
      return blocks;
    }

  if ( time_diff > frac*sec_per_block ) 
   // include first part of hour if it is greater than fraction
    // frac of an hour # RG mod
    {
      block.start_sample = 1;
      block.end_sample = block.start_sample + time_diff/delta;
      block.start_time = begin_time_t;
      blocks.push_back(block);
    }
  else
    {
      // set up block.end_sample so that next part works correctly # RG mod
      block.end_sample = 1 + time_diff/delta;
    }

  // set up details for first block # RG mod
  block.start_sample = block.end_sample +1;  
  block.start_time = end_time_t;
  if ( mode == HOUR_BLOCKS ) 
    {
      end_time.tm_hour++;
    }
  else if ( mode == DAY_BLOCKS )
    {
      end_time.tm_mday++;
    }
  end_time_t = mktime(&end_time);
  time_diff = difftime(end_time_t,block.start_time);
  // # RG mod
  block.end_sample = block.start_sample + time_diff/delta -1;

  // continue breaking into hour blocks until we're left with an 
  // incomplete hour at the end
  while ( block.end_sample < nsamples )
    {
      blocks.push_back(block);
      // # RG mod
      block.start_sample = block.end_sample +1;  
      block.start_time = end_time_t;
      if ( mode == HOUR_BLOCKS )
	{
	  end_time.tm_hour++;
	}
      else if ( mode == DAY_BLOCKS )
	{
	  end_time.tm_mday++;
	}
      end_time_t = mktime(&end_time);
      time_diff = difftime(end_time_t,block.start_time);
      // # RG mod
      block.end_sample = block.start_sample + time_diff/delta -1;

    }
  
  // include the last part block if it is greater than fraction frac 
  // of a block
  block.end_sample = nsamples;
  time_diff = ( block.end_sample - block.start_sample ) * delta;
  if ( time_diff > frac*sec_per_block )
    {
      blocks.push_back(block);
    }
  

  return blocks;

}

// ---------------------------------------------------------------------- //
// function to write out split sac files
// ---------------------------------------------------------------------- //
int write_out_blocks(string infilename,string outdir,
		     vector<Block> blocks,int msec,struct Sac *sac_in)
{
  // break up filename and extract relevant parts
  vector<string> path_parts = tokenize(infilename,"/");
  string filename = path_parts.back();

  vector<string> split_file = tokenize(filename,".");
  string net = split_file.at(6);
  string stn = split_file.at(7);
  string strm = split_file.at(8);
  string comp = split_file.at(9);
  string qual = split_file.at(10);

  struct Sac *sac_out;
  
  // check that blocks is of non-zero length
  if ( blocks.size() == 0 )
    {
      cout << "No data written: file is of insufficient length to split" 
	   << endl;
      return 0;
    }

  // create directory for output
  boost::filesystem::path outpath(outdir);
  if (! boost::filesystem::is_directory(outpath)){
    if ( boost::filesystem::exists(outpath) ){
	cerr << "File with name of output directory already exists" << endl;
	return 1;
      }
    printf("Making directory: %s\n",outdir.c_str());
    boost::filesystem::create_directory(outpath);
  }


  // iterate over blocks
  for (vector<Block>:: iterator it = blocks.begin(); it != blocks.end(); ++it)
    {
      int start = it->start_sample;
      int end = it->end_sample;
      time_t stime = it->start_time;

      // copy across info to new Sac structure
      sac_out = new struct Sac;
      sac_out->header = sac_in->header;
      
      // update header
      sac_out->header.npts = end - start + 1;
      struct tm begin_time = *gmtime(&stime);
      sac_out->header.nzyear = begin_time.tm_year + 1900;
      sac_out->header.nzjday = begin_time.tm_yday + 1;
      sac_out->header.nzhour = begin_time.tm_hour;
      sac_out->header.nzmin = begin_time.tm_min;
      sac_out->header.nzsec = begin_time.tm_sec;
      sac_out->header.nzmsec = msec;
      sac_out->header.b = 0.;

      // copy required section out of array
      sac_out->y = new float[sac_out->header.npts];
      int n = 0;
      for (int i = start - 1; i < end; i++)
	{
	  sac_out->y[n] = sac_in->y[i];
	  n++;
	}
      
      // ** write file to output directory **
      
      // put together file name
      char outstr[100];
      sprintf(outstr,"%04d.%03d.%02d.00.00.%03d0.%s.%s.%s.%s.%s.SAC",
	      sac_out->header.nzyear,sac_out->header.nzjday,
	      sac_out->header.nzhour,sac_out->header.nzmsec,
	      net.c_str(),stn.c_str(),strm.c_str(),
	      comp.c_str(),qual.c_str());

      // concatenate filename with output directory
      boost::filesystem::path outfile = operator/(outpath,outstr);
      string outfilename = outfile.string();
      cout << "Writing: " << outfilename << endl;

      // write out files in output directory
      int stat = write_sacfile(outfilename,sac_out);
      if ( stat ) 
	{
	  fprintf(stderr,"Error writing file: %s", outfilename.c_str());
	  return -1;
	}

      // clean up memory
      delete [] sac_out->y;
      delete [] sac_out;
      
      
    }

  return 0;

}

// ====================================================================== //
// MAIN PROGRAM
// ====================================================================== //
int main(int argc, char **argv)
{

  struct Sac *sac_in;
  vector<Block> blocks;

  // set up parser and variables to store parsed data
  string usage = "split_data [options] input_file.SAC";
  string outdir;
  
  po::options_description options("Allowed options");
  options.add_options()
    ("help,h", "output this help message")
    ("day,d", "split data into day blocks (default)")
    ("hour,u","split data into hour blocks")
    ("output-dir,o",po::value<string>(&outdir)->default_value("split_files"),
     "set output directory")
    ;

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file",po::value<string>(),"input file");

  po::positional_options_description pos;
  pos.add("input-file",1);

  po::options_description all_options("All options, including hidden");
  all_options.add(options).add(hidden);

  po::variables_map vm;

  // parse the command line
  try 
    {
      po::store(po::command_line_parser(argc,argv).options(all_options).
		positional(pos).run(),vm);
      po::notify(vm);

      // validate input
      if ( vm.count("help") )
	{
	  cerr << endl;
	  cerr << usage << endl;
	  cerr << endl;
	  cerr << options << endl;
	  return 1;
	}
      
      if (! vm.count("input-file")) 
	{
	  cerr << "Must supply input SAC filename" << endl;
	  cerr << usage << endl;
	  cerr << endl;
	  cerr << options << endl;
	  return 1;
	}

      if ( vm.count("day") && vm.count("hour") )
	{
	  cerr << "Cannot specify two block lengths for output ... choose one"
	       << endl;
	  cerr << usage << endl;
	  cerr << endl;
	  cerr << options << endl;
	  return 1;
	}

    }
  
  catch(exception& e) 
    {
      cerr << "Error: " << e.what() << endl;
      cerr << usage << endl;
      cerr << endl;
      cerr << options << endl;
      return 1;
    }

  // assign command line options to variables
  string sacfilename = vm["input-file"].as<string>();

  int block_type = DAY_BLOCKS;
  if ( vm.count("day") ) 
    {
      block_type = DAY_BLOCKS;
    }
  else if ( vm.count("hour") )
    {
      block_type = HOUR_BLOCKS;
    }

  // read sacfile and check for errors
  sac_in = read_sacfile(sacfilename);
  if ( sac_in == NULL ) 
    {
      return -1;
    }

  // get reference time as type time_t from SAC header values
  // -- ref_time_t and ref_msec
  time_t ref_time_t = get_time_t(sac_in->header.nzyear,
				 sac_in->header.nzjday,
				 sac_in->header.nzhour,
				 sac_in->header.nzmin,
				 sac_in->header.nzsec);
  int ref_msec = sac_in->header.nzmsec;

  // get begin time as type time_t
  // -- begin_time_t and begin_msec
  time_t begin_time_t = ref_time_t + int(sac_in->header.b);
  int begin_msec = ref_msec + 1000*(sac_in->header.b - int(sac_in->header.b));
  while ( begin_msec >= 1000 )
    {
      begin_time_t++;
      begin_msec -= 1000;
    }
  // define split into blocks
  // variable block_type determines which length blocks to split into
  blocks = split_into_blocks(begin_time_t,sac_in->header.delta,
			     sac_in->header.npts,block_type);


  // write out files
  int stat = write_out_blocks(sacfilename,outdir,blocks,begin_msec,sac_in);
  if ( stat != 0 ) 
    {
      printf("Error writing out files\n");
      return -1;
    }


  return 0;

}


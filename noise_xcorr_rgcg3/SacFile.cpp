// SacFile.cpp
//
// Functions for reading and writing SAC files
//
// Jamie Barron
// January 2010

#include "SacFile.hpp"

using namespace std;


// ************************************************************** //
// function to read in binary SAC file and return a Sac structure //
// which contains the file header and the file data               //
// ************************************************************** //
struct Sac *read_sacfile(const string filename) {

  // generate structure to read sacfile into
  struct Sac *sacf;
  sacf = new struct Sac;

  // size of sac header 
  unsigned int hsize = sizeof(struct SacHeader);

  // check size and existence of file 
  struct stat sb;
  if ( stat(filename.c_str(), &sb) == -1 )
    {
      fprintf(stderr,"File %s does not exist\n",filename.c_str());
      return NULL;
    }
  
  unsigned fsize = (unsigned int) sb.st_size;
  if ( fsize < hsize )
    {
      fprintf(stderr,"File %s is not a valid SAC file\n",filename.c_str());
      return NULL;
    }
  
  // open file stream for reading
  FILE *sac_in;
  sac_in = fopen(filename.c_str(),"rb");
  if ( sac_in == NULL )
    {
      fprintf(stderr,"Error opening binary file for reading: %s\n",
	      filename.c_str());
      return NULL;
    }
  
  // read in header info
  int nread = fread(&sacf->header, hsize, 1, sac_in);
  if ( nread != 1 )
    {
      fprintf(stderr,"Error reading header from file: %s\n", filename.c_str());
      return NULL;
    }

  // check that valid information was read from file
  if (sacf->header.delta <= 0.0 || sacf->header.npts <= 0 ||
      sacf->header.nvhdr < 0 ) 
    {
      fprintf(stderr,"File %s does not have a valid SAC header\n",
	      filename.c_str());
      return NULL;
    }

  if ( (hsize + sacf->header.npts*sizeof(float)) < fsize ) 
    {
      fprintf(stderr,"File %s is corrupt\n",filename.c_str());
      return NULL;
    }

  // allocate memory for the rest of the data
  sacf->y = new float [sacf->header.npts];
  
  // read in remaining data and close file
  nread = fread(sacf->y,sizeof(float),sacf->header.npts,sac_in);
  fclose(sac_in);

  // check that points were read in successfully
  if ( nread != sacf->header.npts )
    {
      fprintf(stderr,"Unable to read data from SAC file %s\n",filename.c_str());
      return NULL;
    }


  return sacf;

}

// 
// function to write out binary SAC file given a Sac structure
// (header and time series)
//
int write_sacfile(const string filename, struct Sac *sacf)
{

  FILE *strm;
  int nwritten;

  // make sure end value is consistent with begin, delta and npts
  sacf->header.e = sacf->header.b + (sacf->header.npts-1)*sacf->header.delta;

  // open stream to write
  strm = fopen(filename.c_str(),"w");
  if (strm == NULL)
    {
      fprintf(stderr,"Error opening binary file for writing: %s\n",
	      filename.c_str());
      return 1;
    }

  // write the file header
  nwritten = fwrite(&sacf->header,sizeof(struct SacHeader),1,strm);
  if ( nwritten != 1 )
    {
      fprintf(stderr,"Error writing to binary file: %s\n",filename.c_str());
      fclose(strm);
      return 1;
    }

  // write the time series to file
  nwritten = fwrite(sacf->y,sacf->header.npts*sizeof(float),1,strm);
  fclose(strm);
  if ( nwritten != 1 )
    {
      fprintf(stderr, "Error writing to binary file: %s\n", filename.c_str());
      return 1;
    }

  return 0;

}




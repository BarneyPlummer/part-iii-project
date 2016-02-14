#pragma once
#include <string>
#include <cstdio>
#include <sys/stat.h>
#include "sacheader.h"

// structure to hold sac data and header 
struct Sac {
  struct SacHeader header;
  float *y;

  Sac()
  {
    header = sac_null;
  }

};

// function prototypes
struct Sac *read_sacfile(const std::string filename);
int write_sacfile(const std::string filename, struct Sac *sac);


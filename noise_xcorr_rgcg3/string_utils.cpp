// **********************************************************************
// string_utils.cpp
//
// String processing utility functions
//
// Author: Jamie Barron
// Date: February 2011
// **********************************************************************

#include "string_utils.hpp"
#include <boost/algorithm/string.hpp>

vector<string> tokenize(string instring,const string sep)
{

  vector<string> tokens;
  boost::split(tokens,instring,boost::is_any_of(sep));
  return tokens;

}

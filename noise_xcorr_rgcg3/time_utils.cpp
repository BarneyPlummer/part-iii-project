// ********************************************************************** //
// time_utils.cpp
//
// Author: Jamie Barron
// Date: Feb 2011
// ********************************************************************** //

#include "time_utils.hpp"

// function to convert Julian day to month/day
void convert_julian_to_month_day(int year,int jday,int *month, int *day)
{
  
  bool leapyear;
  int i;
  int summed_mon_len[12];

  // is year a leapyear
  if ( year % 400 == 0 ) 
    { 
      leapyear = true;
    }
  else if ( year % 100 == 0 )
    { 
      leapyear = false;
    }
  else if ( year % 4 == 0 )
    {
      leapyear = true;
    }
  else
    { 
      leapyear = false;
    }

  // set month length array
  if ( leapyear )
    {
      summed_mon_len[0] = 0;
      summed_mon_len[1] = 31;
      summed_mon_len[2] = 60;
      summed_mon_len[3] = 91;
      summed_mon_len[4] = 121;
      summed_mon_len[5] = 152;
      summed_mon_len[6] = 182;
      summed_mon_len[7] = 213;
      summed_mon_len[8] = 244;
      summed_mon_len[9] = 274;
      summed_mon_len[10] = 305;
      summed_mon_len[11] = 335;
      
      //summed_mon_len = { 0,    31,  60,  91, 121, 152, 
      //		  182, 213, 244, 274, 305, 335 };
    }
  else
    {
      summed_mon_len[0] = 0;
      summed_mon_len[1] = 31;
      summed_mon_len[2] = 59;
      summed_mon_len[3] = 90;
      summed_mon_len[4] = 120;
      summed_mon_len[5] = 151;
      summed_mon_len[6] = 181;
      summed_mon_len[7] = 212;
      summed_mon_len[8] = 243;
      summed_mon_len[9] = 273;
      summed_mon_len[10] = 304;
      summed_mon_len[11] = 334;
      
      //summed_mon_len = { 0,    31,  59,  90, 120, 151, 
      //		 181, 212, 243, 273, 304, 334 };
    }

  
  // calculate month and day
  *month = 12;
  for (i=0; i<12; i++)
    {
      if ( summed_mon_len[i] > jday )
	{
	  *month = i;
	  break;
	}
    }

  *day = jday - summed_mon_len[*month-1];

}

time_t get_time_t(int year, int jday, int hour, int min, int sec)
{

  int month;
  int day;

  // get month and day from year and Julian day
  convert_julian_to_month_day(year,jday,&month,&day);

  // create tm structure and assign to it
  struct tm time;
  time.tm_year = year - 1900;
  time.tm_mon = month - 1;
  time.tm_mday = day;
  time.tm_hour = hour;
  time.tm_min = min;
  time.tm_sec = sec;
  time.tm_isdst = 0;

  // convert to time_t and return
  time_t tt = mktime(&time);
  return tt;

}

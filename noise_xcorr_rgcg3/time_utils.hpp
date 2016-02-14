// time_utils.hpp
#pragma once

// includes
#include <ctime>

// function prototypes
void convert_julian_to_month_day(int jday,int *month, int *day);
time_t get_time_t(int year, int jday, int hour, int min, int sec);

/// _____________________________________________________________________ ///
///                                                                       ///
/// BusyFit 0.3.2 (helperFunctions.hpp) - Busy Function fitting programme ///
/// Copyright (C) 2019 Tobias Westmeier                                   ///
/// _____________________________________________________________________ ///
///                                                                       ///
/// Address:  Tobias Westmeier                                            ///
///           ICRAR M468                                                  ///
///           The University of Western Australia                         ///
///           35 Stirling Highway                                         ///
///           Crawley WA 6009                                             ///
///           Australia                                                   ///
///                                                                       ///
/// E-mail:   tobias.westmeier (at) uwa.edu.au                            ///
/// _____________________________________________________________________ ///
///                                                                       ///
/// This program is free software: you can redistribute it and/or modify  ///
/// it under the terms of the GNU General Public License as published by  ///
/// the Free Software Foundation, either version 3 of the License, or     ///
/// (at your option) any later version.                                   ///
///                                                                       ///
/// This program is distributed in the hope that it will be useful,       ///
/// but WITHOUT ANY WARRANTY; without even the implied warranty of        ///
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          ///
/// GNU General Public License for more details.                          ///
///                                                                       ///
/// You should have received a copy of the GNU General Public License     ///
/// along with this program. If not, see http://www.gnu.org/licenses/.    ///
/// _____________________________________________________________________ ///
///                                                                       ///

#ifndef HELPERFUNCTIONS_HPP
#define HELPERFUNCTIONS_HPP

#define TRIM_BOTH  0
#define TRIM_LEFT  1
#define TRIM_RIGHT 2

#define MATH_CONST_PI 3.141592653589793
#define MATH_CONST_E  2.718281828459045

#define MSG_INFO    0
#define MSG_WARNING 1
#define MSG_ERROR   2

#include <string>
#include <vector>
#include <sstream>
#include <limits>
#include <iomanip>
#include <typeinfo>
#include <cstdlib>
// WARNING: sys/time.h needed for gettimeofday(), but only works in Linux/Unix!
#include <sys/time.h>


// String functions:

void stringToUpper(std::string &strToConvert);
void stringToLower(std::string &strToConvert);

void stringTrim(std::string &strToConvert, int mode = TRIM_BOTH);
int  stringTok(std::string &strToConvert, std::vector<std::string> &tokens, const std::string &delimiters = " ");

int  stringReplace(std::string &strToConvert, const std::string &seachStr, const std::string &replStr);

std::string degToDms(double number);
std::string degToHms(double number);



// Mathematical functions:

int    mathSgn(int value);
int    mathSgn(long value);
int    mathSgn(float value);
int    mathSgn(double value);

int    mathAbs(int value);
long   mathAbs(long value);
float  mathAbs(float value);
double mathAbs(double value);

long   mathRound(float value);
long   mathRound(double value);



// Statistical functions:

double statMean(const std::vector<double> &input);
double statStdDev(const std::vector<double> &input);
double statMedian(const std::vector<double> &input, const bool sort = true);
double statMad(const std::vector<double> &input, const bool sort = true);



// Random number generator:

void   initRand();
double randGauss(double sigma = 1.0);



// Printing error messages:

void printErrorMessage(const char *message, int type = MSG_ERROR);



// String-number conversion templates:

template <typename T> std::string numberToString(T number, int decimals = -1, bool scientific = false)
{
	std::ostringstream convert;
	
	if(decimals >= 0 and (typeid(T) == typeid(float) or typeid(T) == typeid(double)))
	{
		if(decimals > std::numeric_limits<T>::digits10 + 1) decimals = std::numeric_limits<T>::digits10 + 1;
		
		// Choose between fixed-point and scientific notation for float and double types:
		if(scientific == true) convert << std::scientific;
		else                   convert << std::fixed;
		
		// Set the number of digits for float and double types:
		convert << std::setprecision(decimals);
	}
	
	// Let's hope that this makes any sense for a particular type T:
	convert << number;
	
	return convert.str();
}



// Template function to convert std::string to numerical value of type T:
// WARNING: This function does not correctly convert numbers in scientific 
// notation to integer types! For example, "-3.95E+1" will become -39.5 
// when converted to float, but -3 when converted to int!
template <typename T> T stringToNumber(const std::string &strToConvert)
{
	std::stringstream stringStream(strToConvert);
	T result;
	return stringStream >> result ? result : 0;
}



#endif

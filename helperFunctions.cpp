/// ____________________________________________________________________ ///
///                                                                      ///
/// BusyFit 0.3.1 (helperFunctions.cpp) - Busy Function fitting programme///
/// Copyright (C) 2017 Tobias Westmeier                                  ///
/// ____________________________________________________________________ ///
///                                                                      ///
/// Address:  Tobias Westmeier                                           ///
///           ICRAR M468                                                 ///
///           The University of Western Australia                        ///
///           35 Stirling Highway                                        ///
///           Crawley WA 6009                                            ///
///           Australia                                                  ///
///                                                                      ///
/// E-mail:   tobias.westmeier (at) uwa.edu.au                           ///
/// ____________________________________________________________________ ///
///                                                                      ///
/// This program is free software: you can redistribute it and/or modify ///
/// it under the terms of the GNU General Public License as published by ///
/// the Free Software Foundation, either version 3 of the License, or    ///
/// (at your option) any later version.                                  ///
///                                                                      ///
/// This program is distributed in the hope that it will be useful,      ///
/// but WITHOUT ANY WARRANTY; without even the implied warranty of       ///
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         ///
/// GNU General Public License for more details.                         ///
///                                                                      ///
/// You should have received a copy of the GNU General Public License    ///
/// along with this program. If not, see http://www.gnu.org/licenses/.   ///
/// ____________________________________________________________________ ///
///                                                                      ///

#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <ctime>        // Needed for seed of random number generator.
#include <cstdlib>      // Needed for random number generator.

#include "helperFunctions.hpp"



// Convert std::string to upper and lower case:

void stringToUpper(std::string &strToConvert)
{
	if(strToConvert.empty()) return;
	
	std::transform(strToConvert.begin(), strToConvert.end(), strToConvert.begin(), ::toupper);
	
	return;
}

void stringToLower(std::string &strToConvert)
{
	if(strToConvert.empty()) return;
	
	std::transform(strToConvert.begin(), strToConvert.end(), strToConvert.begin(), ::tolower);
	
	return;
}



// Trim std::string by removing whitespace characters:
// mode = 0 (TRIM_BOTH)  ->  trim both ends (default)
//        1 (TRIM_LEFT)  ->  trim leading whitespaces only
//        2 (TRIM_RIGHT) ->  trim trailing whitespaces only

void stringTrim(std::string &strToConvert, int mode)
{
	if(strToConvert.empty()) return;
	
	std::string ws(" \t\n\r\0\v");               // Same whitespace characters as used in PHP's trim()
	size_t foundFirst;
	size_t foundLast;
	
	if(mode < 0 or mode > 2) mode = 0;
	
	if(mode == 0)
	{
		// Trim both ends (default):
		foundFirst = strToConvert.find_first_not_of(ws);
		foundLast  = strToConvert.find_last_not_of(ws);
		
		if(foundFirst != std::string::npos and foundLast != std::string::npos) strToConvert.assign(strToConvert.substr(foundFirst, foundLast - foundFirst + 1));
		else strToConvert.clear();                    // String contains only whitespace
	}
	else if(mode == 1)
	{
		// Trim leading whitespaces:
		foundFirst = strToConvert.find_first_not_of(ws);
		
		if(foundFirst != std::string::npos) strToConvert.assign(strToConvert.substr(foundFirst));
		else strToConvert.clear();                    // String contains only whitespace
	}
	else
	{
		// Trim trailing whitespaces:
		foundLast = strToConvert.find_last_not_of(ws);
		
		if(foundLast != std::string::npos) strToConvert.assign(strToConvert.substr(0, foundLast + 1));
		else strToConvert.clear();                    // String contains only whitespace
	}
	
	return;
}



// Tokenize std::string:

int stringTok(std::string &strToConvert, std::vector<std::string> &tokens, const std::string &delimiters)
{
	if(strToConvert.empty() or delimiters.empty()) return 1;
	
	tokens.clear();
	
	// Ignore initial delimiters:
	std::string::size_type firstPosition = strToConvert.find_first_not_of(delimiters, 0);
	
	if(firstPosition == std::string::npos) return 1;       // String containing delimiters only.
	
	// Find first "true" delimiter:
	std::string::size_type firstDelim = strToConvert.find_first_of(delimiters, firstPosition);
	
	while(firstDelim != std::string::npos or firstPosition != std::string::npos)
	{
		tokens.push_back(strToConvert.substr(firstPosition, firstDelim - firstPosition));
		
		firstPosition = strToConvert.find_first_not_of(delimiters, firstDelim);
		firstDelim    = strToConvert.find_first_of(delimiters, firstPosition);
	}
	
	return 0;
}



// Replace all occurrences of search string with replacement string:

int stringReplace(std::string &strToConvert, const std::string &seachStr, const std::string &replStr)
{
	if(strToConvert.empty() or seachStr.empty()) return 1;
	
	size_t firstPos = 0;
	size_t pos = strToConvert.find(seachStr, firstPos);
	
	if(pos == std::string::npos) return 1;
	
	while(pos != std::string::npos)
	{
		strToConvert = strToConvert.substr(0, pos) + replStr + strToConvert.substr(pos + seachStr.size());
		firstPos     = pos + replStr.size();
		if(firstPos < strToConvert.size()) pos = strToConvert.find(seachStr, firstPos);
		else break;
	}
	
	return 0;
}



// Convert decimal degrees to dd:mm:ss.ss string:

std::string degToDms(double number)
{
	double number2 = fabs(number);
	int    deg;
	int    min;
	double sec;
	
	deg = static_cast<int>(number2);
	min = static_cast<int>(60.0 * (number2 - static_cast<double>(deg)));
	sec = 60.0 * (60.0 * (number2 - static_cast<double>(deg)) - static_cast<double>(min));
	
	std::ostringstream convert;
	
	if(number < 0.0) convert << "-";
	
	convert << std::fixed << std::setw(2) << std::setfill('0') << deg << ":" << std::setw(2) << std::setfill('0') << min << ":" << std::setw(5) << std::setprecision(2) << sec;
	
	return convert.str();
}



// Convert decimal degrees to hh:mm:ss.ss string:

std::string degToHms(double number)
{
	number /= 15.0;
	
	return degToDms(number);
}



// Return sign of value:

int mathSgn(int value)
{
	if(value < 0) return -1;
	else return  1;
}

int mathSgn(long value)
{
	if(value < 0L) return -1;
	else return  1;
}

int mathSgn(float value)
{
	if(value < 0.0F) return -1;
	else return  1;
}

int mathSgn(double value)
{
	if(value < 0.0) return -1;
	else return  1;
}



// Return absolute value of argument:

int mathAbs(int value)
{
	if(value < 0) return -value;
	else return value;
}

long mathAbs(long value)
{
	if(value < 0L) return -value;
	else return value;
}

float mathAbs(float value)
{
	if(value < 0.0F) return -value;
	else return value;
}

double mathAbs(double value)
{
	if(value < 0.0) return -value;
	else return value;
}



// Round value:

long mathRound(float value)
{
	return static_cast<long>(value + 0.5);
}

long mathRound(double value)
{
	return static_cast<long>(value + 0.5);
}



// Initialise random number generator:

void initRand()
{
	timeval t1;
	gettimeofday(&t1, NULL);
	
	srand(t1.tv_usec * t1.tv_sec);          // Initialise random number generator.
	//srand(time(NULL));                    // Don't use this, only operates on 1 s precision!!!
	
	return;
}



// Return Gaussian random deviate:

double randGauss(double sigma)
{
	if(sigma <= 0.0) return 0.0;
	
	double rand1 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
	double rand2 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
	
	return sqrt(-2.0 * sigma * sigma * log(rand1)) * cos(2.0 * MATH_CONST_PI * rand2);
}



// Print error message:

void printErrorMessage(const char *message, int type)
{
	if(type == 2)      std::cerr << "\033[1;31m" << "Error: "   << message << "\033[0m" << '\n' << std::endl;
	else if(type == 1) std::cerr << "\033[1;33m" << "Warning: " << message << "\033[0m" << '\n' << std::endl;
	else               std::cerr << "\033[1;32m" << "Note: "    << message << "\033[0m" << '\n' << std::endl;
	return;
}



// Calculate mean

double statMean(const std::vector<double> &input)
{
	if(input.size() == 0) return std::numeric_limits<double>::quiet_NaN();
	
	double result  = 0.0;
	size_t counter = 0;
	
	for(size_t i = 0; i < input.size(); i++)
	{
		if(std::isfinite(input[i]))
		{
			result += input[i];
			counter++;
		}
	}
	
	if(counter == 0) return std::numeric_limits<double>::quiet_NaN();
	
	return result / static_cast<double>(counter);
}



// Calculate standard deviation

double statStdDev(const std::vector<double> &input)
{
	if(input.size() == 0) return std::numeric_limits<double>::quiet_NaN();
	
	double mean    = statMean(input);
	double result  = 0.0;
	size_t counter = 0;
	
	for(size_t i = 0; i < input.size(); i++)
	{
		if(std::isfinite(input[i]))
		{
			result += (input[i] - mean) * (input[i] - mean);
			counter++;
		}
	}
	
	if(counter < 2) return std::numeric_limits<double>::quiet_NaN();
	
	return sqrt(result / static_cast<double>(counter - 1));
}



// Calculate median

double statMedian(const std::vector<double> &input, const bool sort)
{
	if(input.size() == 0) return std::numeric_limits<double>::quiet_NaN();
	
	std::vector<double> inputCopy;
	
	// Create copy with only non-nan/inf values
	for(size_t i = 0; i < input.size(); i++)
	{
		if(std::isfinite(input[i]))
		{
			inputCopy.push_back(input[i]);
		}
	}
	
	size_t arraySize = inputCopy.size();
	if(arraySize == 0) return std::numeric_limits<double>::quiet_NaN();
	
	// Sort if requested (default behaviour)
	if(sort) std::sort(inputCopy.begin(), inputCopy.end());
	
	if(arraySize % 2 == 0) return (inputCopy[arraySize / 2 - 1] + inputCopy[arraySize / 2]) / 2.0;
	else return inputCopy[arraySize / 2];
}



// Calculate median absolute deviation (MAD)

double statMad(const std::vector<double> &input, const bool sort)
{
	if(input.size() == 0) return std::numeric_limits<double>::quiet_NaN();
	
	std::vector<double> inputCopy;
	
	// Create copy with only non-nan/inf values
	for(size_t i = 0; i < input.size(); i++)
	{
		if(std::isfinite(input[i]))
		{
			inputCopy.push_back(input[i]);
		}
	}
	
	if(inputCopy.size() == 0) return std::numeric_limits<double>::quiet_NaN();
	
	double median = statMedian(inputCopy, sort);
	
	for(size_t i = 0; i < inputCopy.size(); i++)
	{
		inputCopy[i] = fabs(inputCopy[i] - median);
	}
	
	return statMedian(inputCopy);
}

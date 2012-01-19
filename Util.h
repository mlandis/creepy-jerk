/*
 * Util.h
 *
 *  Created on: Jan 3, 2012
 *      Author: mlandis
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>

namespace Util
{

	extern std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems);
	extern std::vector<std::string> split(const std::string &s, char delim);
	extern std::string randString(int len);
	extern std::string getTime(void);

	extern std::string boolToString(bool);
	extern std::string doubleToString(double);
	extern std::string intToString(int);

	extern int stringToInt(std::string);
	extern double stringToDouble(std::string);
	extern bool stringToBool(std::string);

}

#endif /* UTIL_H_ */

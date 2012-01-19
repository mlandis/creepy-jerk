/*
 * Util.cpp
 *
 *  Created on: Jan 3, 2012
 *      Author: mlandis
 */

#include "Util.h"

std::vector<std::string>& Util::split(const std::string &s, char delim, std::vector<std::string> &elems)
{
	std::stringstream ss(s);
	std::string item;
	while(std::getline(ss, item, delim))
	{
		elems.push_back(item);
	}
	return elems;
}


std::vector<std::string> Util::split(const std::string &s, char delim)
{
	std::vector<std::string> elems;
	return split(s, delim, elems);
}

std::string Util::randString(int len)
{

	std::string s;

	static const char alphanum[] =
		"0123456789"
		"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		"abcdefghijklmnopqrstuvwxyz";

	for (int i = 0; i < len; ++i)
	{
		s.push_back( alphanum[rand() % (sizeof(alphanum) - 1)] );
	}

	return s;
}

std::string Util::getTime(void)
{
	std::string s;

	time_t rawtime;
	struct tm * timeinfo;
	char buffer [80];

	time ( &rawtime );
	timeinfo = localtime ( &rawtime );

	// two-digit representation
	// YYMMDDhhmmss
	strftime (buffer,80,"%y%m%d%H%M%S",timeinfo);

	s = buffer;

	return s;
}

bool Util::stringToBool(std::string s)
{
	if (s == "True")
		return true;
	else return false;
}


std::string Util::boolToString(bool tf)
{
	if (tf) return "True";
	else return "\tFalse";
}

std::string Util::intToString(int number)
{
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

std::string Util::doubleToString(double number)
{
	std::stringstream ss;
	ss << number;
	return ss.str();
}

int Util::stringToInt(std::string s)
{
	return atoi(s.c_str());
}

double Util::stringToDouble(std::string s)
{
	return atof(s.c_str());
}

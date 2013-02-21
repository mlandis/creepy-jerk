/*
 * Patron.h
 *
 *  Created on: Mar 14, 2011
 *      Author: mlandis
 */

#ifndef PATRON_H_
#define PATRON_H_

#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include "Parm.h"
#include "Table.h"

class Parm;
class Table;

class Patron {
public:
	Patron(Table* t, std::vector<std::vector<double> > d, double*, int, double stddev, std::string, int);
	Patron(std::vector<std::vector<double> > d, double*, int, double stddev, std::string, int);
	~Patron();

	void stand(void);
	void sit(Table* t);

	Table*										getTable(void) { return table; }
	const std::vector<std::vector<double> >&	getData(void) { return data; }
	double										getData(int taxon, int timepoint) { return data[taxon][timepoint]; }
	gsl_complex									getCfData(int, int, int);
	const std::vector<Parm*>&					getParmVector(void);
	std::string									getName(void) { return name; }
	int											getId(void) { return id; }

	void										print(void);
	void										printCf(void);
	std::string									getPrintStr(void);

private:
	gsl_complex									charFunc(double, double, double);

	std::string									name;
	int											id;
	Table*										table;
	std::vector<std::vector<double> >			data;
	double*										cfData;

	int											numSteps;
	int											numTaxa;
	int											numTimepoints;
	int											oneTaxonSize;
	int											oneTimeSize;

	std::string									printHeaderStr;
	std::string									printBodyStr;
};

#endif /* PATRON_H_ */

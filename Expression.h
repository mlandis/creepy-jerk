/*
 * Expression.h
 *
 *  Created on: Mar 7, 2011
 *      Author: mlandis
 */

#ifndef EXPRESSION_H_
#define EXPRESSION_H_

#include "FileMgr.h"
#include "Patron.h"
#include "Settings.h"

#include <algorithm>
#include <dirent.h>
#include <iostream>
#include <list>
#include <string.h>
#include <sstream>
#include <vector>

class FileMgr;
class Patron;
class Settings;

class Expression {

public:
												Expression(Settings*);
												~Expression();

	void										initializeExprData(void);
	void										initializeData(void);
	void										initializeTaxaNames(void);

	int											getNumTaxa(void)						{ return numTaxa; }
	int											getNumTimepoints(void)					{ return numTimepoints; }
	int											getNumTranscripts(void)					{ return numTranscripts; }
	int											getIndexForTaxon(std::string taxonName);
	std::string									getNameForTaxon(int i)					{ return taxonNames[i]; }

	double										getExpr(int trans, int taxon, int time)	{ return transExpr[trans][taxon][time]; }
	const std::vector<double>&					getExpr(int trans, int taxon)			{ return transExpr[trans][taxon]; }
	const std::vector<std::vector<double> >&	getExpr(int trans)						{ return transExpr[trans]; }
	double										getData(int taxon)						{ return data[taxon]; }

	const std::list<Table*>&					getTableList(void)						{ return tableList; }
	const std::list<Patron*>&					getPatronList(void)						{ return patronList; }

	void										print(void);

private:
	int numTaxa;
	int numTimepoints;
	int numTranscripts;
	double tipStdDev;

	// [trans][taxa][time]
	std::vector<std::vector<std::vector<double> > > transExpr;
	std::list<Patron*> patronList;
	std::list<Table*> tableList;
	std::vector<double> data;

	std::vector<std::string> taxonNames;

	bool useCRP;

	Settings* settingsPtr;
};

#endif /* EXPRESSION_H_ */


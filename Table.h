/*
 * Table.h
 *
 *  Created on: Mar 14, 2011
 *      Author: mlandis
 */

#ifndef TABLE_H_
#define TABLE_H_

#include <iostream>
#include <list>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>

#include "Parm_lambda.h"
#include "Parm_sigma.h"
#include "Parm_tau.h"
#include "Parm.h"
#include "Patron.h"

class Patron;
class Parm;
class Tau;
class Lambda;
class Sigma;

class Table {

public:

	Table(std::list<Table*>* tp, std::vector<Parm*> pv, std::vector<Parm*> bv, int n);
	Table(std::list<Table*>* tp, int);
	Table(Table* t, std::list<Table*>* tp);
	Table(const Table &t);
	~Table();
	Table&						operator=(const Table &t);
	void						setTableListPtr(std::list<Table*>* tp)	{ tableListPtr = tp; }
	void						setParmVector(std::vector<Parm*> pv)	{ parmVector = pv; }
	void						setBranchVector(std::vector<Parm*> bv)	{ branchVector = bv; }
	void						unseatPatron(void);
	void						unseatPatron(Patron* p);
	void						seatPatron(Patron* p);
	void						clone(const Table&);
	void						waitTable(void);
	const std::vector<Parm*>&	getParmVector(void)						{ return parmVector; }
	const std::vector<Parm*>&	getBranchVector(void)					{ return branchVector; }
	const std::list<Patron*>&	getPatronList(void)						{ return patronList; }
	int							getId(void)								{ return id; }

	// Started to implement. LnL values are conserved iff no patrons reseat & no params change
	bool						getUpdateLnL(void)						{ return updateLnL; }
	double						getLnL(void)							{ return lnL; }
	void						setUpdateLnL(bool v)					{ updateLnL = v; }
	void						setLnL(double v)						{ lnL = v; }

	std::string					getPrintStr(void);
	void						print(void);

private:

	std::list<Patron*> patronList;
	std::list<Table*>* tableListPtr;
	std::vector<Parm*> parmVector;
	std::vector<Parm*> branchVector;

	int id;

	bool updateLnL;
	double lnL;

	std::string printHeaderStr;
	std::string printBodyStr;

};


#endif /* TABLE_H_ */

/*
 * Table.cpp
 *
 *  Created on: Mar 14, 2011
 *      Author: mlandis
 */

#include "Table.h"

Table::Table(std::list<Table*>* tp, std::vector<Parm*> pv, std::vector<Parm*> bv, int i)
{
	tableListPtr = tp;
	parmVector = pv;
	branchVector = bv;
	lnL = 0.0;
	updateLnL = true;
	id = i;
}

Table::Table(std::list<Table*>* tp, int i)
{
	tableListPtr = tp;
	std::vector<Parm*> emptyVector;
	parmVector = emptyVector;
	branchVector = emptyVector;
	lnL = 0.0;
	updateLnL = true;
	id = i;
}

Table::Table(const Table &t)
{
	clone(t);
}

Table::Table(Table* t, std::list<Table*>* tp)
{
	clone(*t);
	tableListPtr = tp;
	//std::cout << "Copied table, changed tableListPtr\n";
	print();
}

Table& Table::operator=(const Table &t) {

	if (this != &t)
		clone(t);
	return *this;
}

void Table::clone(const Table &t) {

	tableListPtr = t.tableListPtr;

	parmVector.clear();
	for (std::vector<Parm*>::const_iterator it_p = t.parmVector.begin(); it_p != t.parmVector.end(); it_p++)
	{
		Kappa* derivedPtrK = dynamic_cast<Kappa*> (*it_p);
		Sigma* derivedPtrS = dynamic_cast<Sigma*> (*it_p);
		Alpha* derivedPtrA = dynamic_cast<Alpha*> (*it_p);
		if (derivedPtrK != 0)
		{
			parmVector.push_back(new Kappa(*derivedPtrK, derivedPtrK->getRandomPtr(), derivedPtrK->getName()));
		}
		else if (derivedPtrS != 0)
		{
			parmVector.push_back(new Sigma(*derivedPtrS, derivedPtrS->getRandomPtr(), derivedPtrS->getName()));
		}
		else if (derivedPtrA != 0)
        {
            parmVector.push_back(new Alpha(*derivedPtrA, derivedPtrA->getRandomPtr(), derivedPtrA->getName()));
        }
	}

	branchVector.clear();
	for (std::vector<Parm*>::const_iterator it_p = t.branchVector.begin(); it_p != t.branchVector.end(); it_p++)
	{

	}

	patronList = t.patronList;
	lnL = t.lnL;
	updateLnL = t.updateLnL;
	id = t.id;
}

Table::~Table(void) {

	for (std::vector<Parm*>::iterator it_p = parmVector.begin(); it_p != parmVector.end(); it_p++)
	{
        Alpha* derivedPtrA = dynamic_cast<Alpha*> (*it_p);
		Kappa* derivedPtrK = dynamic_cast<Kappa*> (*it_p);
		Sigma* derivedPtrS = dynamic_cast<Sigma*> (*it_p);
		if (derivedPtrK != 0)
			delete derivedPtrK;
		else if (derivedPtrS != 0)
			delete derivedPtrS;
		else if (derivedPtrA != 0)
		    delete derivedPtrA;
	}
}

void Table::waitTable(void)
{
	if (patronList.size() == 0)
	{
		tableListPtr->remove(this);
		delete this;
	}
}

void Table::unseatPatron(void)
{
	patronList.pop_back();
	/*
	if (patronList.size() == 0)
	{
		tableListPtr->remove(this);
		delete this;
	}
	*/
}

void Table::unseatPatron(Patron* p)
{
	patronList.remove(p);

	/*
	if (patronList.size() == 0)
	{
		tableListPtr->remove(this);
		delete this;
	}
	*/
}

void Table::seatPatron(Patron* p)
{
	patronList.push_back(p);
}

void Table::print(void)
{
	std::cout << "TABLE:\n";
	std::cout << "\tlistPtr: " << &tableListPtr << "\n";
	std::cout << "\tthis:    " << this << "\n";
	std::cout << "\tsize:    " << patronList.size() << "\n";
	std::cout << "\tlnL:     " << lnL << "\n";
	std::cout << "\tavglnL:  " << lnL / patronList.size() << "\n";
	std::cout << "\tPATRONS:" << "\n";


	std::cout << "\t";
	for (std::list<Patron*>::iterator it_i = patronList.begin(); it_i != patronList.end(); it_i++)
	{
		std::cout << "\t" << (*it_i)->getName() << ": (" << &(*it_i);
		std::cout << "\t" << (*it_i)->getData()[0][0] << "  \t" << (*it_i)->getData()[1][0] << ")";
	}
	std::cout << "\n";

	std::cout << "\tPARAMETERS:\n";
	for (unsigned int i = 0; i < parmVector.size(); i++)
	{
		std::cout << "\t\t";
		parmVector[i]->print();
	}

	std::cout << "\tBRANCHES:\n";
	for (unsigned int i = 0; i < branchVector.size(); i++)
	{
		std::cout << "\t\t";
		branchVector[i]->print();
	}
	std::cout << "\n";
}

std::string Table::getPrintStr(void)
{
	std::string printStr = "TABLE_ID\t";

	std::stringstream out;
	out << id;

	printStr += out.str();
	printStr += "\n";
	for (std::list<Patron*>::iterator it_p = patronList.begin(); it_p != patronList.end(); it_p++)
	{
		printStr += (*it_p)->getPrintStr();
		printStr += "\n";
	}

	return printStr;
}


#if 0

/*
void Table::unseatPatron(void)
{
	patronList.pop_back();
	if (patronList.size() == 0)
	{
		tableListPtr->remove(this);
		delete this;
	}
}
*/



/*
void Table::seatPatronFront(Patron* p)
{
	patronList.push_front(p);
}
*/

#endif

/*
 * Patron.cpp
 *
 *  Created on: Mar 14, 2011
 *      Author: mlandis
 */

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

#include "Patron.h"

Patron::Patron(Table* t, std::vector<std::vector<double> > d, double* theta, int ns, double sd, std::string n, int m)
{
	name = n;
	id = m;
	table = NULL;
	sit(t);
	table = t;
	data = d;			// confirm data is stored for fast access (* arith)
	sit(table);

	numSteps = ns;
	numTaxa = data.size();
	numTimepoints = data[0].size();

	oneTimeSize = numSteps * 2;
	oneTaxonSize = numTimepoints * numSteps * 2;

	// NOTE:
	// had problems where elements of this Patron's cfdata were overwritten by another Patron's cfdata
	// e.g. this.cfdata[1][0] would contain another.cfdata[0][0]
	//      another.cfdata[1][,] would be empty
	// Symptoms indicate too little memory was being allocated, but found no evidence of this.
	// Doubled the size of the array which appears to have fixed the problem.
	// Not sure why the extra *2 is needed, but it fixes the problem

	// [numTaxa][numTimepoints][numSteps * (Re,Im)]
	cfData = new double[numTaxa * numTimepoints * numSteps * 2 * 2];

	//print();
	for (int i = 0; i < numTaxa; i++) // taxon
	{
		for (int j = 0; j < numTimepoints; j++) // timepoint
		{
			for (int k = 0; k < numSteps; k++)
			{

				gsl_complex val = charFunc(theta[k], data[i][j], sd);
				REAL(cfData, i * oneTaxonSize + j * oneTimeSize + k*2) = GSL_REAL(val);
				IMAG(cfData, i * oneTaxonSize + j * oneTimeSize + k*2) = GSL_IMAG(val);
				//std::cout << std::setprecision(8) << "\ti:" << i << "\tj:" << j << "\tk:" << k;
				//std::cout << "\t[" << i * oneTaxonSize + j * oneTimeSize + k*2 << "]";
				//std::cout << "\tmu:" << data[i][j] << "\t" << REAL(cfData, i * oneTaxonSize + j * oneTimeSize + k*2);
				//std::cout << " " << IMAG(cfData, i * oneTaxonSize + j * oneTimeSize + k*2) << "\t";
				//std::cout << &cfData[i * oneTaxonSize + j * oneTimeSize + k*2] << "\n";
			}
		}
	}

	//print();
}

Patron::Patron(std::vector<std::vector<double> > d, double* theta, int ns, double sd, std::string n, int m)
{
	name = n;
	id = m;
	table = NULL;
	data = d;			// confirm data is stored for fast access (* arith)

	numSteps = ns;
	numTaxa = data.size();
	numTimepoints = data[0].size();

	oneTimeSize = numSteps * 2;
	oneTaxonSize = numTimepoints * numSteps * 2;

	// NOTE:
	// had problems where elements of this Patron's cfdata were overwritten by another Patron's cfdata
	// e.g. this.cfdata[1][0] would contain another.cfdata[0][0]
	//      another.cfdata[1][,] would be empty
	// Symptoms indicate too little memory was being allocated, but found no evidence of this.
	// Doubled the size of the array which appears to have fixed the problem.
	// Not sure why the extra *2 is needed, but it fixes the problem

	// [numTaxa][numTimepoints][numSteps * (Re,Im)]
	cfData = new double[numTaxa * numTimepoints * numSteps * 2 * 2];

	//print();
	for (int i = 0; i < numTaxa; i++) // taxon
	{
		for (int j = 0; j < numTimepoints; j++) // timepoint
		{
			for (int k = 0; k < numSteps; k++)
			{

				gsl_complex val = charFunc(theta[k], data[i][j], sd);
				REAL(cfData, i * oneTaxonSize + j * oneTimeSize + k*2) = GSL_REAL(val);
				IMAG(cfData, i * oneTaxonSize + j * oneTimeSize + k*2) = GSL_IMAG(val);
				//std::cout << std::setprecision(8) << "\ti:" << i << "\tj:" << j << "\tk:" << k;
				//std::cout << "\t[" << i * oneTaxonSize + j * oneTimeSize + k*2 << "]";
				//std::cout << "\tmu:" << data[i][j] << "\t" << REAL(cfData, i * oneTaxonSize + j * oneTimeSize + k*2);
				//std::cout << " " << IMAG(cfData, i * oneTaxonSize + j * oneTimeSize + k*2) << "\t";
				//std::cout << &cfData[i * oneTaxonSize + j * oneTimeSize + k*2] << "\n";
			}
		}
	}
}

Patron::~Patron() {
	delete [] cfData;
}

void Patron::stand()
{
	if (table != NULL)
		table->unseatPatron(this);
	table = NULL;
}

void Patron::sit(Table* t)
{
	if (table != t)
	{
		t->seatPatron(this);
		//t->print();
		if (table != NULL)
		{
			table->unseatPatron(this);
		}
		table = t;
	}
}

const std::vector<Parm*>& Patron::getParmVector(void)
{
	return table->getParmVector();
}

gsl_complex Patron::getCfData(int taxon, int time, int step)
{
	int i = taxon * oneTaxonSize + time * oneTimeSize + step * 2;
	//std::cout << "TEST:\t" << i << "\t" << REAL(cfData, i) << "\t" << IMAG(cfData, i) << "\n";
	return gsl_complex_rect(REAL(cfData, i), IMAG(cfData, i));
}

void Patron::print(void)
{
	std::cout << "Patron:\tname:" << name << "\tid:" << id << "\tdata:" << &data << "\ttPtr:";
	if (table == NULL)
		std::cout << "NULL";
	else
		std::cout << table;
	std::cout << "\tcfData:" << &cfData << "\n";

	for (unsigned int i = 0; i < data.size(); i++)
	{
		for (unsigned int j = 0; j < data[i].size(); j++)
		{
			std::cout << "\t" << data[i][j];
		}
		std::cout << "\n";
	}

	//printCf();
}

void Patron::printCf(void)
{
	for (int i = 0; i < numTaxa; i++) // taxon
	{
		for (int j = 0; j < numTimepoints; j++) // timepoint
		{
			for (int k = 0; k < numSteps; k++)
			{
				std::cout << std::setprecision(8) << "\ti:" << i << "\tj:" << j << "\tk:" << k;
				std::cout << "\t[" << i * oneTaxonSize + j * oneTimeSize + k*2 << "]";
				std::cout << "\tmu:" << data[i][j] << "\t" << REAL(cfData, i * oneTaxonSize + j * oneTimeSize + k*2);
				std::cout << " " << IMAG(cfData, i * oneTaxonSize + j * oneTimeSize + k*2) << "\t";
				std::cout << &cfData[i * oneTaxonSize + j * oneTimeSize + k*2] << "\n";
			}
		}
	}
}

std::string Patron::getPrintStr(void)
{
	std::string printStr = name + "\t";
	printStr += table->getId() + "\t";

	for (std::vector<Parm*>::const_iterator it_p = table->getParmVector().begin(); it_p != table->getParmVector().end(); it_p++)
	{
		printStr += (*it_p)->getParameterStr();
	}

	return printStr;
}

gsl_complex Patron::charFunc(double T, double m, double s)
{
	gsl_complex val;
	GSL_SET_COMPLEX(&val, (-0.5)*s*s*T*T, m*T);
	return gsl_complex_exp(val);
}

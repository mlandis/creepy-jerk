/*
 * Parm_sigma.cpp
 *
 *  Created on: Jan 26, 2010
 *      Author: mlandis
 */

#include "Parm_sigma.h"

//Var::Var(MbRandom* rp, Tree* tp) {

Var::Var(MbRandom* rp, std::string pn) {

	ranPtr = rp;
	parmName = pn;
	k = ranPtr->exponentialRv(1.0);
	// print();
}

Var::Var(MbRandom* rp, std::string pn, double v)
{
	ranPtr = rp;
	parmName = pn;
	k = v;
	// print();
}

Var::Var(Var &t) {

	clone(t);
}

Var& Var::operator=(const Var &t) {

	if (this != &t)
		clone(t);
	return *this;
}

double Var::change(void) {

	double tuning = log(3.0);
	double oldK = k;
	double newK = oldK * exp(tuning * (ranPtr->uniformRv() - 0.5));
	k = newK;
	return log(newK) - log(oldK);
}

double Var::changeTruncated(double truncateVal)
{

	/*
	double oldK = k;
	double newK = ranPtr->truncatedNormalRv(0, pow(truncateVal,0.5), oldK, 0.5);
	k = newK;
	return ( log(ranPtr->truncatedNormalPdf(0, pow(truncateVal,0.5), oldK, 0.5, newK))
	       - log(ranPtr->truncatedNormalPdf(0, pow(truncateVal,0.5), newK, 0.5, oldK)) );
   */
	double oldK = k;
	double newK = ranPtr->truncatedNormalRv(0, pow(truncateVal,0.5), oldK, 0.5); //.1*truncateVal);
	k = newK;
	return ( log(ranPtr->truncatedNormalPdf(0, pow(truncateVal,0.5), oldK, 0.5, newK))
		   - log(ranPtr->truncatedNormalPdf(0, pow(truncateVal,0.5), newK, 0.5, oldK)) );
}

double Var::change(double newK)
{
	double oldK = k;
	k = newK;
	return log(newK) - log(oldK);
}

double Var::lnProbability(void) {

	return -2.0 * log(1.0 + k);
	//return ranPtr->lnGammaPdf(3.0, (1.0/3.0), k);
}

void Var::print(void) {

	std::cout << parmName << "\t" << std::fixed << std::setprecision(5) << k << std::endl;
}

void Var::clone(const Var &t) {

	ranPtr = t.ranPtr;
	parmName = t.parmName;
//	treePtr = t.treePtr;
	k = t.k;
}

//Sigma::Sigma(MbRandom* rp, Tree* tp, std::string pn) : Parm(rp, pn) {
Sigma::Sigma(MbRandom* rp, std::string pn) : Parm(rp, pn) {

	//vars[0] = new Var(rp, tp);
	vars[0] = new Var(rp, pn);
	vars[1] = new Var(*vars[0]);
}

Sigma::Sigma(MbRandom* rp, std::string pn, double v) : Parm(rp, pn) {

	//vars[0] = new Var(rp, tp);
	vars[0] = new Var(rp, pn, v);
	vars[1] = new Var(*vars[0]);
}

Sigma::Sigma(const Sigma &t, MbRandom* rp, std::string pn) : Parm(rp, pn)
{
	if (this != &t)
	{
		vars[0] = new Var(*(t.getActiveVar()));
		vars[1] = new Var(*(t.getInactiveVar()));
	}
}


Sigma::~Sigma(void) {

	delete vars[0];
	delete vars[1];
}

double Sigma::lnPriorRatio(void) {

	return vars[activeState]->lnProbability() - vars[getInactiveState()]->lnProbability();
}

double Sigma::lnPrior(void)
{
	return vars[activeState]->lnProbability();
}

double Sigma::change(void) {

	numAttemptedChanges++;
	return vars[activeState]->change();
}

double Sigma::getValue(void)
{
	return vars[activeState]->getValue();
}

void Sigma::print(void) {

	vars[activeState]->print();
}

void Sigma::keep(void) {

	*vars[getInactiveState()] = *vars[activeState];
}

void Sigma::restore(void) {

	*vars[activeState] = *vars[getInactiveState()];
}

std::string Sigma::getParameterStr(void) {

	double k = vars[activeState]->getValue();
	char tempCh[50];
	sprintf(tempCh, "%1.6lf\t", k);
	std::string pStr = tempCh;
	return pStr;
}

std::string Sigma::getParameterHeader(void) {

	std::string pStr = parmName + "\t";
	return pStr;
}

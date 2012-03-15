/*
 * Parm_kappa.cpp
 *
 *  Created on: Jan 26, 2010
 *      Author: mlandis
 */

#include "Parm_kappa.h"

//Kurtosis::Kurtosis(MbRandom* rp, Tree* tp) {

Kurtosis::Kurtosis(MbRandom* rp, std::string pn) {

	ranPtr = rp;
	parmName = pn;
	k = ranPtr->exponentialRv(1.0);
	// print();
}

Kurtosis::Kurtosis(MbRandom* rp, std::string pn, double v)
{
	ranPtr = rp;
	parmName = pn;
	k = v;
	// print();
}

Kurtosis::Kurtosis(Kurtosis &t) {

	clone(t);
}

Kurtosis& Kurtosis::operator=(const Kurtosis &t) {

	if (this != &t)
		clone(t);
	return *this;
}

double Kurtosis::change(void) {


	//John's proposal

	double tuning = log(3.0);
	double oldK = k;
	double newK = oldK * exp(tuning * (ranPtr->uniformRv() - 0.5));
	k = newK;
	return log(newK) - log(oldK);


	//Gamma proposal
	/*
	double tuning = 5.0;
	double oldK = k;
	double newK = ranPtr->gammaRv(tuning, tuning/k);
	k = newK;
	return ranPtr->lnGammaPdf(tuning, tuning/k, oldK) - ranPtr->lnGammaPdf(tuning, tuning/oldK, k);
	*/
}

double Kurtosis::changeTruncated(double truncateVal)
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

double Kurtosis::change(double newK)
{
	double oldK = k;
	k = newK;
	return log(newK) - log(oldK);
}

double Kurtosis::lnProbability(void) {

	//return -2.0 * log(1.0 + k);
	double sig = 1;
	return log(2 * sig) - log(sig * sig + k * k);
	//return ranPtr->lnGammaPdf(3.0, (1.0/3.0), k);
}

void Kurtosis::print(void) {

	std::cout << parmName << "\t" << std::fixed << std::setprecision(5) << k << std::endl;
}

void Kurtosis::clone(const Kurtosis &t) {

	ranPtr = t.ranPtr;
	parmName = t.parmName;
//	treePtr = t.treePtr;
	k = t.k;
}

//Kappa::Kappa(MbRandom* rp, Tree* tp, std::string pn) : Parm(rp, pn) {
Kappa::Kappa(MbRandom* rp, std::string pn) : Parm(rp, pn) {

	//kurts[0] = new Kurtosis(rp, tp);
	kurts[0] = new Kurtosis(rp, pn);
	kurts[1] = new Kurtosis(*kurts[0]);
}

Kappa::Kappa(MbRandom* rp, std::string pn, double v) : Parm(rp, pn) {

	//kurts[0] = new Kurtosis(rp, tp);
	kurts[0] = new Kurtosis(rp, pn, v);
	kurts[1] = new Kurtosis(*kurts[0]);
}

Kappa::Kappa(const Kappa &t, MbRandom* rp, std::string pn) : Parm(rp, pn)
{
	if (this != &t)
	{
		kurts[0] = new Kurtosis(*(t.getActiveKurtosis()));
		kurts[1] = new Kurtosis(*(t.getInactiveKurtosis()));
	}
}


Kappa::~Kappa(void) {

	delete kurts[0];
	delete kurts[1];
}

double Kappa::lnPriorRatio(void) {

	return kurts[activeState]->lnProbability() - kurts[getInactiveState()]->lnProbability();
}

double Kappa::lnPrior(void)
{
	return kurts[activeState]->lnProbability();
}

double Kappa::change(void) {

	numAttemptedChanges++;
	return kurts[activeState]->change();
}

double Kappa::getValue(void)
{
	return kurts[activeState]->getValue();
}

void Kappa::print(void) {

	kurts[activeState]->print();
}

void Kappa::keep(void) {

	*kurts[getInactiveState()] = *kurts[activeState];
}

void Kappa::restore(void) {

	*kurts[activeState] = *kurts[getInactiveState()];
}

std::string Kappa::getParameterStr(void) {

	double k = kurts[activeState]->getValue();
	char tempCh[50];
	sprintf(tempCh, "%1.6lf\t", k);
	std::string pStr = tempCh;
	return pStr;
}

std::string Kappa::getParameterHeader(void) {

	std::string pStr = parmName + "\t";
	return pStr;
}

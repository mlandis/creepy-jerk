/*
 * Parm_alpha.cpp
 *
 *  Created on: Jan 26, 2010
 *      Author: mlandis
 */

#include "Parm_alpha.h"

//Stability::Stability(MbRandom* rp, Tree* tp) {

Stability::Stability(MbRandom* rp, std::string pn) {

	ranPtr = rp;
	parmName = pn;
	k = 2 * ranPtr->betaRv(1.0,1.0);
	// print();
}

Stability::Stability(MbRandom* rp, std::string pn, double v)
{
	ranPtr = rp;
	parmName = pn;
	k = v;
	// print();
}

Stability::Stability(Stability &t) {

	clone(t);
}

Stability& Stability::operator=(const Stability &t) {

	if (this != &t)
		clone(t);
	return *this;
}

double Stability::change(void) {


	//John's proposal

	/*
	double tuning = log(3.0);
	double oldK = k;
	double newK = oldK * exp(tuning * (ranPtr->uniformRv() - 0.5));
	k = newK;
	return log(newK) - log(oldK);
*/
	double oldK = k;
	double newK = ranPtr->truncatedNormalRv(0, 2.0, oldK, 0.5); //.1*truncateVal);
	k = newK;
	return ( log(ranPtr->truncatedNormalPdf(0, 2.0, oldK, 0.5, newK))
		   - log(ranPtr->truncatedNormalPdf(0, 2.0, newK, 0.5, oldK)) );

	//Gamma proposal
	/*
	double tuning = 5.0;
	double oldK = k;
	double newK = ranPtr->gammaRv(tuning, tuning/k);
	k = newK;
	return ranPtr->lnGammaPdf(tuning, tuning/k, oldK) - ranPtr->lnGammaPdf(tuning, tuning/oldK, k);
	*/
}

double Stability::changeTruncated(double truncateVal)
{

	/*
	double oldK = k;
	double newK = ranPtr->truncatedNormalRv(0, pow(truncateVal,0.5), oldK, 0.5);
	k = newK;
	return ( log(ranPtr->truncatedNormalPdf(0, pow(truncateVal,0.5), oldK, 0.5, newK))
	       - log(ranPtr->truncatedNormalPdf(0, pow(truncateVal,0.5), newK, 0.5, oldK)) );
   */
	double oldK = k;
	double newK = ranPtr->truncatedNormalRv(0, 2.0, oldK, 0.5); //.1*truncateVal);
	k = newK;
	return ( log(ranPtr->truncatedNormalPdf(0, 2.0, oldK, 0.5, newK))
		   - log(ranPtr->truncatedNormalPdf(0, 2.0, newK, 0.5, oldK)) );
}

double Stability::change(double newK)
{
	double oldK = k;
	k = newK;
	return log(newK) - log(oldK);
}

double Stability::lnProbability(void) {

	return 0.0; // log(0.5);
	//return -2.0 * log(1.0 + k);
	//return ranPtr->lnGammaPdf(3.0, (1.0/3.0), k);
}

void Stability::print(void) {

	std::cout << parmName << "\t" << std::fixed << std::setprecision(5) << k << std::endl;
}

void Stability::clone(const Stability &t) {

	ranPtr = t.ranPtr;
	parmName = t.parmName;
	k = t.k;
}

Alpha::Alpha(MbRandom* rp, std::string pn) : Parm(rp, pn) {

	//stability[0] = new Stability(rp, tp);
	stability[0] = new Stability(rp, pn);
	stability[1] = new Stability(*stability[0]);
}

Alpha::Alpha(MbRandom* rp, std::string pn, double v) : Parm(rp, pn) {

	//stability[0] = new Stability(rp, tp);
	stability[0] = new Stability(rp, pn, v);
	stability[1] = new Stability(*stability[0]);
}

Alpha::Alpha(const Alpha &t, MbRandom* rp, std::string pn) : Parm(rp, pn)
{
	if (this != &t)
	{
		stability[0] = new Stability(*(t.getActiveStability()));
		stability[1] = new Stability(*(t.getInactiveStability()));
	}
}


Alpha::~Alpha(void) {

	delete stability[0];
	delete stability[1];
}

double Alpha::lnPriorRatio(void) {

	return stability[activeState]->lnProbability() - stability[getInactiveState()]->lnProbability();
}

double Alpha::lnPrior(void)
{
	return stability[activeState]->lnProbability();
}

double Alpha::change(void) {

	numAttemptedChanges++;
	return stability[activeState]->change();
}

double Alpha::getValue(void)
{
	return stability[activeState]->getValue();
}

void Alpha::print(void) {

	stability[activeState]->print();
}

void Alpha::keep(void) {

	*stability[getInactiveState()] = *stability[activeState];
}

void Alpha::restore(void) {

	*stability[activeState] = *stability[getInactiveState()];
}

std::string Alpha::getParameterStr(void) {

	double k = stability[activeState]->getValue();
	char tempCh[50];
	sprintf(tempCh, "%1.6lf\t", k);
	std::string pStr = tempCh;
	return pStr;
}

std::string Alpha::getParameterHeader(void) {

	std::string pStr = parmName + "\t";
	return pStr;
}

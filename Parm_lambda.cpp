/*
 * Parm_lambda.cpp
 *
 *  Created on: Jan 26, 2010
 *      Author: mlandis
 */

#include "Parm_lambda.h"

//ExpMean::ExpMean(MbRandom* rp, Tree* tp) {

ExpMean::ExpMean(MbRandom* rp) {

	ranPtr = rp;
	//treePtr = tp;
	k = ranPtr->normalRv(0.0, 4.0);
	// print();
	tuning = 0.5;
}

ExpMean::ExpMean(ExpMean &t) {

	clone(t);
}

ExpMean& ExpMean::operator=(const ExpMean &t) {

	if (this != &t)
		clone(t);
	return *this;
}

double ExpMean::change(void) {

	// update conditional likelihood and transition probability flags

	// TODO: reinstate later

	// Topology* t = treePtr->getActiveTopology();
	// t->updateAllCls(true);
	// t->updateAllTis(true);
	// t->flipAllActiveCls();
	// t->flipAllActiveTis();


	double invertA = 0.05;
	if (ranPtr->uniformRv() < invertA)
	{
		/*
		//std::cerr << "LAMBDA INVERSION\n";
		k = -k;
		return 0;
		//return log(newK) - log(oldK);

		 */

		double oldK = k;
		// double newK = oldK * exp(tuning * (ranPtr->uniformRv() - 0.5));
		double newK = -k;
		k = newK;
		double val1 = ranPtr->normalPdf(oldK, tuning, newK);
		double val2 = ranPtr->normalPdf(newK, tuning, oldK);
		return log(val2) - log(val1);

	}
	else
	{
		//std::cerr << "LAMBDA CHANGE\n";
		// change the transition/transversion rate ratio
		// double tuning = log(5.0);
		double oldK = k;
		// double newK = oldK * exp(tuning * (ranPtr->uniformRv() - 0.5));
		double newK = k + ranPtr->normalRv(0, tuning);
		k = newK;
		double val1 = ranPtr->normalPdf(oldK, tuning, newK);
		double val2 = ranPtr->normalPdf(newK, tuning, oldK);
		return log(val2) - log(val1);
	}
}

double ExpMean::change(double newK)
{
	double oldK = k;
	k = newK;
	double val1 = ranPtr->normalPdf(oldK, tuning, newK);
	double val2 = ranPtr->normalPdf(newK, tuning, oldK);
	return log(val2) - log(val1);
}

double ExpMean::changeSign(void) {
	k = -k;
	return 0;
}

double ExpMean::lnProbability(void) {

	// return -2.0 * log(1.0 + k);
	//return log(ranPtr->normalPdf(0.0, tuning, k));
	return ranPtr->lnNormalPdf(0.0, 10.0, k);
}

void ExpMean::print(void) {

	std::cout << "Lambda = " << std::fixed << std::setprecision(5) << k << std::endl;
}

void ExpMean::clone(const ExpMean &t) {

	ranPtr = t.ranPtr;
//	treePtr = t.treePtr;
	k = t.k;
	tuning = t.tuning;
}

//Lambda::Lambda(MbRandom* rp, Tree* tp, std::string pn) : Parm(rp, pn) {
Lambda::Lambda(MbRandom* rp, std::string pn) : Parm(rp, pn) {

	//expmeans[0] = new ExpMean(rp, tp);
	expmeans[0] = new ExpMean(rp);
	expmeans[1] = new ExpMean(*expmeans[0]);
}

Lambda::Lambda(const Lambda &t, MbRandom* rp, std::string pn) : Parm(rp, pn)
{
	if (this != &t)
	{
		expmeans[0] = new ExpMean(*(t.getActiveExpMean()));
		expmeans[1] = new ExpMean(*(t.getInactiveExpMean()));
	}
}

Lambda::~Lambda(void) {

	delete expmeans[0];
	delete expmeans[1];
}

double Lambda::lnPriorRatio(void) {

	/*
	if (expmeans[activeState]->getRate() == -expmeans[getInactiveState()]->getRate())
	{
		std::cerr << "NOTE: lambda inverted\n";
		std::cerr << expmeans[activeState]->lnProbability() << "\n";
		std::cerr << expmeans[getInactiveState()]->lnProbability() << "\n";
	}
	*/
	return expmeans[activeState]->lnProbability() - expmeans[getInactiveState()]->lnProbability();
}

double Lambda::lnPrior(void)
{
	return expmeans[activeState]->lnProbability();
}

double Lambda::change(void) {

	numAttemptedChanges++;
	return expmeans[activeState]->change();
}

double Lambda::getValue(void)
{
	return expmeans[activeState]->getValue();
}

void Lambda::print(void) {

	expmeans[activeState]->print();
}

void Lambda::keep(void) {

	*expmeans[getInactiveState()] = *expmeans[activeState];
}

void Lambda::restore(void) {

	*expmeans[activeState] = *expmeans[getInactiveState()];
}

std::string Lambda::getParameterStr(void) {

	double k = expmeans[activeState]->getValue();
	char tempCh[50];
	sprintf(tempCh, "%1.6lf\t", k);
	std::string pStr = tempCh;
	return pStr;
}

std::string Lambda::getParameterHeader(void) {

	std::string pStr = parmName + "\t";
	return pStr;
}

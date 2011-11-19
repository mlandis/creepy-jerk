/*
 * Parm_alpha.cpp
 *
 *  Created on: Jan 26, 2010
 *      Author: mlandis
 */

#include "Parm_alpha.h"

//Conc::Conc(MbRandom* rp, Tree* tp) {

Conc::Conc(MbRandom* rp, Model* mp, double a_alpha, double b_alpha) {

	ranPtr = rp;
	modelPtr = mp;
	alpha = 1.0;
	a = a_alpha;
	b = b_alpha;
	x = 1.0;

}

Conc::Conc(Conc &t) {

	clone(t);
}

Conc& Conc::operator=(const Conc &t) {

	if (this != &t)
		clone(t);
	return *this;
}

double Conc::change(void) {

	// update conditional likelihood and transition probability flags

	/*
	// TODO: reinstate later

	// Topology* t = treePtr->getActiveTopology();
	// t->updateAllCls(true);
	// t->updateAllTis(true);
	// t->flipAllActiveCls();
	// t->flipAllActiveTis();

	double invertA = 0.05;
	if (ranPtr->uniformRv() < invertA)
	{
		//std::cerr << "LAMBDA INVERSION\n";
		k = -k;
		return 0;
		//return log(newK) - log(oldK);
	}
	else
	{
		//std::cerr << "LAMBDA CHANGE\n";
		// change the transition/transversion rate ratio
		// double tuning = log(5.0);
		double oldK = k;
		// double newK = oldK * exp(tuning * (ranPtr->uniformRv() - 0.5));
		double newK = k + ranPtr->normalRv(0, 10.0);
		k = newK;
		double val1 = ranPtr->normalPdf(oldK, 10.0, newK);
		double val2 = ranPtr->normalPdf(newK, 10.0, oldK);
		return log(val2) - log(val1);
	}

	*/

	k = modelPtr->getTableListPtr()->size();
	n = modelPtr->getPatronListPtr()->size();
	double u = ranPtr->uniformRv();
	double w = (a + k - 1) / (n * (b - log(x)) + a + k - 1);
	x = ranPtr->betaRv(alpha + 1, n);
	if (u < w)
	{
		alpha = ranPtr->gammaRv(a + k, b - log(x));
	}
	else
	{
		alpha = ranPtr->gammaRv(a + k - 1, b - log(x));
	}

	// Gibbs sampling requires all proposals are accepted
	return 0.0;
}

double Conc::lnProbability(void) {

	// return -2.0 * log(1.0 + k);
	return log(ranPtr->normalPdf(0.0, 10.0, k));
}

void Conc::print(void) {

	std::cout << "Alpha = " << std::fixed << std::setprecision(5) << alpha << " (k=" << k << ",n=" << n << ")" << std::endl;
}

void Conc::clone(const Conc &t) {

	ranPtr = t.ranPtr;
	modelPtr = t.modelPtr;
	alpha = t.alpha;
	a = t.a;
	b = t.b;
	k = t.k;
	n = t.n;
}

Alpha::Alpha(MbRandom* rp, Model* mp, double a, double b, std::string pn) : Parm(rp, pn) {

	//concs[0] = new Conc(rp, tp);
	concs[0] = new Conc(rp, mp, a, b);
	concs[1] = new Conc(*concs[0]);
}

/*
Alpha::Alpha(const Alpha &t, MbRandom* rp, std::string pn) : Parm(rp, pn)
{
	if (this != &t)
	{
		concs[0] = new Conc(*(t.getActiveConc()));
		concs[1] = new Conc(*(t.getInactiveConc()));
	}
}
*/

Alpha::~Alpha(void) {

	delete concs[0];
	delete concs[1];
}

double Alpha::lnPriorRatio(void) {

	/*
	if (concs[activeState]->getRate() == -concs[activeState]->getRate())
	{
		std::cerr << "NOTE: alpha inverted\n";
		std::cerr << concs[activeState]->lnProbability() << "\n";
		std::cerr << concs[getInactiveState()]->lnProbability() << "\n";
	}
	return concs[activeState]->lnProbability() - concs[getInactiveState()]->lnProbability();
	*/

	// Alpha is updated by Gibbs sampling
	return 0.0;
}

double Alpha::change(void) {

	numAttemptedChanges++;
	return concs[activeState]->change();
}

double Alpha::getValue(void)
{
	return concs[activeState]->getValue();
}

void Alpha::print(void) {

	concs[activeState]->print();
}

void Alpha::keep(void) {

	*concs[getInactiveState()] = *concs[activeState];
}

void Alpha::restore(void) {

	*concs[activeState] = *concs[getInactiveState()];
}

std::string Alpha::getParameterStr(void) {

	double k = concs[activeState]->getValue();
	char tempCh[50];
	sprintf(tempCh, "%1.3lf\t", k);
	std::string pStr = tempCh;
	return pStr;
}

std::string Alpha::getParameterHeader(void) {

	std::string pStr = "Alpha\t";
	return pStr;
}

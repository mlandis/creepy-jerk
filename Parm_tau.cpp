/*
 * Parm_tau.cpp
 *
 *  Created on: Jan 26, 2010
 *      Author: mlandis
 */

#include "Parm_tau.h"

//BrLen::BrLen(MbRandom* rp, Tree* tp) {

BrLen::BrLen(MbRandom* rp, std::string pn) {

	ranPtr = rp;
	//treePtr = tp;
	k = ranPtr->exponentialRv(1.0);
	parmName = pn;
	// print();
}

BrLen::BrLen(MbRandom* rp, std::string pn, double v)
{
	ranPtr = rp;
	k = v;
	parmName = pn;
}

BrLen::BrLen(BrLen &t) {

	clone(t);
}

BrLen& BrLen::operator=(const BrLen &t) {

	if (this != &t)
		clone(t);
	return *this;
}

double BrLen::change(void) {

	// update conditional likelihood and transition probability flags

	// TODO: reinstate later

	// Topology* t = treePtr->getActiveTopology();
	// t->updateAllCls(true);
	// t->updateAllTis(true);
	// t->flipAllActiveCls();
	// t->flipAllActiveTis();


	double tuning = log(3.0);
	double oldK = k;
	double newK = oldK * exp(tuning * (ranPtr->uniformRv() - 0.5));
	k = newK;
	return log(newK) - log(oldK);
}

double BrLen::change(double r)
{
	double oldK = k;
	k = r;
	return log(k) - log(oldK);
}

double BrLen::lnProbability(void) {

	return -2.0 * log(1.0 + k);
	//return ranPtr->lnGammaPdf(3.0, 1.0/3.0 ,k);
}

void BrLen::print(void) {

	std::cout << parmName << " " << std::fixed << std::setprecision(5) << k << std::endl;
}

void BrLen::clone(const BrLen &t) {

	ranPtr = t.ranPtr;
	parmName = t.parmName;
//	treePtr = t.treePtr;
	k = t.k;
}

//Tau::Tau(MbRandom* rp, Tree* tp, std::string pn) : Parm(rp, pn) {
Tau::Tau(MbRandom* rp, std::string pn, double v) : Parm(rp, pn) {

	//brlens[0] = new BrLen(rp, tp);
	brlens[0] = new BrLen(rp, pn, v);
	brlens[1] = new BrLen(*brlens[0]);
}

Tau::Tau(MbRandom* rp, std::string pn) : Parm(rp, pn) {

	//brlens[0] = new BrLen(rp, tp);
	brlens[0] = new BrLen(rp, pn);
	brlens[1] = new BrLen(*brlens[0]);
}

Tau::Tau(const Tau &t, MbRandom* rp, std::string pn) : Parm(rp, pn)
{
	if (this != &t)
	{
		brlens[0] = new BrLen(*(t.getActiveBrLen()));
		brlens[1] = new BrLen(*(t.getInactiveBrLen()));
	}
}

Tau::~Tau(void) {

	delete brlens[0];
	delete brlens[1];
}

double Tau::lnPriorRatio(void) {

	return brlens[activeState]->lnProbability() - brlens[getInactiveState()]->lnProbability();
}

double Tau::change(void) {

	numAttemptedChanges++;
	return brlens[activeState]->change();
}

double Tau::getValue(void)
{
	return brlens[activeState]->getValue();
}

void Tau::print(void) {

	brlens[activeState]->print();
}

void Tau::keep(void) {

	*brlens[getInactiveState()] = *brlens[activeState];
}

void Tau::restore(void) {

	*brlens[activeState] = *brlens[getInactiveState()];
}

std::string Tau::getParameterStr(void) {

	double k = brlens[activeState]->getValue();
	char tempCh[50];
	sprintf(tempCh, "%1.3lf\t", k);
	std::string pStr = tempCh;
	return pStr;
}

std::string Tau::getParameterHeader(void) {

	std::string pStr = parmName + "\t";
	return pStr;
}

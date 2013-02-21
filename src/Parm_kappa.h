/*
 * Parm_sigma.h
 *
 *  Created on: Jan 26, 2010
 *      Author: mlandis
 */

// This parameter is the kurtosis (e.g. VG).

#ifndef PARM_KAPPA_H_
#define PARM_KAPPA_H_

#include "MbRandom.h"
#include "Parm.h"

#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <string>

class MbRandom;

class Kurtosis {

public:
	            Kurtosis(MbRandom* rp, std::string pn);
	            Kurtosis(MbRandom* rp, std::string pn, double v);
				Kurtosis(Kurtosis &t);
	Kurtosis	&operator=(const Kurtosis &k);
	double		change(void);
	double		change(double);
	double		changeTruncated(double);
	double		getValue(void)			{ return k; }
	double		lnProbability(void);
	void		print(void);

private:
	void		clone(const Kurtosis &t);
	MbRandom*	ranPtr;
	double		k;
	std::string parmName;

};


class Kappa : public Parm {

public:
	            Kappa(MbRandom* rp, std::string pn);
	            Kappa(MbRandom* rp, std::string pn, double v);
	            Kappa(const Kappa &t, MbRandom *rp, std::string pn);
				~Kappa(void);
	Kurtosis*	getActiveKurtosis(void) const		{return kurts[activeState]; }
	Kurtosis*	getInactiveKurtosis(void) const		{return kurts[getInactiveState()]; }
	double		lnPriorRatio(void);
	double		lnPrior(void);
	double		change(void);
	double		getValue(void);
	void		print(void);
	void		keep(void);
	void		restore(void);
	std::string	getParameterStr(void);
	std::string	getParameterHeader(void);

private:
	Kurtosis		*kurts[2];

};

#endif PARM_KAPPA_H_

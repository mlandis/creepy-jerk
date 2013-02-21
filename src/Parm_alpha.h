/*
 * Parm_sigma.h
 *
 *  Created on: Jan 26, 2010
 *      Author: mlandis
 */

// This parameter is the variance.

#ifndef PARM_ALPHA_H_
#define PARM_ALPHA_H_

#include "MbRandom.h"
#include "Parm.h"

#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <string>

class MbRandom;

class Stability {

public:
	            Stability(MbRandom* rp, std::string pn);
	            Stability(MbRandom* rp, std::string pn, double v);
				Stability(Stability &t);
	Stability			&operator=(const Stability &k);
	double		change(void);
	double		change(double);
	double		changeTruncated(double);
	double		getValue(void)			{ return k; }
	double		lnProbability(void);
	void		print(void);

private:
	void		clone(const Stability &t);
	MbRandom*	ranPtr;
	double		k;
	std::string parmName;

};


class Alpha : public Parm {

public:
	            Alpha(MbRandom* rp, std::string pn);
	            Alpha(MbRandom* rp, std::string pn, double v);
	            Alpha(const Alpha &t, MbRandom *rp, std::string pn);
				~Alpha(void);
	Stability*		getActiveStability(void) const		{return stability[activeState]; }
	Stability*		getInactiveStability(void) const		{return stability[getInactiveState()]; }
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
	Stability		*stability[2];

};

#endif PARM_ALPHA_H_

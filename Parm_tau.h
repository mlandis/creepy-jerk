/*
 * Parm_tau.h
 *
 *  Created on: Jan 26, 2010
 *      Author: mlandis
 */

// TAU is the branch length parameter

#ifndef PARM_TAU_H_
#define PARM_TAU_H_

#include "MbRandom.h"
#include "Parm.h"

#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <string>

class MbRandom;

class BrLen {

public:
				//BrLen(MbRandom* rp, Tree* tp);
	            BrLen(MbRandom* rp, std::string pn);
	            BrLen(MbRandom* rp, std::string pn, double v);
				BrLen(BrLen &t);
	BrLen		&operator=(const BrLen &k);
	double		change(void);
	double		change(double);
	double		getValue(void)			{ return k; }
	double		lnProbability(void);
	void		print(void);

private:
	void		clone(const BrLen &t);
	MbRandom*	ranPtr;
	double		k;
	std::string parmName;

};


class Tau : public Parm {

public:
	            Tau(MbRandom* rp, std::string pn);
	            Tau(MbRandom* rp, std::string pn, double v);
	            Tau(const Tau &t, MbRandom* rp, std::string pn);
				~Tau(void);
	BrLen*		getActiveBrLen(void) const		{return brlens[activeState]; }
	BrLen*		getInactiveBrLen(void) const 	{return brlens[getInactiveState()]; }
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
	BrLen		*brlens[2];

};

#endif /* PARM_TAU_H_ */

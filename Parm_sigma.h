/*
 * Parm_sigma.h
 *
 *  Created on: Jan 26, 2010
 *      Author: mlandis
 */

// This parameter is the variance.

#ifndef PARM_SIGMA_H_
#define PARM_SIGMA_H_

#include "MbRandom.h"
#include "Parm.h"

#include <iomanip>
#include <iostream>
#include <string>

class MbRandom;

class Var {

public:
	            Var(MbRandom* rp, std::string pn);
	            Var(MbRandom* rp, std::string pn, double v);
				Var(Var &t);
	Var			&operator=(const Var &k);
	double		change(void);
	double		change(double);
	double		getValue(void)			{ return k; }
	double		lnProbability(void);
	void		print(void);

private:
	void		clone(const Var &t);
	MbRandom*	ranPtr;
	double		k;
	std::string parmName;

};


class Sigma : public Parm {

public:
	            Sigma(MbRandom* rp, std::string pn);
	            Sigma(MbRandom* rp, std::string pn, double v);
	            Sigma(const Sigma &t, MbRandom *rp, std::string pn);
				~Sigma(void);
	Var*		getActiveVar(void) const		{return vars[activeState]; }
	Var*		getInactiveVar(void) const		{return vars[getInactiveState()]; }
	double		lnPriorRatio(void);
	double		change(void);
	double		getValue(void);
	void		print(void);
	void		keep(void);
	void		restore(void);
	std::string	getParameterStr(void);
	std::string	getParameterHeader(void);

private:
	Var		*vars[2];

};

#endif PARM_SIGMA_H_

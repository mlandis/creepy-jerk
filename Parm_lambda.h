/*
 * Parm_lambda.h
 *
 *  Created on: Jan 26, 2010
 *      Author: mlandis
 */

// Exponential mean

#ifndef PARM_LAMBDA_H_
#define PARM_LAMBDA_H_

#include "MbRandom.h"
#include "Parm.h"

#include <iomanip>
#include <iostream>
#include <string>

class MbRandom;

class ExpMean {

public:
				//ExpMean(MbRandom* rp, Tree* tp);
	            ExpMean(MbRandom* rp);
				ExpMean(ExpMean &t);
	ExpMean		&operator=(const ExpMean &k);
	double		change(void);
	double		change(double);
	double		changeSign(void);
	double		getValue(void)			{ return k; }
	double		lnProbability(void);
	void		print(void);

private:
	void		clone(const ExpMean &t);
	MbRandom*	ranPtr;
	double		k;
	double		tuning;

};


class Lambda : public Parm {

public:
			//	Lambda(MbRandom* rp, Tree* tp, std::string pn);
	            Lambda(MbRandom* rp, std::string pn);
	            Lambda(const Lambda &t, MbRandom* rp, std::string pn);
				~Lambda(void);
	ExpMean*	getActiveExpMean(void) const	{return expmeans[activeState]; }
	ExpMean*	getInactiveExpMean(void) const	{return expmeans[getInactiveState()]; }
	double		lnPriorRatio(void);
	double		change(void);
	double		getValue(void);
	void		print(void);
	void		keep(void);
	void		restore(void);
	std::string	getParameterStr(void);
	std::string	getParameterHeader(void);

private:
	ExpMean		*expmeans[2];

};

#endif /* PARM_LAMBDA_H_ */

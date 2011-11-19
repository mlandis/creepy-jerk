/*
 * Parm_lambda.h
 *
 *  Created on: Jan 26, 2010
 *      Author: mlandis
 */

// CRP conc parameter

#ifndef PARM_ALPHA_H_
#define PARM_ALPHA_H_

#include "MbRandom.h"
#include "Model.h"
#include "Parm.h"

#include <iomanip>
#include <iostream>
#include <string>

class MbRandom;
class Model;

class Conc {

public:
	            Conc(MbRandom* rp, Model* mp, double a, double b);
				Conc(Conc &t);
	Conc		&operator=(const Conc &k);
	void		setK(void);
	double		change(void);
	double		getValue(void)			{ return alpha; }
	double		lnProbability(void);
	void		print(void);

private:
	void		clone(const Conc &t);
	MbRandom*	ranPtr;
	Model*		modelPtr;
	double		alpha;
	double		a;
	double		b;
	double		x;
	int			k;
	int			n;
};


class Alpha : public Parm {

public:
	            Alpha(MbRandom* rp, Model* mp, double a, double b, std::string pn);
				~Alpha(void);
	Conc* getActiveConc(void) const { return concs[activeState]; }
	Conc* getInactiveConc(void) const	{ return concs[getInactiveState()]; }
	double		lnPriorRatio(void);
	double		change(void);
	double		getValue(void);
	void		print(void);
	void		keep(void);
	void		restore(void);
	std::string	getParameterStr(void);
	std::string	getParameterHeader(void);

private:
	Conc *concs[2];

};

#endif /* PARM_ALPHA_H_ */

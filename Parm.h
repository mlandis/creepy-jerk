/*
 * Parm.h
 *
 *  Created on: Jan 24, 2010
 *
 */

#ifndef PARM_H_
#define PARM_H_

#include <iostream>

class MbRandom;

class Parm {

public:
						Parm(MbRandom *rp, std::string pn);
	std::string			getName(void)								{ return parmName; }
	virtual double		lnPriorRatio(void)= 0;
	virtual double		lnPrior(void) = 0;
	virtual double		change(void) = 0;
	virtual void		print(void) = 0;
	virtual void		keep(void) = 0;
	virtual void		restore(void) = 0;
	void				incrementNumAccepted(void)					{ numAcceptedChanges++; }
	int					getActiveState(void) const;
	int					getInactiveState(void) const;
	int					getNumAcceptances(void)						{ return numAcceptedChanges; }
	int					getNumAttempts(void)						{ return numAttemptedChanges; }
	virtual double		getValue(void) = 0;
	virtual std::string	getParameterStr(void) = 0;
	virtual std::string	getParameterHeader(void) = 0;
	MbRandom*			getRandomPtr()								{ return ranPtr; }

protected:
	MbRandom*			ranPtr;
	std::string			parmName;
	int					activeState;
	int					numAttemptedChanges;
	int					numAcceptedChanges;

};

#endif /* PARM_H_ */

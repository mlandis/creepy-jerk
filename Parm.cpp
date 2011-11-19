/*
 * Parm.cpp
 *
 *  Created on: Jan 24, 2010
 *      Author: mlandis
 */

#include "Parm.h"

Parm::Parm(MbRandom* rp, std::string pn) {

	ranPtr = rp;
	parmName = pn;
	activeState = 0;
	numAttemptedChanges = 0;
	numAcceptedChanges = 0;
}

int Parm::getActiveState(void) const {

	if (activeState == 1)	return 1;
	else					return 0;
}

int Parm::getInactiveState(void) const {

	if (activeState == 1)	return 0;
	else					return 1;
}

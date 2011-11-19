/*
 * CondLikes.h
 *
 *  Created on: Mar 30, 2011
 *      Author: mlandis
 */

#ifndef CONDLIKES_H_
#define CONDLIKES_H_

#include "Expression.h"
#include "Patron.h"
#include "Settings.h"


#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <list>
#include <vector>

class Expression;
class Settings;

class CondLikes {

public:
				CondLikes(Expression* ep, Settings* sp);
				~CondLikes(void);

	double*		getClPtr(int space, int trans, int node, int time)			{ return clsPtr[space][trans][node][time]; }
	double*		getClScalerPtr(int space, int trans, int node, int time)	{ return clScalerPtr[space][trans][node][time]; }
	double*		getClPtrUpL(int space, int trans, int node, int time)		{ return clsPtrUpL[space][trans][node][time]; }
	double*		getClScalerPtrUpL(int space, int trans, int node, int time)	{ return clScalerPtrUpL[space][trans][node][time]; }
	double*		getClPtrUpR(int space, int trans, int node, int time)		{ return clsPtrUpR[space][trans][node][time]; }
	double*		getClScalerPtrUpR(int space, int trans, int node, int time)	{ return clScalerPtrUpR[space][trans][node][time]; }
	void		print(void);

	//double*						getCls(int space, int node, int trans, int time, int step)	{ return &cls[space][node][trans][time][step]; }
	//gsl_complex*					getCls(int space, int node, int trans, int time, int step)	{ return &cls[space][node][trans][time][2*step]; }
	//std::vector<double>*			getCls(int space, int node, int trans, int time)			{ return &cls[space][node][trans][time]; }
	//const std::vector<std::vector<std::vector<double> > >&	getCls(int space, int node, int trans, int time);
	//const std::vector<std::vector<std::vector<std::vector<double> > > >&	getCls(int space, int node, int trans, int time, int step);
	//void							setCls(int space, int node, in)

private:

	Expression* expressionPtr;
	Settings* settingsPtr;

	int numNodes;

	int numSpaces;
	int numSteps;
	int numTaxa;
	int numTimepoints;
	int numTranscripts;

	// [space][transcript][node][timepoint][step]
	double*		cls;
	double*****	clsPtr;
	double*		clScalers;
	double*****	clScalerPtr;


	double*		clsUpL;
	double*****	clsPtrUpL;
	double*		clScalersUpL;
	double*****	clScalerPtrUpL;

	double*		clsUpR;
	double*****	clsPtrUpR;
	double*		clScalersUpR;
	double*****	clScalerPtrUpR;



	//std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > > cls;

	//double***** clP;
	//double***** clR;
	//double***** clL;

	//double* clPtrP;
	//double* clPtrL;
	//double* clPtrR;

};

#endif /* CONDLIKES_H_ */

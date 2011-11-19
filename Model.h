/*
 * Model.h
 *
 *  Created on: Mar 7, 2011
 *      Author: mlandis
 */


/*
 * FFT
 * CRP
 *
 *
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <iostream>
#include <iomanip>
#include <list>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include "CondLikes.h"
#include "Expression.h"
#include "MbRandom.h"
#include "Parm.h"
#include "Parm_alpha.h"
#include "Parm_lambda.h"
#include "Parm_sigma.h"
#include "Parm_tau.h"
#include "Parm_tree.h"
#include "Patron.h"
#include "Settings.h"
#include "Table.h"

class Alpha;
class CondLikes;
class Conc;
class Expression;
class MbRandom;
class Settings;
class Topology;


class Model {

public:
	Model(Expression* ep, MbRandom* rp, Settings* sp, Topology* tp);
	virtual ~Model();
	void								updateModel(void);
	void								updateAllFlags(void);
	void								updateAllFlags(int);
	BrLen*								getActiveBrLen(void);
	Var*								getActiveVar(void);
	ExpMean*							getActiveExpMean(void);
	Conc*								getActiveConc(void);
	Alpha*								getAlpha(void)						{ return alpha; }
	int									getTableId(void)					{ return ++tableId; }
	Parm*								pickParmAtRandom(Table* t);
	Parm*								pickBranchAtRandom(Table* t);
	Node*								pickNodeByBrLen(Table* t);
	Table*								pickTableAtRandom(void);
	void								keep(Parm* p);
	void								restore(Parm* p);
	void								printTables(void);
	std::list<Table*>*					getTableListPtr(void)				{ return &tableList; }
	std::list<Patron*>*					getPatronListPtr(void)				{ return &patronList; }
	double								modelLogLikelihood(void);
	double								locusLogLikelihood(Patron*);
	double								tableLogLikelihood(Table*);
	double								tableLogLikelihoodFFT(Table*);
	// jump sampling
	void								sampleJumpsForTree(void);
	void								sampleJumpsForBranch(int, double, double);
	void								sampleJumpsForBranchGivenLambda(int, double);
	void								sampleJumpsForBranchGivenSigma(int, double);
	void								sampleJumpSizesForBranch(int j, double sig_jn);
	void								sampleJumpSizesForTree(void);
	double								getProposalRatioLambda(double oldLambda, double newLambda);
	double								getProposalRatioSigma(double oldSigma, double newSigma);

private:
	void								setCharFunc(Table* t);
	double								scaleIFFT(double, double, int);
	gsl_complex							complexPdf(double, double, double, double); // compound Poisson process w/ skewed Normal
	//  gsl_complex							charFunc(double, double, double);
	gsl_complex							complexPdf(double, double, double, double, double); // jump-diffusion process w/ Normal (non-BM)
	gsl_complex							charFunc(double, double, const std::vector<double>&);
	double								trapInt(double* fn);

	// void								initializeCondLikes(void);
	void								initializeFFT(void);
	void								initializeParms(void);
	void								initializeModelType(void);
	void								initializeTips(void);

	Expression* expressionPtr;
	MbRandom* randomPtr;
	Settings* settingsPtr;
	Topology* topologyPtr;

	std::list<Table*> tableList;
	std::list<Patron*> patronList;
	int tableId;

	int numBranches;
	int numNodes;
	int numTaxa;
	int numTranscripts;
	int numTimepoints;

	// NOTE: will want clever pointing implementation for larger dataset.
	double* like; // NOTE: fixed number of nodes; single var. with additional dimension when topology introduced
	double* tip1;
	double* tip2;
	double* likeVals;

	// Pointer array indexing: [space=2][numNodes][numTrans][numTime][2*numSteps]
	//double***** clsPtr[2];
	//double***** clsPtrUpL[2];
	//double***** clsPtrUpR[2];

	bool useFFT;
	int numSteps;
	int halfSteps;
	double startStep;
	double finalStep;
	double stepSize;
	double ReStepSize;
	double* theta;

	CondLikes* condLikes;

	bool useCRP;
	Alpha* alpha;

	bool fixBranches;
	int modelType;

	std::vector<double> proposalProbs;
};

#endif /* MODEL_H_ */

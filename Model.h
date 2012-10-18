/*
 * Model.h
 *
 *  Created on: Mar 7, 2011
 *      Author: mlandis
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
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>

#include "Expression.h"
#include "MbRandom.h"
#include "Parm.h"
#include "Parm_alpha.h"
#include "Parm_kappa.h"
#include "Parm_sigma.h"
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


// characteristic function definitions for GSL
double poissonNormal(double k, void* params);
double alphaStable(double k, void* params);
double halfSkewNormal(double k, void* params);
double varianceGamma(double k, void* params);
double doubleExponential(double k, void* params);


// class definitions
class Model {

public:
	Model(Expression* ep, MbRandom* rp, Settings* sp, Topology* tp);
	virtual ~Model();
	Parm*								pickParmAtRandom(Table* t);
	Parm*								pickBranchAtRandom(Table* t);
	Node*								pickNodeByBrLen(void);
	Table*								pickTableAtRandom(void);
	void								keep(Parm* p);
	void								restore(Parm* p);
	void								printTables(void);
	std::list<Table*>*					getTableListPtr(void)				{ return &tableList; }
	std::list<Patron*>*					getPatronListPtr(void)				{ return &patronList; }
	double								modelLnLikelihood(int space);
	double								tableLnLikelihood(Table*, int space);
	double								jumpLnLikelihood(Node* p, const std::vector<Parm*>&, int space);
	double								driftLnLikelihood(Node* p, const std::vector<Parm*>&, int space);
	double								besselLnLikelihood(Node*, const std::vector<Parm*>&, int space);
	double								modelDriftOnlyLnLikelihood(int space);
	double								tableDriftOnlyLnLikelihood(Table*, int space);
	double								driftOnlyLnLikelihood(Node* p, const std::vector<Parm*>&, int space);
	double 								proposeJumpSize(Node* p, const std::vector<Parm*>&, int space);
	void								updateGSL(void);



private:

	void								initializeParms(void);
	void								initializeModelType(void);
	void								initializeTips(void);
	void								initializeGSL(void);
	void								freeGSL(void);

	Expression* expressionPtr;
	MbRandom* randomPtr;
	Settings* settingsPtr;
	Topology* topologyPtr;

	std::list<Table*> tableList;
	std::list<Patron*> patronList;
	int tableId;

	// model info
	int numBranches;
	int numNodes;
	int numTaxa;
	int numTranscripts;
	int numTimepoints;

	// model settings
	int modelType;
	double tuningBM;
	bool useJumpKernel;
	int evalType;
	double sigmaJumpProposal;

	// Stepping stone marginal likelihood estimation
	bool useSteppingStone;
	double betaSteppingStone;

	// GSL numerical integration
	gsl_function F;
	gsl_integration_workspace* storeWorkspace;
	gsl_integration_workspace* cycleWorkspace;
	gsl_integration_qawo_table* trigTable;
	int workspaceSize;
	int trigTableSize;
	double integralLength;
	double integralError;

	// MH proposal probabilities
	std::vector<double> proposalProbs;

	// debugging tools
	bool fixBranches;
	bool printStdOut;
};

#endif /* MODEL_H_ */

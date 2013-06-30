/*
 * Mcmc.h
 *
 *  Created on: Dec 16, 2010
 *
 */

#ifndef MCMC_H_
#define MCMC_H_

#include "Expression.h"
//#include "FileHandler.h"
#include "MbRandom.h"
#include "Model.h"
#include "Patron.h"
#include "Settings.h"
#include "Table.h"
#include "Parm_tree.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

class Expression;
class MbRandom;
class Model;
class Patron;
class Settings;
class Table;
class Topology;

class Mcmc
{

public:
	Mcmc(Expression*, Model*, Settings*,  MbRandom*, Topology*);
	virtual ~Mcmc();

private:
	void		runChain(void);
	void		proposeState(void);
	Table*		reseatPatron(void);
	void		proposeParm(Parm*);
	void		proposeChangeJump(void);
	void		proposeRidgeMove(void);

	void		openFiles(std::string outputFilePath);
	void		printChainState(int n, double v);
	void		printChainJumps(int n);
    void        printChainJumpsAsNhx(int n);
	void		printAcceptanceInfo(void);
	std::string printBool(bool);
	std::string printInt(int);
	std::string printDouble(double);
	double		safeExp(double lnX);


	Expression	*expressionPtr;
	MbRandom	*randomPtr;
	Model		*modelPtr;
	Settings	*settingsPtr;
	Topology	*topologyPtr;
	std::list<Table*>* tableListPtr;
	std::list<Patron*>* patronListPtr;

	int			numCycles;
	int			printFreqMH;
	int			printFreqCRP;
	int			printFreqJump;
	int			printFreqStdOut;
    int         snrBurnIn;

	int			numTranscripts;
	int			numTimepoints;
	int			numTaxa;
	int			numNodes;
	int			numBranches;
	int			numParms;

	bool		useSteppingStone;
	double		betaSteppingStone;

	double		alphaCRP;
	double		auxCRP;
	bool		useCRP;

	int			modelType;
	bool		fixBranches;
	bool		useJumpKernel;

	std::vector<double> proposalProbs;
	double		oldLnL, newLnL;
	double		oldKb, newKb;	// lnL from BM
	double		oldKj, newKj;	// lnL from jump kernel
	std::string proposeStr;

	int			patronIndex;

	std::ofstream parmFileStrm, jumpFileStrm, snrFileStrm, tableFileStrm, patronFileStrm;


};

#endif /* MCMC_H_ */

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
//class FileHandler;
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
	void		proposeSeating(void);
	void		proposeParm(Parm*);
	void		proposeBranch(void);
	void		proposeBranchSwap(void);
	void		proposeBranchRescale();
	void		proposeParmRotate(void);
	void		proposeSigmaSwap(void);
	void		proposeParmShift(void);
	void		proposeResizeJumps(void);
	void		proposeAddJump(void);
	void		proposeRemoveJump(void);


	//Parm*		pickParm(void);
	//Parm*		pickBranch(void);



	void		openFiles(std::string outputFilePath);

	void		printChainState(int n, double v);
	void		printPatronState(int n, Patron* p);
	void		printTableState(int n);
	void		printAcceptanceInfo(void);
	std::string printBool(bool);
	double		safeExp(double lnX);
	std::string printInt(int);
	std::string printTableChange(Table*, Table*);

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

	int			numTranscripts;
	int			numTimepoints;
	int			numTaxa;
	int			numNodes;
	int			numBranches;
	int			numParms;

	double		alphaCRP;
	double		auxCRP;
	bool		useCRP;

	int			modelType;
	bool		fixBranches;

	std::vector<double> proposalProbs;
	double		oldLnL, newLnL;
	std::string proposeStr;

	int			patronIndex;

	std::ofstream parmFileStrm, tableFileStrm, patronFileStrm;


};

#endif /* MCMC_H_ */

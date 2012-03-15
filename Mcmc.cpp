/*
 * Mcmc.cpp
 *
 *  Created on: Dec 16, 2010
 *
 */


#define WORKING_JUMP 1
#define DEBUG2 0
#define DEBUG_PROPOSE_JUMP 0
#define DEBUG_PROPOSE_PARM 0
#define JUMP_LIKE_BY_BRANCH 0

#include "Mcmc.h"

Mcmc::Mcmc(Expression *ep, Model *mp, Settings *sp, MbRandom *rp, Topology *tp)
{
	expressionPtr = ep;
	modelPtr = mp;
	randomPtr = rp;
	settingsPtr = sp;
	topologyPtr = tp;
	tableListPtr = modelPtr->getTableListPtr();
	patronListPtr = modelPtr->getPatronListPtr();

	numCycles = settingsPtr->getNumCycles();
	printFreqMH = settingsPtr->getPrintFreqMH();
	printFreqCRP = settingsPtr->getPrintFreqCRP();
	printFreqJump = settingsPtr->getPrintFreqJump();
	printFreqStdOut = settingsPtr->getPrintFreqStdOut();

	// MCMC settings (override)
	//numCycles = 1000000;
	//printFreqMH = 1000;//0000;
	//printFreqJump = 100000;//numCycles+1;
	//printFreqCRP = 50;

	// Data dimension settings
	numTranscripts = expressionPtr->getNumTranscripts();
	numTimepoints = expressionPtr->getNumTimepoints();
	numTaxa = expressionPtr->getNumTaxa();
	numNodes = 2 * numTaxa - 1;
	numBranches = numNodes - 1;
	numParms = 2;

	// Stepping Stone integration settings
	useSteppingStone = settingsPtr->getUseSteppingStone();
	if (useSteppingStone)
		betaSteppingStone = settingsPtr->getBetaSteppingStone();
	else
		betaSteppingStone = 1.0;


	oldLnL = 0.0;
	newLnL = 0.0;
	oldKb = 0.0;
	oldKj = 0.0;
	newKb = 0.0;
	newKj = 0.0;
	proposeStr = "";


	modelType = settingsPtr->getModelType();
	fixBranches = settingsPtr->getFixBranches();
	useJumpKernel = settingsPtr->getUseJumpKernel();


	// PROPOSAL PROBABILITY INITIALIZATION
	double sum = 0.0;

	// propose change to a random parameter
	proposalProbs.push_back( 1.0 );

	// propose change to a random branch
	if (fixBranches)
		proposalProbs.push_back( 0.0 );
	else if (!fixBranches)
		proposalProbs.push_back( 1.0 );

	// propose change to sampled jumps
	if (useJumpKernel)
		proposalProbs.push_back( 5.0 );
	else if (!useJumpKernel)
		proposalProbs.push_back( 0.0 );

	// move along variance ridge
	proposalProbs.push_back(0.0);

	for (int i = 0; i < (int)proposalProbs.size(); i++) sum += proposalProbs[i];
	for (int i = 0; i < (int)proposalProbs.size(); i++) proposalProbs[i] = proposalProbs[i] / sum;

	std::string dataFilePath = settingsPtr->getOutputDirPath() + settingsPtr->getOutputFileName();
	openFiles(dataFilePath);
	runChain();
}

Mcmc::~Mcmc()
{

	// TODO Auto-generated destructor stub
	//delete modelPtr;
	//delete randomPtr;
	//delete settingsPtr;
	//delete expressionPtr;

}


void Mcmc::runChain()
{
	std::cout << "INITIALIZING: MCMC\n";


	// initialize jump samples
	for (int n = 0; n < numNodes; n++)
	{
		Node* p = topologyPtr->getNode(n);
		modelPtr->proposeJumpSize(p, tableListPtr->front()->getParmVector(), 1);
	}
/*
	std::cout << "Topology space 0\n";
	topologyPtr->printJumpSizes(0);
	std::cout << "Topology space 1\n";
	topologyPtr->printJumpSizes(1);
*/

	oldLnL = modelPtr->modelLnLikelihood(1);
	oldKj = topologyPtr->getRoot()->getKj(1);
	oldKb = topologyPtr->getRoot()->getKb(1);

	//std::cout << "INITIALIZING: Kj " << oldKj  << "\n";

	topologyPtr->copyNodeSpaces(0, 1);

	int postBurnN = 0;
	int burnIn = 1000;

	std::cout << "MCMC: Run chain";
	if (useSteppingStone)
		std::cout << " (beta=" << betaSteppingStone << ")";
	std::cout << "\n";

	// get initial likelihood
	for (int n = 0; n <= numCycles; n++)
	{
		double printOldLnL = oldLnL;

		proposeState();

		// print information to the screen
		if (n % printFreqMH == 0)
		{
			printChainState(n, oldLnL);
		}

		if (n % printFreqStdOut == 0)
		{
			std::cout << std::setw(5) << n << " -- ";
			std::cout << std::fixed << std::setprecision(8) << printOldLnL << "\t->\t" << newLnL << "\t" << proposeStr << "\n";
#if DEBUG2
			std::cout << "\n";
#endif
		}


		if (n % printFreqJump == 0)
		{
			printChainJumps(n);
		}

		if (n > burnIn)
		{
			postBurnN++;
		}

	}

	std::cout << "MCMC: Complete.";
	if (useSteppingStone)
		std::cout << " (beta=" << betaSteppingStone << ")";
	std::cout << "\n";
}

void Mcmc::proposeState(void)
{
	double u = randomPtr->uniformRv();
	double sum = 0.0;
	int uIndex = -1;
	for (int i = 0; i < (int)proposalProbs.size(); i++)
	{
		sum += proposalProbs[i];
		if (u < sum)
		{
			uIndex = i;
			break;
		}
	}

	// update random parameter for random table element
	if (uIndex == 0)
	{
		//std::cout << "\tPropose parameter change\n";
		Table* t = modelPtr->pickTableAtRandom();
		Parm* p = modelPtr->pickParmAtRandom(t);
		proposeParm(p);
	}

	else if (uIndex == 1)
	{

		//std::cout << "\tPropose branch change\n";
		Table* t = modelPtr->pickTableAtRandom();
		Parm* p = modelPtr->pickBranchAtRandom(t);
		proposeParm(p);


	}

	else if (uIndex == 2)
	{
	//std::cout << "\tPropose jump size change\n";
		proposeChangeJump();

	}

	else if (uIndex == 3)
	{
		//std::cout << "\tPropose remove jump\n";
		proposeRidgeMove();
	}
}

void Mcmc::proposeParm(Parm* parm)
{
	//parm->print();

	// propose new parameters
	double oldVal = parm->getValue();
	double lnProposalRatio = parm->change();
	double newVal = parm->getValue();

	// calculate the prior ratio
	double lnPriorRatio = parm->lnPriorRatio();

	// calculate the likelihood ratio
	newLnL = modelPtr->modelLnLikelihood(1);
	newKb = topologyPtr->getRoot()->getKb(1);
	newKj = topologyPtr->getRoot()->getKj(1);
	double lnLikelihoodRatio = newLnL - oldLnL;
	if (useSteppingStone)
	{
		//lnLikelihoodRatio *= betaSteppingStone;
		lnLikelihoodRatio = betaSteppingStone * (newKb - oldKb) + (newKj - oldKj);
	}

	// accept or reject MH ratio
	bool acceptState = false;
	double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );
	double u = randomPtr->uniformRv();
	if (u < r)
	{
		acceptState = true;
	}

	#if DEBUG_PROPOSE_PARM
	std::cout << "\t" << parm->getName() << "\t" << oldVal << " -> " << newVal << "\n";
	std::cout << "\tlnProposalRatio:\t" << lnProposalRatio << "\n";
	std::cout << "\tlnPriorRatio:\t\t" << lnPriorRatio << "\n";
	std::cout << "\toldLnL:\t\t\t" << oldLnL << "\n";
	std::cout << "\tnewLnL:\t\t\t" << newLnL << "\n";
	std::cout << "\tlnLikeRatio:\t\t" << lnLikelihoodRatio << "\n";
	std::cout << "\tMH Ratio:\t\t" << r << " (" << lnLikelihoodRatio + lnPriorRatio + lnProposalRatio << ")\n";
	std::cout << "\tu:\t\t\t" << u << "\t\t\t" << "r:\t" << r << "\n";
	#endif



	// update the state of the chain
	if (acceptState == true)
	{
		//topologyPtr->flipAllActiveParms();
		parm->incrementNumAccepted();
		modelPtr->keep(parm);
		topologyPtr->copyNodeSpaces(0, 1);
		oldLnL = newLnL;
		oldKb = newKb;
		oldKj = newKj;
	}
	else
	{
		modelPtr->restore(parm);
		topologyPtr->copyNodeSpaces(1, 0);
	}

	proposeStr = parm->getName() + ",\t\tUpdate : " + printBool(acceptState);
}


void Mcmc::proposeChangeJump(void)
{
	// choose branch
	Table* t = tableListPtr->front();
	int numJumpAccept = 0;
	int numJumpReject = 0;

	// propose new jump displacement for all internal nodes
	for (int i = 0; i < numNodes; i++)
	{
		Node* p = topologyPtr->getDownPassNode(i);

		if (p->getAnc() != NULL)
		{
			// calculate the proposal ratio
			double lnProposalRatio = modelPtr->proposeJumpSize(p, t->getParmVector(), 1);

			// calculate the prior ratio
			double lnPriorRatio = 0.0;

			// calculate the likelihood ratio

#if JUMP_LIKE_BY_BRANCH
			// BY BRANCH
			// save old lnL for node's ancestor
			double oldNodeLnL = p->getAnc()->getK(0);
			double oldNodeKb = p->getAnc()->getKb(0);
			double oldNodeKj = p->getAnc()->getKj(0);

			// set proposed jump for node
			double lnJumpLike = modelPtr->jumpLnLikelihood(p, t->getParmVector(), 1);

			if (p->getLft() == NULL & p->getRht() == NULL)
			{
				p->setK(lnJumpLike, 1);
				p->setKj(lnJumpLike, 1);
			}
			else
			{
				p->setKb(p->getLft()->getKb(1) + p->getRht()->getKb(1), 1); //Check if this is right, update yet?
				p->setKj(lnJumpLike + p->getLft()->getKj(1) + p->getRht()->getKj(1), 1);
				p->setK(lnJumpLike + p->getLft()->getK(1) + p->getRht()->getK(1), 1);
			}

			// save new lnL for node's ancestor
			double newNodeLnL = modelPtr->driftLnLikelihood(p->getAnc(), t->getParmVector(), 1); //check to make sure p->getAnc() is right
			double newNodeKb = p->getAnc()->getKb(1);
			double newNodeKj = p->getAnc()->getKj(1);

			newLnL = oldLnL - oldNodeLnL + newNodeLnL;
			newKb = oldKb - oldNodeKb + newNodeKb;
			newKj = oldKj - oldNodeKj + newNodeKj;

//#if DEBUG
			std::cout << "\tn" << p->getIndex() << "\n";
			std::cout << "\t\toldNode\t" << oldNodeLnL << "\t" << oldNodeKb << "\t" << oldNodeKj << "\n";
			std::cout << "\t\told\t" << oldLnL << "\t" << oldKb << "\t" << oldKj << "\n";
			std::cout << "\t\tnewNode\t" << newNodeLnL << "\t" << newNodeKb << "\t" << newNodeKj << "\n";
			std::cout << "\t\tnew\t" << newLnL << "\t" << newKb << "\t" << newKj << "\n";
//#endif

			double lnLikelihoodRatio = newLnL - oldLnL;
			if (useSteppingStone)
			{
				//lnLikelihoodRatio *= betaSteppingStone;
				lnLikelihoodRatio = betaSteppingStone * (newKb - oldKb) + (newKj - oldKj);

			}


#else
			// BY TREE

			// jump for p
			// drift for model


#if WORKING_JUMP
			//double newK = topologyPtr->getRoot()->getK();

			newKj = modelPtr->jumpLnLikelihood(p, t->getParmVector(), 1);
			oldKj = modelPtr->jumpLnLikelihood(p, t->getParmVector(), 0);
			newKb = modelPtr->modelDriftOnlyLnLikelihood(1);
			//oldKb = modelPtr->modelDriftOnlyLnLikelihood(0);

			//modelPtr->modelLnLikelihood(0);
			//double newKb = topologyPtr->getRoot()->getKb();
			double diffKj0 = topologyPtr->getRoot()->getKj(1) - topologyPtr->getRoot()->getKj(0);
			double diffKb0 = topologyPtr->getRoot()->getKb(1) - topologyPtr->getRoot()->getKb(0);
			double diffKj1 = newKj - oldKj;
			double diffKb1 = newKb - oldKb;
			//	std::cout << "\t\tdiffK:\t" << diffKj0 << "\t" << diffKj1 << "\t" << diffKb0 << "\t" << diffKb1 << "\n";
#else
			newLnL = modelPtr->modelLnLikelihood(1);
			newKj = topologyPtr->getRoot()->getKj(1);
			newKb = topologyPtr->getRoot()->getKb(1);

#endif




//			newLnL += newKb - oldKb;


			double lnLikelihoodRatio = newKb - oldKb + newKj - oldKj;
			if (useSteppingStone)
			{
				lnLikelihoodRatio = betaSteppingStone*(newKb - oldKb) + (newKj - oldKj);
			}
#endif

			// accept or reject MH ratio
			bool acceptState = false;
			double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );
			double u = randomPtr->uniformRv();
			if (u < r)
			{
				acceptState = true;
			}

		#if DEBUG_JUMP
			std::cout << "\tlnProposalRatio:\t" << lnProposalRatio << "\n";
			std::cout << "\tlnPriorRatio:\t\t" << lnPriorRatio << "\n";
			std::cout << "\toldLnL:\t\t\t" << oldLnL << "\n";
			std::cout << "\tnewLnL:\t\t\t" << newLnL << "\n";
			std::cout << "\tlnLikeRatio:\t\t" << lnLikelihoodRatio << "\n";
			std::cout << "\tMH Ratio:\t\t" << r << " (" << lnLikelihoodRatio + lnPriorRatio + lnProposalRatio << ")\n";
			std::cout << "\tu:\t\t\t" << u << "\t\t\t" << "r:\t" << r << "\n";
			#endif

			// update the state of the chain
			if (acceptState == true)
			{
				p->copySpace(0, 1);
				oldKb = newKb;
				oldKj = newKj;
				numJumpAccept++;

			}
			else
			{
				p->copySpace(1, 0);
				numJumpReject++;
			}
		}
	}

	modelPtr->modelLnLikelihood(0);
	oldKj = topologyPtr->getRoot()->getKj(0);
	oldKb = topologyPtr->getRoot()->getKb(0);
	oldLnL = oldKj + oldKb;
	newLnL = oldKj + oldKb;

	proposeStr = "\t\tChngJmp: " + Util::intToString(numJumpAccept) + "/" + Util::intToString(numNodes);
}

void Mcmc::proposeRidgeMove(void)
{

	// get current parameters
	std::vector<Parm*> tempParmVector = tableListPtr->front()->getParmVector();

	Parm* sbmP = tempParmVector[0];	// Sigma-BM
	Parm* ljnP = tempParmVector[1];	// Lambda-JN
	Parm* sjnP = tempParmVector[2];	// Sigma-JN

	Sigma* sbm = dynamic_cast<Sigma*>(sbmP);	// Sigma-BM
	Sigma* ljn = dynamic_cast<Sigma*>(ljnP);	// Lambda-JN
	Sigma* sjn = dynamic_cast<Sigma*>(sjnP);	// Sigma-JN


	// propose new parameters

	double oldSbmVal = sbm->getValue();
	double oldLjnVal = ljn->getValue();
	double oldSjnVal = sjn->getValue();

	double sumVarVal = pow(oldSbmVal, 2) + oldLjnVal * pow(oldSjnVal, 2);
	double lnProposalRatio = sbm->getActiveVar()->changeTruncated( sumVarVal );
	double newSbmVal = sbm->getValue();
	double newLjnVal = ( sumVarVal - pow(newSbmVal,2) ) / pow(oldSjnVal, 2);
	ljn->getActiveVar()->change( newLjnVal );


	// calculate the prior ratio
	double lnPriorRatio = sbmP->lnPriorRatio() + ljnP->lnPriorRatio();

	// calculate the likelihood ratio
	// modelPtr->updateModel();
	newLnL = modelPtr->modelLnLikelihood(1);
	newKb = topologyPtr->getRoot()->getKb(1);
	newKj = topologyPtr->getRoot()->getKj(1);
	double lnLikelihoodRatio = newLnL - oldLnL;
	if (useSteppingStone)
	{
		//lnLikelihoodRatio *= betaSteppingStone;
		lnLikelihoodRatio = betaSteppingStone*(newKb - oldKb) + (newKj - oldKj);
	}

	// accept or reject MH ratio
	bool acceptState = false;
	double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );
	double u = randomPtr->uniformRv();
	if (u < r)
	{
		acceptState = true;
	}

	#if DEBUG2
	std::cout << "\t" << sbmP->getName() << "\t" << oldSbmVal << " -> " << newSbmVal << "\n";
	std::cout << "\t" << ljnP->getName() << "\t" << oldLjnVal << " -> " << newLjnVal << "\n";
	std::cout << "\t" << sjnP->getName() << "\t" << oldSjnVal << " -> " << oldSjnVal << "\n";
	std::cout << "\tlnProposalRatio:\t" << lnProposalRatio << "\n";
	std::cout << "\tlnPriorRatio:\t\t" << lnPriorRatio << "\n";
	std::cout << "\toldLnL:\t\t\t" << oldLnL << "\n";
	std::cout << "\tnewLnL:\t\t\t" << newLnL << "\n";
	std::cout << "\tlnLikeRatio:\t\t" << lnLikelihoodRatio << "\n";
	std::cout << "\tMH Ratio:\t\t" << r << " (" << lnLikelihoodRatio + lnPriorRatio + lnProposalRatio << ")\n";
	std::cout << "\tu:\t\t\t" << u << "\t\t\t" << "r:\t" << r << "\n";
	#endif

	// update the state of the chain
	if (acceptState == true)
	{
		sbmP->incrementNumAccepted();
		ljnP->incrementNumAccepted();
		modelPtr->keep(sbmP);
		modelPtr->keep(ljnP);
		topologyPtr->copyNodeSpaces(0, 1);
		oldLnL = newLnL;
		oldKb = newKb;
		oldKj = newKj;
	}
	else
	{
		modelPtr->restore(sbmP);
		modelPtr->restore(ljnP);
		topologyPtr->copyNodeSpaces(1, 0);
	}

	proposeStr = "\t\tRdgMove: " + printBool(acceptState);
}




void Mcmc::openFiles(std::string fn) {

	std::string tf = fn + ".t";
	std::string pf = fn + ".p";
	std::string jf = fn + ".j";
	std::string patronFile = fn + ".patron";
	std::string tableFile = fn + ".table";

	/*
	treeFileStrm.open( tf.c_str(), std::ios::out );
	if ( !treeFileStrm )
	{
		std::cerr << "ERROR: Problem opening tree output file" << std::endl;
		exit(1);
	}
	*/

	parmFileStrm.open( pf.c_str(), std::ios::out );
	if ( !parmFileStrm )
	{
		std::cerr << "ERROR: Problem opening parameter output file" << std::endl;
		exit(1);
	}

	jumpFileStrm.open( jf.c_str(), std::ios::out );
	if ( !jumpFileStrm )
	{
		std::cerr << "ERROR: Problem opening jump output file" << std::endl;
		exit(1);
	}
	/*
	if (useCRP)
	{
		tableFileStrm.open( tableFile.c_str(), std::ios::out );

		if ( !tableFileStrm )
		{
			std::cerr << "ERROR: Problem opening table output file" << std::endl;
			exit(1);
		}

		patronFileStrm.open( patronFile.c_str(), std::ios::out );
		if ( !patronFileStrm )
		{
			std::cerr << "ERROR: Problem opening patron output file" << std::endl;
			exit(1);
		}
	}
	*/
}

void Mcmc::printChainState(int n, double lnL) {

	if (n == 0)
	{
		std::string pHeaderStr = "";

		for (std::list<Table*>::iterator it_t = tableListPtr->begin(); it_t != tableListPtr->end(); it_t++)
		{
			for (std::vector<Parm*>::const_iterator it_p = (*it_t)->getParmVector().begin(); it_p != (*it_t)->getParmVector().end(); it_p++)
			{
				pHeaderStr += (*it_p)->getParameterHeader();
			}
			if (fixBranches == false)
			{
				for (std::vector<Parm*>::const_iterator it_p = (*it_t)->getBranchVector().begin(); it_p != (*it_t)->getBranchVector().end();it_p++)
				{
					pHeaderStr += (*it_p)->getParameterHeader();
				}
			}
		}

		parmFileStrm << "Cycle\tlnL\tkb\tkj\t" << pHeaderStr << std::endl;
	}

	std::string pStr = "";
	pStr += "\t" + printDouble(oldKb);
	pStr += "\t" + printDouble(oldKj);
	pStr += "\t";

	for (std::list<Table*>::iterator it_t = tableListPtr->begin(); it_t != tableListPtr->end(); it_t++)
	{
		for (std::vector<Parm*>::const_iterator it_p = (*it_t)->getParmVector().begin(); it_p != (*it_t)->getParmVector().end(); it_p++)
		{
			pStr += (*it_p)->getParameterStr();
		}
		if (fixBranches == false)
		{
			for (std::vector<Parm*>::const_iterator it_p = (*it_t)->getBranchVector().begin(); it_p != (*it_t)->getBranchVector().end();it_p++)
			{
				pStr += (*it_p)->getParameterStr();
			}
		}
	}

	parmFileStrm << n << '\t' << std::fixed << std::setprecision(2) << lnL << '\t' << pStr << std::endl;

}

void Mcmc::printChainJumps(int n)
{
	Node* p;
	if (n == 0)
	{
		std::string jHeaderStr = "";

		for (int i = 0; i < numNodes; i++)
		{
			p = topologyPtr->getDownPassNode(i);
			if (p->getAnc() != NULL)
			{
				jHeaderStr += "\t" + printInt(p->getIndex());
			}
		}

		jumpFileStrm << "Cycle" << jHeaderStr << std::endl;
	}

	std::string jStr = "";
	for (int i = 0; i < numNodes; i++)
	{
		p = topologyPtr->getDownPassNode(i);
		if (p->getAnc() != NULL)
		{
			jStr += "\t" + printDouble(p->getSumJumpSize(0));
		}
	}

	jumpFileStrm << n << jStr << std::endl;
}



std::string Mcmc::printBool(bool tf)
{
	if (tf) return "True";
	else return "\tFalse";
}

std::string Mcmc::printInt(int number)
{
   std::stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

std::string Mcmc::printDouble(double number)
{
	std::stringstream ss;
	ss << number;
	return ss.str();
}

double Mcmc::safeExp(double lnX) {

	if (lnX < -300.0)
		return 0.0;
	else if (lnX > 0.0)
		return 1.0;
	else
		return exp(lnX);
}



#if 0

//

/*
for (int i = numRegTables; i < numRegTables + numAuxTables; i++)
{
	if (probIndex == i)
	{
		t = (*it_t);
	}
	else
	{
		std::cout << "remove Aux\n";
		std::vector<Table*>::iterator it_temp = it_t;
		it_temp++;
		auxTableList.remove(*it_t);
		it_t = it_temp;
		std::cout << "ok\n";
	}
	it_t++;
}
*/


/*
std::string Mcmc::trueFalse(bool tf)
{
	if (tf) return "True";
	else return "False";
}
*/


/*
std::cout << "i:" << i << "\n";
if (it_t != auxTableList.end())
{
	it_t++;
	(*it_t)->print();
}
*/

/*
Table* Mcmc::reseatPatron(void)
{
	// 0. Initialize CRP structures
	int numRegTables = tableListPtr->size();
	int numAuxTables = auxCRP * 5 + 1; // constant
	// int numAuxTables = auxCRP * numRegTables + 1; // proportional
	std::vector<double> proposalLnLs;
	std::vector<double> proposalProbs;
	std::vector<Table*> auxTableList;

	// 1. Examine the patrons according to order
	Patron* p = (*patronListPtr->begin());
	patronListPtr->push_back(p);
	patronListPtr->pop_front();

	modelPtr->updateStationaryProbs();
	oldTable = p->getTable();
	Table* pTable = p->getTable();
	p->stand();
	// TODO: Could not get waitTable to work. Seems to be unable to delete itself. Doing it the ugly way.
	// p->getTable()->waitTable();
	if (pTable->getPatronList().size() == 0)
	{
		tableListPtr->remove(pTable);
	}

	// 2. Seat the patron at (1, ..., K) tables and calculate likelihoods
	// std::cout << "numRegTables = " << numRegTables << "\n";
	for (std::vector<Table*>::iterator it_t = tableListPtr->begin(); it_t != tableListPtr->end(); it_t++)
	{
		// sit the patron at the next table
		p->sit(*it_t);

		// update the model
		modelPtr->updateModel();

		// calculate the likelihoods
		proposalLnLs.push_back(modelPtr->modelLnLikelihood());
		// std::cout << "\t" << proposalLnLs.back() << "\n";

		// calculate the likelihood-weighted probabilities
		double val = proposalLnLs.back() * (*it_t)->getPatronList().size() / (numRateElements - 1 + alphaCRP);

		// store the value
		proposalProbs.push_back(val);
		// std::cout << "\t\t" << proposalProbs.back() << "\n";

		// remove patron from table
		p->stand();
	}

	// 3. Seat patrons at (K+1, ..., K+M) auxillary tables
	// std::cout << "numAuxTables = " << numAuxTables << "\n";
	for (int i = 0; i < numAuxTables; i++)
	{
		// draw new auxillary parameters
		std::vector<Parm*> auxParmVector;
		auxParmVector.push_back(new Eta(randomPtr, "Eta"));
		auxParmVector.push_back(modelPtr->getGeoFrequencies());

		// add the table to the auxillary table vector
		auxTableList.push_back(new Table(&auxTableList, auxParmVector));
		Table* auxTable = auxTableList.back();		// TODO: .back() is a const&, .end() is an ::iterator

		// sit the patron at the new table
		p->sit(auxTable);

		// update the model
		modelPtr->updateModel();

		// calculate the likelihoods
		proposalLnLs.push_back(modelPtr->modelLnLikelihood());
		// std::cout << "\t" << proposalLnLs.back() << "\n";

		// calculate probabilities (unnormalized)
		double val = proposalLnLs.back() * (alphaCRP / numAuxTables) / (numRateElements - 1 + alphaCRP);

		// store the value
		proposalProbs.push_back(val);
		// std::cout << "\t\t" << proposalProbs.back() << "\n";

		// remove patron from table
		p->stand();
	}

	// 4. normalize proposed seating probabilities

	// std::cout << "Normalized proposal probs:\n";
	double probNorm = 0.0;
	for (int i = 0; i < (int)proposalProbs.size(); i++)
	{
		if (std::isnan(proposalProbs[i]) == false)
		{
			probNorm += proposalProbs[i];
		}
	}
	for (int i = 0; i < (int)proposalProbs.size(); i++)
	{
		proposalProbs[i] = proposalProbs[i] / probNorm;
		if (proposalProbs[i] < 0.0)
		{
			std::cerr << "ERROR: normalized reseating probabilities contain at least one element < 0.0\n";
		}
		// std::cout << "\t" << proposalProbs[i] << "\n";
	}

	// 5. determine which seating arrangement to keep
	int probIndex = -1;
	double u = randomPtr->uniformRv();
	double probSum = 0.0;
	for (int i = 0; i < (int)proposalProbs.size(); i++)
	{
		if (std::isnan(proposalProbs[i]) == false)
		{
			probSum += proposalProbs[i];
			if (u < probSum)
			{
				probIndex = i;
				break;
			}
		}
		else
		{
			std::cerr << "ERROR: proposalProbs[" << i << "] = NaN\n";
		}
	}
	if (probIndex == -1)
	{
		std::cerr << "ERROR: failed to sample new seating arrangement.\n";
	}
	else if (probIndex > numRegTables + numAuxTables)
	{
		std::cerr << "ERROR: attempted to seat patron outside of tableVector bounds.\n";
	}
	// std::cout << "PICKED: " << probIndex << "\n";


	// Set the unperturbed lnLikelihood for the selected seating arrangement.
	// TODO: Oops, we can set the oldLnL value fine, but the newLnL value should not be updated until the new rate is proposed.
	oldLnL = proposalLnLs[probIndex];

	// commit the new seating arrangement to the proposed MCMC state
	if (probIndex < numRegTables - 1)
	{
		// std::cout << "Reg. table picked\n";
		// find the correct regular table
		std::vector<Table*>::iterator it_t = tableListPtr->begin();
		for (int i = 0; i < probIndex; i++)
		{
			it_t++;
		}
		p->sit(*it_t);
		// (*it_t)->print();
		return *it_t;
	}
	else
	{
		// std::cout << "Aux. table picked\n";
		// find the correct auxillary table
		std::vector<Table*>::iterator it_t = auxTableList.begin();
		Table* t = NULL;

		//
		for (int i = numRegTables; i < probIndex; i++)
		{
			it_t++;
		}
		t = (*it_t);
		auxTableList.remove(*it_t);

		// p->print();
		p->sit(t);
		tableListPtr->push_back(t);
		// p->print();
		// t->print();

		if (t == NULL)
			std::cerr << "ERROR: failed to find auxillary table in auxTableList!\n";

		// clear remaining auxillary tables
		for (std::vector<Table*>::iterator it_a = auxTableList.begin(); it_a != auxTableList.end(); it_a++)
		{
			delete *it_a;
		}

		return t;
	}

	// this should never occur
	return NULL;
}
*/

#endif

#if 0
	/*
	// calculate the likelihood ratio
	double lnLikelihoodRatio = newLnL - oldLnL;

	// calculate likelihood ratio
	double r = safeExp( lnLikelihoodRatio ); // TODO: Is the prior part of Alg. 8? (e.g. non-conj. prior)	// + lnPriorRatio + lnProposalRatio );

	bool acceptState = false;
	u = randomPtr->uniformRv();
	std::cout << "\tu=" << u << "\tr=" << r << "\n";
	if (u < r)
	{
		acceptState = true;
	}

	// update the state of the chain
	if (acceptState == true)
	{
		// update lnL
		oldLnL = newLnL;

		// update seating
		if (probIndex < numRegTables)
		{
			// std::cout << "Reg. table picked\n";
			// find the correct regular table



			newTable = (*tableListPtr)[probIndex];

			if (newTable == NULL)
				std::cerr << "ERROR: failed to find regular table in tableListPtr!\n";

			p->print();
			p->sit(newTable);
			p->print();
			// (*it_t)->print();
			// return *it_t;
			// return lnL
		}
		else
		{
			// std::cout << "Aux. table picked\n";
			// find the correct auxillary table


			//Table* t = NULL;
			//t = (*it_t);
			std::cout << "1\n";

			//newTable = *it_t;
			newTable = auxTableList[probIndex - numRegTables];

			std::cout << "2\n";

			//auxTableList.erase(it_t);
			auxTableList.erase(auxTableList.begin() + probIndex - numRegTables);

			// p->print();
			p->sit(newTable);

			if (newTable == NULL)
				std::cerr << "ERROR: failed to find auxillary table in auxTableList!\n";

			tableListPtr->push_back(newTable);
			// p->print();
			// t->print();

			// clear remaining auxillary tables
			for (std::vector<Table*>::iterator it_a = auxTableList.begin(); it_a != auxTableList.end(); it_a++)
			{
				delete *it_a;
			}

			// return lnL

		}

		// remove old empty regular table
		if (oldTable->getPatronList().size() == 0)
		{
			for (std::vector<Table*>::iterator it_t = tableListPtr->begin(); it_t != tableListPtr->end(); it_t++)
			{
				if (oldTable == *it_t)
				{
					tableListPtr->erase(it_t);
					break;
				}
			}
		}
	}
	else
	{
		// std::cout << "reject seating\n";

		// reseat at old Table
		//p->stand();
		p->sit(oldTable);
		p->print();

		// clean auxillary tables from memory
		for (std::vector<Table*>::iterator it_a = auxTableList.begin(); it_a != auxTableList.end(); it_a++)
		{
			delete *it_a;
		}

		// update any performance statistics
		// ...
		// ...
	}
	*/
#endif

/*
void Mcmc::proposeAddJump(void)
{
	// choose branch
	Table* t = tableListPtr->front();
	Node* p = topologyPtr->getRandomNodeByLength();

	// calculate the proposal ratio
	double lnProposalRatio = modelPtr->proposeAddJump(p, t->getParmVector(), 1);

	// calculate the prior ratio
	double lnPriorRatio = 0.0;

	// calculate the likelihood ratio
	newLnL = modelPtr->modelLnLikelihood();
	double lnLikelihoodRatio = newLnL - oldLnL;

	// accept or reject MH ratio
	bool acceptState = false;
	double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );
	double u = randomPtr->uniformRv();
	if (u < r)
	{
		acceptState = true;
	}

	#if DEBUG2
	std::cout << "\tlnProposalRatio:\t" << lnProposalRatio << "\n";
	std::cout << "\tlnPriorRatio:\t\t" << lnPriorRatio << "\n";
	std::cout << "\toldLnL:\t\t\t" << oldLnL << "\n";
	std::cout << "\tnewLnL:\t\t\t" << newLnL << "\n";
	std::cout << "\tlnLikeRatio:\t\t" << lnLikelihoodRatio << "\n";
	std::cout << "\tMH Ratio:\t\t" << r << " (" << lnLikelihoodRatio + lnPriorRatio + lnProposalRatio << ")\n";
	std::cout << "\tu:\t\t\t" << u << "\t\t\t" << "r:\t" << r << "\n";
	#endif

	// update the state of the chain
	if (acceptState == true)
	{
		p->copySpace(0, 1);
		topologyPtr->incrementNumJumps();
		oldLnL = newLnL;
	}
	else
	{
		p->copySpace(1, 0);
	}

	proposeStr = "n" + printInt(p->getIndex()) + "\t\tAdd jmp: " + printBool(acceptState);
}

void Mcmc::proposeRemoveJump(void)
{
	// choose branch
	Table* t = tableListPtr->front();
	Node* p = topologyPtr->getRandomNodeWithJumps();

	// if the tree contains jumps, remove a jump
	if (p != NULL)
	{
		// calculate the proposal ratio
		double lnProposalRatio = modelPtr->proposeRemoveJump(p, t->getParmVector(), 1);

		// calculate the prior ratio
		double lnPriorRatio = 0.0;

		// calculate the likelihood ratio
		newLnL = modelPtr->modelLnLikelihood();
		double lnLikelihoodRatio = newLnL - oldLnL;

		// accept or reject MH ratio
		bool acceptState = false;
		double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );
		double u = randomPtr->uniformRv();
		if (u < r)
		{
			acceptState = true;
		}

	#if DEBUG2
		std::cout << "\tlnProposalRatio:\t" << lnProposalRatio << "\n";
		std::cout << "\tlnPriorRatio:\t\t" << lnPriorRatio << "\n";
		std::cout << "\toldLnL:\t\t\t" << oldLnL << "\n";
		std::cout << "\tnewLnL:\t\t\t" << newLnL << "\n";
		std::cout << "\tlnLikeRatio:\t\t" << lnLikelihoodRatio << "\n";
		std::cout << "\tMH Ratio:\t\t" << r << " (" << lnLikelihoodRatio + lnPriorRatio + lnProposalRatio << ")\n";
		std::cout << "\tu:\t\t\t" << u << "\t\t\t" << "r:\t" << r << "\n";
		#endif

		// update the state of the chain
		if (acceptState == true)
		{
			p->copySpace(0, 1);
			topologyPtr->decrementNumJumps();
			oldLnL = newLnL;
		}
		else
		{
			p->copySpace(1, 0);
		}

		proposeStr = "n" + printInt(p->getIndex()) + "\t\tRem jmp: " + printBool(acceptState);
	}
	// if the tree contains no jumps, bypass this move
	else
	{
		proposeStr = "\t\tRem jmp: skipped, numJumps=0";
	}
}

*/


#if 0
void Mcmc::proposeBranchRescale(void)
{
	// select the taus for a table
	Table* t = modelPtr->pickTableAtRandom();
	std::vector<Tau*> tauVector;
	std::vector<Parm*> branchVector = t->getBranchVector();

	for (int i = 0; i < numBranches; i++)
	{
		tauVector.push_back(dynamic_cast<Tau*> (branchVector[i]));
	}

	// rescale branches
	// calculate the prior ratio
	// calculate proposal ratio
	double lnProposalRatio = 0.0;
	double lnPriorRatio = 0.0;
	double scale = randomPtr->exponentialRv(1.0);
	double rate;
	for (int i = 0; i < numBranches; i++)
	{
		rate = tauVector[i]->getActiveBrLen()->getValue();
		lnProposalRatio += tauVector[i]->getActiveBrLen()->change(scale * rate);
		lnPriorRatio += tauVector[i]->lnPriorRatio();
	}

	// calculate the likelihood ratio
	modelPtr->updateModel();
	newLnL = modelPtr->modelLnLikelihood();

	double lnLikelihoodRatio = newLnL - oldLnL;

	// accept or reject MH ratio
	double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );

	bool acceptState = false;
	double u = randomPtr->uniformRv();
	if (u < r)
	{
		acceptState = true;
	}

#if DEBUG2
	std::cout << "\tlnProposalRatio:\t" << lnProposalRatio << "\n";
	std::cout << "\tlnPriorRatio:\t\t" << lnPriorRatio << "\n";
	std::cout << "\toldLnL:\t\t\t" << oldLnL << "\n";
	std::cout << "\tnewLnL:\t\t\t" << newLnL << "\n";
	std::cout << "\tlnLikeRatio:\t\t" << lnLikelihoodRatio << "\n";
	std::cout << "\tMH Ratio:\t\t" << r << " (" << lnLikelihoodRatio + lnPriorRatio + lnProposalRatio << ")\n";
	std::cout << "\tu:\t\t\t" << u << "\t\t\t" << "r:\t" << r << "\n";
	#endif

	// update the state of the chain
	if (acceptState == true)
	{
		// TODO: add numswapped
		for (int i = 0; i < numBranches; i++)
		{
			tauVector[i]->incrementNumAccepted();
			modelPtr->keep(tauVector[i]);
		}
		oldLnL = newLnL;
	}
	else
	{
		for (int i = 0; i < numBranches; i++)
		{
			tauVector[i]->incrementNumAccepted();
			modelPtr->restore(tauVector[i]);
		}
	}

	proposeStr = "\t\tRescale: " + printBool(acceptState);

}
#endif


#if 0
void Mcmc::proposeResizeJumps(void)
{

	//std::cout << "\n";
	//parm->print();

	// propose new parameters
	double lnProposalRatio = 0.0;//parm->change();

	// calculate the prior ratio
	double lnPriorRatio = 0.0;//parm->lnPriorRatio();

	// sample new jumps under parameters
	modelPtr->sampleJumpSizesForTree();

	// calculate the likelihood ratio
	// modelPtr->updateModel();
	newLnL = modelPtr->modelLnLikelihood();
	double lnLikelihoodRatio = newLnL - oldLnL;

	// accept or reject MH ratio
	bool acceptState = false;
	double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );
	double u = randomPtr->uniformRv();
	if (u < r)
	{
		acceptState = true;
	}

	/*
	//parm->print();
	std::cout << "\n";
	std::cout << "\tlnProposalRatio:\t" << lnProposalRatio << "\n";
	std::cout << "\tlnPriorRatio:\t\t" << lnPriorRatio << "\n";
	std::cout << "\toldLnL:\t\t\t" << oldLnL << "\n";
	std::cout << "\tnewLnL:\t\t\t" << newLnL << "\n";
	std::cout << "\tlnLikeRatio:\t\t" << lnLikelihoodRatio << "\n";
	std::cout << "\tMH Ratio:\t\t" << r << " (" << lnLikelihoodRatio + lnPriorRatio + lnProposalRatio << ")\n";
	std::cout << "\tu:\t\t\t" << u << "\t\t\t" << "r:\t" << r << "\n";
	*/

	// update the state of the chain
	if (acceptState == true)
	{
		topologyPtr->flipAllActiveParms();
		//parm->incrementNumAccepted();
		//modelPtr->keep(parm);
		oldLnL = newLnL;
	}
	else
	{
		//modelPtr->restore(parm);
	}

	proposeStr = "\t\tJmpReSz: " + printBool(acceptState);

}

#endif

#if 0
void Mcmc::proposeBranchSwap(void)
{

	// select the taus for a table
	Table* t = modelPtr->pickTableAtRandom();
	Parm* p_lambda = t->getParmVector()[0];
	Parm* p_tau1 = modelPtr->pickBranchAtRandom(t);
	Parm* p_tau2;
	do
	{
		p_tau2 = modelPtr->pickBranchAtRandom(t);
	} while(p_tau1 == p_tau2);
	Lambda* l = dynamic_cast<Lambda*> (p_lambda);
	Tau* t1 = dynamic_cast<Tau*> (p_tau1);
	Tau* t2 = dynamic_cast<Tau*> (p_tau2);

/*
	int randomInternalNode = ((numNodes - numTaxa) * randomPtr->uniformRv()) + numTaxa;

	Parm* p_tau1;

//	double t1_rate = 0.0;
//	double t2_rate = 0.0;

	Node* n = topologyPtr->getNode(randomInternalNode);
	if (n->getLft() != NULL && n->getRht() != NULL)
	{

	}*/

	// propose change
	double t1_rate = t1->getActiveBrLen()->getValue();
	double t2_rate = t2->getActiveBrLen()->getValue();
	double lnProposalRatio = t1->getActiveBrLen()->change(t2_rate) + t2->getActiveBrLen()->change(t1_rate) + l->getActiveExpMean()->changeSign();

	// accept or reject
	// calculate the prior ratio
	double lnPriorRatio = p_tau1->lnPriorRatio() + p_tau2->lnPriorRatio() + p_lambda->lnPriorRatio();

	// calculate the likelihood ratio
	modelPtr->updateModel();
	newLnL = modelPtr->modelLnLikelihood();

	double lnLikelihoodRatio = newLnL - oldLnL;

	// accept or reject MH ratio
	double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );

	bool acceptState = false;
	double u = randomPtr->uniformRv();
	if (u < r)
	{
		acceptState = true;
	}

	//std::cout << "\tlnProposalRatio:\t" << lnProposalRatio << "\n";
	//std::cout << "\tlnPriorRatio:\t\t" << lnPriorRatio << "\n";
	//std::cout << "\toldLnL:\t\t\t" << oldLnL << "\n";
	//std::cout << "\tnewLnL:\t\t\t" << newLnL << "\n";
	//std::cout << "\tlnLikeRatio:\t\t" << lnLikelihoodRatio << "\n";
	//std::cout << "\tMH Ratio:\t\t" << r << " (" << lnLikelihoodRatio + lnPriorRatio + lnProposalRatio << ")\n";
	//std::cout << "\tu:\t\t\t" << u << "\t\t\t" << "r:\t" << r << "\n";


	// update the state of the chain
	if (acceptState == true)
	{
		// TODO: add numswapped
		p_lambda->incrementNumAccepted();
		p_tau1->incrementNumAccepted();
		p_tau2->incrementNumAccepted();
		modelPtr->keep(p_lambda);
		modelPtr->keep(p_tau1);
		modelPtr->keep(p_tau2);
		oldLnL = newLnL;
	}
	else
	{
		modelPtr->restore(p_lambda);
		modelPtr->restore(p_tau1);
		modelPtr->restore(p_tau2);
	}

	proposeStr = p_tau1->getName() + " " + p_tau2->getName() + ",\tSwap   : " + printBool(acceptState);

}
void Mcmc::proposeParmRotate(void)
{
	// select (lambda,sigma) for a table
	Table* t = modelPtr->pickTableAtRandom();
	Parm* p_lambda = t->getParmVector()[0];
	Parm* p_sigma = t->getParmVector()[1];
	Lambda* l = dynamic_cast<Lambda*> (p_lambda);
	Sigma* s = dynamic_cast<Sigma*> (p_sigma);

	// propose a rotation by some theta

	double old_lambda = l->getActiveExpMean()->getValue();
	double old_sigma = s->getActiveVar()->getValue();
	double new_lambda = 0.0;
	double new_sigma = 0.0;

	do
	{
		double theta = randomPtr->normalRv(0.0, 3.0);
		new_lambda = old_lambda * cos(theta) - old_sigma * sin(theta);
		new_sigma = old_lambda * sin(theta) + old_sigma * cos(theta);
	}
	while(new_sigma < 0.0);

	// std::cout << "L'=" << new_lambda << " S'=" << new_sigma << "\tL=" << old_lambda << " S=" << old_sigma << "\n";
	double lnProp1 = l->getActiveExpMean()->change(new_lambda);
	double lnProp2 = s->getActiveVar()->change(new_sigma);
	// std::cout << lnProp1 << " " << lnProp2 << "\n";
	double lnProposalRatio = lnProp1 + lnProp2;

	// calculate the prior ratio
	double lnPriorRatio = p_lambda->lnPriorRatio() + p_sigma->lnPriorRatio();

	// calculate the likelihood ratio
	modelPtr->updateModel();
	newLnL = modelPtr->modelLnLikelihood();

	double lnLikelihoodRatio = newLnL - oldLnL;

	// accept or reject MH ratio
	double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );

	bool acceptState = false;
	double u = randomPtr->uniformRv();
	if (u < r)
	{
		acceptState = true;
	}

/*
	std::cout << "\tlnProposalRatio:\t" << lnProposalRatio << "\n";
	std::cout << "\tlnPriorRatio:\t\t" << lnPriorRatio << "\n";
	std::cout << "\toldLnL:\t\t\t" << oldLnL << "\n";
	std::cout << "\tnewLnL:\t\t\t" << newLnL << "\n";
	std::cout << "\tlnLikeRatio:\t\t" << lnLikelihoodRatio << "\n";
	std::cout << "\tMH Ratio:\t\t" << r << " (" << lnLikelihoodRatio + lnPriorRatio + lnProposalRatio << ")\n";
	std::cout << "\tu:\t\t\t" << u << "\t\t\t" << "r:\t" << r << "\n";
*/
	// update the state of the chain
	if (acceptState == true)
	{
		// TODO: add numswapped
		p_lambda->incrementNumAccepted();
		p_sigma->incrementNumAccepted();
		modelPtr->keep(p_lambda);
		modelPtr->keep(p_sigma);
		oldLnL = newLnL;
	}
	else
	{
		modelPtr->restore(p_lambda);
		modelPtr->restore(p_sigma);
	}

	proposeStr = p_lambda->getName() + " " + p_sigma->getName() + ",\tRotate : " + printBool(acceptState);

}
#endif

#if 0
void Mcmc::proposeSigmaSwap(void)
{
	std::string parmNames = "";
	double lnProposalRatio = 0.0;
	double lnPriorRatio = 0.0;

	Table* t = modelPtr->pickTableAtRandom();
	std::vector<Parm*> tempParmVector = t->getParmVector();

	if (modelType == 0 || modelType == 3) // Jump Normal, Drift BM
	{

		Sigma* sbm = dynamic_cast<Sigma*>(tempParmVector[0]);	// Sigma-BM
		Sigma* ljn = dynamic_cast<Sigma*>(tempParmVector[1]);	// Lambda-JN
		Sigma* sjn = dynamic_cast<Sigma*>(tempParmVector[2]);	// Sigma-JN

		double sbm_val = sbm->getValue();
		double ljn_val = ljn->getValue();
		double sjn_val = sjn->getValue();



		double u = randomPtr->uniformRv();
		double tuning = log(3.0);
		double v = exp(tuning * (randomPtr->uniformRv() - 0.5));



		// hold Sigma-JN constant
		if (u < .33)
		{
			modelPtr->sampleJumpsForTree();
			parmNames = sbm->getName() + " " + ljn->getName();
			sbm_val *= v;
			ljn_val *= 1 + (sbm_val * sbm_val * (1 - v * v)) / (ljn_val * sjn_val * sjn_val);
			lnProposalRatio = sbm->getActiveVar()->change(sbm_val) + ljn->getActiveVar()->change(ljn_val);
			lnPriorRatio = sbm->lnPriorRatio() + ljn->lnPriorRatio();
		}

		// hold Lambda-JN constant
		else if (u < .66)
		{
			modelPtr->sampleJumpSizesForTree();
			parmNames = sbm->getName() + " " + sjn->getName();
			sbm_val *= v;
			//sjn_val /= sqrt( 1 - (ljn_val * sjn_val * sjn_val * (v - 1)) / (sbm_val * sbm_val));
			sjn_val = sqrt((sbm_val * sbm_val * ( 1 - v * v)) / (ljn_val * sjn_val * sjn_val) + 1);
			lnProposalRatio = sbm->getActiveVar()->change(sbm_val) + sjn->getActiveVar()->change(sjn_val);
			lnPriorRatio = sbm->lnPriorRatio() + sjn->lnPriorRatio();

		}

		// hold Sigma-BM constant
		else
		{
			modelPtr->sampleJumpsForTree();
			parmNames = ljn->getName() + " " + sjn->getName();
			sjn_val *= v;
			ljn_val /= (v * v);
			lnProposalRatio = sjn->getActiveVar()->change(sjn_val) + ljn->getActiveVar()->change(ljn_val);
			lnPriorRatio = sjn->lnPriorRatio() + ljn->lnPriorRatio();
		}

	}
	else // propose a new state otherwise
	{
		proposeState();
	}

	// std::cout << lnProp1 << " " << lnProp2 << "\n";
//	double lnProposalRatio = lnProp1 + lnProp2;

	// calculate the prior ratio
//	double lnPriorRatio = p_lambda->lnPriorRatio() + p_sigma->lnPriorRatio();

	// calculate the likelihood ratio
	modelPtr->updateModel();
	newLnL = modelPtr->modelLnLikelihood();

	double lnLikelihoodRatio = newLnL - oldLnL;

	// accept or reject MH ratio
	double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );

	bool acceptState = false;
	double u = randomPtr->uniformRv();
	if (u < r)
	{
		acceptState = true;
	}

/*
	std::cout << "\tlnProposalRatio:\t" << lnProposalRatio << "\n";
	std::cout << "\tlnPriorRatio:\t\t" << lnPriorRatio << "\n";
	std::cout << "\toldLnL:\t\t\t" << oldLnL << "\n";
	std::cout << "\tnewLnL:\t\t\t" << newLnL << "\n";
	std::cout << "\tlnLikeRatio:\t\t" << lnLikelihoodRatio << "\n";
	std::cout << "\tMH Ratio:\t\t" << r << " (" << lnLikelihoodRatio + lnPriorRatio + lnProposalRatio << ")\n";
	std::cout << "\tu:\t\t\t" << u << "\t\t\t" << "r:\t" << r << "\n";
*/

	// update the state of the chain
	if (acceptState == true)
	{
		topologyPtr->flipAllActiveParms();
		for (unsigned int i = 0; i < tempParmVector.size(); i++)
		{
			tempParmVector[i]->incrementNumAccepted();
			modelPtr->keep(tempParmVector[i]);
		}
		oldLnL = newLnL;
	}
	else
	{
		for (unsigned int i = 0; i < tempParmVector.size(); i++)
		{
			modelPtr->restore(tempParmVector[i]);
		}
	}

	proposeStr = parmNames + ",\tSigSwap: " + printBool(acceptState);
}

#endif

#if 0
void Mcmc::proposeSeating(void)
{
	// 0. Initialize CRP structures

	oldLnL = modelPtr->modelLnLikelihood();
	int numRegTables = tableListPtr->size();
	int numAuxTables = auxCRP; //auxCRP * 5 + 1;					// constant
	// int numAuxTables = auxCRP * numRegTables + 1;				// proportional
	alphaCRP = modelPtr->getActiveConc()->getValue();				// Run inference on alpha.


	std::vector<double> proposalLnLs;
	std::vector<double> proposalProbs;
	std::list<Table*> auxTableList;

	//std::cout << "START: numRegTables = " << numRegTables << "\n";

	// Reseat all patrons.
	for (std::list<Patron*>::iterator it_p = patronListPtr->begin(); it_p != patronListPtr->end(); it_p++)
	{
		proposalLnLs.clear();
		proposalProbs.clear();
		auxTableList.clear();

		Patron* p = *it_p;
		// (*it_p)->print();

		Table* oldTable = p->getTable();
		//Table* newTable = NULL;
		p->stand();

		// remove table if empty
		if (oldTable->getPatronList().size() == 0)
		{
			tableListPtr->remove(oldTable);
		}
		numRegTables = tableListPtr->size();

		// 2. Seat the patron at (1, ..., K) tables and calculate likelihoods
		//std::cout << "\tnumRegTables = " << numRegTables << "\n";
		for (std::list<Table*>::iterator it_t = tableListPtr->begin(); it_t != tableListPtr->end(); it_t++)
		{
			//std::cout << "+++REGTABLE+++\n";
			//(*it_t)->print();

			// sit the patron at the next table
			p->sit(*it_t);

			// calculate the likelihoods
			proposalLnLs.push_back(modelPtr->locusLogLikelihood(p));

			// remove patron from table
			p->stand();

			// calculate the log likelihood-weighted probabilities
			double val = proposalLnLs.back() + log( (double)((*it_t)->getPatronList().size()) / (numTranscripts + alphaCRP - 1) ); // TEST CRP
			if ((double)(*it_t)->getPatronList().size() == 0)
				std::cerr << "ERROR: # at table [" << &(*it_t) << "] == 0" << "\n";

			// store the value
			proposalProbs.push_back(val);

			//std::cout << "\t\t\t" << proposalLnLs.back() << "\n";
			//std::cout << "\t\t\t" << proposalProbs.back() << "\n";

		}

		// 3. Seat patrons at (K+1, ..., K+M) auxillary tables
		//std::cout << "\tnumAuxTables = " << numAuxTables << "\n";
		for (int i = 0; i < numAuxTables; i++)
		{
			//std::cout << "+++AUXTABLE+++\n";
			// draw new auxillary parameters
			std::vector<Parm*> auxParmVector;

			// TODO: Push parameters based on modelType (consider creating a fn() which returns std::vector<Parm*>
			auxParmVector.push_back(new Sigma(randomPtr, "Lambda"));
			auxParmVector.push_back(new Sigma(randomPtr, "Sigma"));

			std::vector<Parm*> auxBranchVector;
			for (int i = 0; i < numBranches; i++)
			{
				auxBranchVector.push_back(new Tau(randomPtr, "Tau-" + printInt(i)));
			}
			//auxParmVector.push_back(new Tau(randomPtr, "Tau-L"));
			//auxParmVector.push_back(new Tau(randomPtr, "Tau-R"));

			// add the table to the auxillary table vector
			auxTableList.push_back(new Table(&auxTableList, auxParmVector, auxBranchVector, modelPtr->getTableId()));
			Table* auxTable = auxTableList.back();		// TODO: .back() is a const&, .end() is an ::iterator

			//auxTable->print();

			// sit the patron at the new table
			p->sit(auxTable);

			// calculate the likelihoods
			proposalLnLs.push_back(modelPtr->locusLogLikelihood(p));

			// remove patron from table
			p->stand();

			// calculate log probabilities (unnormalized)
			double val = proposalLnLs.back() + log( (double)((alphaCRP / numAuxTables)) / (numTranscripts + alphaCRP - 1) ); // TEST CRP

			// store the value
			proposalProbs.push_back(val);

			//std::cout << "\t\t\t" << proposalLnLs.back() << "\n";
			//std::cout << "\t\t\t" << proposalProbs.back() << "\n";
		}

		// 4. normalize proposed seating probabilities
		//std::cout << "\tNormalized proposal probs:\n";

		// scale log probabilities
		double largest = proposalProbs[0];
		for (int i = 1; i < (int)proposalProbs.size(); i++)
		{
			if (proposalProbs[i] > largest)
			{
				largest = proposalProbs[i]; // proposalProbs[i] = largest;
			}
		}
		//std::cout << "\t\t\t\tlargest:" << largest << "\n";
		for (int i = 0; i < (int)proposalProbs.size(); i++)
		{
			proposalProbs[i] = proposalProbs[i] - largest;
		//	std::cout << "\t\t\t\tlog: " << std::setprecision(8) << proposalProbs[i] << "\n";
			proposalProbs[i] = exp(proposalProbs[i]);
		//	std::cout << "\t\t\t\tdec: " << std::setprecision(8) << proposalProbs[i] << "\n";
		}

		// normalize probabilities
		double probNorm = 0.0;
		for (int i = 0; i < (int)proposalProbs.size(); i++)
		{
			if (std::isnan(proposalProbs[i]) == false)
			{
				probNorm += proposalProbs[i];
			}
		}
		for (int i = 0; i < (int)proposalProbs.size(); i++)
		{
			proposalProbs[i] = proposalProbs[i] / probNorm;
			if (proposalProbs[i] < 0.0)
			{
				std::cerr << "ERROR: normalized reseating probabilities contain at least one element < 0.0\n";
			}
			//std::cout << "\t\t\tp[" << std::setw(5) << i << "]: " << proposalProbs[i] << "\n";
		}

		// 5. determine which seating arrangement to keep
		int probIndex = -1;
		double u = randomPtr->uniformRv();
		double probSum = 0.0;
		for (int i = 0; i < (int)proposalProbs.size(); i++)
		{
			if (std::isnan(proposalProbs[i]) == false)
			{
				probSum += proposalProbs[i];
				if (u < probSum)
				{
					probIndex = i;
					//std::cout << "\t\t\tu       : " << u << "\n";
					break;
				}
			}
			else
			{
				std::cerr << "ERROR: proposalProbs[" << i << "] = NaN\n";
			}
		}

		if (probIndex == -1)
			std::cerr << "ERROR: failed to sample new seating arrangement.\n";
		else if (probIndex > numRegTables + numAuxTables)
			std::cerr << "ERROR: attempted to seat patron outside of tableVector bounds.\n";

		// std::cout << "\tPICKED: " << probIndex << "\n";

		// commit the new seating arrangement to the proposed MCMC state
		if (probIndex < numRegTables)
		{
			// std::cout << "Reg. table picked\n";

			// find the correct regular table
			std::list<Table*>::iterator it_t = tableListPtr->begin();
			for (int i = 0; i < probIndex; i++)
			{
				it_t++;
			}

			p->sit(*it_t);
			// (*it_t)->print();
			// p->print();
		}
		else
		{
			// std::cout << "Aux. table picked\n";

			// find the correct auxillary table
			std::list<Table*>::iterator it_t = auxTableList.begin();

			// over number of aux tables
			for (int i = numRegTables; i < probIndex; i++)
			{
				it_t++;
			}

			Table* newTable = new Table(*(*it_t));
			// std::cout << "\tNEWTABLE " << newTable << " " << &newTable << "\n";
			newTable->setTableListPtr(tableListPtr);
			p->sit(newTable);
			tableListPtr->push_back(newTable);

			// std::cout << "\tnewTable->print()\n";
			// newTable->print();

			// tableListPtr->push_back(new Table(*it_t, tableListPtr));
			// std::cout << "\ttableListPtr->back()->print()\n";
			// tableListPtr->back()->print();
			// p->print();

			// auxTableList.remove(*it_t);

		}

		// clear remaining auxillary tables
		for (std::list<Table*>::iterator it_a = auxTableList.begin(); it_a != auxTableList.end(); it_a++)
		{
			delete *it_a;
		}

		/*
		if (oldTable->getPatronList().size() == 0)
		{
			for (std::list<Table*>::iterator it_t = tableListPtr->begin(); it_t != tableListPtr->end(); it_t++)
			{
				if (oldTable == *it_t)
				{
					tableListPtr->erase(it_t);
					break;
				}
			}
		}
		*/
		// std::cout << "tableListPtr()->size = " << tableListPtr->size() << "\n";
		// std::cout << "\n\n";
	}



	// update the model lnL
	newLnL = modelPtr->modelLnLikelihood();
	// std::cout << "TABLE RESEATING LNL: " << newLnL << "\n";
	numRegTables = tableListPtr->size();
	// std::cout << "END:   numRegTables = " << numRegTables << "\n";

	/*
	for (std::list<Patron*>::iterator it_p = patronListPtr->begin(); it_p != patronListPtr->end(); it_p++)
	{
		(*it_p)->print();
	}

	for (std::list<Table*>::iterator it_t = tableListPtr->begin(); it_t != tableListPtr->end(); it_t++)
	{
		(*it_t)->print();
	}
	*/


	// update alpha
	//modelPtr->getActiveConc()->print();
//	modelPtr->getActiveConc()->change();
//	modelPtr->getAlpha()->keep();
	//modelPtr->getActiveConc()->print();

	// update MCMC Output string
	std::stringstream proposeStrStream;
	proposeStrStream << "Table,\t\tReseat : n_t = " << tableListPtr->size();
	proposeStrStream << "\t[";
	for (std::list<Table*>::iterator it_t = tableListPtr->begin(); it_t != tableListPtr->end(); it_t++)
	{
		proposeStrStream << " " << (*it_t)->getPatronList().size();
	}
	proposeStrStream << " ]";
	proposeStr = proposeStrStream.str();

}
#endif


#if 0
void Mcmc::proposeChangeJump(void)
{
	// choose branch
	Table* t = tableListPtr->front();
	Node* p = topologyPtr->getRandomNodeByLength();

	// calculate the proposal ratio
	//double lnProposalRatio = modelPtr->proposeRemoveJump(p, t->getParmVector(), 1);
	double lnProposalRatio = modelPtr->proposeJumpSize(p, t->getParmVector());

	// calculate the prior ratio
	double lnPriorRatio = 0.0;

	// calculate the likelihood ratio
	newLnL = modelPtr->modelLnLikelihood();
	newKb = topologyPtr->getRoot()->getKb();
	newKj = topologyPtr->getRoot()->getKj();
	double lnLikelihoodRatio = newLnL - oldLnL;
	if (useSteppingStone)
	{
		//lnLikelihoodRatio *= betaSteppingStone;
		lnLikelihoodRatio = betaSteppingStone*(newKb - oldKb) + (newKj - oldKj);
	}

	// accept or reject MH ratio
	bool acceptState = false;
	double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );
	double u = randomPtr->uniformRv();
	if (u < r)
	{
		acceptState = true;
	}

#if DEBUG2
	std::cout << "\tlnProposalRatio:\t" << lnProposalRatio << "\n";
	std::cout << "\tlnPriorRatio:\t\t" << lnPriorRatio << "\n";
	std::cout << "\toldLnL:\t\t\t" << oldLnL << "\n";
	std::cout << "\tnewLnL:\t\t\t" << newLnL << "\n";
	std::cout << "\tlnLikeRatio:\t\t" << lnLikelihoodRatio << "\n";
	std::cout << "\tMH Ratio:\t\t" << r << " (" << lnLikelihoodRatio + lnPriorRatio + lnProposalRatio << ")\n";
	std::cout << "\tu:\t\t\t" << u << "\t\t\t" << "r:\t" << r << "\n";
	#endif

	// update the state of the chain
	if (acceptState == true)
	{
		p->copySpace(0, 1);
		//topologyPtr->decrementNumJumps();
		oldLnL = newLnL;
		oldKb = newKb;
		oldKj = newKj;
	}
	else
	{
		p->copySpace(1, 0);
	}

	proposeStr = "n" + printInt(p->getIndex()) + "\t\tChngJmp: " + printBool(acceptState);
}
#endif
/*
void Mcmc::printTableState(int n)
{
	// HEADER
	// i.e.
	//	CYCLE
	//		TABLE
	//			GENE	PARAMS

	if (n == 0)
	{
		std::string pHeaderStr = "";
		pHeaderStr += "cycle\ttable\tname\t";
		Table* firstTable = modelPtr->getTableListPtr()->front();
		for (std::vector<Parm*>::const_iterator it_p = firstTable->getParmVector().begin(); it_p != firstTable->getParmVector().end(); it_p++)
		{
			pHeaderStr += (*it_p)->getParameterHeader();
		}
	//	pHeaderStr += modelPtr->getAlpha()->getParameterHeader();
	//	pHeaderStr += "\n";
		tableFileStrm << pHeaderStr;
	}


	// BODY

	tableFileStrm << n << "\n";
	for (std::list<Table*>::iterator it_t = tableListPtr->begin(); it_t != tableListPtr->end(); it_t++)
	{
		std::string pStr = "";
		tableFileStrm << "\t" << printInt((*it_t)->getId()) << "\n";
		for (std::vector<Parm*>::const_iterator it_p = (*it_t)->getParmVector().begin(); it_p != (*it_t)->getParmVector().end(); it_p++)
		{
			pStr += (*it_p)->getParameterStr();
		}
		//pStr += modelPtr->getAlpha()->getParameterStr() + "\n";

		for (std::list<Patron*>::const_iterator it_i = (*it_t)->getPatronList().begin(); it_i != (*it_t)->getPatronList().end(); it_i++)
		{
			tableFileStrm << "\t\t" << (*it_i)->getName() << "\t" << pStr;
		}
	}
	tableFileStrm << "\n";

	// parmFileStrm << n << '\t' << std::fixed << std::setprecision(2) << lnL << '\t' << pStr << std::endl;

}

*/

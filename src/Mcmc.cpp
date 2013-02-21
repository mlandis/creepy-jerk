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

	oldLnL = modelPtr->modelLnLikelihood(1);
	oldKj = topologyPtr->getRoot()->getKj(1);
	oldKb = topologyPtr->getRoot()->getKb(1);


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
	std::cout << "\toldKb:\t\t\t" << oldKb << "\n";
	std::cout << "\tnewKb:\t\t\t" << newKb << "\n";
	std::cout << "\toldKj:\t\t\t" << oldKj << "\n";
	std::cout << "\tnewKj:\t\t\t" << newKj << "\n";
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

		#if DEBUG_PROPOSE_JUMP
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

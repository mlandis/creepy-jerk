/*
 * Mcmc.cpp
 *
 *  Created on: Dec 16, 2010
 *
 */

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

	//numCycles = 1000000;
	//numCycles = 10000;
	//numCycles = 0; //1000000;
	//printFreqMH = 1;
	//printFreqCRP = 50;

	numTranscripts = expressionPtr->getNumTranscripts();
	numTimepoints = expressionPtr->getNumTimepoints();
	numTaxa = expressionPtr->getNumTaxa();
	numNodes = 2 * numTaxa - 1;
	numBranches = numNodes - 1;
	numParms = 2;
	alphaCRP = settingsPtr->getAlphaCRP();
	auxCRP = settingsPtr->getAuxCRP();
	useCRP = settingsPtr->getUseCRP();
	// tuning = settingsPtr->getTuning();

	oldLnL = 0.0;
	newLnL = 0.0;
	proposeStr = "";

	modelType = settingsPtr->getModelType();
	fixBranches = settingsPtr->getFixBranches();

	double sum = 0.0;
	proposalProbs.push_back( 1.0 ); // propose change to a random parameter
	if (fixBranches)
		proposalProbs.push_back( 0.0 ); // propose change to a random branch
	else if (!fixBranches)
		proposalProbs.push_back( 1.0 );
	proposalProbs.push_back( 1.0 ); // add branch
	proposalProbs.push_back( 1.0 ); // remove branch
	proposalProbs.push_back( 0.0 ); // swap two neighboring branches and swap sign of a
	proposalProbs.push_back( 0.0 ); // swap two neighboring branches and swap sign of a
	proposalProbs.push_back( 0.0 ); // rescaling the branch lengths, holding total variance, changing number of mutations
	proposalProbs.push_back( 0.0 ); // rotation on a-sigma space
	proposalProbs.push_back( 0.0 );	// propose sigma swap
	proposalProbs.push_back( 0.0 ); // propose jump size resampling

	if (useCRP)
		proposalProbs.push_back( 1.0 ); // propose reseating arrangement
	else if (!useCRP)
		proposalProbs.push_back( 0.0 ); // no reseating

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
	std::cout << "MCMC: Run chain.\n";

	// initialize the chain by keeping the first chain state
	if (useCRP)
	{
		proposeSeating();
		std::cout << "MCMC: Initial seating complete.\n";
		modelPtr->printTables();
	}

	// initialize jump samples
	modelPtr->sampleJumpsForTree();
	topologyPtr->printJumpSamples();
	topologyPtr->flipAllActiveParms();
	topologyPtr->copyNodeSpaces();
	//topologyPtr->copyInactiveJumpsToActiveJumps();

	oldLnL = modelPtr->modelLogLikelihood();

	int sumTables = 0;
	int postBurnN = 0;
	int burnIn = 1000;

	// get initial likelihood
	for (int n = 0; n <= numCycles; n++)
	{
		double printOldLnL = oldLnL;

		proposeState();

		// print information to the screen
		if ( n % printFreqMH == 0 )
		{
			std::cout << std::setw(5) << n << " -- ";
			std::cout << std::fixed << std::setprecision(8) << printOldLnL << " -> " << newLnL << "\t" << proposeStr << "\n";
			//std::cout << "\n";
			/*
			if (acceptState == true)
				std::cout << " -- Accepted change\n";
			else
				std::cout << " -- Rejected change\n";
			*/

			// print to file
			//if (newLnL > 1.0 || n % 10000 == 0)
			//topologyPtr->printJumpSamples();

			printChainState(n, oldLnL);
			if (useCRP)
				printTableState(n);
		}


		if (n % (printFreqCRP) == 0 && useCRP)
		{
			//if (n > burnIn)
			//	std::cout << "*\t*\t*\tmean # tables (" << n << "):" << ((double)sumTables/(double)postBurnN) << " for a=" << alphaCRP << ", n=" << numTranscripts << "\n";
			//std::cout << "\n";
			// modelPtr->printQ();
			//std::cout << "\n";
			modelPtr->printTables();
			modelPtr->getAlpha()->print();
			std::cout << "Simulation Name: " << settingsPtr->getSimName() << "\n";
			std::cout << "Output Name:     " << settingsPtr->getOutputFileName() << "\n";

			// print to file
		}


		if (n > burnIn)
		{
			sumTables += tableListPtr->size();

			postBurnN++;
		}


#if 0
		// oldLnL = newLnL;
		/*
		// update seating via CRP
		// TODO: oldLnL is updated as part of this step
		// TODO: oldTable is updated as part of this step
		newTable = reseatPatron();
		Patron* reseatedPatron = patronListPtr->back();

		// select a parameter to change
		// std::cout << "picked table:\n";
		// table->print();
		// std::cout << "pick parm:\n";
		parm = modelPtr->pickParmAtRandom(newTable);
		// parm->print();

		// propose new parameters
		double lnProposalRatio = parm->change();
		// std::cout << "lnProposalRatio: " << lnProposalRatio << "\n";

		// calculate the prior ratio
		double lnPriorRatio = parm->lnPriorRatio();
		// std::cout << "lnPriorRatio: " << lnPriorRatio << "\n";

		// calculate the likelihood ratio
		modelPtr->updateModel();
		modelPtr->updateStationaryProbs();
		newLnL = modelPtr->modelLogLikelihood();
		double lnLikelihoodRatio = newLnL - oldLnL;
		// std::cout << "lnLikelihoodRatio: " << lnLikelihoodRatio << "\n";

		// calculate likelihood ratio
		double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );
		*/

		// generate a likelihood ratio given a proposed state

		// accept or reject the proposed state as the next state of the chain

		/*
		// print out current chain state
		if ( n == 1 || n % sampleFrequency == 0 || n == numCycles )
			printChainState(n, oldLnL);

		// print summary of chain
		if ( n % summarizeFrequency == 0 )
			treeSummary->printSummary();
		}
		*/

		// treeSummary->print();
		// modelPtr->printAcceptanceInfo();

		//saveCycle();
		//saveAcceptance()

#endif

	}
	//saveChain();
	topologyPtr->printJumpSamples();

	std::cout << "MCMC: Complete.\n";
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
		std::cout << "\tPropose add jump\n";
		//proposeBranchSwap();
		proposeAddJump();

	}

	else if (uIndex == 3)
	{
		std::cout << "\tPropose remove jump\n";
		proposeRemoveJump();
		//proposeBranchRescale();
	}

	else if (uIndex == 4)
	{
		//std::cout << "\tPropose parameter rotation\n";
		proposeParmRotate();
	}

	else if (uIndex == 5)
	{
		// std::cout << "\tPropose sigma swap\n";
		proposeSigmaSwap();
	}

	else if (uIndex == 5)
	{
		//std::cout << "\tPropose seating\n";
		proposeSeating();
	}
	else if (uIndex == 6)
	{
		//std::cout << "\tResize all jumps, hold number of jumps constant\n";

		// NOT IMPLEMENTED
		proposeResizeJumps();
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

	std::cout << "\n";
	std::cout << parm->getName() << "\t" << oldVal << " -> " << newVal << "\n";
//	topologyPtr->printJumpSummary();

	// calculate the likelihood ratio
	// modelPtr->updateModel();
	newLnL = modelPtr->modelLogLikelihood();
	double lnLikelihoodRatio = newLnL - oldLnL;

	// accept or reject MH ratio
	bool acceptState = false;
	double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );
	double u = randomPtr->uniformRv();
	if (u < r) // && newLnL <= 1.0)
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
		parm->incrementNumAccepted();
		modelPtr->keep(parm);
		oldLnL = newLnL;
	}
	else
	{
		modelPtr->restore(parm);
	}

	proposeStr = parm->getName() + ",\t\tUpdate : " + printBool(acceptState);
}

void Mcmc::proposeAddJump(void)
{
	// choose branch
	//Table* t = modelPtr->pickTableAtRandom();
	Table* t = tableListPtr->front();
	double lam_jn = (t->getParmVector()[1])->getValue();
	double sig_jn = (t->getParmVector()[2])->getValue();

	int space = topologyPtr->getRoot()->getActiveParm();
	std::cout << "sanity check\t" << space << "\n";

	// calculate the proposal ratio
	double lnProposalRatio = topologyPtr->proposeAddJump(lam_jn, sig_jn, space);

	// calculate the prior ratio
	double lnPriorRatio = 0.0; // = parm->lnPriorRatio();

	// calculate the likelihood ratio
	// modelPtr->updateModel();
	newLnL = modelPtr->modelLogLikelihood();
	double lnLikelihoodRatio = newLnL - oldLnL;

	// accept or reject MH ratio
	bool acceptState = false;
	double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );
	double u = randomPtr->uniformRv();
	if (u < r) // && newLnL <= 1.0)
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
		topologyPtr->incNumJumps();
		//modelPtr->keepJumps();
		oldLnL = newLnL;
	}
	else
	{
		; // do nothing, I think...
	}

	proposeStr = "\t\t\t\tAdd jmp: " + printBool(acceptState);
}

void Mcmc::proposeRemoveJump(void)
{
	std::cout << "\t\tREM JUMP\n";

	// choose branch
	//Table* t = modelPtr->pickTableAtRandom();
	Table* t = tableListPtr->front();
	double lam_jn = (t->getParmVector()[1])->getValue();
	double sig_jn = (t->getParmVector()[2])->getValue();

	int space = topologyPtr->getRoot()->getActiveParm();

	// calculate the proposal ratio
	double lnProposalRatio = topologyPtr->proposeRemoveJump(lam_jn, sig_jn, space);

	// calculate the prior ratio
	double lnPriorRatio = 0.0; // = parm->lnPriorRatio();

	// std::cout << "\n";
	// std::cout << parm->getName() << "\t" << oldVal << " -> " << newVal << "\n";
	//	topologyPtr->printJumpSummary();

	// calculate the likelihood ratio
	// modelPtr->updateModel();
	newLnL = modelPtr->modelLogLikelihood();
	double lnLikelihoodRatio = newLnL - oldLnL;

	// accept or reject MH ratio
	bool acceptState = false;
	double r = safeExp( lnLikelihoodRatio + lnPriorRatio + lnProposalRatio );
	double u = randomPtr->uniformRv();
	if (u < r) // && newLnL <= 1.0)
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
		topologyPtr->decNumJumps();
		oldLnL = newLnL;
	}
	else
	{
		; // do nothing, I think...
	}

	proposeStr = "\t\t\t\tRem jmp: " + printBool(acceptState);
}

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
	newLnL = modelPtr->modelLogLikelihood();

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
	newLnL = modelPtr->modelLogLikelihood();

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
	newLnL = modelPtr->modelLogLikelihood();

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
	newLnL = modelPtr->modelLogLikelihood();

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
	newLnL = modelPtr->modelLogLikelihood();
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


void Mcmc::proposeSeating(void)
{
	// 0. Initialize CRP structures

	oldLnL = modelPtr->modelLogLikelihood();
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
	newLnL = modelPtr->modelLogLikelihood();
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
	modelPtr->getActiveConc()->change();
	modelPtr->getAlpha()->keep();
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

void Mcmc::openFiles(std::string fn) {

	std::string tf = fn + ".t";
	std::string pf = fn + ".p";
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
}

void Mcmc::printChainState(int n, double lnL) {

	if (n == 0)
	{
		std::string pHeaderStr = "";
		std::string tHeaderStr = "";

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
		parmFileStrm << "Cycle\tlnL\t" << pHeaderStr << std::endl;

		/*
		for (int i = 0; i < modelPtr->getNumParameters(); i++)
		{
			Parm* p = modelPtr->getParameter(i);

			Tree* derivedPtr = dynamic_cast<Tree*> (p);
			if (derivedPtr != 0)
				tHeaderStr += p->getParameterHeader();
			else

			pHeaderStr += p->getParameterHeader();
		}

		//treeFileStrm << tHeaderStr;
		parmFileStrm << "Cycle\tlnL\t" << pHeaderStr << std::endl;
		*/
	}

	std::string pStr = "";
	std::string tStr = "";

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


	/*
	for (int i = 0; i < modelPtr->getNumParameters(); i++)
	{
		Parm* p = modelPtr->getParameter(i);
		Tree* derivedPtr = dynamic_cast<Tree*> (p);
		if (derivedPtr != 0)
			tStr += p->getParameterStr();
		else
			pStr += p->getParameterStr();
	}
	//treeFileStrm << "   tree_" << n << " = " << tStr << ";" << std::endl;
	parmFileStrm << n << '\t' << std::fixed << std::setprecision(2) << lnL << '\t' << pStr << std::endl;
	*/


	if (n == numCycles)
	{
		std::cout << "Mcmc::runChain() complete!\n";
		//treeFileStrm << "end;" << std::endl;
	}

}

void Mcmc::printPatronState(int n, Patron* p)
{

	/*
	 *
	std::string pHeaderStr = "";
	std::string tHeaderStr = "";


	for (std::vector<Parm*>::const_iterator it_p = (*it_t)->getParmVector().begin(); it_p != (*it_t)->getParmVector().end(); it_p++)
	{
		pHeaderStr += (*it_p)->getParameterHeader();
	}

	*/
}

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
		pHeaderStr += modelPtr->getAlpha()->getParameterHeader();
		pHeaderStr += "\n";
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
		pStr += modelPtr->getAlpha()->getParameterStr() + "\n";

		for (std::list<Patron*>::const_iterator it_i = (*it_t)->getPatronList().begin(); it_i != (*it_t)->getPatronList().end(); it_i++)
		{
			tableFileStrm << "\t\t" << (*it_i)->getName() << "\t" << pStr;
		}
	}
	tableFileStrm << "\n";

	// parmFileStrm << n << '\t' << std::fixed << std::setprecision(2) << lnL << '\t' << pStr << std::endl;

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
		proposalLnLs.push_back(modelPtr->modelLogLikelihood());
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
		proposalLnLs.push_back(modelPtr->modelLogLikelihood());
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

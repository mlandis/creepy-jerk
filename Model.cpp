/*
 * Model.cpp
 *
 *  Created on: Mar 7, 2011
 *      Author: mlandis
 */

#include "Model.h"

#define DEBUG 0
#define DEBUG2 0

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

Model::Model(Expression* ep, MbRandom* rp, Settings* sp, Topology* tp)
{

	expressionPtr = ep;
	randomPtr = rp;
	settingsPtr = sp;
	topologyPtr = tp;

	numTaxa = expressionPtr->getNumTaxa();
	numNodes = 2 * numTaxa - 1;
	numBranches = numNodes - 1;	// [numBranches+1 == numNodes] indexes the root which has no branch.
	numTranscripts = expressionPtr->getNumTranscripts();
	numTimepoints = expressionPtr->getNumTimepoints();

	useFFT = settingsPtr->getUseFFT();
	useCRP = settingsPtr->getUseCRP();
	modelType = settingsPtr->getModelType();
	fixBranches = settingsPtr->getFixBranches();

	tableId = 0;


	initializeParms();
	initializeTips();
	//initializeNodes();

	if (useFFT)
	{
		initializeFFT();
		condLikes = new CondLikes(expressionPtr, settingsPtr);
	}


	//condLikes->print();

	//initializeCondLikes();

}

Model::~Model(void) {
	// TODO Auto-generated destructor stub

	for (std::list<Table*>::iterator it_t = tableList.begin(); it_t != tableList.end(); it_t++)
	{
		delete (*it_t);
	}
	for (std::list<Patron*>::iterator it_p = patronList.begin(); it_p != patronList.end(); it_p++)
	{
		delete (*it_p);
	}

	if (useFFT)
	{
		delete [] theta;
		delete [] likeVals;
		delete [] like;
		delete [] tip1;
		delete [] tip2;
		delete condLikes;
	}

	delete alpha;
}


void Model::initializeFFT(void)
{
	std::cout << "INITIALIZING: FFT\n";
	// init FFT variables
	// NOTE: Can do something clever w/r/t fitting numSteps to datapoint.

	numSteps = settingsPtr->getNumSteps();
	halfSteps = numSteps / 2;
	finalStep = settingsPtr->getFinalStep();
	startStep = settingsPtr->getStartStep();
	stepSize = (finalStep - startStep) / (numSteps - 1);
	ReStepSize = PI / finalStep;

#if 0
	// int oneFFTSize = numSteps * 2;
	// int oneLocusSize = numTimepoints * oneFFTSize;
	// int oneNodeSize = numTranscripts * oneLocusSize;
	// int numNodes = 1;
	// int oneClSize = numNodes * oneNodeSize; // NOTE: For general topologies.
	// std::cout << "oneClSize:" << oneClSize << "\n";
	// int oneClSize =
#endif

	// init data frames
	likeVals = new double[numSteps];	// TODO: Why does this need a correction of +1?
	theta = new double[numSteps];
	like = new double[2*numSteps];	// CL
	tip1 = new double[2*numSteps];	// CL
	tip2 = new double[2*numSteps];	// CL

	// FFT workspace
	for (int i = 0; i < numSteps; i++) // oneClSize is multiplied by spaces=2 and divided by (Re,Im)=2
	{
		REAL(like,i) = 0.0;
		IMAG(like,i) = 0.0;
		REAL(tip1,i) = 0.0;
		IMAG(tip1,i) = 0.0;
		REAL(tip2,i) = 0.0;
		IMAG(tip2,i) = 0.0;
		// std::cout << "i:\t" << i << " " << REAL(like,i) << " " << IMAG(like,i) << "\n";
	}

	for (int i = 0; i < numSteps; i++)
	{
		theta[i] = (((i + halfSteps) % numSteps) - halfSteps) * stepSize;
	}

	std::cout << "\n";
}

void Model::initializeParms(void)
{
	std::cout << "INITIALIZING: Parameter Classes (Tables)\n";

	if (modelType == 0) // Poisson-Normal + BM
	{
		proposalProbs.push_back(1.5);			// Sigma-BM
		if (fixBranches)
			proposalProbs.push_back(1.0);		// Lambda-JN
		else if (!fixBranches)
			proposalProbs.push_back(0.0);		// Lambda-JN
		proposalProbs.push_back(1.0);			// Sigma-JN

	}
	else if (modelType == 1) // Poisson-Cauchy + BM
	{
		proposalProbs.push_back(1.0);			// Sigma-BM
		proposalProbs.push_back(1.0);			// Lambda-JC
		proposalProbs.push_back(1.0);			// Gamma-JC

	}
	else if (modelType == 2) // Poisson-Exponential + BM
	{
		proposalProbs.push_back(1.0);			// Sigma-BM
		proposalProbs.push_back(1.0);			// Lambda-JE
		proposalProbs.push_back(1.0);			// a-JE

	}
	else if (modelType == 3) // Sampled Poisson-Normal + BM
	{
		proposalProbs.push_back(1.0);			// Sigma-BM
		proposalProbs.push_back(1.0);			// Lambda-JN
		proposalProbs.push_back(1.0);			// Sigma-JN
	}
	else if (modelType == 4) // Chaix et al
	{
		proposalProbs.push_back(1.0);			// Sigma-N
		proposalProbs.push_back(1.0);			// Lambda-N
	}
	else if (modelType == 5) // Brownian motion
	{
		proposalProbs.push_back(1.0);			// Sigma-BM
	}
	else if (modelType == 6) // Pure Normal jump process
	{
		proposalProbs.push_back(1.0);			// Lambda-JN
		proposalProbs.push_back(1.0);			// Sigma-JN
	}
	else
	{
		std::cerr << "Unrecognized modelType value: " << modelType << "\n";
		exit(1);
	}

	double sum = 0.0;
	for (unsigned int i = 0; i < proposalProbs.size(); i++)
		sum += proposalProbs[i];
	sum = 1 / sum;
	for (unsigned int i = 0; i < proposalProbs.size(); i++)
		proposalProbs[i] *= sum;

	tableList = expressionPtr->getTableList();
	patronList = expressionPtr->getPatronList();

	if (!useCRP)
	{
		// create initial parameters
		std::vector<Parm*> initialParms;

		if (modelType == 0) // Poisson-Normal + BM
		{
			initialParms.push_back(new Sigma(randomPtr, "Sig-BM"));
			if (fixBranches)
				initialParms.push_back(new Sigma(randomPtr, "Lam-JN"));
			else if (!fixBranches)
				initialParms.push_back(new Sigma(randomPtr, "Lam-JN", 1.0));
			initialParms.push_back(new Sigma(randomPtr, "Sig-JN"));

		}
		else if (modelType == 1) // Poisson-Cauchy + BM
		{
			initialParms.push_back(new Sigma(randomPtr, "Sigma-BM"));
			initialParms.push_back(new Sigma(randomPtr, "Lambda-K"));
			initialParms.push_back(new Sigma(randomPtr, "Gamma-K"));

		}
		else if (modelType == 2) // Poisson-Exponential + BM
		{
			initialParms.push_back(new Sigma(randomPtr, "Sigma-BM"));
			initialParms.push_back(new Sigma(randomPtr, "Lambda-K"));
			initialParms.push_back(new Sigma(randomPtr, "a-K"));

		}
		else if (modelType == 3) // Sampled Poisson-Normal + BM
		{
			initialParms.push_back(new Sigma(randomPtr, "Sig-BM"));
			initialParms.push_back(new Sigma(randomPtr, "Lam-JN"));
			initialParms.push_back(new Sigma(randomPtr, "Sig-JN"));
		}
		else if (modelType == 4) // Chaix et al.
		{
			initialParms.push_back(new Sigma(randomPtr, "Sigma-N"));
			initialParms.push_back(new Lambda(randomPtr, "Lambda-N"));
		}
		else if (modelType == 5) // Brownian motion
		{
			initialParms.push_back(new Sigma(randomPtr, "Sig-BM"));
		}
		else if (modelType == 6) // Pure Normal jump
		{
			initialParms.push_back(new Sigma(randomPtr, "Lam-JN"));
			initialParms.push_back(new Sigma(randomPtr, "Sig-JN"));
		}
		else
		{
			std::cerr << "Unrecognized modelType value: " << modelType << "\n";
			exit(1);
		}

		// create initial branches
		std::vector<Parm*> initialBranches;
		for (int i = 0; i < numBranches; i++)
		{
			std::stringstream ss;
			ss << i;
			if (fixBranches)
			{
				double brlen = topologyPtr->getDownPassNode(i)->getV();
				initialBranches.push_back(new Tau(randomPtr, "Tau-" + ss.str(), brlen));
			}
			else if (!fixBranches)
			{
				initialBranches.push_back(new Tau(randomPtr, "Tau-" + ss.str()));
			}
		}

		/*
		if (fixBranches)
		{
			std::cout << initialBranches.size() << "\n";
			for (int j = 0; j < numNodes; j++)
			{
				Node* n = topologyPtr->getDownPassNode(j);

				if (n != NULL && n->getLft() != NULL && n->getRht() != NULL)
				{
					std::cout << "\t\tv:\t" << n->getLft()->getV() << " " << n->getRht()->getV() << "\n";
					std::cout << "\t\tn:\t" << n->getLft()->getIndex() << " " << n->getRht()->getIndex() << "\n";
					// Branch lengths depend on which node's likelihood is being calculated.
					dynamic_cast<BrLen*>(initialBranches[n->getLft()->getIndex()])->change(n->getLft()->getV());
					dynamic_cast<BrLen*>(initialBranches[n->getRht()->getIndex()])->change(n->getRht()->getV());
					std::cout << "ok\n";
				}
			}
		}
		*/

		// Uncomment this code to initialize parameters to values.
		/*
		Lambda* derivedPtrL = dynamic_cast<Lambda*> (initialParms[0]);
		Sigma* derivedPtrS = dynamic_cast<Sigma*> (initialParms[1]);
		Tau* derivedPtrT1 = dynamic_cast<Tau*> (initialParms[2]);
		Tau* derivedPtrT2 = dynamic_cast<Tau*> (initialParms[3]);

		derivedPtrL->getActiveExpMean()->change(-2.0);
		derivedPtrS->getActiveVar()->change(0.1);
		derivedPtrT1->getActiveBrLen()->change(5.0);
		derivedPtrT2->getActiveBrLen()->change(0.5);

		for (unsigned int i = 0; i < initialParms.size(); i++)
		{
			initialParms[i]->keep();
			initialParms[i]->print();
		}
		*/

		tableList.front()->setParmVector(initialParms);
		tableList.front()->setBranchVector(initialBranches);
	}


	double a_alpha = settingsPtr->getACRP();
	double b_alpha = settingsPtr->getBCRP();
	alpha = new Alpha(randomPtr, this, a_alpha, b_alpha, "Alpha");

	std::cout << "\n";
}

void Model::initializeTips(void)
{
	std::cout << "INITIALIZING: Node samples\n";

	std::list<Patron*>::const_iterator it_p = expressionPtr->getPatronList().begin();

	for (int j = 0; j < numTaxa; j++)
	{
		Node* n = topologyPtr->getNode(j);
		//Node* n = topologyPtr->getDownPassNode(j);
		//double x = expressionPtr->getExpr(0,j,0);
		double x = (*it_p)->getData(j, 0);
		std::cout << n->getIndex() << "\t" << x << "\n";
		n->setMu(x);
		n->setSigma(0.001);
		n->setK(0.0);
	}

	/*
	std::cout << "--\n";
	for (int j = 0; j < numNodes; j++)
	{
		Node* n = topologyPtr->getDownPassNode(j);
		//double x = expressionPtr->getExpr(0,j,0);
		std::cout << n->getIndex() << "\t" << n->getMu() << "\t" << n->getSigma() << "\t" << n->getK() << "\n";
	}
	*/
}

void Model::updateModel(void)
{
	// perhaps unnecessary?
}


void Model::sampleJumpsForTree(void)
{
	Table* t = tableList.front();
	double lam_jn = (t->getParmVector()[1])->getValue();
	double sig_jn = (t->getParmVector()[2])->getValue();
	int numJumps = 0;
	Node* p = NULL;
	int space = 0;

	for (int j = 0; j < numNodes; j++)
	{
		sampleJumpsForBranch(j, lam_jn, sig_jn);
		p = topologyPtr->getDownPassNode(j);
		space = p->getActiveParm();
		numJumps += p->getJumpCount(space);
	}

	topologyPtr->setNumJumps(numJumps);
//	topologyPtr->printJumpSamples();

}

void Model::sampleJumpSizesForTree(void)
{
	Table* t = tableList.front();
	//double lam_jn = (t->getParmVector()[1])->getValue();
	double sig_jn = (t->getParmVector()[2])->getValue();

	for (int j = 0; j < numNodes; j++)
	{
		sampleJumpSizesForBranch(j, sig_jn);
	}

//	topologyPtr->printJumpSamples();
}

void Model::sampleJumpsForBranch(int j, double lam_jn, double sig_jn)
{
	Node* n = topologyPtr->getDownPassNode(j);
	if (n->getAnc() != NULL)
	{
		// jump count
		double brLen = n->getV();
		int jumpCount = randomPtr->poissonRv(lam_jn * brLen);
		double lnProbJumpCount = log(randomPtr->poissonProb(lam_jn * brLen, jumpCount));
		n->setJumpCount(jumpCount);
		n->setLnProbJumpCount(lnProbJumpCount);

		// jump sizes
		std::vector<double> jumpSize, lnProbJumpSize;
		double sumJumpSize = 0.0, sumLnProbJumpSize = 0.0;
		double oneJump = 0.0;
		double lnProbOneJump = 0.0;
		for (int i = 0; i < jumpCount; i++)
		{
			oneJump = randomPtr->normalRv(0.0, sig_jn);
			lnProbOneJump = log(randomPtr->normalPdf(0.0, sig_jn, oneJump));
			jumpSize.push_back(oneJump);
			lnProbJumpSize.push_back(lnProbOneJump);

			sumJumpSize += oneJump;
			sumLnProbJumpSize += lnProbOneJump;
		}

		n->setJumpSize(jumpSize);
		n->setLnProbJumpSize(lnProbJumpSize);

		n->setSumJumpSize(sumJumpSize);
		n->setSumLnProbJumpSize(sumLnProbJumpSize);


		// update tip K values
		if (n->getLft() == NULL && n->getRht() == NULL)
		{
			n->setK(sumLnProbJumpSize + lnProbJumpCount);
		}
#if DEBUG2
		//std::cout << n->getIndex() << "\tK:\t" << n->getK()<< " = " << lnProbJumpCount << " + " << sumLnProbJumpSize << "\n";
#endif
	}
}

void Model::sampleJumpSizesForBranch(int j, double sig_jn)
{
	Node* n = topologyPtr->getDownPassNode(j);
	if (n->getAnc() != NULL)
	{
		int inactiveParm = (n->getActiveParm() == 0 ? 1 : 0); // opposite of activeParm

		// jump count is held constant
		int jumpCount = n->getJumpCount(inactiveParm);
		double lnProbJumpCount = n->getLnProbJumpCount(inactiveParm);
		n->setJumpCount(jumpCount);
		n->setLnProbJumpCount(lnProbJumpCount);

		//std::cout << "TEST:\t" << jumpCount << "\t" << lnProbJumpCount << "\n";

		// jump sizes
		std::vector<double> jumpSize, lnProbJumpSize;
		double sumJumpSize = 0.0, sumLnProbJumpSize = 0.0;
		double oneJump = 0.0, lnProbOneJump = 0.0;

		for (int i = 0; i < jumpCount; i++)
		{
			oneJump = randomPtr->normalRv(0.0, sig_jn);
			lnProbOneJump = log(randomPtr->normalPdf(0.0, sig_jn, oneJump));
			//std::cout << "\tlnProbOneJump:\t" << lnProbOneJump << "\n";
			jumpSize.push_back(oneJump);
			lnProbJumpSize.push_back(lnProbOneJump);

			sumJumpSize += oneJump;
			sumLnProbJumpSize += lnProbOneJump;
		}

		n->setJumpSize(jumpSize);
		n->setLnProbJumpSize(lnProbJumpSize);

		n->setSumJumpSize(sumJumpSize);
		n->setSumLnProbJumpSize(sumLnProbJumpSize);

		if (n->getLft() == NULL && n->getRht() == NULL)
		{
			n->setK(sumLnProbJumpSize + lnProbJumpCount);
		}
#if DEBUG2
		//std::cout << n->getIndex() << "\tK:\t" << n->getK()<< " = " << lnProbJumpCount << " + " << sumLnProbJumpSize << "\n";
#endif
	}
}

double Model::getProposalRatioLambda(double oldLambda, double newLambda)
{
	double lnProposalRatio = 1.0;

	int jumpCount = 0;
	int activeParm;

	// get proposal ratio for new samples
	for (int j = 0; j < numNodes; j++)
	{
		Node* n = topologyPtr->getDownPassNode(j);
		activeParm = n->getActiveParm();

		if (n->getAnc() != NULL)
		{
			jumpCount = n->getJumpCount(activeParm);
			lnProposalRatio += log(randomPtr->poissonProb(oldLambda, jumpCount)) - log(randomPtr->poissonProb(newLambda, jumpCount));
		}
	}

	return lnProposalRatio;
}

double Model::getProposalRatioSigma(double oldSigma, double newSigma)
{
	double lnProposalRatio = 1.0;

	std::vector<double> jumpSize;
	std::vector<double>::iterator it_i;
	int activeParm;

	// get proposal ratio for new samples
	for (int j = 0; j < numNodes; j++)
	{
		Node* n = topologyPtr->getDownPassNode(j);
		activeParm = n->getActiveParm();

		if (n->getAnc() != NULL)
		{
			jumpSize = n->getJumpSize(activeParm);
			for (it_i = jumpSize.begin(); it_i < jumpSize.end(); ++it_i)
			{
				lnProposalRatio += log(randomPtr->normalPdf(0.0, oldSigma, (*it_i))) - log(randomPtr->normalPdf(0.0, newSigma, (*it_i)));
			}
		}
	}

	return lnProposalRatio;

}


double Model::modelLogLikelihood(void)
{
	double lnL = 0.0;

	if (!useCRP)
	{
		if (useFFT)
		{
			lnL = tableLogLikelihoodFFT(tableList.front());
		}
		else if (!useFFT)
		{
			lnL = tableLogLikelihood(tableList.front());
		}
	}
	else
	{
		for (std::list<Table*>::iterator it_t = tableList.begin(); it_t != tableList.end(); it_t++)
		{
			double tableLnL = 0.0;
			if (useFFT)
			{
				tableLnL = tableLogLikelihoodFFT(*it_t);
			}
			else if (!useFFT)
			{
				tableLnL = tableLogLikelihood(*it_t);
			}
			//	(*it_t)->setUpdateLnL(false);
			//	(*it_t)->setLnL(tableLnL);
			lnL += tableLnL;
		}
	}

	return lnL;
}

double Model::tableLogLikelihood(Table* t)
{
	double lnL = 0.0;

	const std::vector<Parm*> tempParmVector = t->getParmVector();
	const std::vector<Parm*> tempBranchVector = t->getBranchVector();

	// Collect parameter & branch values from the Table's parameter class
	std::vector<double> parmVals;
	std::vector<double> branchVals;
	for (unsigned int i = 0; i < tempParmVector.size(); i++)
	{
		parmVals.push_back(tempParmVector[i]->getValue());
	}
	for (unsigned int i = 0; i < tempBranchVector.size(); i++)
	{
		branchVals.push_back(tempBranchVector[i]->getValue());
	}

	// Model's inferred sigma_BM
	double sigmaBM = parmVals[0];
	double sigmaBM2 = sigmaBM * sigmaBM;
	double lambdaJN = parmVals[1];
	double sigmaJN = parmVals[2];


#if DEBUG2
	// Calculate likelihoods for single Patron (labeled to a different Parm class)
	std::cout << std::setprecision(4) << std::endl;
	std::cout << "L\t\t\t\tR\t\t\t\tP\n";
	std::cout << "i\tmu\tsigma\tK\tj\tmu\tsigma\tK\tp\tmu\tsigma\tK\n";
#endif


	// update K values at tips according to sampled values
	for (int j = 0; j < numTaxa; j++)
	{
		Node* n = topologyPtr->getNode(j);
		int activeParm = n->getActiveParm();
	//	int inactiveParm = (n->getActiveParm() == 1 ? 0 : 1);
		double x = n->getLnProbJumpCount(activeParm) + n->getSumLnProbJumpSize(activeParm);
	//	double y = n->getLnProbJumpCount(inactiveParm) + n->getSumLnProbJumpSize(inactiveParm);
	//	std::cout << n->getIndex() << ":\t" << x << "\n";
	//	std::cout << n->getIndex() << ":\t" << y << "\n";
		n->setK(x);
	}

	for (int j = 0; j < numNodes; j++)
	{
		Node* n = topologyPtr->getDownPassNode(j);

		if (n->getLft() != NULL && n->getRht() != NULL)
		{
			// Branch lengths depend on which node's likelihood is being calculated.
			double tauL = branchVals[n->getLft()->getIndex()];
			double tauR = branchVals[n->getRht()->getIndex()];

			for (std::list<Patron*>::iterator it_p = patronList.begin(); it_p != patronList.end(); it_p++)
			{
#if DEBUG
				std::cout << "\tPATRON ID:\t" << (*it_p)->getId() << "\n";
				std::cout << "p->getIndex() = " << n->getIndex() << "\tp->getLft()->getIndex() = " << n->getLft()->getIndex() << "\tp->getRht()->getIndex() = " << n->getRht()->getIndex() << "\n";
				std::cout << "Parms [l,s,tL,tR] = " << parmVals[0] << "\t" << parmVals[1] << "\t" << tauL << "\t" << tauR << "\n";

#endif

				int activeParm = n->getActiveParm();

				// get descendants' conditional parameters
				double muL = n->getLft()->getMu() - n->getLft()->getSumJumpSize(activeParm);
				double sigmaL = n->getLft()->getSigma();
				double kL = n->getLft()->getK();

				double muR = n->getRht()->getMu() - n->getRht()->getSumJumpSize(activeParm);
				double sigmaR = n->getRht()->getSigma();
				double kR = n->getRht()->getK();

				// rescale branches according to descendants' variance
				double tL = tauL + (sigmaL * sigmaL) / sigmaBM2;
				double tR = tauR + (sigmaR * sigmaR) / sigmaBM2;

				// update present node's parameters
				double muP = (muL * tR + muR * tL) / (tL + tR);
				double sigmaP = pow(sigmaBM2 * tL * tR / (tL + tR), 0.5);
				n->setMu(muP);
				n->setSigma(sigmaP);

				// calculate likelihood of sampled jumps
				double lnJumpLike = 0.0;
				int jumpCount = n->getJumpCount(activeParm);
				std::vector<double>::iterator jumpSizeIt = n->getJumpSize(activeParm).begin();
				std::vector<double>::iterator jumpSizeEnd = n->getJumpSize(activeParm).end();
				lnJumpLike += randomPtr->poissonProb(lambdaJN, jumpCount);
				for ( ; jumpSizeIt != jumpSizeEnd; jumpSizeIt++)
				{
					lnJumpLike += randomPtr->lnNormalPdf(0.0, sigmaJN, *jumpSizeIt);
				}

				// update normalizing constant
				double kP = -1.0 * pow(muL - muR, 2) / (2 * sigmaBM2 * (tL + tR));
				kP -= log(pow(2 * PI * sigmaBM2 * (tL + tR), 0.5));
				kP += kL + kR;

				//n->setK(kP + n->getSumLnProbJumpSize(activeParm) + n->getLnProbJumpCount(activeParm));
				n->setK(kP + lnJumpLike);

#if DEBUG2
				std::cout << std::setprecision(4);
				std::cout << n->getLft()->getIndex() << "\t" << muL << "\t" << sigmaL << "\t" << kL << "\t";
				std::cout << n->getRht()->getIndex() << "\t" << muR << "\t" << sigmaR << "\t" << kR << "\t";
				std::cout << n->getIndex() << "\t" << muP << "\t" << sigmaP << "\t" << kP << "\n";
#endif

			}
		}
	}

	lnL = topologyPtr->getRoot()->getK();
	return lnL;
}


double Model::tableLogLikelihoodFFT(Table* t)
{
	double lnL = 0.0;

	const std::vector<Parm*> tempParmVector = t->getParmVector();
	const std::vector<Parm*> tempBranchVector = t->getBranchVector();

	// Collect parameter & branch values from the Table's parameter class
	std::vector<double> parmVals;
	std::vector<double> branchVals;
	for (unsigned int i = 0; i < tempParmVector.size(); i++)
	{
		//Parm* tempParm = tempParmVector[i];
		//tempParm->print();
		parmVals.push_back(tempParmVector[i]->getValue());
	}
	for (unsigned int i = 0; i < tempBranchVector.size(); i++)
	{
		branchVals.push_back(tempBranchVector[i]->getValue());
	}

	bool testVals = false;
	if (testVals)
	{
		parmVals.clear();
		parmVals.push_back(0.3);
		parmVals.push_back(1.2);
		parmVals.push_back(0.3);
	}
	bool testBranches = false;
	if (testBranches)
	{
		branchVals.clear();
		branchVals.push_back(0.8);
		branchVals.push_back(1.2);
	}

	// Calculate likelihoods for single Patron (labeled to a different Parm class)
	for (int j = 0; j < numNodes; j++)
	{
		Node* n = topologyPtr->getDownPassNode(j);

		if (n != NULL && n->getLft() != NULL && n->getRht() != NULL)
		{
			// Branch lengths depend on which node's likelihood is being calculated.
			double tauL = branchVals[n->getLft()->getIndex()];
			double tauR = branchVals[n->getRht()->getIndex()];
			int status;

			for (std::list<Patron*>::iterator it_p = patronList.begin(); it_p != patronList.end(); it_p++)
			{
#if DEBUG
				std::cout << "\tPATRON ID:\t" << (*it_p)->getId() << "\n";
				std::cout << "p->getIndex() = " << n->getIndex() << "\tp->getLft()->getIndex() = " << n->getLft()->getIndex() << "\tp->getRht()->getIndex() = " << n->getRht()->getIndex() << "\n";
				std::cout << "Parms [l,s,tL,tR] = " << parmVals[0] << "\t" << parmVals[1] << "\t" << tauL << "\t" << tauR << "\n";

#endif

				double* clL = condLikes->getClPtr(n->getActiveCl(), (*it_p)->getId(), n->getLft()->getIndex(), 0);
				double* clR = condLikes->getClPtr(n->getActiveCl(), (*it_p)->getId(), n->getRht()->getIndex(), 0);
				double* clA = condLikes->getClPtr(n->getActiveCl(), (*it_p)->getId(), n->getIndex(), 0);

				// set data frame to function vals
				for (int i = 0; i < numSteps; i++)
				{
					gsl_complex complex_pdf_L = charFunc(theta[i], tauL, parmVals);
					gsl_complex complex_pdf_R = charFunc(theta[i], tauR, parmVals);

					gsl_complex complex_tip_L = gsl_complex_rect(clL[2*i], clL[2*i+1]);
					gsl_complex complex_tip_R = gsl_complex_rect(clR[2*i], clR[2*i+1]);

					gsl_complex conv_tip_L = gsl_complex_mul(complex_pdf_L, complex_tip_L);
					gsl_complex conv_tip_R = gsl_complex_mul(complex_pdf_R, complex_tip_R);


#if DEBUG
					std::cout << std::setprecision(16);
					if (i == 0)
					{
						std::cout << "\tL\t\t\t\t\t\tR\n";
						std::cout << "n\tcf\t\tcl\t\tcv\t\tcf\t\tcl\t\tcv\n";
					}
					std::cout << "";
					std::cout << "REAL: " << i << ":\t";
					std::cout << GSL_REAL(complex_pdf_L) << "\t" << GSL_REAL(complex_tip_L) << "\t" << GSL_REAL(conv_tip_L) << "\t";
					std::cout << GSL_REAL(complex_pdf_R) << "\t" << GSL_REAL(complex_tip_R) << "\t" << GSL_REAL(conv_tip_R) << "\n";
					std::cout << "IMAG: " << i << ":\t";
					std::cout << GSL_IMAG(complex_pdf_L) << "\t" << GSL_IMAG(complex_tip_L) << "\t" << GSL_IMAG(conv_tip_L) << "\t";
					std::cout << GSL_IMAG(complex_pdf_R) << "\t" << GSL_IMAG(complex_tip_R) << "\t" << GSL_IMAG(conv_tip_R) << "\n";
#endif


					REAL(tip1,i) = GSL_REAL(conv_tip_L);// * 20;
					IMAG(tip1,i) = GSL_IMAG(conv_tip_L);// * pow(10,10);
					REAL(tip2,i) = GSL_REAL(conv_tip_R);// * 20;
					IMAG(tip2,i) = GSL_IMAG(conv_tip_R);// * pow(10,10);

					if (std::isnan(GSL_REAL(conv_tip_L))) std::cerr << "ERROR: conv_tip_L[" << i << "] = NaN\n";
					if (std::isnan(GSL_REAL(conv_tip_R))) std::cerr << "ERROR: conv_tip_R[" << i << "] = NaN\n";

					//std::cout << "\n";
				}


				status = gsl_fft_complex_radix2_inverse(tip1, 1, numSteps);
				if (status != GSL_SUCCESS) std::cerr << "ERROR: tip1 iFFT\n";
				status = gsl_fft_complex_radix2_inverse(tip2, 1, numSteps);
				if (status != GSL_SUCCESS) std::cerr << "ERROR: tip2 iFFT\n";

				//std::cout << "\tpost-iFFT\n";
				for (int i = 0; i < numSteps; i++)
				{
					gsl_complex complex_tip1;
					gsl_complex complex_tip2;

					GSL_SET_COMPLEX(&complex_tip1, REAL(tip1,i), 0.0); // IMAG(tip1,i));
					GSL_SET_COMPLEX(&complex_tip2, REAL(tip2,i), 0.0); //IMAG(tip2,i));

					gsl_complex final_val = gsl_complex_mul(complex_tip1, complex_tip2);

					/*
#if DEBUG

					std::cout << i << ":\t" << GSL_REAL(complex_tip1) << "\t" << GSL_IMAG(complex_tip1) << "\n";
					std::cout << "\t" << GSL_REAL(complex_tip2) << "\t" << GSL_IMAG(complex_tip2) << "\n";
					std::cout << "\t" << GSL_REAL(final_val) << "\t" << GSL_IMAG(final_val) << "\n";
#endif
*/

					REAL(like,i) = GSL_REAL(final_val);
					IMAG(like,i) = 0.0;//GSL_IMAG(final_val);

					//clA[2*i] = GSL_REAL(final_val);
					//clA[2*i+1] = GSL_IMAG(final_val);

					//REAL(like,i) = GSL_REAL(final_val);
					//IMAG(like,i) = GSL_IMAG(final_val);
					//likeVals[i] = GSL_REAL(final_val);

#if DEBUG
					if (i == 0)
					{
						std::cout << "\tL\t\t\t\t\t\tR\n";
						std::cout << "n\tRe(L)\t\tIm(L)\t\tRe(R)\t\tIm(R)\t\tRe(F)\t\tIm(F)\n";
					}
					std::cout << std::setprecision(16);
					std::cout << i << "\t" << GSL_REAL(complex_tip1) << "\t" << GSL_IMAG(complex_tip1);
					std::cout << "\t" << GSL_REAL(complex_tip2) << "\t" << GSL_IMAG(complex_tip2);
					std::cout << "\t" << GSL_REAL(final_val) << "\t" << GSL_IMAG(final_val) << "\n";
#endif

				}

				if (n->getAnc() != NULL)
				{
					gsl_fft_complex_radix2_forward(like, 1, numSteps);
					if (status != GSL_SUCCESS) std::cerr << "ERROR:  fft fwd\n";
				}

				//std::cout << "\tpost-fFFT\n";
				//std::cout << "\tNode: " << j << "\tn->getIndex()" << n->getIndex() << "\n";
				for (int i = 0; i < numSteps; i++)
				{

#if DEBUG
					if (i == 0)
					{
						std::cout << "\tL\t\t\t\t\t\tR\n";
						std::cout << "n\tRe(CL)\tIm(CL)\n";
					}
					std::cout << i << "\t" << REAL(like,i) << "\t" << IMAG(like,i) << "\n";
#endif
					clA[2*i] = REAL(like,i); //* 100;//pow(10,10);
					clA[2*i+1] = IMAG(like,i);//* 100;//;//  * pow(10,10);
					//std::cout << std::setprecision(12) << "\t\tre.like[" << i << "] = " << clA[2*i] << "\tim.like[" << i << "] = " << clA[2*i+1] << "\n";
				}
				//std::cout << "\n";

				// calculate the likelihood at the root
				if (n->getAnc() == NULL)
				{
					//std::cout << "root found\n";
					for (int i = 0; i < numSteps; i++)
					{
						likeVals[i] = clA[2*i];// / pow(10,numNodes-numTaxa);// * pow(10,100);
						//if (likeVals[i] < 0.0) likeVals[i] = 0.0;  // uncommented: 08/12/11

					}

					double val = trapInt(likeVals);

					// NOTE: Prevents negative FFT values from producing NaN in lnL
					if (val <= 0.0 || std::isnan(val)) val = pow(10,-300);

					lnL += log(val); // likelihood of data observed for patron

					// std::cout << std::setprecision(24) << "\t\t\tLIKELIHOOD: " << val << "\n";
					if (std::isnan(log(val)))
					{
						std::cerr << "ERROR: val = NaN\n";
					}

					//std::cout << lnL << "\n";
				}
			}
		}
	}

	// return lnL for all datapoints for a patron at the given table
	return lnL;
}


#if 0

double Model::tableLogLikelihood(Table* t)
{
	double lnL = 0.0;

	// Read in the parameter values.
	double lambda = 0.0;
	double sigma = 0.0;
	double tau1 = 0.0;
	double tau2 = 0.0;

	const std::vector<Parm*> tempParmVector = t->getParmVector();
	bool isRight = false;
	for (std::vector<Parm*>::const_iterator it_p = tempParmVector.begin(); it_p != tempParmVector.end(); it_p++)
	{
		Lambda* derivedPtrL = dynamic_cast<Lambda*> (*it_p);
		Tau* derivedPtrT = dynamic_cast<Tau*> (*it_p);
		Sigma* derivedPtrS = dynamic_cast<Sigma*> (*it_p);
		if (derivedPtrL != 0)
			lambda = derivedPtrL->getActiveExpMean()->getRate();
		else if (derivedPtrS != 0)
			sigma = derivedPtrS->getActiveVar()->getRate();
		else if (derivedPtrT != 0 && isRight == false)
		{
			tau1 = derivedPtrT->getActiveBrLen()->getRate();
			isRight = true;
		}
		else if (derivedPtrT != 0)
		{
			tau2 = derivedPtrT->getActiveBrLen()->getRate();
		}
	}

#if LAMBDA_ONLY
	// TEST 05/06/11 - LAMBDA
	sigma = 0.1;
	tau1 = 2.0;
	tau2 = 0.5;
#elif SIGMA_ONLY
	// TEST 05/06/11 - SIGMA
	lambda = 0.5;
	tau1 = 2.0;
	tau2 = 0.5;
#elif TAU1_ONLY
	// TEST 05/06/11 - TAU1
	lambda = 0.5;
	sigma = 0.1;
	tau2 = 0.5;
#elif TAU2_ONLY
	// TEST 05/06/11 - TAU2
	lambda = 0.5;
	sigma = 0.1;
	tau1 = 2.0;
#endif

	// std::cout << "\t\t" << lambda << " " << sigma << " " << tau1 << " " << tau2 << "\n";

	// Calculate likelihoods for all Patrons of a given table (i.e. labeled to a different Parm classes)
	for (std::list<Patron*>::const_iterator it_p = t->getPatronList().begin(); it_p != t->getPatronList().end(); it_p++)
	{

		int i = 0;

		// set data frame to function vals
		for (i = 0; i < numSteps; i++)
		{
			gsl_complex complex_pdf1 = complexPdf(theta[i], sigma, tau1, lambda);
			gsl_complex complex_pdf2 = complexPdf(theta[i], sigma, tau2, lambda);
			gsl_complex complex_tip1 = (*it_p)->getCfData(0,0,i);
			gsl_complex complex_tip2 = (*it_p)->getCfData(1,0,i);
			gsl_complex conv_tip1 = gsl_complex_mul(complex_pdf1, complex_tip1);
			gsl_complex conv_tip2 = gsl_complex_mul(complex_pdf2, complex_tip2);

			//std::cout << std::setprecision(8) << "\tFn" << i << ":\t" << "tip1\t" << GSL_REAL(complex_pdf1) << " " << GSL_IMAG(complex_pdf1) << "\n";
			//std::cout << std::setprecision(8) << "\tCh" << i << ":\t" << "tip1\t" << GSL_REAL(complex_tip1) << " " << GSL_IMAG(complex_tip1) << "\n";
			//std::cout << std::setprecision(8) << "\tCV" << i << ":\t" << "tip1\t" << GSL_REAL(conv_tip1) << " " << GSL_IMAG(conv_tip1) << "\n";
			//std::cout << std::setprecision(8) << "\tFn" << i << ":\t" << "tip2\t" << GSL_REAL(complex_pdf2) << " " << GSL_IMAG(complex_pdf2) << "\n";
			//std::cout << std::setprecision(8) << "\tCh" << i << ":\t" << "tip2\t" << GSL_REAL(complex_tip2) << " " << GSL_IMAG(complex_tip2) << "\n";
			//std::cout << std::setprecision(8) << "\tCV" << i << ":\t" << "tip2\t" << GSL_REAL(conv_tip2) << " " << GSL_IMAG(conv_tip2) << "\n";

			REAL(tip1,i) = GSL_REAL(conv_tip1);
			IMAG(tip1,i) = GSL_IMAG(conv_tip1);
			REAL(tip2,i) = GSL_REAL(conv_tip2);
			IMAG(tip2,i) = GSL_IMAG(conv_tip2);

			//std::cout << "\n";
		}

		gsl_fft_complex_radix2_inverse(tip1, 1, numSteps);
		gsl_fft_complex_radix2_inverse(tip2, 1, numSteps);

		for (i = 0; i < numSteps; i++)
		{
			gsl_complex complex_tip1;
			GSL_SET_COMPLEX(&complex_tip1, REAL(tip1,i), IMAG(tip1,i));

			gsl_complex complex_tip2;
			GSL_SET_COMPLEX(&complex_tip2, REAL(tip2,i), IMAG(tip2,i));

			gsl_complex final_val = gsl_complex_mul(complex_tip1, complex_tip2);

			REAL(like,i) = GSL_REAL(final_val);
			IMAG(like,i) = GSL_IMAG(final_val);
			likeVals[i] = GSL_REAL(final_val);

			//std::cout << GSL_REAL(final_val) << " " << GSL_IMAG(final_val) << "\n";
			//std::cout << std::setprecision(8) << "\t" << i << ":\t" << "tip1\t" << GSL_REAL(complex_tip1) << " " << GSL_IMAG(complex_tip1) << "\n";
			//std::cout << std::setprecision(8) << "\t\ttip2\t" << GSL_REAL(complex_tip2) << " " << GSL_IMAG(complex_tip2) << "\n";
			//std::cout  << "\t\tfinal\t";// << GSL_REAL(final_val) << "\n"; // << " " << GSL_IMAG(final_val) << "\n";
		}

		double val = trapInt(likeVals);

		// NOTE: Prevents negative FFT values from producing NaN in lnL
		if (val <= 0.0 || std::isnan(val)) val = pow(10,-300);

		// std::cout << "\t\t\tLIKELIHOOD: " << val << "\n";
		//std::cout << "\t\t\tLIKELIHOOD: " << val*val << "\n";

		lnL += log(val); // likelihood of data observed for patron
	}

	// return lnL for all datapoints for a patron at the given table
	return lnL;
}
#endif

/*
double Model::tableLogLikelihood(Table* t)
{
	double lnL = 0.0;

	// Read in the parameter values.
	double p = 0.0;
	double sigmaK = 0.0;
	double sigmaB = 0.0;
	double tauL = 0.0;
	double tauR = 0.0;

	Table* t = p->getTable();
	const std::vector<Parm*> tempParmVector = t->getParmVector();
	const std::vector<Parm*> tempBranchVector = t->getBranchVector();

	for (std::list<Patron*>::iterator it_p = t->getPatronList().begin(); it_p != t->getPatronList().end(); it_p++)
	{
		lnL += locusLogLikelihood(*it_p);
	}

}
*/

double Model::locusLogLikelihood(Patron* p)
{
	double lnL = 0.0;

	// Read in the parameter values.
	double lambdaK = 0.0;
	double sigmaK = 0.0;
	double sigmaB = 0.0;
	double tauL = 0.0;
	double tauR = 0.0;

	Table* t = p->getTable();
	const std::vector<Parm*> tempParmVector = t->getParmVector();
	const std::vector<Parm*> tempBranchVector = t->getBranchVector();

	// Filthy, disgusting hard-coded parameterization. Must allow arbitrary functions & params in future.
	lambdaK = tempParmVector[0]->getValue();
	sigmaK = tempParmVector[1]->getValue();
	sigmaB = tempParmVector[2]->getValue();

	/*
	for (std::vector<Parm*>::const_iterator it_p = tempParmVector.begin(); it_p != tempParmVector.end(); it_p++)
	{
		Lambda* derivedPtrL = dynamic_cast<Lambda*> (*it_p);
		Sigma* derivedPtrS = dynamic_cast<Sigma*> (*it_p);
		if (derivedPtrL != 0)
			lambda = derivedPtrL->getActiveExpMean()->getRate();
		else if (derivedPtrS != 0)
			sigma = derivedPtrS->getActiveVar()->getRate();
	}
	*/

	// Calculate likelihoods for single Patron (labeled to a different Parm class)
	for (int j = 0; j < numNodes; j++)
	{
		Node* n = topologyPtr->getDownPassNode(j);

		if (n != NULL && n->getLft() != NULL && n->getRht() != NULL)
		{

			//std::cout << "\tgene: " << p->getId() << "\tP:" << j << "\t" << n->getName() << "\t" << n->getIndex() << "\n";
			//std::cout << "\tG: " << p->getId() << "\tP:" << n->getIndex() << "\tL:" << n->getLft()->getIndex() << "\tR:" << n->getRht()->getIndex() << "\n";

			int i;

			//std::cout << "\tPATRON ID:\t" << p->getId() << "\n";
			double* clL = condLikes->getClPtr(n->getActiveCl(), p->getId(), n->getLft()->getIndex(), 0);
			double* clR = condLikes->getClPtr(n->getActiveCl(), p->getId(), n->getRht()->getIndex(), 0);
			double* clA = condLikes->getClPtr(n->getActiveCl(), p->getId(), n->getIndex(), 0);

			// get branch lengths (if the branchVector size == 0, initialization occurs)
			if (tempBranchVector.size() > 0)
			{
				//for (int j = 0; j < tempBranchVector.size(); j++)
				//	tempBranchVector[j]->print();
				Tau* derivedPtrT_L = dynamic_cast<Tau*> (tempBranchVector[n->getLft()->getIndex()]);
				Tau* derivedPtrT_R = dynamic_cast<Tau*> (tempBranchVector[n->getRht()->getIndex()]);
				tauL = derivedPtrT_L->getActiveBrLen()->getValue();
				tauR = derivedPtrT_R->getActiveBrLen()->getValue();
			}

			//std::cout << "p->getIndex() = " << n->getIndex() << "\tp->getLft()->getIndex() = " << n->getLft()->getIndex() << "\tp->getRht()->getIndex() = " << n->getRht()->getIndex() << "\n";
			//std::cout << "Parms [l,s,tL,tR] = " << lambda << "\t" << sigma << "\t" << tauL << "\t" << tauR << "\n";

			// set data frame to function vals
			for (i = 0; i < numSteps; i++)
			{
				// NOTE: update 06/28/11
				gsl_complex complex_pdf_L = complexPdf(theta[i], tauL, sigmaB, sigmaK, lambdaK);
				gsl_complex complex_pdf_R = complexPdf(theta[i], tauR, sigmaB, sigmaK, lambdaK);

				//gsl_complex complex_pdf_L = cfJumpDiffusion(theta[i], tauL * sigma, tauL * , lambda);
				//gsl_complex complex_pdf_R = cfJumpDiffusion(theta[i], tauR * sigma, tauR * )

				//gsl_complex complex_tip_L2 = p->getCfData(0,0,i);
				//gsl_complex complex_tip_R2 = p->getCfData(1,0,i);

				//gsl_complex complex_tip_L = gsl_complex_rect(clL->at(2*i), clL->at(2*i+1));
				//gsl_complex complex_tip_R = gsl_complex_rect(clR->at(2*i), clR->at(2*i+1));
				gsl_complex complex_tip_L = gsl_complex_rect(clL[2*i], clL[2*i+1]);
				gsl_complex complex_tip_R = gsl_complex_rect(clR[2*i], clR[2*i+1]);

				gsl_complex conv_tip_L = gsl_complex_mul(complex_pdf_L, complex_tip_L);
				gsl_complex conv_tip_R = gsl_complex_mul(complex_pdf_R, complex_tip_R);

				//std::cout << std::setprecision(16);
				//std::cout << i << ":";
				//std::cout << "\tFn-L\t" << GSL_REAL(complex_pdf_L) << " " << GSL_IMAG(complex_pdf_L) << "\n";
				//std::cout << "\tCh" << i << ":\t" << "tip1\t" << GSL_REAL(complex_tip_L) << " " << GSL_IMAG(complex_tip_L) << "\n";
				//std::cout << std::setprecision(8) << "\tC2" << i << ":\t" << "tip1\t" << GSL_REAL(complex_tip_L2) << " " << GSL_IMAG(complex_tip_L2) << "\n";
				//std::cout << "\tCV-L1\t" << GSL_REAL(conv_tip_L) << " " << GSL_IMAG(conv_tip_L) << "\n";
				//std::cout << "\tFn-R\t" << GSL_REAL(complex_pdf_R) << " " << GSL_IMAG(complex_pdf_R) << "\n";
				//std::cout << std::setprecision(8) << "\tCh" << i << ":\t" << "tip2\t" << GSL_REAL(complex_tip_R) << " " << GSL_IMAG(complex_tip_R) << "\n";
				//std::cout << std::setprecision(8) << "\tC2" << i << ":\t" << "tip2\t" << GSL_REAL(complex_tip_R2) << " " << GSL_IMAG(complex_tip_R2) << "\n";
				//std::cout <<  "\tCV-R\t" << GSL_REAL(conv_tip_R) << " " << GSL_IMAG(conv_tip_R) << "\n";
				//std::cout << "clL:" << clL[2*i] << "\t" << clL[2*i+1] << "\n";
				//std::cout << "clR:" << clR[2*i] << "\t" << clR[2*i+1] << "\n";
				//std::cout << "\n";


				//clL->at(i) = GSL_REAL(conv_tip_L);
				//clR->at(i) = GSL_
				REAL(tip1,i) = GSL_REAL(conv_tip_L);
				IMAG(tip1,i) = GSL_IMAG(conv_tip_L);
				REAL(tip2,i) = GSL_REAL(conv_tip_R);
				IMAG(tip2,i) = GSL_IMAG(conv_tip_R);

				if (std::isnan(GSL_REAL(conv_tip_L))) std::cerr << "ERROR: conv_tip_L[" << i << "] = NaN\n";
				if (std::isnan(GSL_REAL(conv_tip_R))) std::cerr << "ERROR: conv_tip_R[" << i << "] = NaN\n";

				//std::cout << "\n";
			}
			//std::cout << "\tFFT\n";

			gsl_fft_complex_radix2_inverse(tip1, 1, numSteps);
			gsl_fft_complex_radix2_inverse(tip2, 1, numSteps);

			//std::cout << "\tiFFT\n";
			for (i = 0; i < numSteps; i++)
			{
				gsl_complex complex_tip1;
				gsl_complex complex_tip2;

				GSL_SET_COMPLEX(&complex_tip1, REAL(tip1,i), IMAG(tip1,i));
				GSL_SET_COMPLEX(&complex_tip2, REAL(tip2,i), IMAG(tip2,i));

				gsl_complex final_val = gsl_complex_mul(complex_tip1, complex_tip2);

				//std::cout << std::setprecision(16);
				//std::cout << i << ":\t" << GSL_REAL(complex_tip1) << "\t" << GSL_IMAG(complex_tip1) << "\n";
				//std::cout << "\t" << GSL_REAL(complex_tip2) << "\t" << GSL_IMAG(complex_tip2) << "\n";
				//std::cout << "\t" << GSL_REAL(final_val) << "\t" << GSL_IMAG(final_val) << "\n";

				REAL(like,i) = GSL_REAL(final_val);
				IMAG(like,i) = GSL_IMAG(final_val);

				//clA[2*i] = GSL_REAL(final_val);
				//clA[2*i+1] = GSL_IMAG(final_val);

				//REAL(like,i) = GSL_REAL(final_val);
				//IMAG(like,i) = GSL_IMAG(final_val);
				// likeVals[i] = GSL_REAL(final_val);
			}

			if (n->getAnc() != NULL)
			{
				gsl_fft_complex_radix2_forward(like, 1, numSteps);
			}

			//std::cout << "\tNode: " << j << "\tn->getIndex()" << n->getIndex() << "\n";
			for (i = 0; i < numSteps; i++)
			{
				clA[2*i] = REAL(like,i);
				clA[2*i+1] = IMAG(like,i);
				//std::cout << std::setprecision(12) << "\t\tre.like[" << i << "] = " << clA[2*i] << "\tim.like[" << i << "] = " << clA[2*i+1] << "\n";
			}
			//std::cout << "\n";

			// calculate the likelihood at the root
			if (n->getAnc() == NULL)
			{
				//std::cout << "root found\n";
				for (int i = 0; i < numSteps; i++)
				{
					likeVals[i] = clA[2*i];
				}

				double val = trapInt(likeVals);

				// NOTE: Prevents negative FFT values from producing NaN in lnL
				if (val <= 0.0 || std::isnan(val)) val = pow(10,-300);

				lnL += log(val); // likelihood of data observed for patron

				// std::cout << std::setprecision(24) << "\t\t\tLIKELIHOOD: " << val << "\n";
				if (std::isnan(log(val)))
				{
					std::cerr << "ERROR: val = NaN\n";
					std::cerr << "\tlK:    " << lambdaK << "\n";
					std::cerr << "\tsK:    " << sigmaK << "\n";
					std::cerr << "\tsB:    " << sigmaB << "\n";
					std::cerr << "\tt1:    " << tauL << "\n";
					std::cerr << "\tt2:    " << tauR << "\n";
			//		std::cerr << "\ttip1: " << t1 << "\n";
			//		std::cerr << "\ttip2: " << t2 << "\n";

				}

				//std::cout << lnL << "\n";
			}
		}
	}

	// return lnL for all datapoints for a patron at the given table
	return lnL;
}

/*
double Model::locusLogLikelihood2(Patron* p)
{
	double lnL = 0.0;

	// Read in the parameter values.
	double lambda = 0.0;
	double sigma = 0.0;
	double tauL = 0.0;
	double tauR = 0.0;

	Table* t = p->getTable();
	const std::vector<Parm*> tempParmVector = t->getParmVector();
	const std::vector<Parm*> tempBranchVector = t->getBranchVector();

	for (std::vector<Parm*>::const_iterator it_p = tempParmVector.begin(); it_p != tempParmVector.end(); it_p++)
	{
		Lambda* derivedPtrL = dynamic_cast<Lambda*> (*it_p);
		Sigma* derivedPtrS = dynamic_cast<Sigma*> (*it_p);
		if (derivedPtrL != 0)
			lambda = derivedPtrL->getActiveExpMean()->getRate();
		else if (derivedPtrS != 0)
			sigma = derivedPtrS->getActiveVar()->getRate();
	}

	// Calculate likelihoods for single Patron (labeled to a different Parm class)
	for (int j = 0; j < numNodes; j++)
	{
		Node* n = topologyPtr->getDownPassNode(j);

		if (n != NULL && n->getLft() != NULL && n->getRht() != NULL)
		{

			//std::cout << "\tgene: " << p->getId() << "\tP:" << j << "\t" << n->getName() << "\t" << n->getIndex() << "\n";
			//std::cout << "\tG: " << p->getId() << "\tP:" << n->getIndex() << "\tL:" << n->getLft()->getIndex() << "\tR:" << n->getRht()->getIndex() << "\n";

			int i;

			//std::cout << "\tPATRON ID:\t" << p->getId() << "\n";
			double* clL = condLikes->getClPtr(n->getActiveCl(), p->getId(), n->getLft()->getIndex(), 0);
			double* clR = condLikes->getClPtr(n->getActiveCl(), p->getId(), n->getRht()->getIndex(), 0);
			double* clA = condLikes->getClPtr(n->getActiveCl(), p->getId(), n->getIndex(), 0);

			// get branch lengths (if the branchVector size == 0, initialization occurs)
			if (tempBranchVector.size() > 0)
			{
				//for (int j = 0; j < tempBranchVector.size(); j++)
				//	tempBranchVector[j]->print();
				Tau* derivedPtrT_L = dynamic_cast<Tau*> (tempBranchVector[n->getLft()->getIndex()]);
				Tau* derivedPtrT_R = dynamic_cast<Tau*> (tempBranchVector[n->getRht()->getIndex()]);
				tauL = derivedPtrT_L->getActiveBrLen()->getRate();
				tauR = derivedPtrT_R->getActiveBrLen()->getRate();
			}

			//lambda = 2.0;
			//sigma = 0.1;
			//tauL = 1.0;
			//tauR = 1.0;

			//std::cout << "p->getIndex() = " << n->getIndex() << "\tp->getLft()->getIndex() = " << n->getLft()->getIndex() << "\tp->getRht()->getIndex() = " << n->getRht()->getIndex() << "\n";
			//std::cout << "Parms [l,s,tL,tR] = " << lambda << "\t" << sigma << "\t" << tauL << "\t" << tauR << "\n";

			// set data frame to function vals
			for (i = 0; i < numSteps; i++)
			{
				// NOTE: update 06/28/11
				gsl_complex complex_pdf_L = complexPdf(theta[i], sigma, tauL, lambda);
				gsl_complex complex_pdf_R = complexPdf(theta[i], sigma, tauR, lambda);

				//gsl_complex complex_pdf_L = cfJumpDiffusion(theta[i], tauL * sigma, tauL * , lambda);
				//gsl_complex complex_pdf_R = cfJumpDiffusion(theta[i], tauR * sigma, tauR * )

				//gsl_complex complex_tip_L2 = p->getCfData(0,0,i);
				//gsl_complex complex_tip_R2 = p->getCfData(1,0,i);

				//gsl_complex complex_tip_L = gsl_complex_rect(clL->at(2*i), clL->at(2*i+1));
				//gsl_complex complex_tip_R = gsl_complex_rect(clR->at(2*i), clR->at(2*i+1));
				gsl_complex complex_tip_L = gsl_complex_rect(clL[2*i], clL[2*i+1]);
				gsl_complex complex_tip_R = gsl_complex_rect(clR[2*i], clR[2*i+1]);

				gsl_complex conv_tip_L = gsl_complex_mul(complex_pdf_L, complex_tip_L);
				gsl_complex conv_tip_R = gsl_complex_mul(complex_pdf_R, complex_tip_R);

				//std::cout << std::setprecision(16);
				//std::cout << i << ":";
				//std::cout << "\tFn-L\t" << GSL_REAL(complex_pdf_L) << " " << GSL_IMAG(complex_pdf_L) << "\n";
				//std::cout << "\tCh" << i << ":\t" << "tip1\t" << GSL_REAL(complex_tip_L) << " " << GSL_IMAG(complex_tip_L) << "\n";
				//std::cout << std::setprecision(8) << "\tC2" << i << ":\t" << "tip1\t" << GSL_REAL(complex_tip_L2) << " " << GSL_IMAG(complex_tip_L2) << "\n";
				//std::cout << "\tCV-L1\t" << GSL_REAL(conv_tip_L) << " " << GSL_IMAG(conv_tip_L) << "\n";
				//std::cout << "\tFn-R\t" << GSL_REAL(complex_pdf_R) << " " << GSL_IMAG(complex_pdf_R) << "\n";
				//std::cout << std::setprecision(8) << "\tCh" << i << ":\t" << "tip2\t" << GSL_REAL(complex_tip_R) << " " << GSL_IMAG(complex_tip_R) << "\n";
				//std::cout << std::setprecision(8) << "\tC2" << i << ":\t" << "tip2\t" << GSL_REAL(complex_tip_R2) << " " << GSL_IMAG(complex_tip_R2) << "\n";
				//std::cout <<  "\tCV-R\t" << GSL_REAL(conv_tip_R) << " " << GSL_IMAG(conv_tip_R) << "\n";
				//std::cout << "clL:" << clL[2*i] << "\t" << clL[2*i+1] << "\n";
				//std::cout << "clR:" << clR[2*i] << "\t" << clR[2*i+1] << "\n";
				//std::cout << "\n";


				//clL->at(i) = GSL_REAL(conv_tip_L);
				//clR->at(i) = GSL_
				REAL(tip1,i) = GSL_REAL(conv_tip_L);
				IMAG(tip1,i) = GSL_IMAG(conv_tip_L);
				REAL(tip2,i) = GSL_REAL(conv_tip_R);
				IMAG(tip2,i) = GSL_IMAG(conv_tip_R);

				if (std::isnan(GSL_REAL(conv_tip_L))) std::cerr << "ERROR: conv_tip_L[" << i << "] = NaN\n";
				if (std::isnan(GSL_REAL(conv_tip_R))) std::cerr << "ERROR: conv_tip_R[" << i << "] = NaN\n";

				//std::cout << "\n";
			}
			//std::cout << "\tFFT\n";

			gsl_fft_complex_radix2_inverse(tip1, 1, numSteps);
			gsl_fft_complex_radix2_inverse(tip2, 1, numSteps);

			//std::cout << "\tiFFT\n";
			for (i = 0; i < numSteps; i++)
			{
				gsl_complex complex_tip1;
				gsl_complex complex_tip2;

				GSL_SET_COMPLEX(&complex_tip1, REAL(tip1,i), IMAG(tip1,i));
				GSL_SET_COMPLEX(&complex_tip2, REAL(tip2,i), IMAG(tip2,i));

				gsl_complex final_val = gsl_complex_mul(complex_tip1, complex_tip2);

				//std::cout << std::setprecision(16);
				//std::cout << i << ":\t" << GSL_REAL(complex_tip1) << "\t" << GSL_IMAG(complex_tip1) << "\n";
				//std::cout << "\t" << GSL_REAL(complex_tip2) << "\t" << GSL_IMAG(complex_tip2) << "\n";
				//std::cout << "\t" << GSL_REAL(final_val) << "\t" << GSL_IMAG(final_val) << "\n";

				REAL(like,i) = GSL_REAL(final_val);
				IMAG(like,i) = GSL_IMAG(final_val);

				//clA[2*i] = GSL_REAL(final_val);
				//clA[2*i+1] = GSL_IMAG(final_val);

				//REAL(like,i) = GSL_REAL(final_val);
				//IMAG(like,i) = GSL_IMAG(final_val);
				// likeVals[i] = GSL_REAL(final_val);
			}

			if (n->getAnc() != NULL)
			{
				gsl_fft_complex_radix2_forward(like, 1, numSteps);
			}

			//std::cout << "\tNode: " << j << "\tn->getIndex()" << n->getIndex() << "\n";
			for (i = 0; i < numSteps; i++)
			{
				clA[2*i] = REAL(like,i);
				clA[2*i+1] = IMAG(like,i);
				//std::cout << std::setprecision(12) << "\t\tre.like[" << i << "] = " << clA[2*i] << "\tim.like[" << i << "] = " << clA[2*i+1] << "\n";
			}
			//std::cout << "\n";

			// calculate the likelihood at the root
			if (n->getAnc() == NULL)
			{
				//std::cout << "root found\n";
				for (int i = 0; i < numSteps; i++)
				{
					likeVals[i] = clA[2*i];
				}

				double val = trapInt(likeVals);

				// NOTE: Prevents negative FFT values from producing NaN in lnL
				if (val <= 0.0 || std::isnan(val)) val = pow(10,-300);

				lnL += log(val); // likelihood of data observed for patron

				// std::cout << std::setprecision(24) << "\t\t\tLIKELIHOOD: " << val << "\n";
				if (std::isnan(log(val)))
				{
					std::cerr << "ERROR: val = NaN\n";
					std::cerr << "\tl:    " << lambda << "\n";
					std::cerr << "\ts:    " << sigma << "\n";
					std::cerr << "\tt1:   " << tauL << "\n";
					std::cerr << "\tt2:   " << tauR << "\n";
			//		std::cerr << "\ttip1: " << t1 << "\n";
			//		std::cerr << "\ttip2: " << t2 << "\n";

				}

				//std::cout << lnL << "\n";
			}
		}
	}

	// return lnL for all datapoints for a patron at the given table
	return lnL;
}
*/

double Model::scaleIFFT(double T, double s, int a)
{
	return (PI * T) / (s * a);
}


gsl_complex Model::complexPdf(double T, double s, double t, double a)
{

	gsl_complex denom;
	gsl_complex power;

	GSL_SET_COMPLEX(&denom, 1.0, a*T);
	GSL_SET_COMPLEX(&power, -(T*T*s*s/2), a*T);

	/*
	power = gsl_complex_exp(power);
	power = gsl_complex_mul_real(power, t);
	power = gsl_complex_div(power, denom);
	power = gsl_complex_mul_real(gsl_complex_exp(power), exp(-t));
	*/

	power = gsl_complex_exp(power);
	power = gsl_complex_mul_real(power, t);
	power = gsl_complex_div(power, denom);
	power = gsl_complex_sub_real(power,t);
	power = gsl_complex_exp(power);
	//power = gsl_complex_mul_real(gsl_complex_exp(power), exp(-t));

	if (std::isnan(GSL_REAL(power)))
	{
		std::cerr << "ERROR: complexPdf (" << &power << ") = NaN\n";
		std::cerr << "\tT: " << T << "\n";
		std::cerr << "\ts: " << s << "\n";
		std::cerr << "\tt: " << t << "\n";
		std::cerr << "\ta: " << a << "\n";

	}

	// std::cout << "\t\t\tcomplexPdf:\t" << GSL_REAL(power) << " " << GSL_IMAG(power) << "\n";

	return power;
	//return gsl_complex_mul_real(gsl_complex_exp(power), exp(-t));
}

// TODO:	Allow charFunc to represent an arbitrary jump kernel and an arbitrary drift model (e.g. Poisson & Normal)
//			This means the method will allow for an arbitrary number of parameters as required by the supplied fns.


// T	theta
// sb	sigma, Brownian Motion
// sk	sigma, Poisson jump kernel
// lk	lambda, Poisson jump kernel
gsl_complex Model::complexPdf(double T, double t, double sb, double sk, double lk)
{
	double tsk = t * sk;
	double tsb = t * sb;
	gsl_complex val = gsl_complex_rect((-0.5) * tsk * tsk * T * T, 0.0);
	val = gsl_complex_exp(val);
	val = gsl_complex_mul_real(val, lk);
	val = gsl_complex_add(gsl_complex_rect((-0.5) * tsb * tsb * T * T, 0.0), val);

	return gsl_complex_exp(val);

	//return gsl_complex_exp(val);
}

gsl_complex Model::charFunc(double T, double t, const std::vector<double>& parms)
{

	// TODO: Implement modelType selection using inline functions instead of all these control blocks.
	//for (unsigned int i = 0; i < parms.size(); i++) std::cout << parms[i] << "\t";
	//std::cout << "\n";

	gsl_complex driftCf = gsl_complex_rect(0.0, 0.0);
	gsl_complex jumpCf = gsl_complex_rect(0.0, 0.0);

	// Drift process
	if (modelType == 0 || modelType == 5)
		driftCf = gsl_complex_rect((-0.5) * t * parms[0] * parms[0] * T * T, 0.0);

	// Jump process
	if (modelType == 0) // Poisson-Normal + BM
	{
		// parms[0] = Sigma-BM
		// parms[1] = Lambda-K
		// parms[2] = Sigma-K

		jumpCf = gsl_complex_rect((-0.5) * parms[2] * parms[2] * T * T, 0.0);
		jumpCf = gsl_complex_sub_real(gsl_complex_exp(jumpCf), 1.0);
		jumpCf = gsl_complex_mul_real(jumpCf, t * parms[1]);
	}
	else if (modelType == 1) // Poisson-Cauchy + BM
	{
		proposalProbs.push_back(1.0);			// Sigma-BM
		proposalProbs.push_back(1.0);			// Lambda-K
		proposalProbs.push_back(1.0);			// Gamma-K
	}
	else if (modelType == 2) // Poisson-Exponential + BM
	{
		proposalProbs.push_back(1.0);			// Sigma-BM
		proposalProbs.push_back(1.0);			// Lambda-K
		proposalProbs.push_back(1.0);			// a-K

	}
	else if (modelType == 3) // Alpha-stable
	{
		 // a
									// pspspspspeeeu
	}
	else if (modelType == 4) // Chaix et al.
	{
		// parms[0] = Sigma-N
		// parms[1] = Lambda-N

		gsl_complex denom = gsl_complex_rect(1.0, parms[1] * T);
		gsl_complex power = gsl_complex_rect((-0.5) * T * T * parms[0] * parms[0], parms[1] * T);

		//GSL_SET_COMPLEX(&denom, 1.0, a*T);
		//GSL_SET_COMPLEX(&power, -(T*T*s*s/2), a*T);

		/*
		power = gsl_complex_exp(power);
		power = gsl_complex_mul_real(power, t);
		power = gsl_complex_div(power, denom);
		power = gsl_complex_mul_real(gsl_complex_exp(power), exp(-t));
		*/

		power = gsl_complex_exp(power);
		power = gsl_complex_mul_real(power, t);
		power = gsl_complex_div(power, denom);
		power = gsl_complex_sub_real(power,t);
		power = gsl_complex_exp(power);
		//power = gsl_complex_mul_real(gsl_complex_exp(power), exp(-t));

		/*
		if (std::isnan(GSL_REAL(power)))
		{
			std::cerr << "ERROR: complexPdf (" << &power << ") = NaN\n";
			std::cerr << "\tT: " << T << "\n";
			std::cerr << "\ts: " << s << "\n";
			std::cerr << "\tt: " << t << "\n";
			std::cerr << "\ta: " << a << "\n";

		}
		*/
		// std::cout << "\t\t\tcomplexPdf:\t" << GSL_REAL(power) << " " << GSL_IMAG(power) << "\n";

		return power;
	}

	else if (modelType == 5) // Brownian motion, i.e. no jump process
	{
		;
	}

	else if (modelType == 6) // Pure jump process
	{
		// parms[0] = Lambda-K
		// parms[1] = Sigma-K

		jumpCf = gsl_complex_rect((-0.5) * parms[1] * parms[1] * T * T, 0.0);
		jumpCf = gsl_complex_sub_real(gsl_complex_exp(jumpCf), 1.0);
		jumpCf = gsl_complex_mul_real(jumpCf, t * parms[0]);
	}

	else
	{
		std::cerr << "ERROR: Unknown model type: " << modelType << "\n";
		exit(1);
		return gsl_complex_rect(0.0, 0.0);
	}

	return gsl_complex_exp(gsl_complex_add(driftCf, jumpCf));
}


double Model::trapInt(double* fn)
{
	double val = 0.0;
	double x = 0.0;
	// NOTE: This trapezoidal integration is works for FFT indexed real values.

	//std::cout << "ReStepSize: " << ReStepSize << "\n";
	//std::cout << std::setprecision(256) << std::setw(256) << "\n";
	for (int i = 1; i < halfSteps; i++)
	{
		x = (fn[i] + fn[i-1]) / 2;
		val += x;
		//std::cout << "\t\t\t\tval[" << i << "]:\t" << x << "\t" << log(x) << "\t" << val << "\t" << fn[i] << "\t" << fn[i-1] << "\n";
	}
	for (int i = halfSteps + 1; i < numSteps; i++)
	{
		x = (fn[i] + fn[i-1]) / 2;
		val += x;
		//std::cout << "\t\t\t\tval[" << i << "]:\t" << x << "\t" << log(x) << "\t" << val << "\t" << fn[i] << "\t" << fn[i-1] << "\n";
	}
	x = (fn[0] + fn[numSteps-1]) / 2;
	val += x;
	//std::cout << "\t\t\t\tval[" << numSteps << ":\t" << x << "\t" << log(x) << "\t" << val << "\t" << fn[0] << "\t" << fn[numSteps-1] << "\n";

	return val * ReStepSize;
}



Parm* Model::pickParmAtRandom(Table* t)
{

	double u = randomPtr->uniformRv();
	double sum = 0.0;
	for (int i = 0; i < (int)proposalProbs.size(); i++)
	{
		sum += proposalProbs[i];
		if (u < sum)
		{
			return (t->getParmVector())[i];
		}
	}
	return NULL;

}

Parm* Model::pickBranchAtRandom(Table* t)
{

	int i = numBranches * randomPtr->uniformRv();
	return t->getBranchVector()[i];
}

Node* Model::pickNodeByBrLen(Table* t)
{
	double u = randomPtr->uniformRv();
	double sum = 0.0;
	int n = 0;
	Node* p;

	while (sum < u && n < numNodes)
	{
		p = topologyPtr->getNode(n);
		sum += p->getRatio();
		n++;
	}

	if (sum < u)
		return p;
	else
		return NULL;
}

Table* Model::pickTableAtRandom()
{

	// weight probability of selecting a table by the table.#patrons ^ -0.5
	std::vector<double> tableProbs;
	std::list<Table*>::iterator it_t;
	for (it_t = tableList.begin(); it_t != tableList.end(); it_t++)
	{
		tableProbs.push_back( pow( (*it_t)->getPatronList().size(), -0.5) );
	}

	double u = randomPtr->uniformRv();
	double norm = 0.0;
	double sum = 0.0;
	int r_i = -1;
	for (unsigned int i = 0; i < tableProbs.size(); i++)
	{
		norm += tableProbs[i];
	}
	for (unsigned int i = 0; i < tableProbs.size(); i++)
	{
		tableProbs[i] = tableProbs[i] / norm;
	}
	for (unsigned int i = 0; i < tableProbs.size(); i++)
	{
		sum += tableProbs[i];
		if (u < sum)
		{
			r_i = i;
			break;
		}
	}

	it_t = tableList.begin();
	for (int i = 0; i < r_i; i++)
	{
		it_t++;
	}

	return (*it_t);

}


Conc* Model::getActiveConc()
{
	return alpha->getActiveConc();
}

void Model::keep(Parm* p)
{
	p->keep();
}

void Model::restore(Parm* p)
{
	p->restore();
}

void Model::printTables(void)
{
	std::cout << "PRINTING: Tables\n";
	std::cout << "\ttableList.size() = " << tableList.size() << "\n";
	for (std::list<Table*>::iterator it_t = tableList.begin(); it_t != tableList.end(); it_t++)
	{
		(*it_t)->print();
	}
	std::cout << "\n";
	std::cout << "Seating:\t[";
	for (std::list<Table*>::iterator it_t = tableList.begin(); it_t != tableList.end(); it_t++)
	{
		std::cout << " " << (*it_t)->getPatronList().size();
	}
	std::cout << " ]\n";
}

#if 0

// SPOOKY GRAVEYARD OF CODE

clsPtrUpL[n] = new double***[numNodes];
clsPtrUpL[n][0] = new double**[numNodes * numTranscripts];
clsPtrUpL[n][0][0] = new double*[numNodes * numTranscripts * numTimepoints];
clsPtrUpL[n][0][0][0] = new double[numNodes * numTranscripts * numTimepoints * numSteps * 2];

for (int i = 0; i < numNodes; i++)
{
	clsPtrUpL[n][i] = clsPtrUpL[n][i-1] + numTranscripts;
	for (int j = 0; j < numTranscripts; j++)
	{
		clsPtrUpL[n][i][j] = clsPtrUpL[n][i][j-1] + numTimepoints;
		for (int k = 0; k < numTimepoints; k++)
		{
			clsPtrUpL[n][i][j][k] = &tip1[n * oneClSize + i * oneNodeSize + j * oneLocusSize + k * oneFFTSize];
		}
	}
}

clsPtrUpR[n] = new double***[numNodes];
clsPtrUpR[n][0] = new double**[numNodes * numTranscripts];
clsPtrUpR[n][0][0] = new double*[numNodes * numTranscripts * numTimepoints];
clsPtrUpR[n][0][0][0] = new double[numNodes * numTranscripts * numTimepoints * numSteps * 2];

for (int i = 0; i < numNodes; i++)
{
	clsPtrUpR[n][i] = clsPtrUpR[n][i-1] + numTranscripts;
	for (int j = 0; j < numTranscripts; j++)
	{
		clsPtrUpR[n][i][j] = clsPtrUpR[n][i][j-1] + numTimepoints;
		for (int k = 0; k < numTimepoints; k++)
		{
			clsPtrUpR[n][i][j][k] = &tip2[n * oneClSize + i * oneNodeSize + j * oneLocusSize + k * oneFFTSize];
		}
	}
}

#endif

#if 0

	int oneFFTSize = numSteps * 2;
	int oneLocusSize = numTimepoints * oneFFTSize;
	int oneNodeSize = numTranscripts * oneLocusSize;
	//int numNodes = 1;
	int oneClSize = numNodes * oneNodeSize; // NOTE: For general topologies.

	// std::cout << "oneClSize:" << oneClSize << "\n";

	// init data frames
	theta = new double[numSteps];
	like = new double[2*oneClSize];
	tip1 = new double[2*oneClSize];
	tip2 = new double[2*oneClSize];

	for (int i = 0; i < oneClSize; i++) // oneClSize is multiplied by spaces=2 and divided by (Re,Im)=2
	{
		REAL(like,i) = 0.0;
		IMAG(like,i) = 0.0;
		REAL(tip1,i) = 0.0;
		IMAG(tip1,i) = 0.0;
		REAL(tip2,i) = 0.0;
		IMAG(tip2,i) = 0.0;
		// std::cout << "i:\t" << i << " " << REAL(like,i) << " " << IMAG(like,i) << "\n";
	}

	for (int i = 0; i < numSteps; i++)
	{
		theta[i] = (((i + halfSteps) % numSteps) - halfSteps) * stepSize;
	}

	std::cout << "Data frames initialized.\n";

	// initialize pointer information
	for (int n = 0; n < 2; n++)
	{
		clsPtr[n] = new double****[numNodes];
		clsPtr[n][0] = new double***[numNodes * numTranscripts];
		clsPtr[n][0][0] = new double**[numNodes * numTranscripts * numTimepoints];
		clsPtr[n][0][0][0] = new double*[numNodes * numTranscripts * numTimepoints * numSteps * 2];

		//for (int i = 0; i < numNodes; i++)
		//for (int i = 0; i < 1; i++)
		//{
			int i = 0;
			//clsPtr[n][0] = clsPtr[n][i-1] + numTranscripts;

			for (int j = 1; j < numTranscripts; j++)
			{
				std::cout << "j: " << j << "\n";
				clsPtr[n][i][j] = clsPtr[n][i][j-1] + numTimepoints;
				std::cout << "ok\n";
				for (int k = 1; k < numTimepoints; k++)
				{
					std::cout << "k: " << k << "\n";
					clsPtr[n][i][j][k] = clsPtr[n][i][j][k-1] + (numSteps * 2);
					std::cout << "\t" << n * i * j * k << "\n";
				}
			}
		//}

		//std::cout << "test: " << like[46] << "\n";
		for (int i = 0; i < numNodes; i++)
		{
			for (int j = 0; j < numTranscripts; j++)
			{
				for (int k = 0; k < numTimepoints; k++)
				{
					for(int l = 0; l < numSteps; l++)
					{
						//std::cout << "test: " << like[46] << "\n";
						std::cout << "\t" << n << "\t" << i << "\t" << j << "\t" << k << "\t" << l << "\t";
						std::cout << n * oneClSize + i * oneNodeSize + j * oneLocusSize + k * oneFFTSize + l*2 << "\t";
						// clsPtr[n][i][j][k][l] = 0.0;
						std::cout << " = Re:" << like[n * oneClSize + i * oneNodeSize + j * oneLocusSize + k * oneFFTSize + l*2] << "\t";
						std::cout << " = Im:" << like[n * oneClSize + i * oneNodeSize + j * oneLocusSize + k * oneFFTSize + l*2 + 1] << "\t";
						std::cout << "?\n";
						clsPtr[n][i][j][k][l] = &like[n * oneClSize + i * oneNodeSize + j * oneLocusSize + k * oneFFTSize + l*2];
						//std::cout << ".\n";
					}
				}
			}
		}

		/*
		 * for (int i = 0; i < numNodes; i++)
			for (int j = 0; j < numSites; j++)
				clsPtrUpL[n][i][j] = &clsUpL[n * oneClSize + i * oneNodeSize + j * numStates];
		 *
		 */


		clsPtrUpL[n] = new double****[numNodes];
		clsPtrUpL[n][0] = new double***[numNodes * numTranscripts];
		clsPtrUpL[n][0][0] = new double**[numNodes * numTranscripts * numTimepoints];
		clsPtrUpL[n][0][0][0] = new double*[numNodes * numTranscripts * numTimepoints * numSteps * 2];

		for (int i = 1; i < numNodes; i++)
		{
			clsPtrUpL[n][i] = clsPtrUpL[n][i-1] + numTranscripts;
			for (int j = 1; j < numTranscripts; j++)
			{
				clsPtrUpL[n][i][j] = clsPtrUpL[n][i][j-1] + numTimepoints;
				for (int k = 1; k < numTimepoints; k++)
				{
					clsPtrUpL[n][i][j][k] = clsPtrUpL[n][i][j][k-1] + (numSteps * 2);
				}
			}
		}

		for (int i = 0; i < numNodes; i++)
			for (int j = 0; j < numTranscripts; j++)
				for (int k = 0; k < numTimepoints; k++)
					for(int l = 0; l < (2 * numSteps); l++)
						clsPtrUpL[n][i][j][k][l] = &tip1[n * oneClSize + i * oneNodeSize + j * oneLocusSize + k * oneFFTSize + l];


		clsPtrUpR[n] = new double****[numNodes];
		clsPtrUpR[n][0] = new double***[numNodes * numTranscripts];
		clsPtrUpR[n][0][0] = new double**[numNodes * numTranscripts * numTimepoints];
		clsPtrUpR[n][0][0][0] = new double*[numNodes * numTranscripts * numTimepoints * numSteps * 2];

		for (int i = 1; i < numNodes; i++)
		{
			clsPtrUpR[n][i] = clsPtrUpR[n][i-1] + numTranscripts;
			for (int j = 1; j < numTranscripts; j++)
			{
				clsPtrUpR[n][i][j] = clsPtrUpR[n][i][j-1] + numTimepoints;
				for (int k = 1; k < numTimepoints; k++)
				{
					clsPtrUpR[n][i][j][k] = clsPtrUpR[n][i][j][k-1] + (numSteps * 2);
				}
			}
		}

		for (int i = 0; i < numNodes; i++)
			for (int j = 0; j < numTranscripts; j++)
				for (int k = 0; k < numTimepoints; k++)
					for(int l = 0; l < (2 * numSteps); l++)
						clsPtrUpR[n][i][j][k][l] = &tip2[n * oneClSize + i * oneNodeSize + j * oneLocusSize + k * oneFFTSize + l];





	}

	std::cout << "ok\n";

	// initialize tip conditional likelihoods
	// NOTE: done statically because we have no topology
	//for (int i = 0; i < numTaxa; i++)
	for (int i = 0; i < 1; i++)
	{
		std::cout << "i\n";
		for (int j = 0; j < numTranscripts; j++)
		{
			std::cout << "j\n";
			for (int k = 0; k < numTimepoints; k++)
			{
				std::cout << "k\n";
				//double fpkm = expressionPtr->getExpr(j,i,k);
				//std::cout << "\tfpkm:" << fpkm << "\n";

				for(int l = 0; l < (2 * numSteps); l++)
				{
					std::cout << "l\n";
					double* p0 = clsPtr[0][i][j][k][l];
					double* p1 = clsPtr[1][i][j][k][l];
					double* uL0 = clsPtrUpL[0][i][j][k][l];
					double* uL1 = clsPtrUpL[1][i][j][k][l];
					double* uR0 = clsPtrUpR[0][i][j][k][l];
					double* uR1 = clsPtrUpR[1][i][j][k][l];
				}
				// place observed FPKM into FFT bins



				/*
				MbBitfield *bf = alignmentPtr->getCodonMatrixEntry(i, j);
				for (int s = 0; s < alignmentPtr->getNumSenseCodons(); s++)
				{
					if (bf->isBitSet(s) == true)
					{
						p0[s] = 1.0;
						p1[s] = 1.0;
						uL0[s] = 1.0;
						uL1[s] = 1.0;
						uR0[s] = 1.0;
						uR1[s] = 1.0;
					}
				}
				*/
			}
		}
	}
#endif

#if 0

// table likelihood 05/03/11
double Model::tableLogLikelihood(Table* t)
{

	double lnL = 0.0;

	// Read in the parameter values.
	double lambda = 0.0;
	//double sigma = 1.0;
	//double tau1 = 0.0;
	//double tau2 = 0.0;

	const std::vector<Parm*> tempParmVector = t->getParmVector();
	//bool isRight = false;
	for (std::vector<Parm*>::const_iterator it_p = tempParmVector.begin(); it_p != tempParmVector.end(); it_p++)
	{
		Lambda* derivedPtrL = dynamic_cast<Lambda*> (*it_p);
		//Tau* derivedPtrT = dynamic_cast<Tau*> (*it_p);
		// Sigma* derivedPtrS = dynamic_cast<Sigma*> (*it_p);
		if (derivedPtrL != 0)
			lambda = derivedPtrL->getActiveExpMean()->getRate();
		/*
		else if (derivedPtrS != 0)
			sigma = derivedPtrS->getActiveVar()->getRate();
		else if (derivedPtrT != 0 && isRight == false)
		{
			tau1 = derivedPtrT->getActiveBrLen()->getRate();
			isRight = true;
		}
		else if (derivedPtrT != 0)
		{
			tau2 = derivedPtrT->getActiveBrLen()->getRate();
		}
		*/
	}

	// std::cout << "\t\t" << lambda << " " << sigma << " " << tau1 << " " << tau2 << "\n";

	// Calculate likelihoods for all Patrons of a given table (i.e. labeled to a different Parm classes)
	for (std::list<Patron*>::const_iterator it_p = t->getPatronList().begin(); it_p != t->getPatronList().end(); it_p++)
	{
		double t1 = ((*it_p)->getData())[0][0];
		double t2 = ((*it_p)->getData())[1][0];

		double lnLike1 = randomPtr->lnNormalPdf(lambda, 1.0, t1);
		double lnLike2 = randomPtr->lnNormalPdf(lambda, 1.0, t2);

		lnL += lnLike1 + lnLike2;
	}

	return lnL;
}


// locus likelihood 05/03/11
double Model::locusLogLikelihood(Patron* p)
{
	double lnL = 0.0;

	// Read in the parameter values.
	double lambda = 0.0;
	//double sigma = 1.0;
	//double tau1 = 0.0;
	//double tau2 = 0.0;

	const std::vector<Parm*> tempParmVector = p->getTable()->getParmVector();
	//bool isRight = false;
	for (std::vector<Parm*>::const_iterator it_p = tempParmVector.begin(); it_p != tempParmVector.end(); it_p++)
	{
		Lambda* derivedPtrL = dynamic_cast<Lambda*> (*it_p);
		//Tau* derivedPtrT = dynamic_cast<Tau*> (*it_p);
		// Sigma* derivedPtrS = dynamic_cast<Sigma*> (*it_p);
		if (derivedPtrL != 0)
			lambda = derivedPtrL->getActiveExpMean()->getRate();
		/*
		else if (derivedPtrS != 0)
			sigma = derivedPtrS->getActiveVar()->getRate();
		else if (derivedPtrT != 0 && isRight == false)
		{
			tau1 = derivedPtrT->getActiveBrLen()->getRate();
			isRight = true;
		}
		else if (derivedPtrT != 0)
		{
			tau2 = derivedPtrT->getActiveBrLen()->getRate();
		}
		*/
	}

	// std::cout << "\t\t" << lambda << " " << sigma << " " << tau1 << " " << tau2 << "\n";

	// Calculate likelihoods for all Patrons of a given table (i.e. labeled to a different Parm classes)

	double t1 = (p->getData())[0][0];
	double t2 = (p->getData())[1][0];

	double lnLike1 = randomPtr->lnNormalPdf(lambda, 1.0, t1);
	double lnLike2 = randomPtr->lnNormalPdf(lambda, 1.0, t2);



	//double like1 = randomPtr->normalPdf(lambda, 1.0, t1 - t2);

	//std::cout << "like: " << lnLike1 + lnLike2 << "\n";

	lnL += lnLike1 + lnLike2;

	return lnL;
}

#endif

#if 0
double Model::tableLogLikelihood(Table* t)
{
	double lnL = 0.0;

	// Read in the parameter values.
	double lambda = 0.0;
	double sigma = 0.0;
	double tau1 = 0.0;
	double tau2 = 0.0;

	const std::vector<Parm*> tempParmVector = t->getParmVector();
	bool isRight = false;
	for (std::vector<Parm*>::const_iterator it_p = tempParmVector.begin(); it_p != tempParmVector.end(); it_p++)
	{
		Lambda* derivedPtrL = dynamic_cast<Lambda*> (*it_p);
		Tau* derivedPtrT = dynamic_cast<Tau*> (*it_p);
		Sigma* derivedPtrS = dynamic_cast<Sigma*> (*it_p);
		if (derivedPtrL != 0)
			lambda = derivedPtrL->getActiveExpMean()->getRate();
		else if (derivedPtrS != 0)
			sigma = derivedPtrS->getActiveVar()->getRate();
		else if (derivedPtrT != 0 && isRight == false)
		{
			tau1 = derivedPtrT->getActiveBrLen()->getRate();
			isRight = true;
		}
		else if (derivedPtrT != 0)
		{
			tau2 = derivedPtrT->getActiveBrLen()->getRate();
		}
	}

#if LAMBDA_ONLY
	// TEST 05/06/11 - LAMBDA
	sigma = 0.1;
	tau1 = 2.0;
	tau2 = 0.5;
#elif SIGMA_ONLY
	// TEST 05/06/11 - SIGMA
	lambda = 0.5;
	tau1 = 2.0;
	tau2 = 0.5;
#elif TAU1_ONLY
	// TEST 05/06/11 - TAU1
	lambda = 0.5;
	sigma = 0.1;
	tau2 = 0.5;
#elif TAU2_ONLY
	// TEST 05/06/11 - TAU2
	lambda = 0.5;
	sigma = 0.1;
	tau1 = 2.0;
#endif

	// std::cout << "\t\t" << lambda << " " << sigma << " " << tau1 << " " << tau2 << "\n";

	// Calculate likelihoods for all Patrons of a given table (i.e. labeled to a different Parm classes)
	for (std::list<Patron*>::const_iterator it_p = t->getPatronList().begin(); it_p != t->getPatronList().end(); it_p++)
	{

		int i = 0;

		// set data frame to function vals
		for (i = 0; i < numSteps; i++)
		{
			gsl_complex complex_pdf1 = complexPdf(theta[i], sigma, tau1, lambda);
			gsl_complex complex_pdf2 = complexPdf(theta[i], sigma, tau2, lambda);
			gsl_complex complex_tip1 = (*it_p)->getCfData(0,0,i);
			gsl_complex complex_tip2 = (*it_p)->getCfData(1,0,i);
			gsl_complex conv_tip1 = gsl_complex_mul(complex_pdf1, complex_tip1);
			gsl_complex conv_tip2 = gsl_complex_mul(complex_pdf2, complex_tip2);

			//std::cout << std::setprecision(8) << "\tFn" << i << ":\t" << "tip1\t" << GSL_REAL(complex_pdf1) << " " << GSL_IMAG(complex_pdf1) << "\n";
			//std::cout << std::setprecision(8) << "\tCh" << i << ":\t" << "tip1\t" << GSL_REAL(complex_tip1) << " " << GSL_IMAG(complex_tip1) << "\n";
			//std::cout << std::setprecision(8) << "\tCV" << i << ":\t" << "tip1\t" << GSL_REAL(conv_tip1) << " " << GSL_IMAG(conv_tip1) << "\n";
			//std::cout << std::setprecision(8) << "\tFn" << i << ":\t" << "tip2\t" << GSL_REAL(complex_pdf2) << " " << GSL_IMAG(complex_pdf2) << "\n";
			//std::cout << std::setprecision(8) << "\tCh" << i << ":\t" << "tip2\t" << GSL_REAL(complex_tip2) << " " << GSL_IMAG(complex_tip2) << "\n";
			//std::cout << std::setprecision(8) << "\tCV" << i << ":\t" << "tip2\t" << GSL_REAL(conv_tip2) << " " << GSL_IMAG(conv_tip2) << "\n";

			REAL(tip1,i) = GSL_REAL(conv_tip1);
			IMAG(tip1,i) = GSL_IMAG(conv_tip1);
			REAL(tip2,i) = GSL_REAL(conv_tip2);
			IMAG(tip2,i) = GSL_IMAG(conv_tip2);

			//std::cout << "\n";
		}

		gsl_fft_complex_radix2_inverse(tip1, 1, numSteps);
		gsl_fft_complex_radix2_inverse(tip2, 1, numSteps);

		for (i = 0; i < numSteps; i++)
		{
			gsl_complex complex_tip1;
			GSL_SET_COMPLEX(&complex_tip1, REAL(tip1,i), IMAG(tip1,i));

			gsl_complex complex_tip2;
			GSL_SET_COMPLEX(&complex_tip2, REAL(tip2,i), IMAG(tip2,i));

			gsl_complex final_val = gsl_complex_mul(complex_tip1, complex_tip2);

			REAL(like,i) = GSL_REAL(final_val);
			IMAG(like,i) = GSL_IMAG(final_val);
			likeVals[i] = GSL_REAL(final_val);

			//std::cout << GSL_REAL(final_val) << " " << GSL_IMAG(final_val) << "\n";
			//std::cout << std::setprecision(8) << "\t" << i << ":\t" << "tip1\t" << GSL_REAL(complex_tip1) << " " << GSL_IMAG(complex_tip1) << "\n";
			//std::cout << std::setprecision(8) << "\t\ttip2\t" << GSL_REAL(complex_tip2) << " " << GSL_IMAG(complex_tip2) << "\n";
			//std::cout  << "\t\tfinal\t";// << GSL_REAL(final_val) << "\n"; // << " " << GSL_IMAG(final_val) << "\n";
		}

		double val = trapInt(likeVals);

		// NOTE: Prevents negative FFT values from producing NaN in lnL
		if (val <= 0.0 || std::isnan(val)) val = pow(10,-300);

		// std::cout << "\t\t\tLIKELIHOOD: " << val << "\n";
		//std::cout << "\t\t\tLIKELIHOOD: " << val*val << "\n";

		lnL += log(val); // likelihood of data observed for patron
	}

	// return lnL for all datapoints for a patron at the given table
	return lnL;
}

double Model::locusLogLikelihood(Patron* p)
{
	double lnL = 0.0;

	// Read in the parameter values.
	double lambda = 0.0;
	double sigma = 0.0;
	double tau1 = 0.0;
	double tau2 = 0.0;

	Table* t = p->getTable();
	const std::vector<Parm*> tempParmVector = t->getParmVector();
	bool isRight = false;
	for (std::vector<Parm*>::const_iterator it_p = tempParmVector.begin(); it_p != tempParmVector.end(); it_p++)
	{
		Lambda* derivedPtrL = dynamic_cast<Lambda*> (*it_p);
		Tau* derivedPtrT = dynamic_cast<Tau*> (*it_p);
		Sigma* derivedPtrS = dynamic_cast<Sigma*> (*it_p);
		if (derivedPtrL != 0)
			lambda = derivedPtrL->getActiveExpMean()->getRate();
		else if (derivedPtrS != 0)
			sigma = derivedPtrS->getActiveVar()->getRate();
		else if (derivedPtrT != 0 && isRight == false)
		{
			tau1 = derivedPtrT->getActiveBrLen()->getRate();
			isRight = true;
		}
		else if (derivedPtrT != 0)
		{
			tau2 = derivedPtrT->getActiveBrLen()->getRate();
		}
	}

	if (std::isinf(sigma))	std::cerr << "ERROR: sigma = inf\n";
	if (std::isinf(lambda)) std::cerr << "ERROR: lambda = inf\n";
	if (std::isinf(tau1)) std::cerr << "ERROR: tau1 = inf\n";
	if (std::isinf(tau2)) std::cerr << "ERROR: tau2 = inf\n";

	// std::cout << "\t\t" << lambda << " " << sigma << " " << tau1 << " " << tau2 << "\n";

	// Calculate likelihoods for single Patron (labeled to a different Parm class)
	int i;

	double t1 = (p->getData())[0][0];
	double t2 = (p->getData())[1][0];

	// set data frame to function vals
	for (i = 0; i < numSteps; i++)
	{
		gsl_complex complex_pdf1 = complexPdf(theta[i], sigma, tau1, lambda);
		gsl_complex complex_pdf2 = complexPdf(theta[i], sigma, tau2, lambda);
		gsl_complex complex_tip1 = p->getCfData(0,0,i);
		gsl_complex complex_tip2 = p->getCfData(1,0,i);
		gsl_complex conv_tip1 = gsl_complex_mul(complex_pdf1, complex_tip1);
		gsl_complex conv_tip2 = gsl_complex_mul(complex_pdf2, complex_tip2);



		//std::cout << std::setprecision(8) << "\tFn" << i << ":\t" << "tip1\t" << GSL_REAL(complex_pdf1) << " " << GSL_IMAG(complex_pdf1) << "\n";
		//std::cout << std::setprecision(8) << "\tCh" << i << ":\t" << "tip1\t" << GSL_REAL(complex_tip1) << " " << GSL_IMAG(complex_tip1) << "\n";
		//std::cout << std::setprecision(8) << "\tCV" << i << ":\t" << "tip1\t" << GSL_REAL(conv_tip1) << " " << GSL_IMAG(conv_tip1) << "\n";
		//std::cout << std::setprecision(8) << "\tFn" << i << ":\t" << "tip2\t" << GSL_REAL(complex_pdf2) << " " << GSL_IMAG(complex_pdf2) << "\n";
		//std::cout << std::setprecision(8) << "\tCh" << i << ":\t" << "tip2\t" << GSL_REAL(complex_tip2) << " " << GSL_IMAG(complex_tip2) << "\n";
		//std::cout << std::setprecision(8) << "\tCV" << i << ":\t" << "tip2\t" << GSL_REAL(conv_tip2) << " " << GSL_IMAG(conv_tip2) << "\n";

		REAL(tip1,i) = GSL_REAL(conv_tip1);
		IMAG(tip1,i) = GSL_IMAG(conv_tip1);
		REAL(tip2,i) = GSL_REAL(conv_tip2);
		IMAG(tip2,i) = GSL_IMAG(conv_tip2);

		if (std::isnan(likeVals[i]))
		{
			std::cerr << "ERROR: likeVals[" << i << "] = NaN\n";
		}

		if (std::isnan(GSL_REAL(conv_tip1))) std::cerr << "ERROR: conv_tip1[" << i << "] = NaN\n";
		if (std::isnan(GSL_REAL(conv_tip2))) std::cerr << "ERROR: conv_tip1[" << i << "] = NaN\n";


		//std::cout << "\n";
	}

	gsl_fft_complex_radix2_inverse(tip1, 1, numSteps);
	gsl_fft_complex_radix2_inverse(tip2, 1, numSteps);

	for (i = 0; i < numSteps; i++)
	{
		gsl_complex complex_tip1;
		GSL_SET_COMPLEX(&complex_tip1, REAL(tip1,i), IMAG(tip1,i));

		gsl_complex complex_tip2;
		GSL_SET_COMPLEX(&complex_tip2, REAL(tip2,i), IMAG(tip2,i));

		gsl_complex final_val = gsl_complex_mul(complex_tip1, complex_tip2);

		REAL(like,i) = GSL_REAL(final_val);
		IMAG(like,i) = GSL_IMAG(final_val);
		likeVals[i] = GSL_REAL(final_val);

		if (std::isnan(likeVals[i]))
		{
			std::cerr << "ERROR: likeVals[" << i << "] = NaN\n";
		}

		// std::cout << std::setprecision(8) << "\t" << i << ":\t" << "tip1\t" << GSL_REAL(complex_tip1) << " " << GSL_IMAG(complex_tip1) << "\n";
		// std::cout << std::setprecision(8) << "\t\ttip2\t" << GSL_REAL(complex_tip2) << " " << GSL_IMAG(complex_tip2) << "\n";
		// std::cout << std::setprecision(8) << "\t\tfinal\t" << GSL_REAL(final_val) << " " << GSL_IMAG(final_val) << "\n";
	}

	double val = trapInt(likeVals);

	// NOTE: Prevents negative FFT values from producing NaN in lnL
	if (val <= 0.0 || std::isnan(val)) val = pow(10,-300);

	lnL += log(val); // likelihood of data observed for patron

	// std::cout << std::setprecision(24) << "\t\t\tLIKELIHOOD: " << val << "\n";
	if (std::isnan(log(val)))
	{
		std::cerr << "ERROR: val = NaN\n";
		std::cerr << "\tl:    " << lambda << "\n";
 		std::cerr << "\ts:    " << sigma << "\n";
		std::cerr << "\tt1:   " << tau1 << "\n";
		std::cerr << "\tt2:   " << tau2 << "\n";
		std::cerr << "\ttip1: " << t1 << "\n";
		std::cerr << "\ttip2: " << t2 << "\n";

	}


	// return lnL for all datapoints for a patron at the given table
	return lnL;
}
#endif

#if 0

void Model::initializeCondLikes(void)
{
	// [space][node][locus][timepoint][step]
	condLikes = new CondLikes(expressionPtr, settingsPtr);


/*
	// Assigns tip likelihoods equal to the observed expression levels under the characteristic function.
	for (std::list<Patron*>::iterator it_p = patronList.begin(); it_p != patronList.end(); it_p++)
	{
		// transcript index
		int j = (*it_p)->getId();

		//double* cl_real = condLikes->getCls(1,i,j,k,0);
		for (int i = 0; i < numTaxa; i++)
		{
			for (int k = 0; k < numTimepoints; k++)
			{
				for (int l = 0; l < numSteps; l++)
				{
					//(*(cl_real + 2*l)) = GSL_REAL((*it_p)->getCfData(i,k,l));
					//(*(cl_real + 2*l + 1)) = GSL_IMAG((*it_p)->getCfData(i,k,l));

					double re_cf = GSL_REAL((*it_p)->getCfData(i,k,l));
					double im_cf = GSL_IMAG((*it_p)->getCfData(i,k,l));

					double* re_cl = condLikes->getCls(0,i,j,k,2*l);
					double* im_cl = condLikes->getCls(0,i,j,k,2*l+1);
					(*re_cl) = re_cf;
					(*im_cl) = im_cf;

					//double* cl_real2 = condLikes->getCls(1,i,j,k,2*l);
					//(*cl_real2) = GSL_REAL((*it_p)->getCfData(i,k,l));
					//double* cl_imag = condLikes->getCls(1,i,j,k,2*l+1);
					//(*cl_imag) = GSL_IMAG((*it_p)->getCfData(i,k,l));

				}
			}
		}
		std::cout << (*it_p)->getName() << " added.\n";
	}
*/



/*
	for (std::list<Patron*>::iterator it_p = patronList.begin(); it_p != patronList.end(); it_p++)
	{
		for (int i = 0; i < numTaxa; i++)
		{
			for (int j = 0; j < numTranscripts; j++)
			{
				for (int k = 0; k < numTimepoints; k++)
				{
					//double* cl_real = condLikes->getCls(1,i,j,k,0);
					for (int l = 0; l < numSteps; l++)
					{
						//(*(cl_real + 2*l)) = GSL_REAL((*it_p)->getCfData(i,k,l));
						//(*(cl_real + 2*l + 1)) = GSL_IMAG((*it_p)->getCfData(i,k,l));
						double re_cf = GSL_REAL((*it_p)->getCfData(i,k,l))

						double* cl_real2 = condLikes->getCls(1,i,j,k,2*l);
						(*cl_real2) = GSL_REAL((*it_p)->getCfData(i,k,l));
						double* cl_imag = condLikes->getCls(1,i,j,k,2*l+1);
						(*cl_imag) = GSL_IMAG((*it_p)->getCfData(i,k,l));

					}
				}
			}
		}
		//std::cout << (*it_p)->getName() << " added.\n";
	}
*/

/*
	int trans = 0;
	for (std::list<Patron*>::iterator it_p = patronList.begin(); it_p != patronList.end(); it_p++)
	{

		Patron* p = *it_p;
		double* clL = condLikes->getCls(0,0,trans,0,0);
				condLikes->
		for (int i = 0; i < numSteps; i++)
		{
			gsl_complex complex_tip_L1 = p->getCfData(0,0,i);
			gsl_complex complex_tip_R1 = p->getCfData(1,0,i);
			gsl_complex complex_tip_L2 = gsl_complex_rect(clL->at(i), clL->at(i+1));
			gsl_complex complex_tip_R2 = gsl_complex_rect(clR->at(i), clR->at(i+1));
			std::cout << "\t" << GSL_REAL(complex_tip_L1) << "\t" << GSL_IMAG(complex_tip_L1) << "\n";
			std::cout << "\t" << GSL_REAL(complex_tip_R1) << "\t" << GSL_IMAG(complex_tip_R1) << "\n";
			std::cout << "\t" << GSL_REAL(complex_tip_L2) << "\t" << GSL_IMAG(complex_tip_L2) << "\n";
			std::cout << "\t" << GSL_REAL(complex_tip_R2) << "\t" << GSL_IMAG(complex_tip_R2) << "\n";
		}
		trans++;
	}
*/


	// Tip observation initialization using pointer arithmetic. Should vastly improve the speed.
		/*
		for (std::list<Patron*>::iterator it_p = patronList.begin(); it_p != patronList.end(); it_p++)
		{
			int oneTimeSize = 2 * numSteps;
			int oneTransSize = numTimepoints * oneTimeSize;
			int oneNodeSize = numTranscripts * oneTransSize;
			int oneTreeSize = numTaxa * oneNodeSize;
			int oneClSize = 2 * oneTreeSize;

			int locus = (*it_p)->getId();
			double* cl = condLikes->getCls(1, 0, locus, 0, 0);
			double* cf = (*it_p)->getCfData(0,0,0);
			for (int i = 0; i < numTaxa; i++)
			{
				for (int k = 0; k < numTimepoints; k++)
				{
					for (int l = 0; l < numSteps; l++)
					{
						cl[]

						cl[i * j + j * k + k * l + 1] = GSL_REAL(cf[i * k + k * l + l]);
						cl[i * j + j * k + k * l + l + 1] = GSL_IMAG(cf[i * k + k * l + l]);
						//(*(cl_real + 2*l)) = GSL_REAL((*it_p)->getCfData(i,k,l));
						//(*(cl_real + 2*l + 1)) = GSL_IMAG((*it_p)->getCfData(i,k,l));


						//double* cl_real2 = condLikes->getCls(1,i,j,k,2*l);
						//(*cl_real2) = GSL_REAL((*it_p)->getCfData(i,k,l));
						//double* cl_imag = condLikes->getCls(1,i,j,k,2*l+1);
						//(*cl_imag) = GSL_IMAG((*it_p)->getCfData(i,k,l));


						std::cout << cl) << "\t" << *(cl_imag) << "\n";
						std::cout << *(cl_real + 2*l) << "\t" << *(cl_real + 2*l + 1) << "\n";
						std::cout << "\n";

					}
				}
			}
			//std::cout << (*it_p)->getName() << " added.\n";
		}
		*/

	condLikes->print();

	/*
		double* x = condLikes->getCls(1,0,0,0,0);
		std::cout << "\t" << (*x) << "\n";
		(*x) = 1.0;
		std::cout << "\t" << (*x) << "\n";
		x = condLikes->getCls(1,0,0,0,0);
		std::cout << "\t" << (*x) << "\n";
	*/

	std::cout << "\n";
}
#endif

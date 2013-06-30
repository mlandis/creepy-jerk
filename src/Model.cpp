/*
 * Model.cpp
 *
 *  Created on: Mar 7, 2011
 *      Author: mlandis
 */

#include "Model.h"

#define DEBUG 0
#define DEBUG2 0
#define DEBUG3 0

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
//#define PI 3.1415926535

#define DEBUG_MSG_QAWO 0
#define DEBUG_QAWO 0
#define DEBUG_BESSEL 0

#define NORM_INV_GAUSS_BESSEL 0
#define ALPHA_STABLE_NUMINT 1
#define JUMP_NORM_NUMINT 2
#define JUMP_NORM_PDF 3
#define VAR_GAMMA_BESSEL 4
#define BM_ONLY 5
#define VAR_GAMMA_NUMINT 6
#define BM_JUMP_PDF 7
#define DBL_EXP_NUMINT 8
#define SKEW_NORMAL_NUMINT 9

#define PDF 0
#define NUMINT 1
#define BESSEL 2

double poissonNormal(double k, void* params)
{
	// params[0] == sig_bm
	// params[1] == lam_jn
	// params[2] == sig_jn
	// params[3] == t
	// params[4] == tuningBM
	// params[5] == x

	std::vector<double> p =	*(std::vector<double>*)params;
	double sig_bm = p[0];
	double lam_jn = p[1];
	double sig_jn = p[2];
	double t = p[3];
	double tuning_BM = p[4];
	double x = p[5];

	double left = -pow(sig_bm,2) * t * (1 - tuning_BM) * k * k / (2 * x * x);
	double right = t * lam_jn * ( exp(-pow(sig_jn,2) * k * k / (2 * x * x) ) - 1 );
	return exp(left + right);
}


double alphaStable(double k, void* params)
{
	// params[0] == sig_bm
	// params[1] == alpha
	// params[2] == c
	// params[3] == t
	// params[4] == tuningBM
	// params[5] == x

	std::vector<double> p =	*(std::vector<double>*)params;
	double sig_bm = p[0];
	double alpha = p[1];
	double c = p[2];
	double t = p[3];
	double tuning_BM = p[4];
	double x = p[5];

	///double left = -pow(sig_bm,2) * t * (1 - tuning_BM) * k * k / (2 * x * x);
	double right = -pow(fabs(c * k / x), alpha) * t;

	return exp(right);
}

// returns half of a skew normal
double halfSkewNormal(double k, void* params)
{
	// params[0] == sig_bm
	// params[1] == lam_jn
	// params[2] == sig_jn
	// params[3] == alpha
	// params[4] == t
	// params[5] == tuningBM
	// params[6] == x

	std::vector<double> p =	*(std::vector<double>*)params;
	double sig_bm = p[0];
	double alpha = p[1];
	double c = p[2];
	double t = p[3];
	double tuning_BM = p[4];
	double x = p[5];

	double left = -pow(sig_bm,2) * t * (1 - tuning_BM) * k * k / (2 * x * x);
	double right = -pow(fabs(c * k / x), alpha) * t;

	return exp(left + right);
}

double varianceGamma(double k, void* params)
{
	// params[0] == sig_bm
	// params[1] == kap_vg
	// params[2] == sig_vg
	// params[3] == t
	// params[4] == tuningBM
	// params[5] == x

	std::vector<double> p =	*(std::vector<double>*)params;
	double sig_bm = p[0];
	double kap_vg = p[1];
	double sig_vg = p[2];
	double t = p[3];
	double tuning_BM = p[4];
	double x = p[5];

	double left = -pow(sig_bm,2) * t * (1 - tuning_BM) * k * k / (2 * x * x);
	double right = pow(1.0 + sig_vg * sig_vg * k * k / (2 * x * x * kap_vg) , -kap_vg * t);

	return exp(left + right);
}

double doubleExponential(double k, void* params)
{
	// params[0] == sig_bm
	// params[1] == kap_vg
	// params[2] == sig_vg
	// params[3] == t
	// params[4] == tuningBM
	// params[5] == x

	std::vector<double> p =	*(std::vector<double>*)params;
	double sig_bm = p[0];
	double lam_de = p[1];
	double bee_de = p[2];
	double t = p[3];
	double tuning_BM = p[4];
	double x = p[5];

	double left = -pow(sig_bm,2) * t * (1 - tuning_BM) * k * k / (2 * x * x);
	double right = t * lam_de * pow(1 + k * k * bee_de * bee_de, -1) - 1;

	return exp(left + right);
}

Model::Model(Expression* ep, MbRandom* rp, Settings* sp, Topology* tp)
{

	expressionPtr = ep;
	randomPtr = rp;
	settingsPtr = sp;
	topologyPtr = tp;

	numTaxa = expressionPtr->getNumTaxa();
	numNodes = 2 * numTaxa - 1;
	numBranches = numNodes - 1;	// [numBranches+1 == numNodes] indexes the root, which has no branch.
	numTranscripts = expressionPtr->getNumTranscripts();
	numTimepoints = expressionPtr->getNumTimepoints();

	modelType = settingsPtr->getModelType();
	fixBranches = settingsPtr->getFixBranches();

	tableId = 0;

	tuningBM = settingsPtr->getTuningBM();
	useJumpKernel = settingsPtr->getUseJumpKernel();
	sigmaJumpProposal = settingsPtr->getSigmaJumpProposal();
	printStdOut = settingsPtr->getPrintStdOut();

	modelShiftList = settingsPtr->getModelShiftList();

	useSteppingStone = settingsPtr->getUseSteppingStone();
	betaSteppingStone = settingsPtr->getBetaSteppingStone();

	initializeParms();
	initializeTips();
	initializeGSL();
}

Model::~Model(void) {

	for (std::list<Table*>::iterator it_t = tableList.begin(); it_t != tableList.end(); it_t++)
	{
		delete (*it_t);
	}
	freeGSL();
}


void Model::initializeParms(void)
{
	if (printStdOut) std::cout << "INITIALIZING: Parameter Classes (Tables)\n";


	// create initial parameters
	std::vector<Parm*> initialParms;
	initialParms.clear();

	if (modelType == 0) // undef
	{
		;
	}
	else if (modelType == ALPHA_STABLE_NUMINT) // NumInt Alpha-stable + BM
	{
		initialParms.push_back(new Sigma(randomPtr, "sigma-AS"));
		initialParms.push_back(new Alpha(randomPtr, "alpha-AS"));
		initialParms.push_back(new Sigma(randomPtr, "beta-AS"));
		proposalProbs.push_back(1.0);
		proposalProbs.push_back(1.0);
		proposalProbs.push_back(1.0);
		useJumpKernel = true;
		evalType = NUMINT;
		F.function = &alphaStable;
		tuningBM = 1.0;							// does not require conditioning on brLen

	}
	else if (modelType == JUMP_NORM_NUMINT) // NumInt Poisson-Normal + BM
	{
		initialParms.push_back(new Sigma(randomPtr, "sigma-JN"));
		initialParms.push_back(new Sigma(randomPtr, "lambda-JN"));
		initialParms.push_back(new Sigma(randomPtr, "delta-JN"));
		proposalProbs.push_back(1.0);
		proposalProbs.push_back(1.0);
		proposalProbs.push_back(1.0);
		useJumpKernel = true;
		evalType = NUMINT;
		F.function = &poissonNormal;

	}
	else if (modelType == JUMP_NORM_PDF) // Pdf Sampled Poisson-Normal + BM
	{
		initialParms.push_back(new Sigma(randomPtr, "sigma-JN"));
		initialParms.push_back(new Sigma(randomPtr, "lambda-JN"));
		initialParms.push_back(new Sigma(randomPtr, "delta-JN"));
		proposalProbs.push_back(1.0);
		proposalProbs.push_back(1.0);
		proposalProbs.push_back(1.0);
		useJumpKernel = true;
		evalType = PDF;
		//tuningBM = 1.0;							// does not require conditioning on brLen
	}
	else if (modelType == VAR_GAMMA_BESSEL) // Bessel Variance Gamma
	{
		initialParms.push_back(new Sigma(randomPtr, "sigma-VG"));
		initialParms.push_back(new Kappa(randomPtr, "kappa-VG"));
		initialParms.push_back(new Sigma(randomPtr, "tau-VG"));
		proposalProbs.push_back(1.0);
		proposalProbs.push_back(1.0);
		proposalProbs.push_back(1.0);
		useJumpKernel = true;
		evalType = BESSEL;
		tuningBM = 1.0;							// does not require conditioning on brLen
	}
	else if (modelType == BM_ONLY) // Brownian motion only
	{
		initialParms.push_back(new Sigma(randomPtr, "sigma-BM"));
		proposalProbs.push_back(1.0);
		useJumpKernel = false;
	}
	else if (modelType == VAR_GAMMA_NUMINT) // NumInt Variance Gamma
	{
		initialParms.push_back(new Sigma(randomPtr, "sigma-VG"));
		initialParms.push_back(new Sigma(randomPtr, "kappa-VG"));
		initialParms.push_back(new Sigma(randomPtr, "tau-VG"));
		proposalProbs.push_back(1.0);
		proposalProbs.push_back(1.0);
		proposalProbs.push_back(1.0);
		useJumpKernel = true;
		evalType = NUMINT;
		F.function = &varianceGamma;
	}
	else if (modelType == BM_JUMP_PDF) // Pdf Brownian motion (w/ fake jump)
	{
		initialParms.push_back(new Sigma(randomPtr, "sigma-BM"));
		proposalProbs.push_back(1.0);
		useJumpKernel = true;
		evalType = PDF;
	}
	else if (modelType == DBL_EXP_NUMINT)
	{
		initialParms.push_back(new Sigma(randomPtr, "Sig-BM"));
		initialParms.push_back(new Sigma(randomPtr, "Lam-DE"));
		initialParms.push_back(new Sigma(randomPtr, "Sig-DE"));
		proposalProbs.push_back(1.0);
		proposalProbs.push_back(1.0);
		proposalProbs.push_back(1.0);
		useJumpKernel = true;
		evalType = NUMINT;
		F.function = &doubleExponential;
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

	tableList.front()->setParmVector(initialParms);
	//tableList.front()->setBranchVector(initialBranches);

	if (printStdOut) std::cout << "\n";
}

void Model::initializeTips(void)
{
	if (printStdOut) std::cout << "INITIALIZING: Node samples\n";


	for (int j = 0; j < numTaxa; j++)
	{
		Node* p = topologyPtr->getNode(j);
		//Node* n = topologyPtr->getDownPassNode(j); // wrong
		//double x = expressionPtr->getExpr(0,j,0); // old Patron-based Expression
		double x = expressionPtr->getData(j);
		std::cout << "\t" << p->getName() << "\t" << p->getIndex() << "\t" << x << "\n";
		p->setMu(x, 1);
		p->setSigma(0.0, 1);
		if (useJumpKernel)
		{
			double sumJumpSize = randomPtr->normalRv(0, 0.5);
			if (p->getV() == 0.0)
				sumJumpSize = 0.0;
			p->setSumJumpSize(sumJumpSize, 1);
		}
		p->setK(0.0, 1);
		p->setKb(0.0, 1);
		p->setKj(0.0, 1);
		p->copySpace(0, 1);
	}

	if (printStdOut) std::cout << "\n";
}

void Model::initializeGSL(void)
{

	std::cout << "INITIALIZING: GSL memory\n";
	workspaceSize = settingsPtr->getWorkspaceSize();
	trigTableSize = settingsPtr->getTrigTableSize();
	integralLength = settingsPtr->getIntegralLength();
	integralError = settingsPtr->getIntegralError();

	storeWorkspace = gsl_integration_workspace_alloc (workspaceSize);
	cycleWorkspace = gsl_integration_workspace_alloc (workspaceSize);
	trigTable = gsl_integration_qawo_table_alloc (0.0, integralLength, GSL_INTEG_COSINE, trigTableSize);
	gsl_integration_qawo_table_set(trigTable, 1.0, trigTableSize, GSL_INTEG_COSINE);

}

void Model::freeGSL(void)
{
	gsl_integration_workspace_free(storeWorkspace);
	gsl_integration_workspace_free(cycleWorkspace);
	gsl_integration_qawo_table_free(trigTable);
}

double Model::proposeJumpSize(Node* p, const std::vector<Parm*>& parmVector, int space)
{
	double lnProposalRatio = 0.0;

	// normal jump kernel
	if (modelType != BM_ONLY)
	{
		double jumpSize = p->getSumJumpSize(space);
		double jumpSizeNew = randomPtr->normalRv(jumpSize, sigmaJumpProposal);
#if DEBUG2
		std::cout << "\tjump at n" << p->getIndex() << ":\t" << jumpSize << " ->\t" << jumpSizeNew << "\n";
#endif
		if (p->getV() == 0.0)
			jumpSizeNew = 0.0;
		p->setSumJumpSize(jumpSizeNew, space);
		return lnProposalRatio;
	}

	// should never reach this
	if (modelType != BM_ONLY)
		std::cerr << "ERROR: proposeJumpSize() bypassed all conditions\n";

	return lnProposalRatio;
}

void Model::updateGSL(void)
{

}

double Model::jumpLnLikelihood(Node* p, const std::vector<Parm*>& parmVector, int space)
{


	Table* t = tableList.front();
	double k = 0.0;
	double newK = 0.0;

	// for time series, do nothing
	if (p->getV() == 0.0) {
		return 0.0;
	}

	if (modelType == ALPHA_STABLE_NUMINT || modelType == JUMP_NORM_NUMINT || modelType == VAR_GAMMA_NUMINT || modelType == DBL_EXP_NUMINT)
	{
		double x = p->getSumJumpSize(space);


#if DEBUG_QAWO
		// "exact" prob
		double sig_bm = (t->getParmVector()[0])->getValue();
		double lam_jn = (t->getParmVector()[1])->getValue();
		double sig_jn = (t->getParmVector()[2])->getValue();
        double sig_bm_2 = sig_bm * sig_bm;
        double sig_jn_2 = sig_jn * sig_jn;

		double v = p->getV();
		double w = v * (1.0 - tuningBM);
        double w_sqrt = pow(w,0.5);

		// REMOVE LATER
		//return randomPtr->lnNormalPdf(0.0, sig_bm, x);

		// sum over all possible numbers of jumps (truncated to n<100)
		int n = 0;
		do
		{
			if (lam_jn == 0.0)
			{
				double y = randomPtr->lnNormalPdf(0.0, pow(w, 0.5) * sig_bm, x);
				return y;
			}
			newK = randomPtr->poissonProb(lam_jn * v, n);
			newK *= randomPtr->normalPdf(0.0, pow(n * pow(sig_jn,2) + w * pow(sig_bm,2), 0.5), x);
			expected += newK;
			n++;
		}while(n < 100);
#endif

		//tuningBM = 0.9995;

		// pass GSL_function F with params
		std::vector<double> params;
		params.push_back((t->getParmVector()[0])->getValue()); // sig_bm
		params.push_back((t->getParmVector()[1])->getValue()); // lam_jn OR alf_as
		params.push_back((t->getParmVector()[2])->getValue()); // sig_jn OR cee_as
		params.push_back(p->getV());
		params.push_back(tuningBM);
		params.push_back(fabs(x));
		F.params = &params;

		// integrate
		double a = 0.0; // integrate from a to +Inf
		double result;
		double error;

#if DEBUG_QAWO
		 printf ("x %.8f\tsig_bm %.8f\tlam_jn %.8f\tsig_jn %.8f\tt %.8f\n", x, sig_bm, lam_jn, sig_jn, v);
#endif

		gsl_integration_qawf(&F, a, integralError, workspaceSize, storeWorkspace, cycleWorkspace, trigTable, &result, &error);
		// gsl_integration_qawo(&F, a, 0, 1e-7, workspaceSize, storeWorkspace, trigTable, &result, &error);
		k = log(result/(PI * fabs(x)));

#if DEBUG_MSG_QAWO
		printf ("\tNumInt\n");
		printf ("\t\tk               = % .18f\n", k);
		printf ("\t\tresult          = % .18f\n", result/(PI*fabs(x)));
		// printf ("exact result	  = % .18f\n", expected);
		printf ("\t\testimated error = % .18f\n", error);
		//printf ("actual error	  = % .18f\n", result/(PI*fabs(x)) - expected);
		printf ("\t\tintervals       =  %d\n\n", storeWorkspace->size);
		// printf ("sig_bm %.8f \tlam_jn %.8f\tsig_jn %.8f\tomega %.8f\tt %.8f\n", sig_bm, lam_jn, sig_jn, omega, v);
#endif

	}

	// X ~ BM + Poisson * Normal jumps
	else if (modelType == JUMP_NORM_PDF)
	{
		double sig_bm = (t->getParmVector()[0])->getValue();
		double lam_jn = (t->getParmVector()[1])->getValue();
		double sig_jn = (t->getParmVector()[2])->getValue();
        double sig_bm_2 = sig_bm * sig_bm;
        double sig_jn_2 = sig_jn * sig_jn;

		double v = p->getV();
		double w = v * (1.0 - tuningBM);
        double w_sqrt = pow(w,0.5);
		double x = p->getSumJumpSize(space);

        double th = 1e-20;
        double oldK = 0.0;
        bool pastMax = false;
        
		// sum over all possible numbers of jumps (truncated to n<100)
		int n = 0;
        //std::cout << p->getIndex() << "\n";
        //std::cout << v << " " << x << " " << sig_bm << " " << sig_jn << " " << lam_jn << "\n";
		do
		{
			if (lam_jn == 0.0)
			{
				std::cout << "Model::jumpLnLikelihood():\tlam_jn == 0.0 ... but how!?\n";
//				double y = randomPtr->lnNormalPdf(0.0, pow(w, 0.5) * sig_bm, x);
                double y = randomPtr->lnNormalPdf(0.0, w_sqrt * sig_bm, x);
#if DEBUG2
				std::cout << "\t\tn" << p->getIndex() << "\tv:\t" << v << "\tx:\t" << x << "\tlnL:\t" << y <<  "\n";
#endif
				return y;
			}
            oldK = newK;
			newK = randomPtr->poissonProb(lam_jn * v, n);
			newK *= randomPtr->normalPdf(0.0, pow(n * sig_jn_2 + w * sig_bm_2, 0.5), x);
            // newK *= randomPtr->normalPdf(0.0, pow(n * pow(sig_jn,2) + w * pow(sig_bm,2), 0.5), x);
			//newK *= randomPtr->normalPdf(0.0, pow(n * pow(sig_jn,2),0.5),x);// + w * pow(sig_bm,2), 0.5), x);
			k += newK;
			n++;
            //std::cout << "\t" << k << "\n";
		}while(n < 1000 && ((newK - oldK) >= 0
                        || -(newK - oldK) >= th)
        );
        //std::cout << "n " << n << "\n";
        //std::cout << k << " ";
		k = log(k);
        //std::cout << k << "\n\n";
#if DEBUG2
			std::cout << "\t\tn:" << p->getIndex() <<  "\tv:\t" << v << "\tw:\t" << w << "\tx:\t" << x << "\tlnL:\t" << k <<  "\n";
			std::cout << "\t\tn*sig_jn^2:\t" << n*pow(sig_jn,2) << "\tw*sig_bm^2:\t" << w * pow(sig_bm,2) << "\n";
			std::cout << "\t\t\tsig_bm:\t" << sig_bm << "\tlam_jn:\t" << lam_jn << "\tsig_jn:\t" << sig_jn << "\n";
			//std::cout << "\n";
#endif
	}

	else if (modelType == BM_JUMP_PDF)
	{
		double sig_bm = (t->getParmVector()[0])->getValue();
		double v = p->getV();
		double w = v * (1.0 - tuningBM);
		double x = p->getSumJumpSize(space);

		k = randomPtr->lnNormalPdf(0.0, pow(w,0.5) * sig_bm, x);
	}
	if (modelType == VAR_GAMMA_BESSEL)
	{
		k = besselLnLikelihood(p, parmVector, space);
	}

	return k;
}

double Model::besselLnLikelihood(Node* p, const std::vector<Parm*>& parmVector, int space)
{

	double k = 0.0;

	if (modelType == VAR_GAMMA_BESSEL)
	{
		double v = p->getV();
		double sig_bm = parmVector[0]->getValue();
		double kap_vg = parmVector[1]->getValue();
		double sig_vg = parmVector[2]->getValue();
		double xx = fabs(p->getSumJumpSize(space));

		double x = xx;
		//double x = xx - kap_vg * v * log(1 - 0.5 * sig_vg * sig_vg / kap_vg);

		/*
		double h_bessel_arg1 = fabs(v * kap_vg - 0.5);
		double h_bessel_arg2 = pow(2 * kap_vg, 0.5) * x / sig_vg;
		double h_bessel = gsl_sf_bessel_lnKnu(h_bessel_arg1, h_bessel_arg2);

		double h_top = (0.75 - 0.5 * v * kap_vg) * log(2) + (0.25 + 0.5 * v * kap_vg) * log(kap_vg) + (0.5 - v * kap_vg) * log(sig_vg / x);
		double h_bottom = 0.5 * log(PI * sig_vg * sig_vg) + gsl_sf_lngamma(v * kap_vg);
		k = h_top + h_bessel - h_bottom;
		*/

		double h_bessel_arg1 = fabs(v / kap_vg - 0.5);
		double h_bessel_arg2 = pow(2 / kap_vg, 0.5) * x / sig_vg;
		double h_bessel = gsl_sf_bessel_lnKnu(h_bessel_arg1, h_bessel_arg2);

		double h_top = (0.75 - 0.5 * v / kap_vg) * log(2) + (0.25 + 0.5 * v / kap_vg) * -log(kap_vg) + (0.5 - v / kap_vg) * log(sig_vg / x);
		double h_bottom = 0.5 * log(PI * sig_vg * sig_vg) + gsl_sf_lngamma(v / kap_vg);
		k = h_top + h_bessel - h_bottom;

#if DEBUG_BESSEL
		std::cout << "\tBessel\n";
		std::cout << "\t\tv\t" << v << "\tkap_vg\t" << kap_vg << "\tsig_vg\t" << sig_vg << "\tsig_bm\t" << sig_bm << "\tx\t" << x << "\n";
		std::cout << "\t\thba1\t" << h_bessel_arg1 << "\thba2\t" << h_bessel_arg2 <<"\n";
		std::cout << "\t\tk\t" << k <<  " = " << h_top << " + " << h_bessel << " - " << h_bottom << "\n";
#endif
	}

	else if (modelType == NORM_INV_GAUSS_BESSEL)
	{
		double v = p->getV();
		double sig_bm = parmVector[0]->getValue();
		double kap_vg = parmVector[1]->getValue();
		double sig_vg = parmVector[2]->getValue();
		double xx = fabs(p->getSumJumpSize(space));

		double x = xx;
		//double x = xx - kap_vg * v * log(1 - 0.5 * sig_vg * sig_vg / kap_vg);

		/*
		double h_bessel_arg1 = fabs(v * kap_vg - 0.5);
		double h_bessel_arg2 = pow(2 * kap_vg, 0.5) * x / sig_vg;
		double h_bessel = gsl_sf_bessel_lnKnu(h_bessel_arg1, h_bessel_arg2);

		double h_top = (0.75 - 0.5 * v * kap_vg) * log(2) + (0.25 + 0.5 * v * kap_vg) * log(kap_vg) + (0.5 - v * kap_vg) * log(sig_vg / x);
		double h_bottom = 0.5 * log(PI * sig_vg * sig_vg) + gsl_sf_lngamma(v * kap_vg);
		k = h_top + h_bessel - h_bottom;
		*/

		double h_bessel_arg1 = fabs(v / kap_vg - 0.5);
		double h_bessel_arg2 = pow(2 / kap_vg, 0.5) * x / sig_vg;
		double h_bessel = gsl_sf_bessel_lnKnu(h_bessel_arg1, h_bessel_arg2);

		double h_top = (0.75 - 0.5 * v / kap_vg) * log(2) + (0.25 + 0.5 * v / kap_vg) * -log(kap_vg) + (0.5 - v / kap_vg) * log(sig_vg / x);
		double h_bottom = 0.5 * log(PI * sig_vg * sig_vg) + gsl_sf_lngamma(v / kap_vg);
		k = h_top + h_bessel - h_bottom;

#if DEBUG_BESSEL
		std::cout << "\tBessel\n";
		std::cout << "\t\tv\t" << v << "\tkap_vg\t" << kap_vg << "\tsig_vg\t" << sig_vg << "\tsig_bm\t" << sig_bm << "\tx\t" << x << "\n";
		std::cout << "\t\thba1\t" << h_bessel_arg1 << "\thba2\t" << h_bessel_arg2 <<"\n";
		std::cout << "\t\tk\t" << k <<  " = " << h_top << " + " << h_bessel << " - " << h_bottom << "\n";
#endif
	}

	return k;
}

double Model::modelLnLikelihood(int space)
{
	double lnL = 0.0;

	lnL = tableLnLikelihood(tableList.front(), space);

	return lnL;
}

double Model::tableLnLikelihood(Table* t, int space)
{
	double lnL = 0.0;
	const std::vector<Parm*> parmVector = t->getParmVector();


	// update lnL values at tips according to sampled values
	if (useJumpKernel)
	{
		for (int j = 0; j < numTaxa; j++)
		{
			Node* p = topologyPtr->getNode(j);

			double lnJumpLike = jumpLnLikelihood(p, parmVector, space);

			p->setK(lnJumpLike, space);
			p->setKj(lnJumpLike, space);
		}
	}

	//  drift is accounted for in besselLnLikelihood()

	// Pruning-wise, calculate the likelihood for each branch in the tree
	// NOTE: driftLnLikelihood() calls jumpLnLikelihood()
	for (int j = 0; j < numNodes; j++)
	{
		Node* p = topologyPtr->getDownPassNode(j);

		if (p->getLft() != NULL && p->getRht() != NULL)
		{
			driftLnLikelihood(p, parmVector, space);
		}
	}

	lnL = topologyPtr->getRoot()->getK(space);

	return lnL;
}


double Model::modelDriftOnlyLnLikelihood(int space)
{
	double lnL = 0.0;

	lnL = tableDriftOnlyLnLikelihood(tableList.front(), space);

	return lnL;
}

double Model::tableDriftOnlyLnLikelihood(Table* t, int space)
{
	double lnL = 0.0;
	const std::vector<Parm*> parmVector = t->getParmVector();

	/*
	// update lnL values at tips according to sampled values
	if (useJumpKernel)
	{
		for (int j = 0; j < numTaxa; j++)
		{
			Node* p = topologyPtr->getNode(j);

			//double lnJumpLike = jumpLnLikelihood(p, parmVector, space);

			p->setK(lnJumpLike, space);
			p->setKj(lnJumpLike, space);
		}
	}
	*/

	// Pruning-wise, calculate the likelihood for each branch in the tree
	// NOTE: driftLnLikelihood() calls jumpLnLikelihood()
	for (int j = 0; j < numNodes; j++)
	{
		Node* p = topologyPtr->getDownPassNode(j);

		if (p->getLft() != NULL && p->getRht() != NULL)
		{
			driftOnlyLnLikelihood(p, parmVector, space);
		}
	}

	lnL = topologyPtr->getRoot()->getKb(space);
	return lnL;
}


double Model::driftOnlyLnLikelihood(Node* p, const std::vector<Parm*>& parmVector, int space)
{
	double tauL = p->getLft()->getV();
	double tauR = p->getRht()->getV();
	double sigmaBM = parmVector[0]->getValue();
	double sigmaBM2 = sigmaBM * sigmaBM;
	double lnL = 0.0;



	#if DEBUG
					std::cout << "p->getIndex() = " << p->getIndex() << "\tp->getLft()->getIndex() = " << p->getLft()->getIndex() << "\tp->getRht()->getIndex() = " << p->getRht()->getIndex() << "\n";
	#endif

	// get descendants' conditional parameters
	double muL = p->getLft()->getMu(space);
	double sigmaL = p->getLft()->getSigma(space);
	double kL = p->getLft()->getK(space);

	double muR = p->getRht()->getMu(space);
	double sigmaR = p->getRht()->getSigma(space);
	double kR = p->getRht()->getK(space);


	// rescale branches according to descendants' variance
	double tL = tauL + (sigmaL * sigmaL) / sigmaBM2;
	double tR = tauR + (sigmaR * sigmaR) / sigmaBM2;

	// rescale branches and displace node values if using a jump kernel
	if (useJumpKernel)
	{
		muL -= p->getLft()->getSumJumpSize(space);
		muR -= p->getRht()->getSumJumpSize(space);
		tL = tL * tuningBM;
		tR = tR * tuningBM;
	}

	// update present node's parameters
	double muP = (muL * tR + muR * tL) / (tL + tR);
	double sigmaP = pow(sigmaBM2 * tL * tR / (tL + tR), 0.5);
	p->setMu(muP, space);
	p->setSigma(sigmaP, space);

	// calculate likelihood of sampled jumps
	double lnJumpLike = 0.0;

	/*
	// ignore the root since it does not sample any jumps
	if (useJumpKernel && p->getAnc() != NULL)
	{
		lnJumpLike = jumpLnLikelihood(p, parmVector, space);
	}
	*/

	// update normalizing constant
	double kP = -1.0 * pow(muL - muR, 2) / (2 * sigmaBM2 * (tL + tR));
	kP -= log(pow(2 * PI * sigmaBM2 * (tL + tR), 0.5));

	p->setKb(kP + p->getRht()->getKb(space) + p->getLft()->getKb(space), space);
	//p->setKj(lnJumpLike + p->getRht()->getKj(space) + p->getLft()->getKj(space), space);

	//kP += kL + kR;
	lnL = kP + p->getRht()->getKb(space) + p->getLft()->getKb(space);// + lnJumpLike;
	//p->setK(kP, space);
#if DEBUG3
		std::cout << std::setprecision(4);
		std::cout << "\t\t*\tnL:\t" << p->getLft()->getIndex() << "\tmuL:\t" << muL << "\tsigL:\t" << sigmaL << "\tkL:\t" << kL << "\n";
		std::cout << "\t\t\tnR:\t" << p->getRht()->getIndex() << "\tmuR:\t" << muR << "\tsigR:\t" << sigmaR << "\tkR:\t" << kR << "\n";
		std::cout << "\t\t\tnP:\t" << p->getIndex() << "\tmuP:\t" << muP << "\tsigP:\t" << sigmaP << "\tkP:\t" << kP << "\n";
#endif

	//}
	return lnL;
}



double Model::driftLnLikelihood(Node* p, const std::vector<Parm*>& parmVector, int space)
{
	double tauL = p->getLft()->getV();
	double tauR = p->getRht()->getV();
	double sigmaBM = parmVector[0]->getValue();
	double sigmaBM2 = sigmaBM * sigmaBM;
	double lnL = 0.0;

	#if DEBUG
					std::cout << "p->getIndex() = " << p->getIndex() << "\tp->getLft()->getIndex() = " << p->getLft()->getIndex() << "\tp->getRht()->getIndex() = " << p->getRht()->getIndex() << "\n";
	#endif

	// get descendants' conditional parameters
	double muL = p->getLft()->getMu(space);
	double sigmaL = p->getLft()->getSigma(space);
	double kL = p->getLft()->getK(space);

	double muR = p->getRht()->getMu(space);
	double sigmaR = p->getRht()->getSigma(space);
	double kR = p->getRht()->getK(space);

	// rescale branches according to descendants' variance
	double tL = tauL + (sigmaL * sigmaL) / sigmaBM2;
	double tR = tauR + (sigmaR * sigmaR) / sigmaBM2;

	// rescale branches and displace node values if using a jump kernel
	if (useJumpKernel)
	{
		muL -= p->getLft()->getSumJumpSize(space);
		muR -= p->getRht()->getSumJumpSize(space);
		tL = tL * tuningBM;
		tR = tR * tuningBM;
	}

	// update present node's parameters
	double muP = (muL * tR + muR * tL) / (tL + tR);
	double sigmaP = pow(sigmaBM2 * tL * tR / (tL + tR), 0.5);
	p->setMu(muP, space);
	p->setSigma(sigmaP, space);

	// calculate likelihood of sampled jumps
	double lnJumpLike = 0.0;

	// ignore the root since it does not sample any jumps
	if (useJumpKernel && p->getAnc() != NULL)
	{
		lnJumpLike = jumpLnLikelihood(p, parmVector, space);
	}

	// update normalizing constant
	double kP = -1.0 * pow(muL - muR, 2) / (2 * sigmaBM2 * (tL + tR));
	kP -= log(pow(2 * PI * sigmaBM2 * (tL + tR), 0.5));

	p->setKb(kP + p->getRht()->getKb(space) + p->getLft()->getKb(space), space);
	p->setKj(lnJumpLike + p->getRht()->getKj(space) + p->getLft()->getKj(space), space);

	kP += kL + kR;
	lnL = kP + lnJumpLike;
	p->setK(kP + lnJumpLike, space);
#if DEBUG3
		std::cout << std::setprecision(4);
		std::cout << "\t\t*\tnL:\t" << p->getLft()->getIndex() << "\tmuL:\t" << muL << "\tsigL:\t" << sigmaL << "\tkL:\t" << kL << "\n";
		std::cout << "\t\t\tnR:\t" << p->getRht()->getIndex() << "\tmuR:\t" << muR << "\tsigR:\t" << sigmaR << "\tkR:\t" << kR << "\n";
		std::cout << "\t\t\tnP:\t" << p->getIndex() << "\tmuP:\t" << muP << "\tsigP:\t" << sigmaP << "\tkP:\t" << kP << "\n";
#endif

//	}
	return lnL;
}

double Model::getProcessKurtosis(void)
{
    Table* t = tableList.front();
    std::vector<Parm*> p = t->getParmVector();
    double ret = 0.0;

    if (modelType == NORM_INV_GAUSS_BESSEL)
    {
        ;
    }
    else if (modelType == ALPHA_STABLE_NUMINT)
    {
        ret = -1.0;
    }
    else if (modelType == JUMP_NORM_NUMINT)
    {
        double p0 = p[0]->getValue();
        double p1 = p[1]->getValue();
        double p2 = p[2]->getValue();
        ret = (3 * p1 * pow(p2,4)) / (p0*p0 + p1*p2*p2);
    }
    else if (modelType == JUMP_NORM_PDF)
    {
        double p0 = p[0]->getValue();
        double p1 = p[1]->getValue();
        double p2 = p[2]->getValue();
        ret = (3 * p1 * pow(p2,4)) / (p0*p0 + p1*p2*p2);
    }
    else if (modelType == VAR_GAMMA_BESSEL)
    {
        double p0 = p[0]->getValue();
        double p1 = p[1]->getValue();
        double p2 = p[2]->getValue();
        ret = (3 * p1 * pow(p2,4)) / (p0*p0 + p2*p2);

    }
    else if (modelType == BM_ONLY)
    {
        ;
    }
    else if (modelType == VAR_GAMMA_NUMINT)
    {
        double p0 = p[0]->getValue();
        double p1 = p[1]->getValue();
        double p2 = p[2]->getValue();
        ret = (3 * p1 * pow(p2,4)) / (p0*p0 + p2*p2);
    }
    else if (modelType == BM_JUMP_PDF)
    {
        ;
    }
    else if (modelType == DBL_EXP_NUMINT)
    {
        ;
    }
    else if (modelType == SKEW_NORMAL_NUMINT)
    {
        ;
    }
    return ret;
}

double Model::getProcessVariance(void)
{
    Table* t = tableList.front();
    std::vector<Parm*> p = t->getParmVector();
    double ret = 0.0;
    
    if (modelType == NORM_INV_GAUSS_BESSEL)
    {
        ;
    }
    else if (modelType == ALPHA_STABLE_NUMINT)
    {
        ret = -1.0;
    }
    else if (modelType == JUMP_NORM_NUMINT)
    {
        double p0 = p[0]->getValue();
        double p1 = p[1]->getValue();
        double p2 = p[2]->getValue();
        ret = p0*p0 + p1*p2*p2;
    }
    else if (modelType == JUMP_NORM_PDF)
    {
        double p0 = p[0]->getValue();
        double p1 = p[1]->getValue();
        double p2 = p[2]->getValue();
        ret = p0*p0 + p1*p2*p2;
    }
    else if (modelType == VAR_GAMMA_BESSEL)
    {
        double p0 = p[0]->getValue();
        double p1 = p[1]->getValue();
        double p2 = p[2]->getValue();
        ret = p0*p0 + p2*p2;
    }
    else if (modelType == BM_ONLY)
    {
        double p0 = p[0]->getValue();
        ret = p0*p0;
    }
    else if (modelType == VAR_GAMMA_NUMINT)
    {
        double p0 = p[0]->getValue();
        double p1 = p[1]->getValue();
        double p2 = p[2]->getValue();
        ret = p0*p0 + p2*p2;
    }
    else if (modelType == BM_JUMP_PDF)
    {
        double p0 = p[0]->getValue();
        ret = p0*p0;
    }
    else if (modelType == DBL_EXP_NUMINT)
    {
        ;
    }
    else if (modelType == SKEW_NORMAL_NUMINT)
    {
        ;
    }
    return ret;
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

Node* Model::pickNodeByBrLen(void)
{
	double u = randomPtr->uniformRv();
	double sum = 0.0;
	int n = 0;
	Node* p = NULL;

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

	// weight probability of selecting a table by the table. #patrons ^ -0.5
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


void Model::keep(Parm* p)
{
	p->keep();
}

void Model::restore(Parm* p)
{
	p->restore();
}


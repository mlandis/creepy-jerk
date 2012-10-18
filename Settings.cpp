
/*
 * Settings.cpp
 *
 *  Created on: Dec 16, 2010
 *
 */

#include "Settings.h"

#define NORM_INV_GAUSS 0
#define ALPHA_STABLE_NUMINT 1
#define JUMP_NORM_NUMINT 2
#define JUMP_NORM_PDF 3
#define VAR_GAMMA_BESSEL 4
#define BM_ONLY 5
#define VAR_GAMMA_NUMINT 6
#define BM_JUMP_PDF 7
#define DBL_EXP_NUMINT 8
#define SKEW_NORMAL_NUMINT 9


Settings::Settings(int argc, char** argv)
{

	// Default initialization for Settings

	// Model settings
	seed = 52;
	modelType = 3;
	fixBranches = true;
	tuningBM = 0.9995;
	useJumpKernel = true;
	if (modelType == BM_ONLY)
		useJumpKernel = false;
	sigmaJumpProposal = 1.0;

	// Stepping Stone settings
	useSteppingStone = false;
	betaSteppingStone = 1.000000;
	printStdOut = true;

	// FFT settings -- currently unused
	useFFT = false;
	tipStdDev = 0.001;
	numSteps = pow(2, 16);	// 2^19 is the upper limit
	startStep = -10000.0;	//-300.0;
	finalStep = 10000.0;	//300.0;

	// MCMC settings
	numCycles = pow(10,1)*2;
	printFreqMH = 1;
	printFreqJump = printFreqMH;
	printFreqStdOut = printFreqMH;
	printFreqCRP = printFreqMH;
	printStdOut = true;

	// CRP settings
	alphaCRP = 0.75;
	auxCRP = 3.0;
	aCRP = 1.0;
	bCRP = 1.0;
	useCRP = false;

	// GSL settings (NOTE: different models are numerically unstable for different GSL settings)
	workspaceSize = 30;
	integralLength = 5000.0;
	trigTableSize = 50;
	integralError = 1e-6;

	// File names
	inputDirPath = "/Users/mlandis/data/expr_phylo/input/";
	outputDirPath = "/Users/mlandis/data/expr_phylo/output/";
	outputFileName = "";


	simName = "cj." + Util::getTime();
	exprFileName = simName + ".data.txt";
	treeFileName = simName + ".tree.txt";
	taxaFileName = simName + ".taxa.txt";


	simName = "gogo2";
	taxaFileName = "primates.eastman.isler_pruned.taxa.txt";
	treeFileName = "primates.eastman.isler_pruned.tree.txt";
	exprFileName = "primates.eastman.isler_pruned.mass.data.txt";


	// Assign settings arguments
	setArguments(argc, argv);

	// Assign settings string
	std::stringstream ss;
	ss << ".model" << modelType;
	ss << ".seed" << seed;
	if (useSteppingStone)
	{
		ss << ".betaSS" << betaSteppingStone;
	}

	if (outputFileName == "")
	{
		outputFileName = simName + ss.str() + "." + exprFileName;
	}

	//if (printStdOut)
	print();
}

Settings::Settings(void)
{
	// do nothing
}

void Settings::setArguments(int argc, char** argv)
{
	std::map<std::string, std::string> argMap;
	std::map<std::string, std::string>::iterator argMapIt;
	std::pair<std::map<std::string, std::string>::iterator, bool> ret;

	std::vector<std::string> argVector;

	for (int i = 0; i < argc; i++)
	{
		std::stringstream ss(argv[i]);
		argVector = Util::split(ss.str(), '=');
		if (argVector.size() == 2)
		{
			//std::cout << argVector[0] << "\t" << argVector[1] << "\n";
			ret = argMap.insert(std::pair<std::string, std::string>(argVector[0], argVector[1]));
		}
	}


	std::string argName = "";
	std::string argVal = "";

	for (argMapIt = argMap.begin(); argMapIt != argMap.end(); argMapIt++)
	{

		argName = (*argMapIt).first;
		argVal = (*argMapIt).second;

		//if (printStdOut) std::cout << "Setting " << argName << " = " << argVal << "\n";

		if (argName == "" || argVal == "") ;

		// INPUT SETTINGS
		else if (argName == "-inputDirPath")
			inputDirPath = argVal;
		else if (argName == "-settingsFilePath")
			settingsFilePath = argVal;
		else if (argName == "-exprFileName")
			exprFileName = argVal;
		else if (argName == "-treeFileName")
			treeFileName = argVal;
		else if (argName == "-taxaFileName")
			taxaFileName = argVal;


		// OUTPUT SETTINGS
		else if (argName == "-outputDirPath")
			outputDirPath = argVal;
		else if (argName == "-outputFileName")
			outputFileName = argVal;
		else if (argName == "-simName")
			simName = argVal;
		else if (argName == "-printStdOut")
			printStdOut = Util::stringToBool(argVal);
		else if (argName == "-printFreqStdOut")
			printFreqStdOut = Util::stringToInt(argVal);


		// MODEL SETTINGS
		else if (argName == "-seed")
			seed = Util::stringToInt(argVal);
		else if (argName == "-modelType")
		{
			modelType = Util::stringToInt(argVal);
			if (modelType == BM_ONLY)
				useJumpKernel = false;
			else
				useJumpKernel = true;
			std::cout << "modelType = " << modelType << "\n";
		}
		else if (argName == "-fixBranches")
			fixBranches = Util::stringToBool(argVal);
		else if (argName == "-tuningBM")
			tuningBM = Util::stringToDouble(argVal);
		else if (argName == "-useJumpKernel")
			useJumpKernel = Util::stringToBool(argVal);
		else if (argName == "-sigmaJumpProposal")
			sigmaJumpProposal = Util::stringToDouble(argVal);


		// MODEL TESTING SETTINGS
		else if (argName == "-useSteppingStone")
			useSteppingStone = Util::stringToBool(argVal);
		else if (argName == "-betaSteppingStone")
			betaSteppingStone = Util::stringToDouble(argVal);

		// MCMC SETTINGS
		else if (argName == "-numCycles")
			numCycles = Util::stringToInt(argVal);
		else if (argName == "-printFreqMH")
			printFreqMH = Util::stringToInt(argVal);
		else if (argName == "-printFreqJump")
			printFreqJump = Util::stringToInt(argVal);

	}
}


void Settings::print(void)
{
	std::cout << "Settings\n";

	std::cout << "\tInput\n";
	std::cout << "\t\tinputDirPath      = " << inputDirPath << "\n";
	std::cout << "\t\texprFileName      = " << exprFileName << "\n";
	std::cout << "\t\ttreeFileName      = " << treeFileName << "\n";
	std::cout << "\t\ttaxaFileName      = " << taxaFileName << "\n";
	std::cout << "\t\tsettingsFilePath  = " << settingsFilePath << "\n";

	std::cout << "\tOutput\n";
	std::cout << "\t\tsimName           = " << simName << "\n";
	std::cout << "\t\toutputDirPath     = " << outputDirPath << "\n";
	std::cout << "\t\toutputFileName    = " << outputFileName << "\n";
	std::cout << "\t\tprintStdOut       = " << Util::boolToString(printStdOut) << "\n";
	std::cout << "\t\tprintFreqStdOut   = " << Util::intToString(printFreqStdOut) << "\n";

	std::cout << "\tModel\n";
	std::cout << "\t\tseed              = " << seed << "\n";
	std::cout << "\t\tmodelType         = " << modelType << "\n";
	std::cout << "\t\ttuningBM          = " << tuningBM << "\n";
	std::cout << "\t\tuseJumpKernel     = " << Util::boolToString(useJumpKernel) << "\n";
	std::cout << "\t\tsigmaJumpProposal = " << sigmaJumpProposal << "\n";

	std::cout << "\tMCMC\n";
	std::cout << "\t\tnumCycles         = " << numCycles << "\n";
	std::cout << "\t\tprintFreqMH       = " << printFreqMH << "\n";
	std::cout << "\t\tprintFreqJump     = " << printFreqJump << "\n";
	std::cout << "\t\tprintFreqStdOut   = " << printFreqStdOut << "\n";

	std::cout << "\tModel Testing\n";
	std::cout << "\t\tuseSteppingStone  = " << Util::boolToString(useSteppingStone) << "\n";
	std::cout << "\t\tbetaSteppingStone = " << betaSteppingStone << "\n";

	std::cout << "\n";
}


/*
 * Settings.cpp
 *
 *  Created on: Dec 16, 2010
 *
 */

#include "Settings.h"

Settings::Settings(int argc, char** argv)
{

	// Default initialization for Settings

	// Model settings
	seed = 2;
	modelType = 3;
	fixBranches = true;
	tuningBM = 0.9995;
	useJumpKernel = false;

	// Stepping Stone settings
	useSteppingStone = false;
	betaSteppingStone = 1.0;
	numSteppingStone = 0;
	printStdOut = true;

	// FFT settings
	useFFT = false;
	tipStdDev = 0.001;
	numSteps = pow(2, 16);	// 2^19 is the upper limit
	startStep = -10000.0;	//-300.0;
	finalStep = 10000.0;	//300.0;

	// MCMC settings
	numCycles = 100;//00000;
	printFreqCRP = 0;
	printFreqMH = 1;//000;
	printFreqJump = numCycles+1;

	// CRP settings
	alphaCRP = 0.75;
	auxCRP = 3.0;
	aCRP = 1.0;
	bCRP = 1.0;
	useCRP = false;

	// File names
	inputDirPath = "/Users/mlandis/data/expr_phylo/input/";
	outputDirPath = "/Users/mlandis/data/expr_phylo/output/";
	outputFileName = "";

	simName = "";
	exprFileName = simName + ".data.txt";
	treeFileName = simName + ".tree.txt";
	taxaFileName = simName + ".taxa.txt";
//	exprFileName = "fake_data_9_9_11/tip_50_Sig_BM_0.9_Lam.JN_0.06_Sig.JN_1.06.data.txt";
//	exprFileName = "fake_data_9_9_11/tip_50_Sig_BM_1.5_Lam.JN_0.34_Sig.JN_1.8.data.txt";
//	exprFileName = "fake_data_9_9_11/tip_50_Sig_BM_0.97_Lam.JN_1.14_Sig.JN_0.37.data.txt";
	exprFileName = "disentangle.data.txt";
	treeFileName = "Sig.BM_.5_Lam.JN_0_Sig.JN_0_Sig.Tip_.001.tree.txt";
	//exprFileName = "2tip_Sig.BM_2.0_Lam.JN_0.0_Sig.JN_0.0_Sig.Tip_0.0.data.txt";
	//treeFileName = "2tip_Sig.BM_2.0_Lam.JN_0.0_Sig.JN_0.0_Sig.Tip_0.0.tree.txt";
	//taxaFileName = "2tip_Sig.BM_2.0_Lam.JN_0.0_Sig.JN_0.0_Sig.Tip_0.0.taxa.txt";
//	exprFileName = "Sig.BM_2.0_Lam.JN_0_Sig.JN_0_Sig.Tip_0.0_4.data.txt";
//	exprFileName = "Sig.BM_.5_Lam.JN_1.2_Sig.JN_.9_Sig.Tip_0.0_4.data.txt";
//	treeFileName = "Sig.BM_.5_Lam.JN_0_Sig.JN_0_Sig.Tip_.001.tree.txt";
//	taxaFileName = "Sig.BM_0.5_Lam.JN_0.0_Sig.JN_0.0_Sig.Tip_0.0.taxa.txt";
	taxaFileName = "data1.taxa.txt";

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

		// MODEL SETTINGS
		else if (argName == "-seed")
			seed = Util::stringToInt(argVal);
		else if (argName == "-modelType")
			modelType = Util::stringToInt(argVal);
		else if (argName == "-fixBranches")
			fixBranches = Util::stringToBool(argVal);
		else if (argName == "-tuningBM")
			tuningBM = Util::stringToDouble(argVal);
		else if (argName == "-useJumpKernel")
			useJumpKernel = Util::stringToBool(argVal);



		// MODEL TESTING SETTINGS
		else if (argName == "-useSteppingStone")
			useSteppingStone = Util::stringToBool(argVal);
		else if (argName == "-numSteppingStone")
			numSteppingStone = Util::stringToInt(argVal);
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

	std::cout << "\tModel\n";
	std::cout << "\t\tseed              = " << seed << "\n";
	std::cout << "\t\tmodelType         = " << modelType << "\n";
	std::cout << "\t\ttuningBM          = " << tuningBM << "\n";
	std::cout << "\t\tuseJumpKernel     = " << Util::boolToString(useJumpKernel) << "\n";

	std::cout << "\tMCMC\n";
	std::cout << "\t\tnumCycles         = " << numCycles << "\n";
	std::cout << "\t\tprintFreqMH       = " << printFreqMH << "\n";
	std::cout << "\t\tprintFreqJump     = " << printFreqJump << "\n";

	std::cout << "\tModel Testing\n";
	std::cout << "\t\tuseSteppingStone  = " << Util::boolToString(useSteppingStone) << "\n";
	std::cout << "\t\tnumSteppingStone  = " << numSteppingStone << "\n";
	std::cout << "\t\tbetaSteppingStone = " << betaSteppingStone << "\n";

	std::cout << "\n";
}

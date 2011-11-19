
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

	// FFT settings
	useFFT = false;
	tipStdDev = 0.001;
	numSteps = pow(2, 16);	// 2^19 is the upper limit
	startStep = -10000.0;	//-300.0;
	finalStep = 10000.0;	//300.0;

	// MCMC settings
	numCycles = 2000000;//00000;
	printFreqCRP = 100;
	printFreqMH = 1;//000;

	// CRP settings
	alphaCRP = 0.75;
	auxCRP = 3.0;
	aCRP = 1.0;
	bCRP = 1.0;
	useCRP = false;

	// Assign settings string
	std::stringstream ss;
	//if (fixBranches) ss << ".fixbr"; else ss << ".estbr";
	//if (useCRP) ss << ".useCRP"; else ss << ".noCRP";
	if (useFFT)
	{
		ss << ".useFFT";
		ss << ".s" << numSteps;
	}
	else
	{
		ss << ".noFFT";
	}

	ss << ".model" << modelType;
	ss << ".seed" << seed;


	// File names
	inputDirPath = "/Users/mlandis/data/expr_phylo/input/";
	outputDirPath = "/Users/mlandis/data/expr_phylo/output/";
	simName = "tip_50";
	std::string outName = ".A";
	exprFileName = simName + ".data.txt";
	treeFileName = simName + ".tree.txt";
	taxaFileName = simName + ".taxa.txt";
//	exprFileName = "fake_data_9_9_11/tip_50_Sig_BM_0.9_Lam.JN_0.06_Sig.JN_1.06.data.txt";
	exprFileName = "fake_data_9_9_11/tip_50_Sig_BM_0.97_Lam.JN_1.14_Sig.JN_0.37.data.txt";
//	exprFileName = "fake_data_9_9_11/tip_50_Sig_BM_1.5_Lam.JN_0.34_Sig.JN_1.8.data.txt";
	treeFileName = "fake_data_9_9_11/tip_50.tree.txt";
	taxaFileName = "data1.taxa.txt";
//	taxaFileName = "data1.tip_10.taxa.txt";
	outputFileName = simName + outName + ss.str();

}

Settings::Settings(void)
{
	// do nothing
}

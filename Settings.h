/*
 * Settings.h
 *
 *  Created on: Dec 16, 2010
 *
 *	This class handles all user-defined settings for the program. Specifically:
 *
 *		Filepaths (input, output, tree, logs, etc.)
 *		Model settings (parameters, etc)
 *		MCMC settings (number of cycles, etc)
 *
 */

#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>

#include "FileMgr.h"
#include "Util.h"

#ifndef SETTINGS_H_
#define SETTINGS_H_

class Settings
{
public:
	Settings(void);
	Settings(int argc, char** argv);
	std::string getInputDirPath(void)			{ return inputDirPath; }
	std::string getOutputDirPath(void)			{ return outputDirPath; }
	std::string getSettingsFilePath(void)		{ return settingsFilePath; }
	std::string getExprFileName(void)			{ return exprFileName; }
	std::string getTreeFileName(void)			{ return treeFileName; }
	std::string getTaxaFileName(void)			{ return taxaFileName; }
	std::string getOutputFileName(void)			{ return outputFileName; }
	std::string getSimName(void)				{ return simName; }

	int			getSeed(void)					{ return seed; }
	bool		getFixBranches(void)			{ return fixBranches; }
	double		getTuningBM(void)				{ return tuningBM; }
	bool		getUseJumpKernel(void)			{ return useJumpKernel; }

	bool		getUseSteppingStone(void)		{ return useSteppingStone; }
	int			getNumSteppingStone(void)		{ return numSteppingStone; }
	double		getBetaSteppingStone(void)		{ return betaSteppingStone; }

	double		getAlphaCRP(void)				{ return alphaCRP; }
	double		getAuxCRP(void)					{ return auxCRP; }
	double		getACRP(void)					{ return aCRP; }
	double		getBCRP(void)					{ return bCRP; }
	bool		getUseCRP(void)					{ return useCRP; }

	int			getNumCycles(void)				{ return numCycles; }
	int			getPrintFreqMH(void)			{ return printFreqMH; }
	int			getPrintFreqCRP(void)			{ return printFreqCRP; }
	int			getPrintFreqJump(void)			{ return printFreqJump; }
	bool		getPrintStdOut(void)			{ return printStdOut; }

	bool		getUseFFT(void)					{ return useFFT; }
	int			getNumSteps(void)				{ return numSteps; }
	double		getFinalStep(void)				{ return finalStep; }
	double		getStartStep(void)				{ return startStep; }
	double		getTipStdDev(void)				{ return tipStdDev; }

	int			getNumTaxa(void)				{ return numTaxa; }
	int			getNumTimepoints(void)			{ return numTimepoints; }

	int			getModelType(void)				{ return modelType; }
	void		print(void);

private:

	void		setArguments(int argc, char** argv);

	// I/O
	std::string inputDirPath;
	std::string outputDirPath;
	std::string settingsFilePath;
	std::string exprFileName;
	std::string treeFileName;
	std::string taxaFileName;
	std::string outputFileName;
	std::string simName;

	// Model
	int seed;
	int modelType;
	bool fixBranches;
	double tuningBM;
	bool useJumpKernel;

	// Stepping Stone
	bool useSteppingStone;
	int numSteppingStone;
	double betaSteppingStone;

	// CRP
	double alphaCRP;
	double aCRP;
	double bCRP;
	double auxCRP;
	bool useCRP;

	// MCMC
	int numCycles;
	int printFreqMH;
	int printFreqCRP;
	int printFreqJump;
	bool printStdOut;

	// FFT
	bool useFFT;
	int numSteps;
	double finalStep;
	double startStep;
	double tipStdDev;

	// Expression
	int numTaxa;
	int numTimepoints;
};

#endif /* SETTINGS_H_ */

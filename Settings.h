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
#include <sstream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>

#include "FileMgr.h"

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

	double		getAlphaCRP(void)				{ return alphaCRP; }
	double		getAuxCRP(void)					{ return auxCRP; }
	double		getACRP(void)					{ return aCRP; }
	double		getBCRP(void)					{ return bCRP; }
	bool		getUseCRP(void)					{ return useCRP; }

	int			getNumCycles(void)				{ return numCycles; }
	int			getPrintFreqMH(void)			{ return printFreqMH; }
	int			getPrintFreqCRP(void)			{ return printFreqCRP; }

	bool		getUseFFT(void)					{ return useFFT; }
	int			getNumSteps(void)				{ return numSteps; }
	double		getFinalStep(void)				{ return finalStep; }
	double		getStartStep(void)				{ return startStep; }
	double		getTipStdDev(void)				{ return tipStdDev; }

	int			getNumTaxa(void)				{ return numTaxa; }
	int			getNumTimepoints(void)			{ return numTimepoints; }

	int			getModelType(void)				{ return modelType; }


private:

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

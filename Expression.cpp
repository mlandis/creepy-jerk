/*
 * Expression.cpp
 *
 *  Created on: Mar 7, 2011
 *      Author: mlandis
 */

#include "Expression.h"

// strategies on allocating large scale data containers

Expression::Expression(Settings* sp) {

	// load in data from settingsPtr
	settingsPtr = sp;
	numTimepoints = settingsPtr->getNumTimepoints();
	useCRP = settingsPtr->getUseCRP();
	tipStdDev = settingsPtr->getTipStdDev();

	initializeTaxaNames();
	initializeExprData();

}

Expression::~Expression(void)
{

}


void Expression::initializeTaxaNames(void)
{

	if (settingsPtr->getPrintStdOut()) std::cout << "INITIALIZING: Taxa names\n";

	// read data from file
	std::string taxaFileName = settingsPtr->getInputDirPath() + settingsPtr->getTaxaFileName();
	FileMgr* taxaFile = new FileMgr(taxaFileName);
	std::ifstream taxaStream;
	if (taxaFile->openFile(taxaStream) == false)
	{
		std::cerr << "Cannot open file \"" + taxaFile->getFileName() + "\"\n";
		exit(1);
	}

	int taxonCount = 0;
	std::string taxonLine = "";
	while (taxaStream.good())
	{
		std::getline(taxaStream, taxonLine);
		taxonNames.push_back(taxonLine);
		taxonCount++;
	}

	numTaxa = taxonCount;

	if (settingsPtr->getPrintStdOut())
		std::cout << "\tRead in " << numTaxa << " taxa.\n\n";

}

void Expression::initializeExprData(void)
{

	if (settingsPtr->getPrintStdOut())
		std::cout << "INITIALIZING: Expression data\n";


	// Initialize FFT settings for tip data
	int numSteps = settingsPtr->getNumSteps();
	int halfSteps = numSteps / 2;
	double finalStep = settingsPtr->getFinalStep();
	double startStep = settingsPtr->getStartStep();
	double stepSize = (finalStep - startStep) / (numSteps - 1);
	double* theta = new double[numSteps];
	for (int i = 0; i < numSteps; i++)
	{
		theta[i] = (((i + halfSteps) % numSteps) - halfSteps) * stepSize;
	}

	// Prepare to assign Expression data to Patrons
	std::vector<Patron*> patronVector;

	// Prepare to assign all Patrons to Tables (via CRP or not)
	Table* firstTable;
	if (useCRP)
	{
		firstTable = new Table(NULL,-1); // TEST 05/05/11
	}
	else if (!useCRP)
	{
		firstTable = new Table(&tableList, 0);
		tableList.push_back(firstTable);
	}

	// Read data from file
	std::string exprFileName = settingsPtr->getInputDirPath() + settingsPtr->getExprFileName();
	FileMgr* exprFile = new FileMgr(exprFileName);
	std::ifstream exprStream;
	if (exprFile->openFile(exprStream) == false)
	{
		std::cerr << "Cannot open file \"" + exprFile->getFileName() + "\"\n";
		exit(1);
	}

	// Initialize file parsing variables
	std::string field = "";
	std::string line = "";
	std::string geneName = "";
	std::vector<std::vector<double> > tempValues;

	int lineCount = 0;
	int taxaIndex = 0;
	int transIndex = 0;
	int timeIndex = 0;

	// std::cout << "Reading " << exprFile->getFileName() << "\n";
	while (exprStream.good())
	{
		std::getline(exprStream, line);
		std::istringstream exprLine(line);

		while (exprLine.good())
		{
			exprLine >> field;
			timeIndex = 0;

			// gene name
			if (exprLine.peek() == EOF && lineCount % (numTaxa + 1) == 0) // NOTE: +1 factor for geneName + numTaxa
			{
				taxaIndex = 0;
				geneName = field;
				tempValues.clear();
				tempValues.resize(numTaxa);
				//std::cout << "\t" << geneName << "\n";
			}

			// expression data
			else
			{
				tempValues[taxaIndex].push_back(atof(field.c_str()));
				timeIndex++;
				//std::cout << "\t\t" << atof(field.c_str()) << "\n";

				// if there are no new timepoints, this taxon locus data has been read
				if (exprLine.peek() == EOF)
				{
					taxaIndex++;

					// if all taxa have been read, instantiate a new Patron
					if (taxaIndex == numTaxa)
					{
						patronList.push_back(new Patron(firstTable, tempValues, theta, numSteps, tipStdDev, geneName, transIndex));
						transIndex++;
						//std::cout << "\tADDED.\n";
					}
				}

			}
		}
		lineCount++;
	}

	// push the final patron from the dataset onto the initialization table
	//patronVector.push_back(new Patron(firstTable, tempValues, theta, numSteps, geneName, index_trans));
	//std::cout << "Adding: " << geneName << " " << index_trans << "\n";
	std::random_shuffle(patronVector.begin(), patronVector.end());
	// patronList.push_back(new Patron(firstTable, tempValues, theta, numSteps, geneName)); // TEST 04/28/11

	for (unsigned int i = 0; i < patronVector.size(); i++)
	{
		patronList.push_back(patronVector[i]);
	}

	// Read in the number of transcripts and timepoints as determined by the input data.
	numTranscripts = patronList.size();
	//std::cout << timeIndex << "\n"; exit(1);
	numTimepoints = timeIndex;

	if (settingsPtr->getPrintStdOut())
		std::cout << "\tRead in " << numTranscripts << " transcripts.\n\n";

	delete [] theta;
	delete exprFile;
}

int Expression::getIndexForTaxon(std::string queryName)
{
	for (unsigned int i = 0; i < taxonNames.size(); i++)
	{
		if (taxonNames[i] == queryName)
		{
			return i;
		}
	}
	return -1;
}

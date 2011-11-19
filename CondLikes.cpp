/*
 * CondLikes.cpp
 *
 *  Created on: Mar 30, 2011
 *      Author: mlandis
 */

#include "CondLikes.h"

CondLikes::CondLikes(Expression* ep, Settings* sp) {

	std::cout << "INITIALIZING: Conditional likelihoods\n";

	expressionPtr = ep;
	settingsPtr = sp;

	// assign dimension sizes
	// [space][transcript][node][timepoint][step]
	numSpaces = 2;
	numTranscripts = expressionPtr->getNumTranscripts();
	numTaxa = expressionPtr->getNumTaxa();
	numNodes = 2 * numTaxa - 1;
	numTimepoints = expressionPtr->getNumTimepoints();
	numSteps = settingsPtr->getNumSteps();

	// allocate memory
	int oneTimeSize = numSteps * 2;
	int oneNodeSize = numTimepoints * oneTimeSize;
	int oneTransSize = numNodes * oneNodeSize;
	int oneClSize = numTranscripts * oneTransSize;

	std::cout << "\tnumTaxa:      " << numTaxa << "\n";
	std::cout << "\toneClSize:    " << oneClSize << "\n";
	std::cout << "\toneTransSize: " << oneTransSize << "\n";
	std::cout << "\toneNodeSize:  " << oneNodeSize << "\n";
	std::cout << "\toneTimeSize:  " << oneTimeSize << "\n";

	cls = new double[2 * oneClSize];
	for (int i = 0; i < 2 * oneClSize; i++)
		cls[i] = 0.0;
	clsUpL = new double[2 * oneClSize];
	for (int i = 0; i < 2 * oneClSize; i++)
		clsUpL[i] = 0.0;
	clsUpR = new double[2 * oneClSize];
	for (int i = 0; i < 2 * oneClSize; i++)
		clsUpR[i] = 0.0;



	// initialize pointer information

	clsPtr = new double****[numSpaces];
	for (int n = 0; n < numSpaces; n++)
	{
		clsPtr[n] = new double***[numTranscripts];
		for (int i = 0; i < numTranscripts; i++)
		{
			clsPtr[n][i] = new double**[numNodes];
			for (int j = 0; j < numNodes; j++)
			{
				clsPtr[n][i][j] = new double*[numTimepoints];
				for (int k = 0; k < numTimepoints; k++)
				{
					int m = n * oneClSize + i * oneTransSize + j * oneNodeSize + k * oneTimeSize;
					//std::cout << n << " " << i << " " << j << " " << k << " : " << m << "\n";
					clsPtr[n][i][j][k] = &cls[m];
					//std::cout << "cls[" << m << "]: " << cls[m] << "\n";
				}
			}
		}
	}

	clsPtrUpL = new double****[numSpaces];
	for (int n = 0; n < numSpaces; n++)
	{
		clsPtrUpL[n] = new double***[numTranscripts];
		for (int i = 0; i < numTranscripts; i++)
		{
			clsPtrUpL[n][i] = new double**[numNodes];
			for (int j = 0; j < numNodes; j++)
			{
				clsPtrUpL[n][i][j] = new double*[numTimepoints];
				for (int k = 0; k < numTimepoints; k++)
				{
					int m = n * oneClSize + i * oneTransSize + j * oneNodeSize + k * oneTimeSize;
					clsPtrUpL[n][i][j][k] = &clsUpL[m];
				}
			}
		}
	}

	clsPtrUpR = new double****[numSpaces];
	for (int n = 0; n < numSpaces; n++)
	{
		clsPtrUpR[n] = new double***[numTranscripts];
		for (int i = 0; i < numTranscripts; i++)
		{
			clsPtrUpR[n][i] = new double**[numNodes];
			for (int j = 0; j < numNodes; j++)
			{
				clsPtrUpR[n][i][j] = new double*[numTimepoints];
				for (int k = 0; k < numTimepoints; k++)
				{
					int m = n * oneClSize + i * oneTransSize + j * oneNodeSize + k * oneTimeSize;
					clsPtrUpR[n][i][j][k] = &clsUpR[m];
				}
			}
		}
	}

	// initialize tip conditional likelihoods

	std::list<Patron*>::const_iterator it_p = expressionPtr->getPatronList().begin();
	gsl_complex tipVal;

	for (int i = 0; i < numTranscripts; i++)
	{
		for (int j = 0; j < numTaxa; j++)
		{
			for (int k = 0; k < numTimepoints; k++)
			{
				double* p0 = clsPtr[0][i][j][k];
				double* p1 = clsPtr[1][i][j][k];
				double* uL0 = clsPtrUpL[0][i][j][k];
				double* uL1 = clsPtrUpL[1][i][j][k];
				double* uR0 = clsPtrUpR[0][i][j][k];
				double* uR1 = clsPtrUpR[1][i][j][k];

				for (int l = 0; l < 2 * numSteps; l += 2)
				{

					tipVal = (*it_p)->getCfData(j, k, (l/2));

					// std::cout << i << " " << j << " " << k << " " << l << " " << ": " << GSL_REAL(tipVal) << " " << GSL_IMAG(tipVal) << "\n";

					p0[l] 		= GSL_REAL(tipVal);
					p0[l + 1]	= GSL_IMAG(tipVal);
					p1[l] 		= GSL_REAL(tipVal);
					p1[l + 1]	= GSL_IMAG(tipVal);
					uL0[l] 		= GSL_REAL(tipVal);
					uL0[l + 1]	= GSL_IMAG(tipVal);
					uL1[l] 		= GSL_REAL(tipVal);
					uL1[l + 1]	= GSL_IMAG(tipVal);
					uR0[l] 		= GSL_REAL(tipVal);
					uR0[l + 1]	= GSL_IMAG(tipVal);
					uR1[l] 		= GSL_REAL(tipVal);
					uR1[l + 1]	= GSL_IMAG(tipVal);

				}
			}
		}
		//std::cout << (*it_p)->getName() << " added.\n";

		it_p++;
	}

	// allocate memory for the scalers

#if 0
	clScalers = new double[2 * numTranscripts * numNodes * numTimepoints];
	for (int i = 0; i < 2 * numTranscripts * numNodes * numTimepoints; i++)
	{
		clScalers[i] = 0.0;
	}
	for (int n = 0; n < 2; n++)
	{
		clScalerPtr[n] = new double***[numTranscripts * numNodes];
		clScalerPtr[n][0] = new double**[numTranscripts * numNodes * numTimepoints];
		for (int i=1; i<numNodes; i++)
			clScalerPtr[n][i] = clScalerPtr[n][i-1] + numSites;
		for (int i=0; i<numNodes; i++)
			for (int j=0; j<numSites; j++)
				clScalerPtr[n][i][j] = &clScalers[n*(numNodes*numSites) + i*(numSites) + j];
	}
#endif


	// allocate memory for the scalers
	clScalers = new double[2 * oneClSize];
	for (int i = 0; i < 2 * oneClSize; i++)
		clScalers[i] = 0.0;
	clScalersUpL = new double[2 * oneClSize];
	for (int i = 0; i < 2 * oneClSize; i++)
		clScalersUpL[i] = 0.0;
	clScalersUpR = new double[2 * oneClSize];
	for (int i = 0; i < 2 * oneClSize; i++)
		clScalersUpR[i] = 0.0;

	// assign pointers to scalers
	clScalerPtr = new double****[numSpaces];
	for (int n = 0; n < numSpaces; n++)
	{
		clScalerPtr[n] = new double***[numTranscripts];
		for (int i = 0; i < numTranscripts; i++)
		{
			clScalerPtr[n][i] = new double**[numNodes];
			for (int j = 0; j < numNodes; j++)
			{
				clScalerPtr[n][i][j] = new double*[numTimepoints];
				for (int k = 0; k < numTimepoints; k++)
				{
					// need one clScaler per timepoint
					int m = n * oneClSize + i * oneTransSize + j * oneNodeSize + k; // + k * oneTimeSize
					//std::cout << n << " " << i << " " << j << " " << k << " : " << m << "\n";
					clScalerPtr[n][i][j][k] = &clScalers[m];
					//std::cout << "cls[" << m << "]: " << cls[m] << "\n";
				}
			}
		}
	}

	clScalerPtrUpL = new double****[numSpaces];
	for (int n = 0; n < numSpaces; n++)
	{
		clScalerPtrUpL[n] = new double***[numTranscripts];
		for (int i = 0; i < numTranscripts; i++)
		{
			clScalerPtrUpL[n][i] = new double**[numNodes];
			for (int j = 0; j < numNodes; j++)
			{
				clScalerPtrUpL[n][i][j] = new double*[numTimepoints];
				for (int k = 0; k < numTimepoints; k++)
				{
					// need one clScaler per timepoint
					int m = n * oneClSize + i * oneTransSize + j * oneNodeSize + k; // + k * oneTimeSize
					//std::cout << n << " " << i << " " << j << " " << k << " : " << m << "\n";
					clScalerPtrUpL[n][i][j][k] = &clScalersUpL[m];
					//std::cout << "cls[" << m << "]: " << cls[m] << "\n";
				}
			}
		}
	}

	clScalerPtrUpR = new double****[numSpaces];
	for (int n = 0; n < numSpaces; n++)
	{
		clScalerPtrUpR[n] = new double***[numTranscripts];
		for (int i = 0; i < numTranscripts; i++)
		{
			clScalerPtrUpR[n][i] = new double**[numNodes];
			for (int j = 0; j < numNodes; j++)
			{
				clScalerPtrUpR[n][i][j] = new double*[numTimepoints];
				for (int k = 0; k < numTimepoints; k++)
				{
					// need one clScaler per timepoint
					int m = n * oneClSize + i * oneTransSize + j * oneNodeSize + k; // + k * oneTimeSize
					//std::cout << n << " " << i << " " << j << " " << k << " : " << m << "\n";
					clScalerPtrUpR[n][i][j][k] = &clScalersUpR[m];
					//std::cout << "cls[" << m << "]: " << cls[m] << "\n";
				}
			}
		}
	}

	//print();
	std::cout << "\n";
}

CondLikes::~CondLikes(void) {

	delete [] cls;
	for (int n=0; n<2; n++)
	{
		delete [] clsPtr[n][0];
		delete [] clsPtr[n];
	}

	delete [] clsUpL;
	for (int n=0; n<2; n++)
	{
		delete [] clsPtrUpL[n][0];
		delete [] clsPtrUpL[n];
	}


	delete [] clsUpR;
	for (int n=0; n<2; n++)
	{
		delete [] clsPtrUpR[n][0];
		delete [] clsPtrUpR[n];
	}

	delete [] clScalers;
	for (int n=0; n<2; n++)
	{
		delete [] clScalerPtr[n][0];
		delete [] clScalerPtr[n];
	}

	delete [] clScalersUpL;
	for (int n=0; n<2; n++)
	{
		delete [] clScalerPtrUpL[n][0];
		delete [] clScalerPtrUpL[n];
	}

	delete [] clScalersUpR;
	for (int n=0; n<2; n++)
	{
		delete [] clScalerPtrUpR[n][0];
		delete [] clScalerPtrUpR[n];
	}
}

void CondLikes::print(void) {

	for (int j = 0; j < numNodes; j++)
	{
		std::cout << "NODE[" << j << "]:\n";
		for (int i = 0; i < numTranscripts; i++)
		{
			for (int k = 0; k < numTimepoints; k++)
			{
				std::cout << std::setw(5) << i << ",";
				std::cout << std::setw(5) << k << " -- ";
				double* p0 = getClPtr(0, i, j, k);
				for (int s = 0; s < 2 * numSteps; s += 2)
				{
					std::cout << std::fixed << std::setprecision(6) << "(" << p0[s] << "," << p0[s+1] << ") ";
				}
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}
}


// GRAVEYARD.


/*
	clsPtrUpL[n] = new double**[numNodes];
	clsPtrUpL[n][0] = new double*[numNodes * numSites];
	for (int i=1; i<numNodes; i++)
		clsPtrUpL[n][i] = clsPtrUpL[n][i-1] + numSites;
	for (int i=0; i<numNodes; i++)
		for (int j=0; j<numSites; j++)
			clsPtrUpL[n][i][j] = &clsUpL[n*oneClSize + i*oneNodeSize + j*numStates];

	clsPtrUpR[n] = new double**[numNodes];
	clsPtrUpR[n][0] = new double*[numNodes * numSites];
	for (int i=1; i<numNodes; i++)
		clsPtrUpR[n][i] = clsPtrUpR[n][i-1] + numSites;
	for (int i=0; i<numNodes; i++) b
		for (int j=0; j<numSites; j++)
			clsPtrUpR[n][i][j] = &clsUpR[n*oneClSize + i*oneNodeSize + j*numStates];
	 */




#if 0
		for (int n = 0; n < 2; n++)
		{

		/*
		clsPtr[n] = new double***[numTranscripts];
		clsPtr[n][0] = new double**[numTranscripts * numNodes * numTimepoints];
		clsPtr[n][0][0] = new double*[numNodes * numTimepoints];
		for (int i = 1; i < numTranscripts; i++)
		{
			clsPtr[n][i] = clsPtr[n][i-1] + (numNodes * numTimepoints);
			//clsPtr[n][i][0] = new double*[numNodes * numTimepoints];
			for (int j = 1; j < numNodes; j++)
			{
				std::cout << n << " " << i << " " << j << "\n";
				clsPtr[n][i][j] = clsPtr[n][i][j-1] + numTimepoints;
				//clsPtr[n][i][j][0] = new double[numTimepoints];
				//for (int k = 1; k < numTimepoints; k++)
				//{
				//	clsPtr[n][i][j][k] = clsPtr[n][i][j][k-1] + num
				//}
			}
		}
		std::cout << "ok2\n";
		for (int i = 0; i < numTranscripts; i++)
		{
			for (int j = 0; j < numNodes; j++)
			{
				for (int k = 0; k < numTimepoints; k++)
				{
					std::cout << i << " " << j << " " << k << " " << "\n";
					std::cout << oneClSize << " " << oneTransSize << " " << oneNodeSize << " " << oneTimeSize << "\n";
					std::cout << n*oneClSize + i*oneTransSize + j*oneNodeSize + k*oneTimeSize << "\n";
					clsPtr[n][i][j][k] = &cls[n*oneClSize + i*oneTransSize + j*oneNodeSize + k*oneTimeSize];
				}
			}
		}

		*/

		std::cout << "1\n";

		clsPtrUpL[n] = new double***[numNodes];
		clsPtrUpL[n][0] = new double**[numTranscripts * numNodes * numTimepoints];
		for (int i = 1; i < numTranscripts; i++)
		{
			clsPtrUpL[n][i] = clsPtrUpL[n][i-1] + (numNodes * numTimepoints);
			clsPtrUpL[n][i][0] = new double*[numNodes * numTimepoints];
			for (int j = 1; j < numNodes; j++)
			{
				clsPtrUpL[n][i][j] = clsPtrUpL[n][i][j-1] + numTimepoints;
				//clsPtr[n][i][j][0] = new double
			}
		}
		for (int i = 0; i < numTranscripts; i++)
		{
			for (int j = 0; j < numNodes; j++)
			{
				for (int k = 0; k < numTimepoints; k++)
				{
					for (int l = 0; l < numSteps; l++)
					{
						clsPtrUpL[n][i][j][k] = &clsUpL[n*oneClSize + i*oneTransSize + j*oneNodeSize + k*oneTimeSize];
					}
				}
			}
		}

		std::cout << "2\n";

		clsPtrUpR[n] = new double***[numNodes];
		clsPtrUpR[n][0] = new double**[numTranscripts * numNodes * numTimepoints];
		for (int i = 1; i < numTranscripts; i++)
		{
			clsPtrUpR[n][i] = clsPtrUpR[n][i-1] + (numNodes * numTimepoints);
			clsPtrUpR[n][i][0] = new double*[numNodes * numTimepoints];
			for (int j = 1; j < numNodes; j++)
			{
				clsPtrUpR[n][i][j] = clsPtrUpR[n][i][j-1] + numTimepoints;
				//clsPtr[n][i][j][0] = new double
			}
		}
		for (int i = 0; i < numTranscripts; i++)
		{
			for (int j = 0; j < numNodes; j++)
			{
				for (int k = 0; k < numTimepoints; k++)
				{
					for (int l = 0; l < numSteps; l++)
					{
						clsPtrUpR[n][i][j][k] = &clsUpR[n*oneClSize + i*oneTransSize + j*oneNodeSize + k*oneTimeSize];
					}
				}
			}
		}

		std::cout << "3\n";
#endif

#if 0

	double***** clsPtr; // [space][trans][node][time][step]
	double* cls;

	int oneTimeSize = 2 * numSteps;
	int oneNodeSize = oneTimeSize * numTime;
	int oneTransSize = oneNodeSize * numNodes;
	int oneClSize = oneTransSize * numTrans;

	std::cout << "oneClSize:    " << oneClSize << "\n";
	std::cout << "oneTransSize: " << oneTransSize << "\n";
	std::cout << "oneNodeSize:  " << oneNodeSize << "\n";
	std::cout << "oneTimeSize:  " << oneTimeSize << "\n";

	// set cls
	std::cout << "Initializing values.\n";
	cls = new double[2 * oneClSize];
	for (int i = 0; i < 2 * oneClSize; i++)
	{
		cls[i] = i * 1.0;
		std::cout << " " << std::setw(8) << std::setprecision(8) << cls[i];
		if (i % 10 == 9) std::cout << "\n";
	}
	std::cout << "\n\n";

	std::cout << "Setting pointers.\n";

	//clsPtr = new double[numSpaces][numTrans][numNodes][numTime][numSteps];


	clsPtr = new double****[numSpaces];
	for (int n = 0; n < numSpaces; n++)
	{
		clsPtr[n] = new double***[numTrans];
		for (int i = 0; i < numTrans; i++)
		{
			clsPtr[n][i] = new double**[numNodes];
			for (int j = 0; j < numNodes; j++)
			{
				clsPtr[n][i][j] = new double*[numTime];
				for (int k = 0; k < numTime; k++)
				{
					int m = n * oneClSize + i * oneTransSize + j * oneNodeSize + k * oneTimeSize;
					std::cout << n << " " << i << " " << j << " " << k << " : " << m << "\n";
					clsPtr[n][i][j][k] = NULL;
					std::cout << "cls[" << m << "]: " << cls[m] << "\n";
					clsPtr[n][i][j][k] = &cls[m];
					std::cout << "cls[" << m << "]: " << cls[m] << "\n";
				}
			}
		}
	}


#endif
		/*
		clsPtrUpL[n] = new double**[numNodes];
		clsPtrUpL[n][0] = new double*[numNodes * numSites];
		for (int i=1; i<numNodes; i++)
			clsPtrUpL[n][i] = clsPtrUpL[n][i-1] + numSites;
		for (int i=0; i<numNodes; i++)
			for (int j=0; j<numSites; j++)
				clsPtrUpL[n][i][j] = &clsUpL[n*oneClSize + i*oneNodeSize + j*numStates];

		clsPtrUpR[n] = new double**[numNodes];
		clsPtrUpR[n][0] = new double*[numNodes * numSites];
		for (int i=1; i<numNodes; i++)
			clsPtrUpR[n][i] = clsPtrUpR[n][i-1] + numSites;
		for (int i=0; i<numNodes; i++) b
			for (int j=0; j<numSites; j++)
				clsPtrUpR[n][i][j] = &clsUpR[n*oneClSize + i*oneNodeSize + j*numStates];
				*/
	//}


	/*
	for (int i=0; i<numTaxa; i++)
	{
		for (int j=0; j<numSites; j++)
		{
			double* p0 = clsPtr[0][i][j];
			double* p1 = clsPtr[1][i][j];
			double* uL0 = clsPtrUpL[0][i][j];
			double* uL1 = clsPtrUpL[1][i][j];
			double* uR0 = clsPtrUpR[0][i][j];
			double* uR1 = clsPtrUpR[1][i][j];

			MbBitfield *bf = alignmentPtr->getCodonMatrixEntry(i, j);
			for (int k=0; k<numOmegaCategories; k++)
			{
				for (int s=0; s<alignmentPtr->getNumSenseCodons(); s++)
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
				p0 += alignmentPtr->getNumSenseCodons();
				p1 += alignmentPtr->getNumSenseCodons();
				uL0 += alignmentPtr->getNumSenseCodons();
				uL1 += alignmentPtr->getNumSenseCodons();
				uR0 += alignmentPtr->getNumSenseCodons();
				uR1 += alignmentPtr->getNumSenseCodons();
			}
		}
	}
*/


	/*
		// allocate memory for the scalers
		clScalers = new double[2 * numTranscripts * numNodes * numTimepoints];
		for (int i = 0; i < 2 * numTranscripts * numNodes * numTimepoints; i++)
		{
			clScalers[i] = 0.0;
		}
		for (int n = 0; n < 2; n++)
		{
			clScalerPtr[n] = new double***[numTranscripts * numNodes];
			clScalerPtr[n][0] = new double**[numTranscripts * numNodes * numTimepoints];
			for (int i=1; i<numNodes; i++)
				clScalerPtr[n][i] = clScalerPtr[n][i-1] + numSites;
			for (int i=0; i<numNodes; i++)
				for (int j=0; j<numSites; j++)
					clScalerPtr[n][i][j] = &clScalers[n*(numNodes*numSites) + i*(numSites) + j];
		}

		clScalersUpL = new double[2 * numNodes * numSites];
		for (int i=0; i<2*numNodes*numSites; i++)
			clScalersUpL[i] = 0.0;
		for (int n=0; n<2; n++)
		{
			clScalerPtrUpL[n] = new double**[numNodes];
			clScalerPtrUpL[n][0] = new double*[numNodes * numSites];
			for (int i=1; i<numNodes; i++)
				clScalerPtrUpL[n][i] = clScalerPtrUpL[n][i-1] + numSites;
			for (int i=0; i<numNodes; i++)
				for (int j=0; j<numSites; j++)
					clScalerPtrUpL[n][i][j] = &clScalersUpL[n*(numNodes*numSites) + i*(numSites) + j];
		}

		clScalersUpR = new double[2 * numNodes * numSites];
		for (int i=0; i<2*numNodes*numSites; i++)
			clScalersUpR[i] = 0.0;
		for (int n=0; n<2; n++)
		{
			clScalerPtrUpR[n] = new double**[numNodes];
			clScalerPtrUpR[n][0] = new double*[numNodes * numSites];
			for (int i=1; i<numNodes; i++)
				clScalerPtrUpR[n][i] = clScalerPtrUpR[n][i-1] + numSites;
			for (int i=0; i<numNodes; i++)
				for (int j=0; j<numSites; j++)
					clScalerPtrUpR[n][i][j] = &clScalersUpR[n*(numNodes*numSites) + i*(numSites) + j];
		}
		*/

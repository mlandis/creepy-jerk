/*
 * main.cpp
 *
 *  Created on: Mar 7, 2011
 *      Author: mlandis
 *
 */

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#define RESEAT 0
#define LAMBDA_ONLY 0
#define SIGMA_ONLY 0
#define TAU1_ONLY 0
#define TAU2_ONLY 1

#include <iostream>
#include <math.h>

#include "Expression.h"
#include "FileMgr.h"
#include "MbRandom.h"
#include "Model.h"
#include "Mcmc.h"
#include "Settings.h"
#include "Parm_tree.h"

int main(int argc, char *argv[])
{

	std::cout << "test2\n";

	Settings* mySettings = new Settings(argc, argv);

	MbRandom* myRandom = new MbRandom(mySettings->getSeed());

	Expression* myExpression = new Expression(mySettings);

	Topology* myTopology = new Topology(myRandom, myExpression, mySettings, 1.0);

	Model* myModel = new Model(myExpression, myRandom, mySettings, myTopology);

	Mcmc* myMcmc = new Mcmc(myExpression, myModel, mySettings, myRandom, myTopology);

	delete myMcmc;
	delete myModel;
	delete myExpression;
	delete myRandom;
	delete mySettings;

	std::cout << "COMPLETE.\n";

	return 0;
}

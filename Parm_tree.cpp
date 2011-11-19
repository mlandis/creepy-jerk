/*
 * Parm_tree.cpp
 *
 *  Created on: Jan 24, 2010
 *      Author: mlandis
 */

#include "Parm_tree.h"

Node::Node(void) {

	lft			= NULL;
	rht			= NULL;
	anc			= NULL;
	index		= 0;
	dnPassIndex	= 0;
	q			= 0;
	v			= 1.0;
	ratio		= 0.0;
	name		= "";
	activeCl	= 0;
	activeTi	= 0;
	updateCl	= false;
	updateTi	= false;
	flag		= false;
	part		= new MbBitfield;

	// node values
	activeParm = 0;
	updateParm = false;
	mu = 0.0;
	sigma = 0.0;
	k = 0.0;

	// sampled jump values
	jumpCount[0] = 0;
	//lnProbJumpCount[0] = 0.0;
	//jumpSize[0].push_back(0.0);
	//lnProbJumpSize[0].push_back(0.0);
	sumJumpSize[0] = 0.0;
	//sumLnProbJumpSize[0] = 0.0;

}

void Node::addJump(double j, int space)
{
	std::cout << "OK333\n";

	jumpCount[space] = jumpCount[space] + 1;
	std::cout << "OK3\n";
	//lnProbJumpCount[space] = lnPrJumpCount;
	std::cout << "OK3\n";

	jumpSize[space].push_back(j);
	//lnProbJumpSize[space].push_back(lnPrJumpSize);

	//std::cout << "OK4\n";
	sumJumpSize[space] += j;
	//sumLnProbJumpSize[space] += lnPrJumpSize;
	//std::cout << "OK5\n";
}

void Node::removeJump(int j, int space)
{
	double x = jumpSize[space][j];
	//double p = lnProbJumpSize[space][j];

	jumpCount[space] = jumpCount[space] - 1;
	//lnProbJumpCount[space] = lnPrJumpCount;

	jumpSize[space].erase(jumpSize[space].begin() + j);
	//lnProbJumpSize[space].erase(lnProbJumpSize[space].begin() + j);
	sumJumpSize[space] -= x;
	//sumLnProbJumpSize[space] -= p;

}

void Node::copySpace(void)
{
	int inactiveParm = activeParm == 0 ? 1 : 0;
	jumpCount[activeParm] = jumpCount[inactiveParm];
	jumpSize[activeParm] = jumpSize[inactiveParm];
	sumJumpSize[activeParm] = sumJumpSize[inactiveParm];
}

Node::~Node(void) {

	delete part;
}


Topology::Topology(MbRandom *rp, Expression *ep, double vp) {

	std::cout << "INITIALIZING: Topology\n";

	// set base state for some pointers
	nodes = NULL;
	downPassSequence = NULL;

	// initialize some basic tree variables
	ranPtr = rp;
	numTaxa = ep->getNumTaxa();
	numNodes = 2 * numTaxa - 2; // -2 accounts for explicitly declared root
	brlenLambda = vp;
	treeLength = 0.0;
	numJumps = 0;

	// build a random tree
	buildRandomTree(ep);
	initializeTaxonBipartitions();
	setBranchRatios();
	getTaxonBipartitions();
	print();

	std::cout << "\n";
}

Topology::Topology(MbRandom *rp, Expression *ep, Settings *sp, double vp) {

	std::cout << "INITIALIZING: Topology\n";

	// set base state for some pointers
	nodes = NULL;
	downPassSequence = NULL;

	// initialize some basic tree variables
	ranPtr = rp;
	settingsPtr = sp;
	numTaxa = ep->getNumTaxa();
	numNodes = 2 * numTaxa - 1; // -2 accounts for explicitly declared root
	brlenLambda = vp;
	treeLength = 0.0;
	numJumps = 0;

	// get Newick string
	std::ifstream treeFileStream;
	std::string treeFileName = settingsPtr->getInputDirPath() + settingsPtr->getTreeFileName();
	FileMgr treeFileMgr = FileMgr(treeFileName);
	if (treeFileMgr.openFile(treeFileStream) == false)
	{
		std::cerr << "ERROR: Could not open file: " << treeFileName << "\n";
	}
	std::string treeStr = "";
	getline(treeFileStream, treeStr);
	treeFileMgr.closeFile(treeFileStream);

	// build tree from string
	buildTreeFromNewickString(ep, treeStr);
	initializeTaxonBipartitions();
	setBranchRatios();
	getTaxonBipartitions();
	print();

	std::cout << "\n";
}

Topology::Topology(Topology &t) {

	numTaxa = 0;
	numNodes = 0;
	nodes = NULL;
	downPassSequence = NULL;
	clone(t);
}

Topology::~Topology(void) {

	delete [] nodes;
	delete [] downPassSequence;
}

Topology& Topology::operator=(const Topology &t) {

	if (this != &t)
		clone(t);
	return *this;
}

void Topology::buildRandomTree(Expression *ep) {

	if (numTaxa < 3) {
		std::cerr << "ERROR: Too few taxa!" << std::endl;
		exit(1);
	}

	// allocate the nodes & down pass sequence vector
	nodes = new Node[numNodes];
	downPassSequence = new Node*[numNodes];

	// set the node indices
	for (int i = 0; i < numNodes; i++) {
		nodes[i].setIndex(i);
	}

	for (int i = 0; i < numNodes; i++) {
	//	nodes[i].setPathSize(ep->getNumCodonSites());
		// nodes[i].initPathSize(ep->getNumCodonSites()); // removed 060111
	}

	// set the taxon names
	for (int i = 0; i < numTaxa; i++)
		;// nodes[i].setName(ep->getNameForTaxon(i)); // removed 060111

	// build the three species tree
	Node **availableNodes = new Node*[numNodes];
	availableNodes[0] = NULL;
	int numAvailableNodes = 0;
	int nextIntNode = numTaxa;
	Node *p = &nodes[0];
	root = p;
	Node *q = &nodes[nextIntNode++];
	availableNodes[numAvailableNodes++] = q;
	p->setLft(q);
	q->setAnc(p);
	p = &nodes[1];
	availableNodes[numAvailableNodes++] = p;
	q->setLft(p);
	q->setAnc(q);
	p = &nodes[2];
	availableNodes[numAvailableNodes++] = p;
	q->setRht(p);
	p->setAnc(q);

	// add the remaining taxa to the tree
	int nextTipNode = 3;
	for (int i = 3; i < numTaxa; i++) {
		// pick a branch to attach the new branch to
		int whichNode = (int)(ranPtr->uniformRv()*numAvailableNodes);
		p = availableNodes[whichNode];
		Node *pAnc = p->getAnc();
		Node *newTip = &nodes[nextTipNode++];
		Node *newInt = &nodes[nextIntNode++];

		// p is to the left of pAnc
		if (p == pAnc->getLft()) {
			pAnc->setLft(newInt);
			newInt->setAnc(pAnc);
			newInt->setLft(p);
			p->setAnc(newInt);
			newInt->setRht(newTip);
			newTip->setAnc(newInt);
		}

		// p is to the right of pAnc
		else {
			pAnc->setRht(newInt);
			newInt->setAnc(pAnc);
			newInt->setRht(p);
			p->setAnc(newInt);
			newInt->setLft(newTip);
			newTip->setAnc(newInt);
		}

	}

	delete [] availableNodes;

	// initialize the down pass sequence for the tree
	getDownPassSequence();

	// set the branch lengths
	for (int n = 0; n < numNodes; n++) {
		p = &nodes[n];
		if (p->getAnc() != NULL)
		{
			p->setV(ranPtr->exponentialRv(brlenLambda));
			treeLength += nodes[n].getV();
		}
	}
}

void Topology::buildTreeFromNewickString(Expression *ep, std::string &ts) {

	// parse the tree string and put each token into a vector of strings
	std::vector<std::string> parsedNewick;
	std::string temp = "";
	bool readingBrlen = false;
	int nt = 0;

	for (int i = 0; i < (int)ts.size(); i++) {
		char c = ts[i];

		// ignore character (whitespace)
		if (c == ' ')
			continue;

		// the character is punctuation
		if (c == '(' || c == ')' || c == ',' || c == ':' || c == ';') {
			temp = c;
			parsedNewick.push_back(temp);
			if (c == ':')
				readingBrlen = true;
			else
				readingBrlen = false;
		}

		// the character is part of a taxon name
		else {
			int j = i;
			std::string taxonName = "";
			while (ts[j] != '(' && ts[j] != ')' && ts[j] != ',' && ts[j] != ':' && ts[j] != ';') {
				taxonName += ts[j];
				j++;
			}
			parsedNewick.push_back(taxonName);
			i = j -1;
			if (readingBrlen == false)
				nt++;
			readingBrlen = false;
		}

		if (c == ';')
			break;
	}


	// check that the number of taxa in the tree description is the same as the number of taxa in the alignment
	if (nt != numTaxa) {
		std::cerr << "ERROR: The tree file is not the right size" << std::endl;
		std::cout << "nt = " << nt << " numTaxa = " << numTaxa << std::endl;
		exit(1);
	}

	// show the tokens

	// for (std::vector<std::string>::iterator p = parsedNewick.begin(); p != parsedNewick.end(); p++)
	//	std::cout << "token = \"" << (*p) << "\"" << "\n";

	// allocate the nodes
	nodes = new Node[numNodes];
	downPassSequence = new Node*[numNodes];

	// set the node indices
	for (int i = 0; i < numNodes; i++) {
		nodes[i].setIndex(i);
	}

	// initialize the branch lengths
	for (int i = 0; i < numNodes; i++) {
		nodes[i].setV(ranPtr->exponentialRv(brlenLambda));
		treeLength += nodes[i].getV();
	}

	// set the taxon names
	for (int i = 0; i < numTaxa; i++)
	{
		nodes[i].setName(ep->getNameForTaxon(i));
	}

	// build up the tree using the information stored in the parsed vector
	//int nextInteriorNode = numTaxa; // CHANGED: 060611 MJL
	int nextInteriorNode = numNodes-1;
	Node *p = NULL;

	for (std::vector<std::string>::iterator t = parsedNewick.begin(); t != parsedNewick.end(); t++) {
		// std::cout << "Token: " << (*t) << "\n";
		// add a new interior node
		if ((*t) == "(") {
			if (p == NULL) {
				//p = &nodes[nextInteriorNode++]; // CHANGED: 060611 MJL
				p = &nodes[nextInteriorNode--];
			}
			else {
				//Node *q = &nodes[nextInteriorNode++];
				Node *q = &nodes[nextInteriorNode--]; // CHANGED: 060611 MJL
				if (p->getLft() == NULL) {
					p->setLft(q);
					q->setAnc(p);
					p = q;
				}
				else if (p->getRht() == NULL) {
					p->setRht(q);
					q->setAnc(p);
					p = q;
				}
				else if (p->getAnc() == NULL) {
					p->setAnc(q);
					q->setLft(p);
					p = q;
					p->setFlag(true);
				}
				else {
					std::cout << "ERROR: Problem reading the Newick-formatted tree (1)" << std::endl;
					exit(1);
				}
			}
			readingBrlen = false;
		}

		// we hit a right parenthesis, so we should go down the tree -- unless!
		// we should go up
		else if ((*t) == ")") {
			if (p->getFlag() == false && p->getAnc() != NULL)
				p = p->getAnc();
			else if (p->getFlag() == true && p->getLft() != NULL)
				p = p->getLft();
			else {
				std::cout << "ERROR: Problem reading the Newick-formatted tree (2)" << std::endl;
				exit(1);
			}
			readingBrlen = false;
		}

		// we hit a comma, so we should go down the tree -- unless!
		// we should go up
		else if ((*t) == ",") {
			if (p->getFlag() == false && p->getAnc() != NULL)
				p = p->getAnc();
			else if (p->getFlag() == true && p->getLft() != NULL)
				p = p->getLft();
			else {
				std::cout << "ERROR: Problem reading the Newick-formatted tree (3)" << std::endl;
				exit(1);
			}
			readingBrlen = false;
		}

		else if ((*t) == ":") {
			readingBrlen = true;
		}

		// we are at the end of the tree description, nothing to do
		else if ((*t) == ";") {	}

		else {
			// read in a taxon name and add the node to the tree
			if (readingBrlen == false) {

				std::string theName = (*t);
				int theIndex = ep->getIndexForTaxon((*t));
				//int theIndex = 0;
				if (theIndex == -1) {
					std::cerr << "ERROR: Could not find taxon " << (*t) << " in the Expression" << std::endl;
					exit(1);
				}
				Node *q = &nodes[theIndex];
				q->setName(theName);

				if (p->getLft() == NULL) {
					p->setLft(q);
					q->setAnc(p);
					p = q;
				}
				else if (p->getRht() == NULL) {
					p->setRht(q);
					q->setAnc(p);
					p = q;
				}
				else if (p->getAnc() == NULL) {
					p->setAnc(q);
					q->setLft(p);
					p = q;
					p->setFlag(true);
				}
				else {
					std::cout << "ERROR: Problem reading the Newick-formatted tree (4)" << std::endl;
					exit(1);
				}
			}

			// reading a branch length
			else {
				double x = 0.0;
				std::istringstream buf(*t);
				buf >> x;
				if (x < 0.00001)
					x = 0.0001;
				if (p->getFlag() == false)
					p->setV(x);
				else
					p->getLft()->setV(x);
				readingBrlen = false;
				treeLength += x;
			}
		}
	}

	// set the pointer to the root
	Node *q = p;
	while (q->getAnc() != NULL)
		q = q->getAnc();
	root = q;

	// initialize the down pass sequence for the tree
	getDownPassSequence();

	// set the root branch length to zero
	root->setV(0.0);
}

void Topology::clone(const Topology &t) {

	numTaxa		= t.numTaxa;
	numNodes	= t.numNodes;
	brlenLambda	= t.brlenLambda;
	treeLength  = t.treeLength;
	numJumps    = t.numJumps;
	ranPtr		= t.ranPtr;

	if (nodes == NULL)
		nodes = new Node[t.numNodes];
	if (downPassSequence == NULL)
		downPassSequence = new Node*[t.numNodes];

	// copy from p-nodes to the q-nodes
	for (int n = 0; n < numNodes; n++) {
		Node *p = &t.nodes[n];
		Node *q = &nodes[n];

		if (p->getLft() == NULL)
			q->setLft(NULL);
		else
			q->setLft(&nodes[p->getLft()->getIndex()]);

		if (p->getRht() == NULL)
			q->setRht(NULL);
		else
			q->setRht(&nodes[p->getRht()->getIndex()]);

		if (p->getAnc() == NULL)
			q->setAnc(NULL);
		else
			q->setAnc(&nodes[p->getAnc()->getIndex()]);

		q->setIndex(p->getIndex());
		q->setV(p->getV());
		q->setName(p->getName());
		q->setActiveCl(p->getActiveCl());
		q->setActiveTi(p->getActiveTi());
		q->setUpdateCl(p->getUpdateCl());
		q->setUpdateTi(p->getUpdateTi());
		q->setFlag(p->getFlag());
		q->setActiveParm(p->getActiveParm());
		q->setUpdateParm(p->getUpdateParm());
		int space = p->getActiveParm();
		q->setJumpCount(p->getJumpCount(space));
		q->setJumpSize(p->getJumpSize(space));
		q->setSumJumpSize(p->getSumJumpSize(space));

		if (p == t.root)
			root = q;

		p = t.downPassSequence[n];
		downPassSequence[n] = &nodes[p->getIndex()];
	}
}

double Topology::change(void) {

	// update conditional likelihood and transition probability flags
	updateAllCls(true);
	updateAllTis(true);
	flipAllActiveCls();
	flipAllActiveTis();

	// pick a branch at random
	Node* p = NULL;
	do {
		p = &nodes[(int)(ranPtr->uniformRv() * numNodes)];
	}while(p->getAnc() == NULL);

	// update the flags
	// while (q != NULL) {
	// 	q->flipActiveTi();
	//	q->setUpdateTi(true);
	//	q->flipActiveCl();
	//	q->setUpdateCl(true);
	//	q = q->getAnc();
	// };

	double oldV = p->getV();
	double tuning = log(4.0);
	double newV = oldV * exp(tuning * (ranPtr->uniformRv() - 0.5));
	p->setV(newV);

	return log(newV) - log(oldV);
}

void Topology::setBranchRatios(void)
{
	for (int n = 0; n < numNodes; n++)
	{
		Node *p = &nodes[n];
		if (p->getAnc() != NULL)
		{
			p->setRatio(p->getV() / treeLength);
		}
	}
}

double Topology::proposeAddJump(double lambda, double sigma, int space)
{
	int u = ranPtr->uniformRv();
	double sumLength = 0.0;
	//int branchJumpCount = 0;
	Node* p = NULL;

	// select a branch proportional to its length
	for (int n = 0; n < numNodes && sumLength < u; n++)
	{
		p = getNode(n);
		sumLength += p->getRatio();
	}

	double x = ranPtr->normalRv(0.0, sigma);
	//double lnJumpCountProb = log(ranPtr->poissonProb(lambda, branchJumpCount + 1));
	//double lnJumpSizeProb = ranPtr->lnNormalPdf(0.0, sigma, x);

	std::cout << "OK2\n";
	p->addJump(x, space);

	return 1.0 / ((numJumps + 1) * ranPtr->lnNormalPdf(0.0, sigma, x) * p->getRatio());
}

double Topology::proposeRemoveJump(double lambda, double sigma, int space)
{
	int u = ranPtr->uniformRv() * numJumps;
	int sumJumps = 0;
	int branchJumpCount = 0;
	Node* p = NULL;

	// select a branch proportionally to the number of jumps it possesses
	for (int n = 0; n < numNodes && u < sumJumps; n++)
	{
		p = getNode(n);
		branchJumpCount = p->getJumpCount(space);
		sumJumps += branchJumpCount;
	}

	//double lnJumpCountProb = log(ranPtr->poissonProb(lambda, branchJumpCount - 1));
	int n = ranPtr->uniformRv() * branchJumpCount;
	double x = p->getJumpSize(space)[n];
	p->removeJump(n, space);

	return (numJumps + 1) * ranPtr->lnNormalPdf(0.0, sigma, x) * p->getRatio();
}

void Topology::copyNodeSpaces(void)
{
	Node* p = NULL;

	for (int n = 0; n < numNodes; n++)
	{
		p = getNode(n);
		p->copySpace();
	}
}

int Topology::dex(Node *p) {

	if (p != NULL)
		return p->getIndex();
	return -1;
}

double Topology::lnProbability() {

	double lnP = 0.0;
	for (int n = 0; n < numNodes; n++) {
		Node *p = &nodes[n];
		if (p->getAnc() != NULL)
			lnP += ranPtr->lnExponentialPdf(brlenLambda, p->getV());
	}
	return lnP;
}

void Topology::updateAllCls(bool tf) {

	for (int n = 0; n < numNodes; n++) {
		Node *p = &nodes[n];
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
			p->setUpdateCl(tf);
	}
}

void Topology::updateAllTis(bool tf) {

	for (int n = 0; n < numNodes; n++) {
		Node *p = &nodes[n];
		if (p->getAnc() != NULL)
			p->setUpdateTi(tf);
	}
}

void Topology::updateAllParms(bool tf) {

	for (int n = 0; n < numNodes; n++) {
		Node *p = &nodes[n];
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
			p->setUpdateParm(tf);
	}
}

void Topology::flipAllActiveCls() {

	for (int n = 0; n < numNodes; n++) {
		Node *p = &nodes[n];
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL)
			//// TODO: disabled temporarily
			//// p->flipActiveCl();
			;
	}
}

void Topology::flipAllActiveTis() {

	for (int n = 0; n < numNodes; n++) {
		Node *p = &nodes[n];
		if (p->getAnc() != NULL)
			p->flipActiveTi();
	}
}

void Topology::flipAllActiveParms() {

	for (int n = 0; n < numNodes; n++) {
		Node *p = &nodes[n];
		if (p->getAnc() != NULL)
			p->flipActiveParm();
	}
}

void Topology::getDownPassSequence(void) {

	int i = 0;
	passDn(root, &i);
}

void Topology::passDn(Node *p, int *x) {

	if (p != NULL) {
		passDn(p->getLft(), x);
		passDn(p->getRht(), x);
		p->setDnPassIndex((*x));
		downPassSequence[(*x)++] = p;
	}
}

void Topology::print(void) {

	showNodes(root, 0);
}

void Topology::showNodes(Node *p, int indent) {

	if (p != NULL) {
		for (int i = 0; i < indent; i++)
			std::cout << " ";

		// (a,b,c) L puppy <== Root
		std::cout << std::setw(2) << dex(p) << " (";
		std::cout << std::setw(2) << dex(p->getLft()) << ",";
		std::cout << std::setw(2) << dex(p->getRht()) << ",";
		std::cout << std::setw(2) << dex(p->getAnc()) << ") ";
		std::cout << p->getPartition() << " ";
		std::cout << std::fixed << std::setprecision(5) << p->getV() << " ";
		std::cout << p->getName() << " ";
		if (p == root)
			std::cout << " <== Root ";
		std::cout << std::endl;

		showNodes(p->getLft(), indent+2);
		showNodes(p->getRht(), indent+2);
	}
}

std::string Topology::getNewick() {

	std::stringstream ss;
	writeTree(root->getLft(), ss);
	std::string newick = ss.str();
	return newick;
}

void Topology::markPathDownFromNode(Node* p) {

	for (int i = 0; i < numNodes; i++)
		nodes[i].setFlag(false);
	Node* q = p;
	while(q != NULL) {
		q->setFlag(true);
		q = q->getAnc();
	}
}

void Topology::writeTree(Node *p, std::stringstream &ss) {

	if (p != NULL) {
		if (p->getLft() == NULL)
			ss << p->getIndex() + 1 << ":" << std::fixed << std::setprecision(5) << p->getV();
		else {
			if (p->getAnc() != NULL)
				ss << "(";
			writeTree(p->getLft(), ss);
			ss << ",";
			writeTree(p->getRht(), ss);
			if (p->getAnc() != NULL) {
				if (p->getAnc()->getAnc() == NULL)
					ss << "," << p->getAnc()->getIndex() + 1 << ":" << std::fixed << std::setprecision(5) << p->getV();
				if (p->getAnc()->getAnc() != NULL)
					ss << "):" << std::fixed << std::setprecision(5) << p->getV();
				else
					ss << ")";
			}
			else {
				if (p->getAnc() == NULL)
					ss << ")";
				else
					ss << "):" << std::fixed << std::setprecision(5) << p->getV();
			}
		}
	}
}

void Topology::initializeTaxonBipartitions(void) {

	for (int n = 0; n < numNodes; n++) {
		Node *p = &nodes[n];
		MbBitfield *tipBf = new MbBitfield(numTaxa);
		p->getPartition() = *tipBf;
		delete tipBf;
		if (p->getLft() == NULL || p->getRht() == NULL || p->getAnc() == NULL)
			p->getPartition().setBit(p->getIndex());
	}
}

void Topology::getTaxonBipartitions(void) {

	for (int n = 0; n < numNodes; n++) {
		Node *p = getDownPassNode(n);
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL) {
			MbBitfield &lp = p->getLft()->getPartition();
			MbBitfield &rp = p->getRht()->getPartition();
			MbBitfield &pp = p->getPartition();
			pp = (lp | rp);
		}
	}

	for (int n = 0; n < numNodes; n++) {
		Node *p = getDownPassNode(n);
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL) {
			MbBitfield &pp = p->getPartition();
			if (pp.isBitSet(0) == true)
				pp.flipBits();
		}
	}
}

void Topology::printTaxonBipartitions(void) {

	for (int n = 0; n < numNodes; n++) {
		Node* p = &nodes[n];
		if (p->getLft() != NULL && p->getRht() != NULL && p->getAnc() != NULL) {
			MbBitfield &pp = p->getPartition();
			std::cout << p->getIndex() << " -- " << pp << std::endl;
		}
	}
}

void Topology::printJumpSamples(void) {

	std::cout << std::setw(10) << "i\tn\tlnpr_n\tx\tlnpr_x\tx_i\n";

	std::vector<double>::iterator it_x, it_p;
	for (int n = 0; n < numNodes; n++) {
		Node* p = &nodes[n];
		std::vector<double> jumpSize, lnProbJumpSize;
		if (p->getAnc() != NULL) {
			int activeParm = p->getActiveParm();
			jumpSize = p->getJumpSize(activeParm);
			lnProbJumpSize = p->getLnProbJumpSize(activeParm);
			std::cout << p->getIndex() << "\t" << p->getJumpCount(activeParm) << "\t" << p->getLnProbJumpCount(activeParm) << "\t";
			std::cout << p->getSumJumpSize(activeParm) << "\t" << p->getSumLnProbJumpSize(activeParm) << "";


			if (p->getJumpCount(activeParm) != (int)(p->getJumpSize(activeParm).size()))
				std::cerr << "ERROR: jumpCount != jumpSize.size()\n";

			if (p->getJumpCount(activeParm) != 0)
			{
				std::cout << "\t[";

				it_x = p->getJumpSize(activeParm).begin();
				it_p = p->getLnProbJumpSize(activeParm).begin();

				while (it_x != p->getJumpSize(activeParm).end() && it_p != p->getLnProbJumpSize(activeParm).end())
				{
					std::cout << " " << *it_x << "(" << *it_p << ")";
					it_x++;
					it_p++;
				}
				std::cout << " ]";
			}

			std::cout << std::endl;
		}
	}
}

void Topology::printJumpSummary(void) {

	int jumpCount = 0;
	int sumJumpCount = 0;
	double lnProbJumpCount = 0.0;
	double sumLnProbJumpCount = 0.0;

	double jumpSize = 0.0;
//	double meanJumpSize = 0.0;
//	double varJumpSize = 0.0;
	double lnProbSumJumpSize = 0.0;
	double sumLnProbSumJumpSize = 0.0;
	double sumJumpSize = 0.0;

	std::vector<double>::iterator it_x, it_p;

	// for each node
	for (int n = 0; n < numNodes; n++)
	{
		Node* p = &nodes[n];
		//std::vector<double> jumpSize;
		//std::vector<double> lnProbJumpSize;

		// exclude the root
		if (p->getAnc() != NULL)
		{
			int activeParm = p->getActiveParm();


			jumpCount = p->getJumpCount(activeParm);
			lnProbJumpCount = p->getLnProbJumpCount(activeParm);
			jumpSize = p->getSumJumpSize(activeParm);
			lnProbSumJumpSize = p->getSumLnProbJumpSize(activeParm);

			sumJumpCount += jumpCount;
			sumLnProbJumpCount += lnProbJumpCount;
			sumJumpSize += jumpSize;
			sumLnProbSumJumpSize += lnProbSumJumpSize;

			if (p->getJumpCount(activeParm) != (int)(p->getJumpSize(activeParm).size()))
			{
				std::cerr << "ERROR: jumpCount != jumpSize.size()\n";
			}
		}
	}


	// mean
	//meanJumpSize = sumJumpSize / sumJumpCount;

	// var

	std::cout << "n " << sumJumpCount << "\tpr_n " << sumLnProbJumpCount << "\tx " << sumJumpSize << "\tpr_x " << lnProbSumJumpSize << "\n";
}

void Topology::copyInactiveJumpsToActiveJumps(void)
{
	for (int n = 0; n < numNodes; n++)
	{
		Node* p = &nodes[n];
		int activeParm = p->getActiveParm();
		int inactiveParm = (activeParm == 0 ? 1 : 0);

		p->setJumpCount(p->getJumpCount(inactiveParm));
		p->setJumpSize(p->getJumpSize(inactiveParm));
		p->setLnProbJumpCount(p->getLnProbJumpCount(inactiveParm));
		p->setSumJumpSize(p->getSumJumpSize(inactiveParm));
		p->setSumLnProbJumpSize(p->getSumLnProbJumpSize(inactiveParm));
	}
}

Tree::Tree(MbRandom *rp, std::string pn, Expression *ep, double vp) : Parm(rp, pn) {

	trees[0] = new Topology(rp, ep, vp);
	trees[1] = new Topology(*trees[0]);
}

Tree::Tree(MbRandom *rp, std::string pn, Expression *ep, Settings *sp, double vp) : Parm(rp, pn) {

	trees[0] = new Topology(rp, ep, sp, vp);
	trees[1] = new Topology(*trees[0]);

}

Tree::~Tree(void) {

	delete trees[0];
	delete trees[1];
}

double Tree::change(void) {

	numAttemptedChanges++;
	return trees[activeState]->change();
}

double Tree::lnPriorRatio(void) {

	return trees[activeState]->lnProbability() - trees[getInactiveState()]->lnProbability();
}

double Tree::getValue(void)
{
	return 0.0;
}

void Tree::print(void) {

	trees[activeState]->print();
}

void Tree::keep(void) {

	*trees[getInactiveState()] = *trees[activeState];
}

void Tree::restore(void) {

	*trees[activeState] = *trees[getInactiveState()];
}

std::string Tree::getParameterStr(void) {

	std::string pStr = trees[activeState]->getNewick();
	return pStr;
}

std::string Tree::getParameterHeader(void) {

	Topology* t = trees[activeState];

	std::string pStr ="";
	pStr += "#Nexus\n\n";
	pStr += "begin trees;\n";
	pStr += "   translate\n";

	int i = 1;
	for (int n = 0; n < t->getNumNodes(); n++) {
		Node* p = t->getDownPassNode(n);
		if (p->getLft() == NULL && p->getRht() == NULL && p->getAnc() != NULL) {
			char temp[100];
			sprintf(temp, "   %d %s,\n", i++, p->getName().c_str());
			pStr += temp;
		}
		else if (p->getAnc() == NULL)
		{
			char temp[100];
			sprintf(temp, "   %d %s;\n", i++, p->getName().c_str());
			pStr += temp;
		}
	}

	return pStr;
}

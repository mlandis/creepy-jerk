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
	v			= 1.0;
	ratio		= 0.0;
	name		= "";
	activeCl	= 0;
	activeTi	= 0;
	updateCl	= false;
	updateTi	= false;
	flag		= false;

	// node values
	activeParm = 0;
	updateParm = false;
	mu[0] = 0.0;
	mu[1] = 0.0;
	sigma[0] = 0.0;
	sigma[1] = 0.0;
	k[0] = 0.0;
	k[1] = 0.0;
	kb[0] = 0.0;
	kb[1] = 0.0;
	kj[0] = 0.0;
	kj[1] = 0.0;

	// sampled jump values
	jumpCount[0] = 0;
	jumpCount[1] = 0;
	jumpSize[0].clear();
	jumpSize[1].clear();
	sumJumpSize[0] = 0.0;
	sumJumpSize[1] = 0.0;

}

void Node::addJump(double j, int space)
{
	jumpCount[space] = jumpCount[space] + 1;
	jumpSize[space].push_back(j);
	sumJumpSize[space] += j;
}

void Node::removeJump(int j, int space)
{
	double x = jumpSize[space][j];
	jumpCount[space] = jumpCount[space] - 1;
	jumpSize[space].erase(jumpSize[space].begin() + j);
	sumJumpSize[space] -= x;
}

void Node::copySpace(int i, int j)
{
	jumpCount[i] = jumpCount[j];
	jumpSize[i] = jumpSize[j];
	sumJumpSize[i] = sumJumpSize[j];
	kb[i] = kb[j];
	kj[i] = kj[j];
	k[i] = k[j];
	mu[i] = mu[j];
	sigma[i] = sigma[j];
}

Node::~Node(void) {

	//delete part;
}


Topology::Topology(MbRandom *rp, Expression *ep, double vp) {

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

//	if (sp->getPrintStdOut())
	{
		std::cout << "INITIALIZING: Topology\n";
	}

	// build a random tree
	buildRandomTree(ep);
	//initializeTaxonBipartitions();
	setBranchRatios();
	//getTaxonBipartitions();
//	if (sp->getPrintStdOut())
	{
		print();
		std::cout << "\n";
	}
}

Topology::Topology(MbRandom *rp, Expression *ep, Settings *sp, double vp) {



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

	if (sp->getPrintStdOut())
	{
		std::cout << "TOPOLOGY: Initializing from Newick string\n";
	}

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

	std::cout << "TOPOLOGY: Newick string\n" << treeStr << "\n";

	// build tree from string
	buildTreeFromNewickString(ep, treeStr);
	//initializeTaxonBipartitions();
	setBranchRatios();
	//getTaxonBipartitions();

	if (sp->getPrintStdOut())
	{
		std::cout << "TOPOLOGY: Initialization successful\n";
		print();
		std::cout << "\n";
	}
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
	//for (std::vector<std::string>::iterator p = parsedNewick.begin(); p != parsedNewick.end(); p++)
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
				/*
				if (x < 0.00001)
					x = 0.0001;
				*/
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
		q->setJumpCount(p->getJumpCount(space), space);
		q->setJumpSize(p->getJumpSize(space), space);
		q->setSumJumpSize(p->getSumJumpSize(space), space);

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
			p->setRatio( p->getV() / treeLength);
		}
	}
}


/*
double Topology::proposeAddJump(double lambda, double sigma, int space)
{
	int u = ranPtr->uniformRv();
	double sumLength = 0.0;
	Node* p = NULL;

	// select a branch proportional to its length
	for (int n = 0; n < numNodes && sumLength < u; n++)
	{
		p = getNode(n);
		sumLength += p->getRatio();
	}

	double x = ranPtr->normalRv(0.0, sigma);

	// add jump to node
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
*/

void Topology::copyNodeSpaces(int i, int j)
{
	Node* p = NULL;

	for (int n = 0; n < numNodes; n++)
	{
		p = getNode(n);
		p->copySpace(i, j);
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

Node* Topology::getRandomNode(void) {

	int u = ranPtr->uniformRv() * numNodes;
	return &nodes[u];
}

Node* Topology::getRandomNodeWithJumps(void) {

	int u = ranPtr->uniformRv() * numJumps;
	int jumpSum = 0;
	int i = 0;
	for ( ; i < numNodes; i++)
	{
		jumpSum += nodes[i].getJumpCount(0);
		if (jumpSum > u)
		{
			//std::cout << "remove jump " << i << "\n";
			return &nodes[i];
		}
	}
	std::cout << "ERROR: attempted to remove a jump when jumpCount == 0\n";
	return NULL;
}

Node* Topology::getRandomNodeNotRoot(void) {

	if (numNodes < 1) return NULL;

	Node* p = NULL;
	do
	{
		p = getRandomNode();
	}while(p->getAnc() == NULL);

	return p;
}

Node* Topology::getRandomNodeByLength(void) {

//	std::cout << "treeLength\t" << treeLength << "\n";

	double lengthSum = 0.0;
	for (int j = 0; j < numNodes; j++)
	{
		lengthSum += nodes[j].getRatio();
	}


	double u = ranPtr->uniformRv() * lengthSum;
	double l = 0.0;
	int i = 0;
	for (; i < numNodes; i++)
	{
		l += nodes[i].getRatio();
		//std::cout << l << "\t" << u << "\t" << i << "\n";
		if (l > u)
		{
			return &nodes[i];
		}
	}

	return NULL;
}

Node* Topology::getRandomNodeNotRootByLength(void) {

	if (numNodes < 1) return NULL;

	Node* p = NULL;
	do
	{
		p = getRandomNode();
	}while(p->getAnc() == NULL);

	return p;
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
	std::cout << "treeLength = " << treeLength << "\n";
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
		//std::cout << p->getPartition() << " ";
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

/*
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
*/


/*
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

*/


void Topology::printJumpSamples(int space) {

	std::cout << std::setw(10) << "i\tv\tn\tx\tx_i\n";

	std::vector<double>::iterator it_x;
	for (int n = 0; n < numNodes; n++) {
		Node* p = &nodes[n];
		std::vector<double> jumpSize;
		if (p->getAnc() != NULL) {
			jumpSize = p->getJumpSize(space);
			std::cout << p->getIndex() << "\t" << p->getV() << "\t" << p->getJumpCount(space) << "\t";
			std::cout << p->getSumJumpSize(space) << "\t";


			if (p->getJumpCount(space) != (int)(p->getJumpSize(space).size()))
				std::cerr << "ERROR: jumpCount != jumpSize.size()\n";

			if (p->getJumpCount(space) != 0)
			{
				std::cout << "\t[";

				it_x = p->getJumpSize(space).begin();

				while (it_x != p->getJumpSize(space).end())
				{
					std::cout << " " << *it_x;
					it_x++;
				}
				std::cout << " ]";
			}

			std::cout << std::endl;
		}
	}

	std::cout << "numJumps:\t" << numJumps << "\n";
}

void Topology::printJumpSizes(int space) {

	std::cout << std::setw(10) << "index\tlength\tjump\n";

	std::vector<double>::iterator it_x;
	for (int n = 0; n < numNodes; n++)
	{
		Node* p = &nodes[n];
		std::vector<double> jumpSize;
		if (p->getAnc() != NULL)
		{
			std::cout << p->getIndex() << "\t" << p->getV() << "\t" << p->getSumJumpSize(space) << "\n";
		}
	}
}

void Topology::printJumpSummary(void) {

	int jumpCount = 0;
	int sumJumpCount = 0;
	double lnProbJumpCount = 0.0;
	double sumLnProbJumpCount = 0.0;

	double jumpSize = 0.0;
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

	std::cout << "n " << sumJumpCount << "\tpr_n " << sumLnProbJumpCount << "\tx " << sumJumpSize << "\tpr_x " << lnProbSumJumpSize << "\n";
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

double Tree::lnPrior(void)
{
	return trees[activeState]->lnProbability();
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

/*
// write tree string
std::string treeStr = "";
treeStr = addNodeNhxToString(treePtr->getRoot(), treeStr);
//std::cout << "nhxStr\n" << treeStr << "\n";
nhxStrm << "tree TREE1 = " << treeStr << "\n";
nhxStrm << "End;\n";


*/


std::string Topology::getNhxString(void)
{
    std::string s = "";
    return addNodeNhxToString(root, s);
}

std::string Topology::addNodeNhxToString(Node* p, std::string s)
{
    
    if (p != NULL)
    {
        // define divergence events
        Node* lft = p->getLft();
        Node* rht = p->getRht();
        if (lft != NULL && rht != NULL)
        {
            s += "(";
            s = addNodeNhxToString(lft,s);
            s += ",";
            s = addNodeNhxToString(rht,s);
            s += ")";            
        }
        
        // define node & branch values
        std::stringstream ss;
        
        // only label tips
        if (lft == NULL && rht == NULL)
            ss << p->getName();
        
        ss << "[&j=";
        ss << std::setprecision(16) << std::fixed << p->getSumJumpSize(0);
        ss << "]";
        ss << ":" << p->getV();
        s += ss.str();
    }
    
    // string complete
    if (p == root)
        s += ";";
    
    return s;
}

std::string Topology::getNhxStringForSnr(double& boundary)
{
    std::string s = "";
    std::cout << "BOUNDARY: " << boundary << "\n";
    return addNodeNhxForSnrToString(root, s, boundary);
}

std::string Topology::addNodeNhxForSnrToString(Node* p, std::string s, double& boundary)
{
    
    if (p != NULL)
    {
        
        // define divergence events
        Node* lft = p->getLft();
        Node* rht = p->getRht();
        if (lft != NULL && rht != NULL)
        {
            s += "(";
            s = addNodeNhxForSnrToString(lft,s,boundary);
            s += ",";
            s = addNodeNhxForSnrToString(rht,s,boundary);
            s += ")";
        }
        
        
        // get signal to noise ratio for branch
        std::vector<double> storedJumps = p->getStoredJumps();

        double mean = 0.0;
        for (int i = 0; i < storedJumps.size(); i++)
            mean += storedJumps[i];
        mean /= storedJumps.size();
        
        double sd = 0.0;
        for (int i = 0; i < storedJumps.size(); i++)
            sd += pow(mean - storedJumps[i],2);
        sd /= storedJumps.size();
        sd = pow(sd,0.5);
        
        double snr = 0.0;
        if (p->getV() > 0.0 && sd > 0.0)
            snr = mean / (sd * pow(p->getV(),0.5));
        
        //std::cout << "** " << p->getName() << " " <<  p->getIndex() << " " << p->getV() << " " << mean << " " << sd  << " " << snr << " " << " " << boundary << " " << storedJumps.size() << "\n";
        
        // get largest magnitude snr to create symmetric color bar about 0.0
        if (fabs(snr) > boundary)
            boundary = fabs(snr);

        //std::cout << s << "\n";
        //std::cout << p->getName() << "\t" << p->getIndex() << "\t" << p->getV() << "\t" << mean << "\t" << sd << "\t" << snr << "\t" << boundary << "\n";
        
        // define node & branch values
        std::stringstream ss;
        
        // only label tips
        if (lft == NULL && rht == NULL)
            ss << p->getName();
        
        ss << "[&j=";
        ss << std::setprecision(16) << std::fixed << snr;
        ss << "]";
        ss << ":" << p->getV();
        s += ss.str();
        
    }
    
    // string complete
    if (p == root)
        s += ";";
    
    return s;
}

void Topology::storeAllJumps(void)
{
    for (int i = 0; i < numNodes; i++)
        nodes[i].pushStoredJump();
}

/*
 * Parm_tree.h
 *
 *  Created on: Jan 24, 2010
 *      Author: mikee
 */

#ifndef PARM_TREE_H_
#define PARM_TREE_H_

#include "Expression.h"
#include "FileMgr.h"
#include "MbBitfield.h" // Don't need? 06/01/11
#include "MbRandom.h"
#include "Parm.h"

#include <algorithm>
#include <stdio.h>
#include <sstream>
#include <string>
#include <vector>

class Expression;
class FileMgr;
class MbRandom;
class Parm;

class Node {

public:
							Node(void);
							~Node(void);
	void					flipActiveCl(void)		{ activeCl == 0 ? activeCl = 1 : activeCl = 0; }
	void					flipActiveTi(void)		{ activeTi == 0 ? activeTi = 1 : activeTi = 0; }
	void					flipActiveParm(void)	{ activeParm == 0 ? activeParm = 1 : activeParm = 0; }
	Node*					getLft(void)			{ return lft; }
	Node*					getRht(void)			{ return rht; }
	Node*					getAnc(void)			{ return anc; }
	double					getV(void)				{ return v; }
	std::string				getName(void)			{ return name; }
	int						getIndex(void)			{ return index; }
	int						getDnPassIndex(void)	{ return dnPassIndex; }
	int						getActiveCl(void)		{ return activeCl; }
	int						getActiveTi(void)		{ return activeTi; }
	bool					getUpdateCl(void)		{ return updateCl; }
	bool					getUpdateTi(void)		{ return updateTi; }
	bool					getFlag(void)			{ return flag; }
	void					setLft(Node *p)			{ lft = p; }
	void					setRht(Node *p)			{ rht = p; }
	void					setAnc(Node *p)			{ anc = p; }
	void					setV(double x)			{ v = x; }
	void					setRatio(double x)		{ ratio = x; }
	void					setQ(int x)				{ q = x; }
	void					setName(std::string s)	{ name = s; }
	void					setIndex(int x)			{ index = x; }
	void					setDnPassIndex(int x)	{ dnPassIndex = x; }
	void					setActiveCl(int x)		{ activeCl = x; }
	void					setActiveTi(int x)		{ activeTi = x; }
	void					setUpdateCl(bool tf)	{ updateCl = tf; }
	void					setUpdateTi(bool tf)	{ updateTi = tf; }
	void					setFlag(bool tf)		{ flag = tf; }
	MbBitfield&				getPartition(void)		{ return *part; }


	//
	void					setMu(double x)			{ mu = x; }
	void					setSigma(double x)		{ sigma = x; }
	void					setK(double x)			{ k = x; }
	void					setKb(double x)			{ kb = x; }
	void					setKj(double x)			{ kj = x; }
	void					addK(double x)			{ k += x; }
	void					addMu(double x)			{ mu += x; }
	double					getMu(void)				{ return mu; }
	double					getSigma(void)			{ return sigma; }
	double					getK(void)				{ return k; }
	double					getKb(void)				{ return kb; }
	double					getKj(void)				{ return kj; }
	double					getRatio(void)			{ return ratio; }
	void					like(void);

	void					addJump(double j, int space);
	void					removeJump(int j, int space);
	void					copySpace(int i, int j);

	int						getActiveParm(void)		{ return activeParm; }
	bool					getUpdateParm(void)		{ return updateParm; }
	int						getJumpCount(int x)		{ return jumpCount[x]; }
	std::vector<double>&	getJumpSize(int x)		{ return jumpSize[x]; }
	int						getJumpCount(void)		{ return jumpCount[activeParm]; }
	std::vector<double>&	getJumpSize(void)		{ return jumpSize[activeParm]; }
	double					getSumJumpSize(void)	{ return sumJumpSize[activeParm]; }
	double					getLnProbJumpCount(int x)	{ return lnProbJumpCount[x]; }
	std::vector<double>&	getLnProbJumpSize(int x)	{ return lnProbJumpSize[x]; }
	double					getSumJumpSize(int x)	{ return sumJumpSize[x]; }
	double					getSumLnProbJumpSize(int x)	{ return sumLnProbJumpSize[x]; }


	void					setActiveParm(int x)	{ activeParm = x; }
	void					setUpdateParm(bool tf)	{ updateParm = tf;}
	void					setJumpCount(int x)	{ jumpCount[1] = x; }
	void					setLnProbJumpCount(double x) { lnProbJumpCount[1] = x; }
	void					setJumpSize(std::vector<double>& x)	{ jumpSize[1] = x; }
	void					setLnProbJumpSize(std::vector<double>& x) { lnProbJumpSize[1] = x; }
	void					setSumJumpSize(double x) { sumJumpSize[1] = x; }
	void					setSumLnProbJumpSize(double x) { sumLnProbJumpSize[1] = x; }




private:
	Node*					lft;
	Node*					rht;
	Node*					anc;
	int						index;
	int						dnPassIndex;
	int						q;						// number of changes
	double					v;						// branch length
	double					ratio;
	std::string				name;
	int						activeCl;
	int						activeTi;
	bool					updateCl;
	bool					updateTi;
	bool					flag;
	MbBitfield*				part;

	// sampling (non-FFT) jump diffusion
	int						activeParm;
	bool					updateParm;
	double					mu;
	double					sigma;
	double					k;						// likelihood at node
	double					kb;						// Brownian motion likelihood component
	double					kj;						// jump kernel likelihood component
	int						jumpCount[2];
	double					lnProbJumpCount[2];
	std::vector<double>		jumpSize[2];
	std::vector<double>		lnProbJumpSize[2];
	double					sumJumpSize[2];
	double					sumLnProbJumpSize[2];

};


class Topology {

public:
					Topology(MbRandom *rp, Expression *ep, double vp);
					// Topology(MbRandom *rp, Expression *ep, std::string treeStr, double vp);
					Topology(MbRandom *rp, Expression *ep, Settings *sp, double vp);
					Topology(Topology &t);
					~Topology(void);
	Topology		&operator=(const Topology &t);
	double			change(void);
	Node*			getDownPassNode(int i)			{ return downPassSequence[i]; }
	Node*			getNode(int i)					{ return &nodes[i]; }
	Node*			getRoot(void)					{ return root; }
	void			getDownPassSequence(void);
	int				getNumNodes(void)				{ return numNodes; }
	std::string		getNewick(void);
	double			getTreeLength(void)				{ return treeLength; }
	int				getNumJumps(void)				{ return numJumps; }
	Node*			getRandomNode(void);
	Node*			getRandomNodeWithJumps(void);
	Node*			getRandomNodeNotRoot(void);
	Node*			getRandomNodeByLength(void);
	Node*			getRandomNodeNotRootByLength(void);
	void			setNumJumps(int x, int space)				{ numJumps = x; }
	void			incrementNumJumps(void)				{ numJumps++; }
	void			decrementNumJumps(void)				{ numJumps--; }

	void			setBranchRatios(void);
	//double			proposeAddJump(double lambda, double sigma, int space);
	//double			proposeRemoveJump(double lambda, double sigma, int space);
	void			copyNodeSpaces(int i, int j);

	double			lnProbability(void);
	void			print(void);
	void			updateAllCls(bool tf);
	void			updateAllTis(bool tf);
	void			updateAllParms(bool tf);
	void			flipAllActiveCls(void);
	void			flipAllActiveTis(void);
	void			flipAllActiveParms(void);
	void			copyInactiveJumpsToActiveJumps(void);
	void			initializeTaxonBipartitions(void);
	void			getTaxonBipartitions(void);
	void			printTaxonBipartitions(void);
	void			printJumpSamples(int space);
	void			printJumpSummary(void);
	void			printJumpSizes(int space);
	void			markPathDownFromNode(Node* p);

private:
	void			buildRandomTree(Expression *ep);
	void			buildTreeFromNewickString(Expression *ep, std::string &ts);
	void			clone(const Topology &t);
	int				dex(Node *p);
	void			passDn(Node *p, int *x);
	void			showNodes(Node *p, int indent);
	void			writeTree(Node *p, std::stringstream &ss);
	int				numTaxa;
	int				numNodes;
	Node			*nodes;
	MbRandom		*ranPtr;
	Settings		*settingsPtr;
	Node			*root;
	Node			**downPassSequence;
	double			brlenLambda;
	double			treeLength;
	int				numJumps;

};


class Tree : public Parm  {

public:
					Tree(MbRandom *rp, std::string pn, Expression *ep, double vp);
					Tree(MbRandom *rp, std::string pn, Expression *ep, Settings *sp, double vp);
					~Tree(void);
	Topology*		getActiveTopology(void)			{ return trees[activeState]; }
	double			lnPriorRatio(void);
	double			lnPrior(void);
	double			change(void);
	double			getValue(void);
	void			print(void);
	void			keep(void);
	void			restore(void);
	std::string		getParameterStr(void);
	std::string		getParameterHeader(void);

private:
	Topology		*trees[2];

};

#endif /* PARM_TREE_H_ */
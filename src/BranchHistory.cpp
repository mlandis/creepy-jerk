#include "BranchHistory.h"

BranchHistory::BranchHistory(Node* n, int na) {
    
    nodePtr = n;
    numAreas = na;

    ancestralState.resize(numAreas);
    for (int i=0; i<ancestralState.size(); i++)
        ancestralState[i] = false;

    ancestralStateCounts.resize(numAreas);
    for (int i=0; i<ancestralStateCounts.size(); i++)
        ancestralStateCounts[i] = 0;
}

BranchHistory::~BranchHistory(void) {
    
    removeAllModelEvents();
}

BranchHistory& BranchHistory::operator=(const BranchHistory& h) {
    
    if (this != &h)
        {
        numAreas = h.numAreas;
        nodePtr = h.nodePtr;
        ancestralState = h.ancestralState;
        
        for (std::set<ModelEvent*,comp_history>::iterator it = branchChanges.begin(); it != branchChanges.end(); it++)
            delete (*it);
        branchChanges.clear();

        for (std::set<ModelEvent*,comp_history>::iterator it = h.branchChanges.begin(); it != h.branchChanges.end(); it++)
            {
            ModelEvent* newChange = new ModelEvent( *(*it) );
            branchChanges.insert(newChange);
            }
        }
    return (*this);
}

/*

// count model and rate switching events

int BranchHistory::getNumLosses(void) {
    
    std::vector<bool> s = getAncestralState();
    int numLosses = 0;
    for (std::set<ModelEvent*,comp_history>::iterator it = branchChanges.begin(); it != branchChanges.end(); it++)
        {
        int a = (*it)->getArea();
        if (s[a] == true)
            numLosses++;
        s[a].flip();
        }
    return numLosses;
}
*/

void BranchHistory::print(void) {
    
    std::cout << "Changes along branch " << nodePtr->getIndex() << ":" << std::endl;    
    int n = 0;
    std::vector<bool> s = getAncestralState();
 
    std::cout << std::setw(6) << n << "       " << std::fixed << std::setprecision(6) << 0.0 << " -- ";
    for (int i=0; i<s.size(); i++)
        std::cout << s[i];
    std::cout << std::endl;

    for (std::set<ModelEvent*,comp_history>::iterator it = branchChanges.begin(); it != branchChanges.end(); it++)
        {
        n++;
        int a = (*it)->getModelType();
        double pos = (*it)->getPosition();
        s[a].flip();
        std::cout << std::setw(6) << n << std::setw(6) << a << " " << std::fixed << std::setprecision(6) << pos << " -- ";
        for (int i=0; i<s.size(); i++)
            std::cout << s[i];
        std::cout << std::endl;
        }
}

void BranchHistory::removeChanges(std::set<int>& areaSet) {
    
    std::set<ModelEvent*,comp_history>::iterator it = branchChanges.begin();
    while (it != branchChanges.end())
        {
        //std::set<ModelEvent*,comp_history>::iterator tmp = it++;
        
        std::set<int>::iterator p = areaSet.find( (*it)->getModelType() );
        if (p != areaSet.end())
            {
            //std::cout << "removing area change in area " << (*it)->getArea() << " from branch " << nodePtr->getIndex() << std::endl;
            delete (*it);
            branchChanges.erase(it);
            it = branchChanges.begin();
            }
        else
            it++;
        }
}

void BranchHistory::removeAllModelEvents(void) {

    for (std::set<ModelEvent*>::iterator it = branchChanges.begin(); it != branchChanges.end(); it++)
        delete (*it);
    branchChanges.clear();
}

void BranchHistory::setAncestralArea(int a, int s) {

    if (s == 0)
        ancestralState[a] = false;
    else
        ancestralState[a] = true;
}

void BranchHistory::updateAncestralStateCounts(void) {
    for (int i = 0; i < numAreas; i++)
        if (ancestralState[i] == true)
            ancestralStateCounts[i] += 1;
}

#ifndef BranchHistory_H
#define BranchHistory_H

#include <string>
#include <set>
#include <vector>
#include "ModelEvent.h"
#include "Parm_tree.h"
#include <iomanip>
#include <iostream>

class Node;

class BranchHistory {
    
    public:
                                            BranchHistory(Node* n, int na);  
                                           ~BranchHistory(void);
        BranchHistory&                      operator=(const BranchHistory& h);
        void                                addModelEvent(ModelEvent* a) { branchChanges.insert(a); }
        std::vector<bool>                   getAncestralState(void) { return ancestralState; }
        std::set<ModelEvent*,comp_history>& getBranchChanges(void) { return branchChanges; }
        bool                                getAncestralStateForArea(int a) { return ancestralState[a]; }
        int                                 getNumChanges(void) { return (int)branchChanges.size(); }
        int                                 getNumGains(void);
        int                                 getNumLosses(void);
        void                                print(void);
        void                                removeAllModelEvents(void);
        void                                removeChanges(std::set<int>& areaSet);
        void                                setAncestralArea(int a, int s);
        void                                updateAncestralStateCounts(void);
        int                                 getAncestralStateCountForArea(int a) { return ancestralStateCounts[a]; }
        
    private:
        int                                 numAreas;
        std::set<ModelEvent*,comp_history>  branchChanges;
        Node*                               nodePtr;
        std::vector<bool>                   ancestralState;
        std::vector<int>                    ancestralStateCounts;
};

#endif

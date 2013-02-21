#ifndef MODELEVENT_H
#define MODELEVENT_H

#include <iomanip>
#include <iostream>
#include <vector>

class ModelEvent {
    
    public:
                                    ModelEvent(int a, double p);
                                    ModelEvent(ModelEvent& a);
                                   ~ModelEvent(void);
        bool                        operator<(const ModelEvent& a) const;
        int                         getModelType(void) { return modelType; }
        double                      getPosition(void) { return position; }
        std::vector<std::vector<double> > getParameterMultipliers(void) { return parameterMultipliers; }
        void                        print(void);
    
    private:
        int                         modelType;
        double                      position;
        std::vector<std::vector<double> > parameterMultipliers;
};

class comp_history {
    
    public:
        bool                        operator()(ModelEvent* m1, ModelEvent* m2) const { return (*m1 < *m2); }
};

#endif

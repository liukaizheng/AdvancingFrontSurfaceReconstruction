#ifndef GLPKSOLVER_H
#define GLPKSOLVER_H

#include "LinearVariable.h"
#include "LinearConstraint.h"
#include <vector>

class GLPKSolver
{
public:
    GLPKSolver(){}
    ~GLPKSolver();

    void CreateVariables(int n, LinearVariable::Type type);
    void CreateConstraints(int n);
    void AddConstraints(LinearConstraint* c);
    void CrateObjective(bool isMin);
    bool Solve();
    LinearVariable* GetVariable(int i);
    LinearConstraint* GetConstraint(int i);
    LinearConstraint* GetObjective();
    const double* GetResult() const;

protected:
    int BoundType(double lb, double ub);

private:
    std::vector<LinearVariable*> mVariableList;
    std::vector<LinearConstraint*> mConstraintList;
    LinearConstraint* mObjective;
    bool mIsMin;
    std::vector<double> mResult;
};

#endif

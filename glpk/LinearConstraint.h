#ifndef LINEARCONSTRAINT_H
#define LINEARCONSTRAINT_H

#include <unordered_map>

class LinearConstraint
{
public:
    LinearConstraint(double lb = -std::numeric_limits<double>::infinity(), double ub = std::numeric_limits<double>::infinity()):mLb(lb), mUb(ub) {}
    ~LinearConstraint() {}

    void Insert(int i, double val);
    void GetBounds(double& lb, double& ub) const;
    void SetUpBound(const double& ub);
    void SetLowBound(const double& lb);
    std::unordered_map<int, double>& GetCoeffs();

private:
    double mLb;
    double mUb;
    std::unordered_map<int, double> mCoeffs; 
    std::string mName;
};

#endif

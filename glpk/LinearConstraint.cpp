#include "LinearConstraint.h"

void LinearConstraint::Insert(int i, double val)
{
    mCoeffs.emplace(std::make_pair(i, val));
}

void LinearConstraint::GetBounds(double& lb, double& ub) const
{
    lb = mLb;
    ub = mUb;
}

std::unordered_map<int, double>& LinearConstraint::GetCoeffs()
{
    return mCoeffs;
}

void LinearConstraint::SetUpBound(const double& ub)
{
    mUb = ub;
}
void LinearConstraint::SetLowBound(const double& lb)
{
    mLb = lb;
}

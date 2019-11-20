#include "LinearVariable.h"
#include <sstream>
LinearVariable::LinearVariable(int idx, Type type, double lb, double ub)
    :mIdx(idx), mType(type), mLb(lb), mUb(ub) 
{
    if(type == BINARY) 
    { 
        mLb = 0; 
        mUb = 1;
    }
    std::stringstream ss;
    ss << "x" << idx;
    mName = ss.str();
}

int LinearVariable::GetIndex() const 
{
    return mIdx;
}

LinearVariable::Type LinearVariable::GetType() const
{
    return mType;
}

void LinearVariable::GetBounds(double& lb, double &ub) const
{
    lb = mLb;
    ub = mUb;
}

void LinearVariable::SetBounds(const double& lb, const double& ub)
{
    mLb = lb;
    mUb = ub;
}

void LinearVariable::SetUpBound(const double& ub)
{
    mUb = ub;
}

void LinearVariable::SetLowBound(const double& lb)
{
    mLb = lb;
}

void LinearVariable::SetIndex(int idx)
{
    mIdx = idx;
}

void LinearVariable::SetType(const LinearVariable::Type& type) 
{
    mType = type;
    if(type == BINARY)
    {
        mLb = 0.0;
        mUb = 1.0;
    }
}

std::string LinearVariable::GetName() const
{
    return mName;
}

void LinearVariable::SetName(const std::string& name)
{
    mName = name;
}

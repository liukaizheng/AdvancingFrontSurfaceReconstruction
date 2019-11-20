#ifndef LINEARVARIABLE_H
#define LINEARVARIABLE_H

#include <limits>
#include <string>

class LinearVariable
{
public:
    enum Type
    {
        BINARY,
        INTERGER,
        CONTINUOUS
    };
    LinearVariable(int idx, Type type, double lb = -std::numeric_limits<double>::infinity(), double ub = std::numeric_limits<double>::infinity());
    ~LinearVariable() {}

    void GetBounds(double& lb, double &ub) const;
    void SetBounds(const double& lb, const double& ub);
    void SetUpBound(const double& ub);
    void SetLowBound(const double& lb);
    int GetIndex() const;
    void SetIndex(int idx); 
    Type GetType() const;
    void SetType(const Type& type);
    std::string GetName() const;
    void SetName(const std::string& name);

private:
    int mIdx;
    Type mType;
    double mLb;
    double mUb;
    std::string mName;
};
#endif

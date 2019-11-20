#include "GLPKSolver.h"
#include "glpk.h"
#include <sstream>

GLPKSolver::~GLPKSolver()
{
    for(unsigned i = 0; i < mVariableList.size(); i++)
        delete mVariableList[i];

    for(unsigned i = 0; i < mConstraintList.size(); i++)
        delete mConstraintList[i];
}

void GLPKSolver::CreateVariables(int n, LinearVariable::Type type)
{
    mVariableList.reserve(n);
    for(int i = 0; i < n; i++)
    {
        mVariableList.emplace_back(new LinearVariable(i, type));
    }
}
void GLPKSolver::CreateConstraints(int n)
{
    mConstraintList.reserve(n);
    for(int i = 0; i < n; i++)
    {
        mConstraintList.emplace_back(new LinearConstraint);
    }
}

void GLPKSolver::CrateObjective(bool isMin)
{
    mIsMin = isMin;
    mObjective = new LinearConstraint;
}

LinearVariable* GLPKSolver::GetVariable(int i)
{
    return mVariableList[i];
}

LinearConstraint* GLPKSolver::GetConstraint(int i)
{
    return mConstraintList[i];
}

void GLPKSolver::AddConstraints(LinearConstraint* c)
{
    mConstraintList.emplace_back(c);
}

LinearConstraint* GLPKSolver::GetObjective()
{
    return mObjective;
}

const double* GLPKSolver::GetResult() const
{
    return mResult.data();
}

int GLPKSolver::BoundType(double lb, double ub)
{

    if (lb <= -std::numeric_limits<double>::infinity() && ub >= std::numeric_limits<double>::infinity())
        return GLP_FR;		// free (unbounded) variable

    else if (lb > -std::numeric_limits<double>::infinity() && ub >= std::numeric_limits<double>::infinity())
        return GLP_LO;		// variable with lower bound

    else if (lb <= -std::numeric_limits<double>::infinity() && ub < std::numeric_limits<double>::infinity())
        return GLP_UP;		// variable with upper bound

    else {// lb > -std::numeric_limits<double>::infinity() && ub < std::numeric_limits<double>::infinity()
        if (lb == ub)
            return GLP_FX;	// fixed variable
        else
            return GLP_DB;  // double-bounded variable
    }
}

bool GLPKSolver::Solve()
{
    glp_prob* lp = glp_create_prob();
    if(!lp) return false;

    int num_integer_variables = 0;

    int num_variables = static_cast<int>(mVariableList.size());
    glp_add_cols(lp, num_variables);

    for(int i = 0; i < num_variables; i++)
    {
        const LinearVariable* var = mVariableList[i];
        glp_set_col_name(lp, i + 1, var->GetName().data());

        if(var->GetType() == LinearVariable::INTERGER)
        {
            glp_set_col_kind(lp, i + 1, GLP_IV);
            num_integer_variables++;
        }
        else if(var->GetType() == LinearVariable::BINARY)
        {
            glp_set_col_kind(lp, i + 1, GLP_BV);
            num_integer_variables++;
        }
        else
            glp_set_col_kind(lp, i + 1, GLP_CV);

        double lb, ub;
        var->GetBounds(lb, ub);
        int type = BoundType(lb, ub);
        glp_set_col_bnds(lp, i + 1, type, lb, ub);
    }

    int num_constraints = static_cast<int>(mConstraintList.size());
    glp_add_rows(lp, num_constraints);
    for(int i = 0; i < num_constraints; i++)
    {
        LinearConstraint* c = mConstraintList[i];
        auto coeffs = c->GetCoeffs();
        std::vector<int> indices(coeffs.size() + 1, 0);
        std::vector<double> coefficients(coeffs.size() + 1, 0.0);
        int idx = 1;

        for(auto it = coeffs.begin(); it != coeffs.end(); it++)
        {
            int var_idx = it->first;
            double coeff = it->second;
            indices[idx] = var_idx + 1;
            coefficients[idx] = coeff;
            idx++;
        }

        glp_set_mat_row(lp, i + 1, static_cast<int>(coeffs.size()), indices.data(), coefficients.data());
        double lb, ub;
        c->GetBounds(lb, ub);
        int type = BoundType(lb, ub);
        glp_set_row_bnds(lp, i + 1, type, lb, ub);

        std::stringstream ss;
        ss << "c" << i;
        glp_set_row_name(lp, i + 1, ss.str().data());
    }
    auto objCoeffs = mObjective->GetCoeffs();

    for(auto it = objCoeffs.begin(); it != objCoeffs.end(); it++)
    {
        int var_idx = it->first;
        double coeff = it->second;
        glp_set_obj_coef(lp, var_idx + 1, coeff);
    }

    glp_set_obj_dir(lp, mIsMin ? GLP_MIN : GLP_MAX);

    int status = -1;
    if(num_integer_variables == 0)
    {
        glp_smcp parm;
        glp_init_smcp(&parm);
        parm.msg_lev = GLP_MSG_ALL;
        status = glp_simplex(lp, &parm);
    }
    else
    {
        glp_iocp parm;
        glp_init_iocp(&parm);
        parm.presolve = GLP_ON;
        status = glp_intopt(lp, &parm);
    }

    switch(status)
    {
        case 0:
            {
                mResult.resize(num_variables);
                for(int i = 0; i < num_variables; i++)
                {
                    if(num_integer_variables == 0)
                        mResult[i] = glp_get_col_prim(lp, i + 1);
                    else
                        mResult[i] = glp_mip_col_val(lp, i + 1);
                }
                break;
            }

    }
    glp_delete_prob(lp);

    return status == 0;
}

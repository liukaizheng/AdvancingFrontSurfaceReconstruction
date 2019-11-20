#include "DelaunayRadius.h"
#include "PreTypeDefine.h"

template <typename DerivedV, typename DerivedT, typename DerivedC, typename DerivedR, typename TypeuF2F>
void DelaunayRadius(
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedT>& T,
        const Eigen::MatrixBase<DerivedC>& C,
        const Eigen::MatrixBase<DerivedR>& TR,
        const std::vector<std::vector<TypeuF2F>>& uF2F,
        const std::vector<bool>& computable,
        Eigen::PlainObjectBase<DerivedR>& R)
{
    typedef typename DerivedT::Scalar Index;
    typedef Eigen::Matrix<typename DerivedV::Scalar, -1, 1> VectorXd;
    const Index numTets = static_cast<Index>(T.rows());
    R.resize(TR.size());

    for(Index i = 0; i < R.size(); i++)
    {
#ifndef NDebug
        assert(uF2F[i].size() > 0 && uF2F[i].size() < 3);
#endif
        if(uF2F[i].size() == 1)
        {
            R[i] = TR[i];
        }
        else
        {
            const auto& cf = uF2F[i][0];
            const auto& nf = uF2F[i][1];
            const auto ct = cf % numTets;
            const auto nt = nf % numTets;
            if(!computable[ct] || !computable[nt])
            {
                R[i] = TR[i];
            }
            else
            {
                const VectorXd cc = C.row(ct);
                const VectorXd nc = C.row(nt);

                const VectorXd p = V.row(T(ct, (cf / numTets + 1) % 4));

                const VectorXd v0 = cc - nc, v1 = cc - p, v2 = p - nc; 
                const auto d1 = v0.dot(v1), d2 = v0.dot(v2);
                if(d1 > 0 && d2 > 0)
                {
                    R[i] = TR[i];
                }
                else
                {
                    if(d1 < 0)
                        R[i] = (cc - p).squaredNorm();
                    else
                        R[i] = (nc - p).squaredNorm();
                }

            }
        }
    }
}

template 
void DelaunayRadius<RowMatrixXd, RowMatrixXi, RowMatrixXd, Eigen::VectorXd, int>(
        const Eigen::MatrixBase<RowMatrixXd>& V,
        const Eigen::MatrixBase<RowMatrixXi>& T,
        const Eigen::MatrixBase<RowMatrixXd>& C,
        const Eigen::MatrixBase<Eigen::VectorXd>& TR,
        const std::vector<std::vector<int>>& uF2F,
        const std::vector<bool>& computable,
        Eigen::PlainObjectBase<Eigen::VectorXd>& R);

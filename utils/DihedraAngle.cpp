#include "DihedraAngle.h"
#include "PreTypeDefine.h"

template <typename DerivedV, typename DerivedF>
bool CosDihedralAngle(
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedF>& F,
        const typename DerivedF::Scalar& f1,
        const typename DerivedF::Scalar& f2,
        const bool reverseNormal, 
        typename DerivedV::Scalar& cosA)
{
    typedef Eigen::Matrix<typename DerivedV::Scalar, 1, 3> RowVector;
    RowVector v11 = V.row(F(f1, 1)) - V.row(F(f1, 0));
    RowVector v12 = V.row(F(f1, 2)) - V.row(F(f1, 0));
    auto v1 = v11.cross(v12);
    if(v1.squaredNorm() < std::numeric_limits<typename DerivedV::Scalar>::epsilon())
        return false;

    RowVector v21 = V.row(F(f2, 1)) - V.row(F(f2, 0));
    RowVector v22 = V.row(F(f2, 2)) - V.row(F(f2, 0));
    auto v2 = v21.cross(v22);
    if(v2.squaredNorm() < std::numeric_limits<typename DerivedV::Scalar>::epsilon())
        return false;

    v1.normalize();
    v2.normalize();

    cosA = v1.dot(v2);
    if(reverseNormal)
        cosA = -cosA;
    return true;
}

template bool CosDihedralAngle(const Eigen::MatrixBase<RowMatrixXd>& V,
                               const Eigen::MatrixBase<RowMatrixXi>& F,
                               const int& f1, const int& f2,
                               const bool reverseNormal, double& cosA);

#include "TriangleCircleRadius.h"
#include "PreTypeDefine.h"
#include <limits>

template <typename DerivedV, typename DerivedF, typename DerivedC>
void TriangleCircleCenter(
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedF>& F,
        Eigen::PlainObjectBase<DerivedC>& C,
        std::vector<bool>& computable)
{
    typedef typename DerivedV::Scalar Scalar;
    C.resize(F.rows(), 3);
    computable.resize(F.rows(), true);
    for(Eigen::Index i = 0; i < F.rows(); i++)
    {

        Scalar a1, b1, c1, d1;
        Scalar a2, b2, c2, d2;
        Scalar a3, b3, c3, d3;

	    const Scalar &x1 = V(F(i, 0), 0), &y1 = V(F(i, 0), 1), &z1 = V(F(i, 0), 2);
	    const Scalar &x2 = V(F(i, 1), 0), &y2 = V(F(i, 1), 1), &z2 = V(F(i, 1), 2);
	    const Scalar &x3 = V(F(i, 2), 0), &y3 = V(F(i, 2), 1), &z3 = V(F(i, 2), 2);

        a1 = (y1*z2 - y2*z1 - y1*z3 + y3*z1 + y2*z3 - y3*z2);
        b1 = -(x1*z2 - x2*z1 - x1*z3 + x3*z1 + x2*z3 - x3*z2);
        c1 = (x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
        d1 = -(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);

        a2 = 2 * (x2 - x1);
        b2 = 2 * (y2 - y1);
        c2 = 2 * (z2 - z1);
        d2 = x1 * x1 + y1 * y1 + z1 * z1 - x2 * x2 - y2 * y2 - z2 * z2;

        a3 = 2 * (x3 - x1);
        b3 = 2 * (y3 - y1);
        c3 = 2 * (z3 - z1);
        d3 = x1 * x1 + y1 * y1 + z1 * z1 - x3 * x3 - y3 * y3 - z3 * z3;

        Scalar det = a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1;
        if(std::abs(det) > std::numeric_limits<Scalar>::epsilon())
        {
            C(i, 0) = -(b1*c2*d3 - b1*c3*d2 - b2*c1*d3 + b2*c3*d1 + b3*c1*d2 - b3*c2*d1) / det;
            C(i, 1) =  (a1*c2*d3 - a1*c3*d2 - a2*c1*d3 + a2*c3*d1 + a3*c1*d2 - a3*c2*d1) / det;
            C(i, 2) = -(a1*b2*d3 - a1*b3*d2 - a2*b1*d3 + a2*b3*d1 + a3*b1*d2 - a3*b2*d1) / det;
        }
        else 
        {
            computable[i] = false;
        }
    }
}

template <typename DerivedV, typename DerivedF, typename DerivedR>
void SquaredTriangleCircleRadius(
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedF>& F,
        Eigen::PlainObjectBase<DerivedR>& R)
{
    std::vector<bool> computable;
    DerivedV C;
    TriangleCircleCenter(V, F, C, computable);
    R.resize(F.rows());
    for(Eigen::Index i = 0; i < F.rows(); i++)
    {
        if(computable[i])
        {
            R(i) = (V.row(F(i, 0)) - C.row(i)).squaredNorm();
        }
        else
            R(i) = static_cast<typename DerivedR::Scalar>(0);
    }
}

template void TriangleCircleCenter<RowMatrixXd, RowMatrixXi, RowMatrixXd>(
    const Eigen::MatrixBase<RowMatrixXd>& V,
    const Eigen::MatrixBase<RowMatrixXi>& F,
    Eigen::PlainObjectBase<RowMatrixXd>& C, std::vector<bool>& computable);

template void
SquaredTriangleCircleRadius<RowMatrixXd, RowMatrixXi, Eigen::VectorXd>(
    const Eigen::MatrixBase<RowMatrixXd>& V,
    const Eigen::MatrixBase<RowMatrixXi>& F,
    Eigen::PlainObjectBase<Eigen::VectorXd>& R);

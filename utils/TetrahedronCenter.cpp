#include "TetrahedronCenter.h"
#include "PreTypeDefine.h"

template <typename DerivedV, typename DerivedT, typename DerivedC>
void TetrahedronCenter(const Eigen::MatrixBase<DerivedV>& V,
                       const Eigen::MatrixBase<DerivedT>& T,
                       Eigen::PlainObjectBase<DerivedC>& C,
                       std::vector<bool>& computable) {
    typedef typename DerivedV::Scalar Scalar;
    C.resize(T.rows(), 3);
    computable.resize(T.rows(), true);

    for (Eigen::Index i = 0; i < T.rows(); i++) {
        const Scalar &x0 = V(T(i, 0), 0), &y0 = V(T(i, 0), 1),
                     &z0 = V(T(i, 0), 2);
        const Scalar &x1 = V(T(i, 1), 0), &y1 = V(T(i, 1), 1),
                     &z1 = V(T(i, 1), 2);
        const Scalar &x2 = V(T(i, 2), 0), &y2 = V(T(i, 2), 1),
                     &z2 = V(T(i, 2), 2);
        const Scalar &x3 = V(T(i, 3), 0), &y3 = V(T(i, 3), 1),
                     &z3 = V(T(i, 3), 2);

        const Scalar d =
            (x1 * y2 * z0 - x1 * y3 * z0 - x0 * y2 * z1 + x0 * y3 * z1 -
             x1 * y0 * z2 + x0 * y1 * z2 - x0 * y3 * z2 + x1 * y3 * z2 +
             x3 * (y1 * z0 - y2 * z0 - y0 * z1 + y2 * z1 + y0 * z2 - y1 * z2) +
             x1 * y0 * z3 - x0 * y1 * z3 + x0 * y2 * z3 - x1 * y2 * z3 +
             x2 * (y3 * z0 + y0 * z1 - y3 * z1 - y0 * z3 + y1 * (-z0 + z3))) *
            2;
        if (std::abs(d) < std::numeric_limits<Scalar>::epsilon()) {
            computable[i] = false;
            continue;
        }
        C(i, 0) =
            ((-y0 + y3) *
                 ((x0 * x0 - x1 * x1 + y0 * y0 - y1 * y1 + z0 * z0 - z1 * z1) *
                      (z0 - z2) -
                  (z0 - z1) * (x0 * x0 - x2 * x2 + y0 * y0 - y2 * y2 + z0 * z0 -
                               z2 * z2)) +
             ((y0 - y2) *
                  (x0 * x0 - x1 * x1 + y0 * y0 - y1 * y1 + z0 * z0 - z1 * z1) -
              (y0 - y1) *
                  (x0 * x0 - x2 * x2 + y0 * y0 - y2 * y2 + z0 * z0 - z2 * z2)) *
                 (z0 - z3) +
             (y2 * (z0 - z1) + y0 * (z1 - z2) + y1 * (-z0 + z2)) *
                 (x0 * x0 - x3 * x3 + y0 * y0 - y3 * y3 + z0 * z0 - z3 * z3)) / d;
        C(i, 1) =
            ((x0 - x3) *
                 ((x0 * x0 - x1 * x1 + y0 * y0 - y1 * y1 + z0 * z0 - z1 * z1) *
                      (z0 - z2) -
                  (z0 - z1) * (x0 * x0 - x2 * x2 + y0 * y0 - y2 * y2 + z0 * z0 -
                               z2 * z2)) +
             ((-x0 + x2) *
                  (x0 * x0 - x1 * x1 + y0 * y0 - y1 * y1 + z0 * z0 - z1 * z1) +
              (x0 - x1) *
                  (x0 * x0 - x2 * x2 + y0 * y0 - y2 * y2 + z0 * z0 - z2 * z2)) *
                 (z0 - z3) +
             (x2 * (-z0 + z1) + x1 * (z0 - z2) + x0 * (-z1 + z2)) *
                 (x0 * x0 - x3 * x3 + y0 * y0 - y3 * y3 + z0 * z0 - z3 * z3)) / d;
        C(i, 2) =
            ((-y0 + y3) * ((-x0 + x2) * (x0 * x0 - x1 * x1 + y0 * y0 - y1 * y1 +
                                         z0 * z0 - z1 * z1) +
                           (x0 - x1) * (x0 * x0 - x2 * x2 + y0 * y0 - y2 * y2 +
                                        z0 * z0 - z2 * z2)) +
             (x0 - x3) * ((-y0 + y2) * (x0 * x0 - x1 * x1 + y0 * y0 - y1 * y1 +
                                        z0 * z0 - z1 * z1) +
                          (y0 - y1) * (x0 * x0 - x2 * x2 + y0 * y0 - y2 * y2 +
                                       z0 * z0 - z2 * z2)) +
             (x2 * (y0 - y1) + x0 * (y1 - y2) + x1 * (-y0 + y2)) *
                 (x0 * x0 - x3 * x3 + y0 * y0 - y3 * y3 + z0 * z0 - z3 * z3)) / d;
    }
}

template void TetrahedronCenter<RowMatrixXd, RowMatrixXi, RowMatrixXd>(
    const Eigen::MatrixBase<RowMatrixXd>& V,
    const Eigen::MatrixBase<RowMatrixXi>& T,
    Eigen::PlainObjectBase<RowMatrixXd>& C, std::vector<bool>& computable);

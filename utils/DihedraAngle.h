#ifndef DIHEDRAANGLE_H
#define DIHEDRAANGLE_H

#include <Eigen/Dense>

template <typename DerivedV, typename DerivedF>
bool CosDihedralAngle(
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedF>& F,
        const typename DerivedF::Scalar& f1,
        const typename DerivedF::Scalar& f2,
        const bool reverseNormal, 
        typename DerivedV::Scalar& cosA);

#endif

#ifndef TriangleCircleRadius_h
#define TriangleCircleRadius_h

#include <Eigen/Dense>
#include <vector>

template <typename DerivedV, typename DerivedF, typename DerivedC>
void TriangleCircleCenter(
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedF>& F,
        Eigen::PlainObjectBase<DerivedC>& C,
        std::vector<bool>& computable);

template <typename DerivedV, typename DerivedF, typename DerivedR>
void SquaredTriangleCircleRadius(
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedF>& F,
        Eigen::PlainObjectBase<DerivedR>& R);

#endif

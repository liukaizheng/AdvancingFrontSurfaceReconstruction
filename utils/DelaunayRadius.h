#ifndef DELAUNAYRADIUS_H
#define DELAUNAYRADIUS_H

#include <Eigen/Dense>
#include <vector>

template <typename DerivedV, typename DerivedT, typename DerivedC, typename DerivedR, typename TypeuF2F>
void DelaunayRadius(
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedT>& T,
        const Eigen::MatrixBase<DerivedC>& C,
        const Eigen::MatrixBase<DerivedR>& TR,
        const std::vector<std::vector<TypeuF2F>>& uF2F,
        const std::vector<bool>& computable,
        Eigen::PlainObjectBase<DerivedR>& R);

#endif

#ifndef FACETOORIENTEDFACE_H
#define FACETOORIENTEDFACE_H

#include <Eigen/Dense>
#include <vector>

template <typename DerivedI>
void FaceToOrientedFace(
    const Eigen::MatrixBase<DerivedI>& F2uF,
    std::vector<std::vector<typename DerivedI::Scalar>>& uF2F);

#endif

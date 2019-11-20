#include "FaceToOrientedFace.h"
#include "PreTypeDefine.h"

template <typename DerivedI>
void FaceToOrientedFace(
    const Eigen::MatrixBase<DerivedI>& F2uF,
    std::vector<std::vector<typename DerivedI::Scalar>>& uF2F)
{
    uF2F.clear();
    auto numUniqueFaces = static_cast<std::size_t>(F2uF.maxCoeff()) + 1;
    uF2F.resize(numUniqueFaces);
    for(Eigen::Index i = 0; i < F2uF.size(); i++)
    {
        uF2F[F2uF(i)].emplace_back(static_cast<typename DerivedI::Scalar>(i));
    }
}

template void FaceToOrientedFace<Eigen::VectorXi>(
    const Eigen::MatrixBase<Eigen::VectorXi>& F2uF,
    std::vector<std::vector<int>>& uF2F);

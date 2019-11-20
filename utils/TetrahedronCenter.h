#ifndef TETRAHEDRONCENTER_H
#define TETRAHEDRONCENTER_H

#include <Eigen/Dense>
#include <vector>

template <typename DerivedV, typename DerivedT, typename DerivedC>
void TetrahedronCenter(
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedT>& T,
        Eigen::PlainObjectBase<DerivedC>& C,
        std::vector<bool>& computable);
#endif

#ifndef WRITEOBJ_H
#define WRITEOBJ_H

#include <Eigen/Dense>

template <typename DerivedV, typename DerivedF>
bool WriteObj(
        const std::string& name,
        const Eigen::MatrixBase<DerivedV>& V,
        const Eigen::MatrixBase<DerivedF>& F);
#endif

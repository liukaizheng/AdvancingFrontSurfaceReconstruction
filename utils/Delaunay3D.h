#ifndef DELAUNAY3D_H
#define DELAUNAY3D_H

#include <Eigen/Dense>

template <typename DerivedV, typename DerivedF>
void Delaunay3D(const Eigen::MatrixBase<DerivedV>& V,
                       Eigen::PlainObjectBase<DerivedF>& F);

#endif

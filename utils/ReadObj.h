#ifndef READOBJ_H
#define READOBJ_H

#include <vector>
#include <Eigen/Dense>

/*
 *V: vertex
 *TC: texture coordinate
 *N: vertex normal
 *F: face
 */
template <typename Scalar, typename Index>
bool ReadObj(
        const std::string& fileName,
        std::vector<std::vector<Scalar>>& V,
        std::vector<std::vector<Scalar>>& TC,
        std::vector<std::vector<Scalar>>& N,
        std::vector<std::vector<Index>>& F,
        std::vector<std::vector<Index>>& FTC,
        std::vector<std::vector<Index>>& FN);

/*
 *V: vertex
 *F: face
 */
template <typename DerivedV, typename DerivedF>
bool ReadObj(
        const std::string& fileName,
        Eigen::PlainObjectBase<DerivedV>& V, 
        Eigen::PlainObjectBase<DerivedF>& F);

#endif

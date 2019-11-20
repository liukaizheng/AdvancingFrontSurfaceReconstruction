#ifndef SORTMATRIX_H
#define SORTMATRIX_H
#include <Eigen/Dense>

template <typename DerivedA>
void SortCols(
        const Eigen::MatrixBase<DerivedA>& A,
        bool ascending,
        Eigen::PlainObjectBase<DerivedA>& C);

/*
 *A: input matrix,  m x n
 *C: output matrix, m x n
 *IC:index vector, m, C(i, j) = A(IC(i), j)
 */
template <typename DerivedA, typename DerivedI>
void SortRows(
        const Eigen::MatrixBase<DerivedA>& A,
        Eigen::PlainObjectBase<DerivedA>& C,
        bool ascending,
        Eigen::PlainObjectBase<DerivedI>& IC);

/*
 *A: input matrix,  m1 x n
 *C: output matrix, m2 x n (m2 <= m1)
 *IA:index vector,  m1,  A(i) --> C(IA(i))
 *IC:index vector,  m2,  C(i) --> A(IC(i))
 */
template <typename DerivedA, typename DerivedI>
void UniqueRows(
        const Eigen::MatrixBase<DerivedA>& A,
        Eigen::PlainObjectBase<DerivedA>& C,
        Eigen::PlainObjectBase<DerivedI>& IA,
        Eigen::PlainObjectBase<DerivedI>& IC);


template <typename DerivedA, typename DerivedI>
void QuickUniqueRows(
        const Eigen::MatrixBase <DerivedA>& A,
        Eigen::PlainObjectBase<DerivedA>& C,
        Eigen::PlainObjectBase<DerivedI>& IA,
        Eigen::PlainObjectBase<DerivedI>& IC);

#endif

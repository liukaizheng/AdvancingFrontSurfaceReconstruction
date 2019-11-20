#include "SortMatrix.h"
#include "PreTypeDefine.h"
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <ctime>
#include <iostream>

template <typename DerivedA>
void SortCols(
        const Eigen::MatrixBase<DerivedA>& A,
        bool ascending,
        Eigen::PlainObjectBase<DerivedA>& C)
{
    typedef typename DerivedA::Scalar Scalar;
    C.resize(A.rows(), A.cols());
    std::vector<Scalar> ss(A.cols());
    for(Eigen::Index i = 0; i < A.rows(); i++) {
        for(Eigen::Index j = 0; j < A.cols(); j++)
            ss[j] = A(i, j);
        std::sort(ss.begin(), ss.end(), [&](const Scalar a, const Scalar b)->bool {
            if(ascending) return a < b;
            else return a > b;
        });
        for(Eigen::Index j = 0; j < A.cols(); j++)
            C(i, j) = ss[j];
    }
}


template <typename DerivedA, typename DerivedI>
void SortRows(
        const Eigen::MatrixBase<DerivedA>& A,
        Eigen::PlainObjectBase<DerivedA>& C,
        bool ascending,
        Eigen::PlainObjectBase<DerivedI>& IC)
{    
    typedef Eigen::Index Index;
    const Index m = A.rows();
    IC.resize(m);
    for(Index i = 0; i < m; i++)
        IC(i) = i;
    if(ascending) {
        auto ascending = [&](const Index &i, const Index &j) ->bool {
            for(Index c = 0; c < A.cols(); c++) {
                if(A(i, c) < A(j, c)) return true;
                else if(A(i, c) > A(j, c)) return false;
            }
            return false;
        };
        std::sort(IC.data(), IC.data() + m, ascending);
    }
    else {
        auto descending = [&](const Index &i, const Index &j) ->bool {
            for(Index c = 0; c < A.cols(); c++) {
                if(A(i, c) < A(j, c)) return false;
                else if(A(i, c) > A(j, c)) return true;
            }
            return false;
        };
        std::sort(IC.data(), IC.data() + m, descending);
    }

    C.resize(m, A.cols());
    for(Index i = 0; i < m; i++)
        for(Index j = 0; j < A.cols(); j++)
            C(i, j) = A(IC(i), j);
}

template <typename DerivedA, typename DerivedI>
void UniqueRows(
        const Eigen::MatrixBase<DerivedA>& A,
        Eigen::PlainObjectBase<DerivedA>& C,
        Eigen::PlainObjectBase<DerivedI>& IA,
        Eigen::PlainObjectBase<DerivedI>& IC)
{
    typedef Eigen::Index Index;
    const Index m = A.rows();
    DerivedA cA, rA;    //cA: sorted by col; rA: sorted by row
    DerivedI IM;
    SortCols(A, true, cA);
    SortRows(cA, rA, true, IM);
    std::vector<Index> vIC(m);
    for (Index i = 0; i < m; i++) vIC[i] = i;
    auto equalIndex = [&](const Index& i, const Index& j) -> bool {
        for (Index c = 0; c < A.cols(); c++) {
            if (rA(i, c) != rA(j, c)) return false;
        }
        return true;
    };
    vIC.erase(std::unique(vIC.begin(), vIC.end(), equalIndex), vIC.end());

    IA.resize(m);
    IC.resize(vIC.size());
    Index count = 0;
    IC(0) = IM(0);
    for (Index i = 0; i < m; i++) {
        if (count < vIC.size() - 1 && i == vIC[count + 1]) {
            IA(IM(i)) = ++count;
            IC(count) = IM(i);
        } else
            IA(IM(i)) = count;
    }
    C.resize(vIC.size(), A.cols());
    for (Index i = 0; i < vIC.size(); i++) {
        C.row(i) = A.row(IC(i));
    }
}

template <typename Scalar>
struct _Point
{
    _Point() {mP.clear();}
    std::vector<Scalar> mP;
};

template <typename Scalar>
struct std::less<_Point<Scalar>*> 
{
    bool operator()(const _Point<Scalar>* p1, const _Point<Scalar>* p2)
    {
        for(unsigned i = 0; i < p1->mP.size(); i++)
        {
            if(p1->mP[i] < p2->mP[i])
                return true;
            else if(p1->mP[i] > p2->mP[i])
                return false;
        }
        return false;
    }
};

template <typename DerivedA, typename DerivedI>
void QuickUniqueRows(
        const Eigen::MatrixBase <DerivedA>& A,
        Eigen::PlainObjectBase<DerivedA>& C,
        Eigen::PlainObjectBase<DerivedI>& IA,
        Eigen::PlainObjectBase<DerivedI>& IC)
{
    typedef typename DerivedA::Scalar Scalar;
    typedef typename DerivedI::Scalar Index;
    typedef _Point<Scalar> Point;

    DerivedA cA;    //cA: sorted by col; rA: sorted by row
    SortCols(A, true, cA);

    Point* pool = new Point[A.rows()];
    std::map<Point*, Index> pMap;
    Index numPoints = static_cast<Index>(cA.rows());

    auto start = clock();
    std::unordered_map<Point*, std::vector<Index>> spMap;  //same point 
    for(Index i = 0; i < numPoints; i++)
    {
        auto& p = pool[i].mP;
        p.resize(A.cols());
        for(Eigen::Index j = 0; j < cA.cols(); j++)
            p[j] = cA(i, j);
        auto pb = pMap.emplace(std::make_pair(&pool[i], i)); //point and bool value
        if(!pb.second)
        {
            spMap[pb.first->first].emplace_back(i);
        }
    }
    auto finish = clock();
    std::cout << "step1: " << finish - start << "\n";
    start = clock();
    C.resize(pMap.size(), A.cols());
    IA.resize(A.rows());
    IC.resize(pMap.size());
    int idx = 0;
    for(auto it = pMap.begin(); it != pMap.end(); it++)
    {
        auto& p = it->first;
        for(Eigen::Index i = 0; i < A.cols(); i++)
        {
            C(idx, i) = p->mP[i];
        }
        IA(it->second) = idx;
        IC(idx) = it->second;

        for(unsigned i = 0; i < spMap[p].size(); i++)
        {
            IA(spMap[p][i]) = idx;
        }
        idx++;
    }
    finish = clock();
    std::cout << "step2: " << finish - start << "\n";

    delete[] pool;
}

template void SortCols<RowMatrixXi>(const Eigen::MatrixBase<RowMatrixXi>& A,
                                    bool ascending,
                                    Eigen::PlainObjectBase<RowMatrixXi>& C);

template void SortRows<RowMatrixXi, Eigen::VectorXi>(
    const Eigen::MatrixBase<RowMatrixXi>& A,
    Eigen::PlainObjectBase<RowMatrixXi>& C, bool ascending,
    Eigen::PlainObjectBase<Eigen::VectorXi>& IC);

template void UniqueRows<RowMatrixXi, Eigen::VectorXi>(
    const Eigen::MatrixBase<RowMatrixXi>& A,
    Eigen::PlainObjectBase<RowMatrixXi>& C,
    Eigen::PlainObjectBase<Eigen::VectorXi>& IA,
    Eigen::PlainObjectBase<Eigen::VectorXi>& IC);

template void UniqueRows<RowMatrixXd, Eigen::VectorXi>(
    const Eigen::MatrixBase<RowMatrixXd>& A,
    Eigen::PlainObjectBase<RowMatrixXd>& C,
    Eigen::PlainObjectBase<Eigen::VectorXi>& IA,
    Eigen::PlainObjectBase<Eigen::VectorXi>& IC);

template void QuickUniqueRows<RowMatrixXi, Eigen::VectorXi>(
    const Eigen::MatrixBase<RowMatrixXi>& A,
    Eigen::PlainObjectBase<RowMatrixXi>& C,
    Eigen::PlainObjectBase<Eigen::VectorXi>& IA,
    Eigen::PlainObjectBase<Eigen::VectorXi>& IC);

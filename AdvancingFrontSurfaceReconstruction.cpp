#include "AdvancingFrontSurfaceReconstruction.h"
#include "Delaunay3D.h"
#include "PreTypeDefine.h"
#include "Plausibility.h"
#include "SortMatrix.h"
#include "FaceToOrientedFace.h"
#include "TetrahedronCenter.h"
#include "TriangleCircleRadius.h"
#include "DelaunayRadius.h"
#include <ctime>
#include <algorithm>
#include <iterator>
#include <Eigen/Sparse>
#include <queue>

#include "WriteObj.h"
#include <iostream>


typedef Eigen::SparseMatrix<double> SpMat;

#define PrintMSG(msg) \
    std::cout << #msg << ": " << msg << "\n";

AdvancingFrontSurfaceReconstruction::AdvancingFrontSurfaceReconstruction(const std::vector<double>& V, AdvancingFrontSurfaceReconstruction::Mode mode, const double K, const double cosA, const double cosB)
    :mV(V), mMode(mode), mK(K), mCosA(cosA), mCosB(cosB), STANDBY_CANDIDATE(2), NEVER_VALID_CANDIDATE(3)
{
    RowMatrixXd _V;
    RowMatrixXi T;
    _V.resize(V.size() / 3, 3);
    std::copy(V.begin(), V.end(), _V.data());


    Delaunay3D(_V, T);
    mT.resize(T.size());
    std::copy(T.data(), T.data() + mT.size(), mT.begin());

    auto start = clock();
    RowMatrixXi rF; //raw faces
    rF.resize(T.rows() * 4, 3);
    //first face
    rF.topRows(T.rows()) = T.rightCols(3);
    //second face
    rF.block(T.rows(), 0, T.rows(), 1) = T.col(2);
    rF.block(T.rows(), 1, T.rows(), 1) = T.col(3);
    rF.block(T.rows(), 2, T.rows(), 1) = T.col(0);
    //third face
    rF.block(T.rows() * 2, 0, T.rows(), 1) = T.col(3);
    rF.block(T.rows() * 2, 1, T.rows(), 1) = T.col(0);
    rF.block(T.rows() * 2, 2, T.rows(), 1) = T.col(1);
    //fourth face
    rF.bottomRows(T.rows()) = T.leftCols(3);

    RowMatrixXi urF;   //unique raw faces
    Eigen::VectorXi F2uF;
    {
        Eigen::VectorXi _uF2F;
        UniqueRows(rF, urF, F2uF, _uF2F);
        std::cout << "extract faces: " << clock() - start << "\n";
    }
    mNumFaces = static_cast<int>(urF.rows());
    mF.resize(urF.size());
    std::copy(urF.data(), urF.data() + mF.size(), mF.begin());
    mF2uF.resize(F2uF.size());
    std::copy(F2uF.data(), F2uF.data() + mF2uF.size(), mF2uF.begin());
    FaceToOrientedFace(F2uF, mUF2F);

    RowMatrixXi E(urF.rows() * 3, 2);
    E.topRows(urF.rows()) = urF.rightCols(2);
    E.block(urF.rows(), 0, urF.rows(), 1) = urF.col(2);
    E.block(urF.rows(), 1, urF.rows(), 1) = urF.col(0);
    E.bottomRows(urF.rows()) = urF.leftCols(2);
    RowMatrixXi uE;
    Eigen::VectorXi E2uE;
    {
        start = clock();
        Eigen::VectorXi _uE2E;
        UniqueRows(E, uE, E2uE, _uE2E);
        std::cout << "extract edges:" << clock() - start << "\n";
        mE.resize(uE.size());
        std::copy(uE.data(), uE.data() + mE.size(), mE.begin());
        mV2E.resize(mV.size() / 3);
        for(Eigen::Index i = 0; i < uE.rows(); i++)
        {
            mV2E[uE(i, 0)].emplace_back(i);
            mV2E[uE(i, 1)].emplace_back(i);
        }

        mVE2E.resize(mV.size() / 3);

        for(int i = 0; i < mF.size() / 3; i++)
        {
            int* fv = &mF[i*3];
            for(int j = 0; j < 3; j++)
            {
                int ei[2] = {(j + 1) % 3, (j + 2) % 3};
                int e[2] = {E2uE(ei[0] * mNumFaces + i), E2uE(ei[1] * mNumFaces + i)};
                mVE2E[fv[j]][e[0]].emplace_back(e[1]);
                mVE2E[fv[j]][e[1]].emplace_back(e[0]);
            }
        }
    }
    mE2uE.resize(E2uE.size());
    std::copy(E2uE.data(), E2uE.data() + mE2uE.size(), mE2uE.begin());
    FaceToOrientedFace(E2uE, mUE2E);
	mInvalidFaces.resize(mUE2E.size());

	
	std::vector<bool> computable;
	RowMatrixXd C;
	TetrahedronCenter(_V, T, C, computable);
    Eigen::VectorXd TR, R;
    SquaredTriangleCircleRadius(_V, urF, TR);
    DelaunayRadius(_V, T, C, TR, mUF2F, computable, R);
    mR.resize(R.size());
    //std::copy(R.data(), R.data() + mR.size(), mR.begin());
	std::copy(TR.data(), TR.data() + mR.size(), mR.begin());

    mPP = static_cast<Plausibility<double, int>*>(malloc(sizeof(Plausibility<double, int>) * mUE2E.size()));

    init();
}

AdvancingFrontSurfaceReconstruction::~AdvancingFrontSurfaceReconstruction()
{
    free(mPP);
}

void AdvancingFrontSurfaceReconstruction::run()
{
    const auto recomputeValue = [this](int& f, int& ei)->int {
        const auto uf = std::abs(f) - 1;
        const auto& ue = mE2uE[uf + ei * mNumFaces];
        const auto it = mBEF.find(ue);
        if(it == mBEF.end())
            return -1;

        mInvalidFaces[ue].emplace(uf);

        const int uif = it->second.first;
        const int iei = it->second.second;

        int cf, cei;
        double cp;
        CandidateTriangle(uif, iei, cf, cei, cp);
        mPP[ue].f = cf;
        mPP[ue].ei = cei;
        mPP[ue].p = cp;
        return ue;
    };
    while(!mQ.empty())
    {
        auto qit = mQ.begin();
        std::set<Plausibility<double, int>*>::iterator::difference_type ind = 0; 
        while(!mQ.empty())
        {
            qit = mQ.begin();
            Plausibility<double, int>* pb = *qit;
            PrintMSG(mS.size());
#ifndef NDEBUG
            std::stringstream ss;
            ss << "models/m_" << mS.size() << ".obj";
            WriteObj(mS, ss.str());
#endif
            mQ.erase(qit);
            if(pb->p < STANDBY_CANDIDATE)
            {
                int ret = StitchTriangle(pb->f, pb->ei, pb->p);
                if(ret == 1)
                {
                    /*qit = mQ.begin();
                    ind = 0;*/
                }
                else if(ret == 0)
                {
                    const int ue = recomputeValue(pb->f, pb->ei);
                    if(ue >= 0)
                    {
                        if(mPP[ue].f != mNumFaces + 1)
                        {
                            mQ.emplace(&mPP[ue]);
                        }
                    }
                }
                else if(ret == 2)
                {
                    const int uf = std::abs(pb->f) - 1;
                    int* fv = &mF[uf * 3];
                    mVGP[fv[0]].emplace_back(pb);
                    mVGP[fv[1]].emplace_back(pb);
                    mVGP[fv[2]].emplace_back(pb);
                }
            }
            else 
            {
                /*mQ.erase(qit);
                qit = mQ.begin();
                ind = 0;*/
            }
        }
        if(!mQ.empty())
        {
            qit = mQ.begin();
            Plausibility<double, int>* pb = *qit;
            mQ.erase(qit);
            const int ue = recomputeValue(pb->f, pb->ei);
            if(mPP[ue].f != mNumFaces + 1)
            {
                mQ.emplace(&mPP[ue]);
            }
            qit = mQ.begin();
            ind = 0;
        }
    }
}

void AdvancingFrontSurfaceReconstruction::GetFace(std::vector<int>& faces)
{
    faces.resize(mS.size() * 3);
    int ind = 0;

    for(auto it = mS.begin(); it != mS.end(); it++)
    {
        int* ifv = &faces[ind * 3];
        int* ofv = &mF[it->first * 3];
        ifv[0] = ofv[0];
        ifv[1] = ofv[1];
        ifv[2] = ofv[2];

        if(it->second)
            std::swap(ifv[0], ifv[2]);
        ind++;
    }
}

void AdvancingFrontSurfaceReconstruction::GetBoundaryEdge(std::vector<std::vector<int>>& edges)
{
    std::unordered_map<int, std::vector<int>> V2E;
    //V2E.reserve(mBEF.size());
    for(auto it = mBEF.begin(); it != mBEF.end(); it++)
    {
        int uf = std::abs(it->second.first) - 1;
        int ei = it->second.second;

        const int* fvs = &mF[uf * 3];
        int v[2] = {fvs[(ei + 1) % 3], fvs[(ei + 2) % 3]};
        const int ue = mE2uE[ei * 3 + uf];
        V2E[v[0]].emplace_back(it->first);
        V2E[v[1]].emplace_back(it->first);
    }

    std::unordered_set<int> visited;
    //visited.reserve(mBEF.size());

    std::vector<std::unordered_set<int>> BE;

    for(auto it = mBEF.begin(); it != mBEF.end(); it++)
    {
        const auto eit = visited.find(it->first);
        if(eit != visited.end()) continue;

        std::queue<int> Q;
        Q.emplace(it->first);
        std::unordered_set<int> edge_component;
        while(!Q.empty())
        {
            const int ue = Q.front();
            Q.pop();
            edge_component.emplace(ue);
            visited.emplace(ue);
            const int* ev = &mE[ue * 2];
            for(unsigned i = 0; i < 2; i++)
            {
                for(unsigned j = 0; j < V2E[ev[i]].size(); j++)
                {
                    const int tue = V2E[ev[i]][j];
                    const auto teit = edge_component.find(tue);
                    if(teit == edge_component.end())
                        Q.emplace(tue);
                }
            }
        }
        BE.emplace_back(edge_component);
    }

    edges.resize(BE.size());
    for(unsigned i = 0; i < BE.size(); i++)
    {
        auto se = *BE[i].begin() + 1;
        int e = se;
        do {
            edges[i].emplace_back(e);
            auto ue = std::abs(e) - 1;
            int* ev = &mE[ue * 2];
            int cv = e < 0 ? ev[0] : ev[1];
            for(unsigned j = 0; j < V2E[cv].size(); j++)
            {
                const auto nue = V2E[cv][j];
                if(nue == ue) continue;

                int idx = -1;
                const int* nev = &mE[nue * 2];
                if(nev[0] == cv)
                    idx = 1;
                else if(nev[1] == cv)
                    idx = 0;

                if(idx >= 0)
                {
                    e = nue + 1;
                    if(idx == 0)
                        e = - e;
                    break;
                }
            }

        }while(e != se);
    }
}

void AdvancingFrontSurfaceReconstruction::init()
{
    if(mMode == PAGE_RANK)
        InitialMatrix();

    unsigned minFace = mNumFaces;
    double minRadius = std::numeric_limits<double>::infinity();

    std::unordered_set<int> invalidFaces;
    bool unfinished =true;
	int ind = 0;
    while(unfinished)
    {
		PrintMSG(ind++);
        for(unsigned i = 0; i < mR.size(); i++)
        {
            if(mR[i] < minRadius && mR[i] > std::numeric_limits<double>::epsilon())
            {
                auto it = invalidFaces.find(i);
                if(it == invalidFaces.end())
                {
                    minFace = i;
                    minRadius = mR[i];
                }
            }
        }
        if(minFace != mNumFaces)
        {
            std::vector<Plausibility<double, int>*> candidates;
            const int f = minFace + 1;
            for(unsigned i = 0; i < 3; i++)
            {
                int cf, cei;
                double cp;
                CandidateTriangle(f, i, cf, cei, cp);
                if(cp < -1.0)
                {
                    auto pb = &mPP[mE2uE[minFace + i * mNumFaces]];
                    pb->p = cp;
                    pb->f = cf;
                    pb->ei = cei;
                    candidates.emplace_back(pb);
                }
            }
            if(candidates.size() == 3)
            {
                for(unsigned i = 0; i < 3; i++)
                    mQ.emplace(candidates[i]);
                unfinished = false;
            }
            else
            {
                invalidFaces.emplace(minFace);
                minFace = mNumFaces;
                minRadius = std::numeric_limits<double>::infinity();
            }
        }
        else
            unfinished = false;
    }


    if(minFace != mNumFaces)
    {
        int fe[3] = {mE2uE[minFace], mE2uE[minFace + mNumFaces], mE2uE[minFace + 2 * mNumFaces]};
        int* fv = &mF[minFace * 3];
        const int f = minFace + 1;

        mS.clear();
        mS.emplace(std::make_pair(minFace, false));

        mSV.clear();
        mSV.resize(mV.size() / 3, false);
        mSV[fv[0]] = true;
        mSV[fv[1]] = true;
        mSV[fv[2]] = true;

        mSE.clear();
        mSE.emplace(fe[0]);
        mSE.emplace(fe[1]);
        mSE.emplace(fe[2]);

        mBEF.clear();
        mBEF.emplace(std::make_pair(fe[0], std::make_pair(f, 0)));
        mBEF.emplace(std::make_pair(fe[1], std::make_pair(f, 1)));
        mBEF.emplace(std::make_pair(fe[2], std::make_pair(f, 2)));

        mBVE.clear();
        mBVE.resize(mV.size() / 3);
        mBVE[fv[0]].emplace_back(fe[1]);
        mBVE[fv[0]].emplace_back(fe[2]);
        mBVE[fv[1]].emplace_back(fe[2]);
        mBVE[fv[1]].emplace_back(fe[0]);
        mBVE[fv[2]].emplace_back(fe[0]);
        mBVE[fv[2]].emplace_back(fe[1]);

        mVGP.clear();
        mVGP.resize(mV.size() / 3);
    }

}

void AdvancingFrontSurfaceReconstruction::InitialMatrix()
{
    SpMat prMat(mNumFaces, mNumFaces);
    std::vector<Eigen::Triplet<double>> triplets;
    for(int i = 0; i < mNumFaces; i++)
    {

        std::set<int> nfs;
        for(int j = 0; j < 3; j++)
        {
            const auto ue = mE2uE[i + j * mNumFaces];
            const auto& nes = mUE2E[ue];
            if(nes.size() > 0)
            {
                for(int k = 0; k < nes.size(); k++) 
                {
                    const auto nf = nes[k] % mNumFaces;
                    if(nf != i)
                    {
                        nfs.emplace(nf);
                    }
                }
            }
        }

        for(auto it = nfs.begin();it != nfs.end(); it++)
        {
            triplets.emplace_back(Eigen::Triplet<double>(*it, i, 1.0 / nfs.size()));
        }
    }
    prMat.setFromTriplets(triplets.begin(), triplets.end());

    double diff = std::numeric_limits<double>::infinity();
    auto v = Eigen::VectorXd::Constant(mNumFaces, 1, 1.0 / mNumFaces).eval();
    while(diff > 0.0001)
    {
        auto v1 = prMat * v;
        diff = (v1 - v).cwiseAbs().sum();
        v = v1;
    }
    mPRVec.resize(mNumFaces);
    std::copy(v.data(), v.data() + mNumFaces, mPRVec.begin());
}

inline bool AdvancingFrontSurfaceReconstruction::DihedralAngle(const int& f1, const int& f2, double& value)
{
    const int uf1 = std::abs(f1) - 1;
    const int uf2 = std::abs(f2) - 1;
    const int& f10 = mF[3*uf1];
    const int& f11 = mF[3*uf1 + 1];
    const int& f12 = mF[3*uf1 + 2];
    const int& f20 = mF[3*uf2];
    const int& f21 = mF[3*uf2 + 1];
    const int& f22 = mF[3*uf2 + 2];
    Eigen::Vector3d v11(mV[3 * f11] - mV[3 * f10], mV[3 * f11 + 1] - mV[3*f10 + 1], mV[3*f11 + 2] - mV[3*f10 + 2]);
    Eigen::Vector3d v12(mV[3 * f12] - mV[3 * f10], mV[3 * f12 + 1] - mV[3*f10 + 1], mV[3*f12 + 2] - mV[3*f10 + 2]);
    Eigen::Vector3d v1 = v11.cross(v12);
    if(v1.squaredNorm() < std::numeric_limits<double>::epsilon())
        return false;

    Eigen::Vector3d v21(mV[3 * f21] - mV[3 * f20], mV[3 * f21 + 1] - mV[3*f20 + 1], mV[3*f21 + 2] - mV[3*f20 + 2]);
    Eigen::Vector3d v22(mV[3 * f22] - mV[3 * f20], mV[3 * f22 + 1] - mV[3*f20 + 1], mV[3*f22 + 2] - mV[3*f20 + 2]);
    Eigen::Vector3d v2 = v21.cross(v22);
    if(v2.squaredNorm() < std::numeric_limits<double>::epsilon())
        return false;

    v1.normalize();
    v2.normalize();
    value = v1.dot(v2);
    if((f1 < 0) ^ (f2 < 0))
        value = -value;
    return true;
}

inline void AdvancingFrontSurfaceReconstruction::CandidateTriangle(const int& f, const int& ei, int& cf, int& cei, double& value)
{
    const int uf = std::abs(f) - 1;
    cf = mNumFaces + 1;
    double cr = std::numeric_limits<double>::infinity(), ca;    //radius and angle
    const auto ue = mE2uE[uf + ei * mNumFaces];
    const auto& adjEdges = mUE2E[ue];
    int v11 = mF[3 * uf + (ei + 1) % 3];
    int v12 = mF[3 * uf + (ei + 2) % 3];
    if(f < 0)
        std::swap(v11, v12);

    for(unsigned i = 0; i < adjEdges.size(); i++)
    {
        int af = adjEdges[i] % mNumFaces;
        if(std::abs(af) == uf) continue;
        const auto it = mInvalidFaces[ue].find(af);
        if(it != mInvalidFaces[ue].end()) continue;

        const int aei = adjEdges[i] / mNumFaces;
        const int& v21 = mF[3*af + (aei + 1) % 3]; 
        const int& v22 = mF[3*af + (aei + 2) % 3]; 
        double da;
        if(v11 == v21)
        {
            af = -(af + 1);
#ifndef NDEBUG
            assert(v12 == v22);
#endif
        }
        else
        {
            af = af + 1;
#ifndef NDEBUG
            assert(v11 == v22 && v12 == v21);
#endif
        }

        const int uaf =  std::abs(af) - 1;
        if(DihedralAngle(f, af, da))
        {
            if(da > mCosA)
            {
                if(mMode == DEFAULT)
                {
                    if(mR[uaf] < cr)
                    {
                        cf = af;
                        cei = aei;
                        cr = mR[uaf];
                        ca = da;
                    }
                }
                else if(mMode == PAGE_RANK)
                {
                    if(mR[uaf] > std::numeric_limits<double>::epsilon())
                    {
                        /*double val = -1.0 / (mPRVec[uaf] * mR[uaf]);*/
                        double val = -1.0 /mR[uaf] * std::pow(1.0 / mPRVec[uaf], 2.0);
                        if(val < cr)
                        {
                            cf = af;
                            cei = aei;
                            cr = val;
                        }
                    }
                }
            }
        }
    }

    const auto ucf = std::abs(cf) - 1;
    if(ucf == mNumFaces || mR[ucf] < std::numeric_limits<double>::epsilon())
    {
        value = NEVER_VALID_CANDIDATE;
        return;
    }

    if(mMode == PAGE_RANK)
    {
        value = cr;
        return;
    }

    if(ca > mCosB)
    {
        value = -(1.0 + 1.0 / cr);
    }
    else
    {
        if( mR[ucf] < mK * mR[uf] || mR[uf] < std::numeric_limits<double>::epsilon())
        {
            value = -ca;
        }
        else
        {
            value = NEVER_VALID_CANDIDATE;      //reasonable??
        }
    }
}

inline AdvancingFrontSurfaceReconstruction::Type AdvancingFrontSurfaceReconstruction::Validate(const int& f, const int& ei)
{
    const auto uf = std::abs(f) - 1;
    if(!mSV[mF[uf * 3 + ei]])
        return EXTENSION;

    const auto& ue1 = mE2uE[uf + ((ei + 1) % 3) * mNumFaces];
    const auto& ue2 = mE2uE[uf + ((ei + 2) % 3) * mNumFaces];
    const auto bit1 = mBEF.find(ue1);
    const auto bit2 = mBEF.find(ue2);
    if(bit1 != mBEF.end() && bit2 != mBEF.end())
        return HOLE_FILLING;

    const auto sit1 = mSE.find(ue1);
    const auto sit2 = mSE.find(ue2);

    if(bit1 != mBEF.end() && sit2 == mSE.end())
        return EARING_FILLING_1;

    if(bit2 != mBEF.end() && sit1 == mSE.end())
        return EARING_FILLING_2;

    if(sit1 == mSE.end() && sit2 == mSE.end() && !mBVE[mF[uf * 3 + ei]].empty())
        return GLUING;

    return NOT_VALID;
}

inline int AdvancingFrontSurfaceReconstruction::StitchTriangle(const int& f, const int& ei, const double& p)
{
    const auto uf = std::abs(f) - 1;
    const int* fv = &mF[3 * uf];

    const auto addFace = [&, this]() {
        mS.emplace(std::make_pair(uf, f < 0));
        mSV[fv[0]] = true;
        mSV[fv[1]] = true;
        mSV[fv[2]] = true;
        mSE.emplace(mE2uE[uf]);
        mSE.emplace(mE2uE[uf + mNumFaces]);
        mSE.emplace(mE2uE[uf + 2 * mNumFaces]);
    };

    const auto addGluingCandidate = [&, this](const int& i) {
        const int& v = fv[i];
        if(mVGP[v].empty()) return;

        for(auto it = mVGP[v].begin(); it != mVGP[v].end(); it++) 
        {
            mQ.emplace(*it);
        }
    };

    const auto removeEdge = [&, this](const int& i) {
        const auto& ue = mE2uE[uf + i * mNumFaces];
        auto eit = mBEF.find(ue);
#ifndef NDEBUG
        assert(eit != mBEF.end());
#endif
        if(eit != mBEF.end())
        {
            mBEF.erase(eit);
            int ev[2] = {fv[(i+1)%3], fv[(i+2)%3]};
            for(unsigned j = 0; j < 2; j++)
            {
                auto vit = std::find(mBVE[ev[j]].begin(), mBVE[ev[j]].end(), ue);
#ifndef NDEBUG
                assert(vit != mBVE[ev[j]].end());
#endif
                mBVE[ev[j]].erase(vit);
            }
        }
        auto qit = mQ.find(&mPP[ue]);
        if(qit != mQ.end())
        {
            mQ.erase(qit);
        }
    };

    const auto addEdge = [&, this](const int& i) {
        const auto ue = mE2uE[uf + i * mNumFaces];
        const auto eit = mBEF.find(ue);
#ifndef NDEBUG
        assert(eit == mBEF.end());
#endif
        mBEF.emplace(std::make_pair(ue, std::make_pair(f, i)));
        mBVE[fv[(i+1)%3]].emplace_back(ue);
        mBVE[fv[(i+2)%3]].emplace_back(ue);
    };

    const auto addCandidate = [&, this](const int& i) {
        int cf, cei;
        double cp;
        CandidateTriangle(f, i, cf, cei, cp);
        const auto& ue = mE2uE[uf + i * mNumFaces];
        mPP[ue].p = cp;
        mPP[ue].f = cf;
        mPP[ue].ei = cei;
        mQ.emplace(&mPP[ue]);
    };

#ifndef NDEBUG
    const auto _ue = mE2uE[uf + ei * mNumFaces];
    const auto _it = mBEF.find(_ue);
    if(_it != mBEF.end())
        WriteCandidate(_it->second.first, _it->second.second);
    /*PrintMSG(fv[0]);
    PrintMSG(fv[1]);
    PrintMSG(fv[2]);
    PrintMSG(_ue);
    std::cout << "\n";*/
#endif
    Type type = Validate(f, ei);
    int eis[2] = {(ei + 1) % 3, (ei + 2) % 3};
    if(type == EXTENSION)
    {
        addFace();
        removeEdge(ei);
        addEdge(eis[0]);
        addEdge(eis[1]);
        addCandidate(eis[0]);
        addCandidate(eis[1]);
        addGluingCandidate(0);
        addGluingCandidate(1);
        addGluingCandidate(2);
        return 1;
    }
    else if(type == HOLE_FILLING)
    {
        addFace();
        removeEdge(ei);
        removeEdge(eis[0]);
        removeEdge(eis[1]);
        mVGP[fv[0]].clear();
        mVGP[fv[1]].clear();
        mVGP[fv[2]].clear();
        return 1;
    }
    else if(type == EARING_FILLING_1 || type == EARING_FILLING_2)
    {
        addFace();
        int ie = eis[0], re = eis[1];   //ie: edge to be inserted, re: edge to be removed
        if(type == EARING_FILLING_1)
            std::swap(ie, re);
        removeEdge(ei);
        removeEdge(re);
        addEdge(ie);
        addCandidate(ie);
        mVGP[fv[ie]].clear();
        addGluingCandidate(re);
        addGluingCandidate(ei);
        return 1;
    }
    else if(type == GLUING)
    {
        const auto& v0 = fv[ei];
#ifndef NDEBUG
        assert(mBVE[v0].size() == 2);
#endif
        const int& ue = mE2uE[uf + ei * mNumFaces];
        const auto face_pair = mBEF[ue];

        addFace();
        removeEdge(ei);
        addEdge(eis[0]);
        addEdge(eis[1]);

        int abv1[2];
        int abv2[2];

        for(unsigned i = 0; i < 2; i++)
        {
            const auto& be = mBVE[v0][i];
            const auto it = mBEF.find(be);
#ifndef NDEBUG
            assert(it != mBEF.end());
#endif
            const int& abf = it->second.first;
            const int& abei = it->second.second;
            const int uabf = std::abs(abf) - 1;
            const int* abfv = &mF[uabf * 3];
            abv1[i] = abfv[(abei + 1) % 3];
            abv2[i] = abfv[(abei + 2) % 3];
            if(abf < 0)
                std::swap(abv1[i], abv2[i]);
        }

        /*std::vector<std::tuple<int, int, double>> candidates;
        std::vector<unsigned> indices;*/

        std::vector<int> ge;            //gluing triangle edges
        std::vector<int> indices;       
        for(unsigned i = 0; i < 2; i++)
        {
            const auto cue = mE2uE[uf + eis[i] * mNumFaces];
            //const auto cit = mGE.find(cue);
            ge.emplace_back(cue);
            Plausibility<double, int>* pb = &mPP[cue];
            //if(cit == mGE.end())
            //{
                int cf, cei;
                double cp;
                CandidateTriangle(f, eis[i], cf, cei, cp);
                pb->f = cf;
                pb->ei = cei;
                pb->p = cp;
                //mGE.emplace(cue);
            //}
            if(pb->p < STANDBY_CANDIDATE)
            {
                const int* cfv = &mF[(std::abs(pb->f) - 1) * 3];
                int cvei = (pb->ei + 1) % 3;
                if(cfv[cvei] == v0)
                    cvei = (cvei + 1) % 3;
                int cv[2] = {cfv[(cvei + 1) % 3], cfv[(cvei + 2) % 3]};
                if(pb->f < 0)
                    std::swap(cv[0], cv[1]);

                if((abv1[0] == cv[1] && abv2[0] == cv[0]) ||
                   (abv1[1] == cv[1] && abv2[1] == cv[0]))
                {
                    Type t = Validate(pb->f, pb->ei);
                    if(t == EXTENSION || t == HOLE_FILLING || t == EARING_FILLING_1 || t == EARING_FILLING_2)
                    {
                        if(pb->p <= p)
                            indices.emplace_back(i);
                    }
                }
            }
        }
        if(indices.empty())
        {
            const auto sit = mS.find(uf);
#ifndef NDEBUG
            assert(sit != mS.end());
#endif
            mS.erase(sit);
            for(unsigned i = 1; i < 3; i++)
            {
                const int rei = (ei + i) % 3;
                const int& ae = mE2uE[uf + rei * mNumFaces];
                const auto eit = mBEF.find(ae);
#ifndef NDEBUG
                assert(eit != mBEF.end());
#endif
                mBEF.erase(eit);
                int rv[2] = {fv[(rei + 1) % 3], fv[(rei + 2) % 3]};
                for(unsigned j = 0; j < 2; j++)
                {
                    const auto vit = std::find(mBVE[rv[j]].begin(), mBVE[rv[j]].end(), ae);
#ifndef NDEBUG
                    assert(vit != mBVE[rv[j]].end());
#endif
                    mBVE[rv[j]].erase(vit);
                }

                const auto seit = mSE.find(ae);
#ifndef NDEBUG
                assert(seit != mSE.end());
#endif
                mSE.erase(seit);
            }

            mBEF.emplace(std::make_pair(ue, face_pair));
            mBVE[fv[(ei + 1) % 3]].emplace_back(ue);
            mBVE[fv[(ei + 2) % 3]].emplace_back(ue);
            //return 2;
            return 0;
        }

        int ind = -1;
        if(indices.size() == 1)
        {
            ind = indices[0] == 0 ? 1 : 0;
        }
        else if(indices.size() == 2)
        {
            if(mPP[ge[indices[1]]].p < mPP[ge[indices[0]]].p)
            {
                indices[0] = 1;
                indices[1] = 0;
            }
        }

        addGluingCandidate(0);
        addGluingCandidate(1);
        addGluingCandidate(2);

        for(unsigned i = 0; i < indices.size(); i++)
        {
            int ret = StitchTriangle(mPP[ge[indices[i]]].f, mPP[ge[indices[i]]].ei, mPP[ge[indices[i]]].p);
#ifndef NDEBUG
            if(i == 0)
                assert(ret == 1);
#endif
            if(ret != 1)
            {
                ind = indices[i];
            }
        }
        if(ind >= 0)
        {

            if(mPP[ge[ind]].f != mNumFaces + 1) 
            {
                mQ.emplace(&mPP[ge[ind]]);
                std::string here = "here";
                PrintMSG(here);
            }
        }
        return 1;
    }
    return 0;
}

#ifndef NDEBUG
void AdvancingFrontSurfaceReconstruction::WriteObj(const std::unordered_map<int, bool>& S, const std::string name)
{
    unsigned ind = 0;
    RowMatrixXi F(S.size(), 3);
    for(auto it = S.begin(); it != S.end(); it++) 
    {
        const int* fv =  &mF[it->first * 3];
        F(ind, 0) = fv[0];
        F(ind, 1) = fv[1];
        F(ind, 2) = fv[2];
        if(it->second)
            F.row(ind) = F.row(ind).reverse().eval();
        ind++;
    }

    RowMatrixXd V(mV.size() / 3 ,3);
    std::copy(mV.begin(), mV.end(), V.data());
    ::WriteObj(name, V, F);
}

void AdvancingFrontSurfaceReconstruction::WriteCandidate(const int& f, const int& ei)
{

    const int uf = std::abs(f) - 1;
    double cv = mCosA;
    const auto& adjEdges = mUE2E[mE2uE[uf + ei * mNumFaces]];
    int v11 = mF[3 * uf + (ei + 1) % 3];
    int v12 = mF[3 * uf + (ei + 2) % 3];
    if(f < 0)
        std::swap(v11, v12);

    std::unordered_map<int, bool> S;
    for(unsigned i = 0; i < adjEdges.size(); i++)
    {
        int af = adjEdges[i] % mNumFaces;
        const int aei = adjEdges[i] / mNumFaces;
        const int& v21 = mF[3*af + (aei + 1) % 3]; 
        const int& v22 = mF[3*af + (aei + 2) % 3]; 
        if(v11 == v21)
        {
            af = -af;
        }
        S.emplace(std::make_pair(std::abs(af), af < 0)); 
    }
    std::stringstream ss;
    ss << "models/m_" << mS.size() <<"c.obj";
    WriteObj(S, ss.str());
}

int AdvancingFrontSurfaceReconstruction::GetBoundaryEdge(int v1, int v2)
{
    for(auto it = mBEF.begin(); it != mBEF.end(); it++)
    {
        int uf = std::abs(it->second.first) - 1;
        int ei = it->second.second;

        const int* fvs = &mF[uf * 3];
        int v[2] = {fvs[(ei + 1) % 3], fvs[(ei + 2) % 3]};
        if ((v[0] == v1 && v[1] == v2) || (v[0] == v2 && v[1] == v1))
            return it->first;
    }
    return - 1;
}
#endif

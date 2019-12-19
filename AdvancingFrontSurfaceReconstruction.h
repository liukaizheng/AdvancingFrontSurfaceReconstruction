#ifndef ADVANCINGFRONTSURFACERECONSTRUCTION_H
#define ADVANCINGFRONTSURFACERECONSTRUCTION_H

#include "Plausibility.h"
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <cmath>

/*template <typename Scalar, typename Index>
struct Plausibility;*/


class AdvancingFrontSurfaceReconstruction
{
public:
    enum Type {
        EXTENSION,
        EARING_FILLING_1,
        EARING_FILLING_2,
        HOLE_FILLING,
        GLUING,
        NOT_VALID
    };
    AdvancingFrontSurfaceReconstruction(const std::vector<double>& V, const double K = 5.0, const double cosA = std::cos(M_PI / 6.0 * 5.0), const double cosB = std::cos(M_PI / 6.0));
    virtual ~AdvancingFrontSurfaceReconstruction();

    void run();
    void GetFace(std::vector<int>& face);
    void GetBoundaryEdge(std::vector<std::vector<int>>& edges);

#ifndef NDEBUG
    void WriteObj(const std::unordered_map<int, bool>& S, const std::string name);
    void WriteCandidate(const int& f, const int& ei);
    int GetBoundaryEdge(int v1, int v2);
#endif
protected:
    void init();
    bool DihedralAngle(const int& f1, const int& f2, double& value);
    void CandidateTriangle(const int& f, const int& ei, int& cf, int& cei, double& value, bool consideringRadius = true);
    Type Validate(const int& f, const int& ei);
    int StitchTriangle(const int& f, const int& ei, const double& p);

    const std::vector<double>& mV;          //vertices
    double mK;
    double mCosA;
    double mCosB;
    const double STANDBY_CANDIDATE;
    const double NEVER_VALID_CANDIDATE;

    std::vector<int> mT;                                    //tetrahedrons
    std::vector<int> mF;                                    //faces
    std::vector<int> mE;                                    //edges
    int mNumFaces;
    std::vector<int> mF2uF;                                 //every oriented face maps to its unique face
    std::vector<std::vector<int>> mUF2F;                    //neighbors of every face
    std::vector<int> mE2uE;                                 //every oriented edge mpas to its unique edge
    std::vector<std::vector<int>> mUE2E;                    //neighbors of every edge
    std::vector<double> mR;                                 //delaunay radius
    double mAR;                                             //average radius

    Plausibility<double, int>* mPP;                         //Plausibility pool
    std::unordered_map<int, bool> mS;                       //surface
    std::vector<bool> mSV;                                  //vertices on surface     
    std::unordered_set<int> mSE;                            //edges on surface
    std::unordered_map<int, std::pair<int, int>> mBEF;      //map from boundary edge to its parent triangle
    std::vector<std::vector<int>> mBVE;                     //boundary edges ownd by boundary vertices
    std::vector<std::vector<Plausibility<double, int>*> > mVGP;      
    std::set<Plausibility<double, int>*> mQ;                //Plausibility queue
    std::vector<std::unordered_set<int>> mInvalidFaces;     //invalid faces which cannot be selected


    std::vector<std::vector<int>> mV2E;
    std::vector<std::unordered_map<int, std::vector<int>>> mVE2E;
};

#endif

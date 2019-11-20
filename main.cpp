#include "AdvancingFrontSurfaceReconstruction.h"
#include "ReadObj.h"
#include "WriteObj.h"
#include "SortMatrix.h"
#include <algorithm>
#include <iostream>
#include <ctime>

#define PrintMSG(msg) \
    std::cout << #msg << ": " << msg << "\n";


int main()
{
    /*auto FF = Eigen::Matrix<int, -1, -1, 1>::Random(1000000, 3).eval();
    auto C = FF;
    Eigen::VectorXi IA, IC;
    auto start = clock();
    QuickUniqueRows(FF, C, IA, IC);
    auto finish = clock();
    PrintMSG((finish - start));
    start = clock();
    UniqueRows(FF, C, IA, IC);
    finish = clock();
    PrintMSG((finish - start));
*/
    Eigen::Matrix<double, -1, -1, 1> _V, V;
    Eigen::Matrix<int, -1, -1, 1> F;
    ReadObj("bunny.obj", _V, F);
    {
        Eigen::VectorXi _IA, _IC;
        UniqueRows(_V, V, _IA, _IC);
    }
	V = _V;
    std::vector<double> vdata(V.size());
    std::copy(V.data(), V.data() + vdata.size(), vdata.begin());
    AdvancingFrontSurfaceReconstruction afsr(vdata, AdvancingFrontSurfaceReconstruction::DEFAULT, 5.0, std::cos(M_PI / 6.0 * 5.0), std::cos(M_PI / 6.0));
    afsr.run();
    std::vector<int> faces;
    afsr.GetFace(faces);
    /*std::vector<std::vector<int>> boundary_edge;
    afsr.GetBoundaryEdge(boundary_edge);*/

#ifndef NDEBUG
    /*PrintMSG(afsr.GetBoundaryEdge(58, 68));*/
#endif

    F.resize(faces.size() / 3, 3);
    std::copy(faces.begin(), faces.end(), F.data());
    WriteObj("123.obj", V, F);
}

#include "Delaunay3D.h"
#include "PreTypeDefine.h"
#include "tetgen.h"

template <typename DerivedV, typename DerivedF>
void Delaunay3D(const Eigen::MatrixBase<DerivedV>& V,
                       Eigen::PlainObjectBase<DerivedF>& F)
{
    typedef Eigen::Map<Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> MapXdr;
    typedef Eigen::Map<Eigen::Matrix<typename DerivedF::Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> MapXir;
    
    tetgenio _in;
    _in.numberofpoints = V.rows();
    _in.pointlist = static_cast<typename DerivedV::Scalar*>(calloc(V.size(), sizeof(typename DerivedV::Scalar)));
    MapXdr pointList(_in.pointlist, V.rows(), V.cols());
    pointList = V;

    std::string str = "";
    char* switches = const_cast<char *>(str.c_str());
    tetgenio out;
    tetrahedralize(switches, &_in, &out);

    F = MapXir(out.tetrahedronlist, out.numberoftetrahedra, 4);
}


template void Delaunay3D<RowMatrixXd, RowMatrixXi>(const Eigen::MatrixBase<RowMatrixXd>&, Eigen::PlainObjectBase<RowMatrixXi>&);

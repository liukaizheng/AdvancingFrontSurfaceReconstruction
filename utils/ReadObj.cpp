#include "ReadObj.h"
#include "PreTypeDefine.h"
#include <cstdio>
#include <iostream>

template <typename Scalar, typename Index>
bool ReadObj(
        const std::string& fileName,
        std::vector<std::vector<Scalar>>& V,
        std::vector<std::vector<Scalar>>& TC,
        std::vector<std::vector<Scalar>>& N,
        std::vector<std::vector<Index>>& F,
        std::vector<std::vector<Index>>& FTC,
        std::vector<std::vector<Index>>& FN)
{    
    FILE *obj_file = fopen(fileName.c_str(), "r");
    if(obj_file == NULL) {
        std::cout << "maybe error in filename\n";
        return false;
    }
    V.clear();
    TC.clear();
    N.clear();
    F.clear();
    FTC.clear();
    FN.clear();

    std::string v("v");
    std::string vn("vn");
    std::string vt("vt");
    std::string f("f");
#define MAX_LINE 2048
    char line[MAX_LINE];
    std::vector<Scalar> tv(3);
    int count = 0;
    while(fgets(line, MAX_LINE, obj_file)) {
        char type[MAX_LINE];
        if(sscanf(line, "%s", type) == 1) {
            char *c = &line[strlen(type)];
            if(type == v) {
                count = sscanf(c, "%lf %lf %lf\n", &tv[0], &tv[1], &tv[2]);
                if(count < 3 ) {
                    std::cout << "ERROR:vertex should have at least 3 coordinates\n";
                    fclose(obj_file);
                    return false;
                }
                V.push_back(tv);
            }
            else if(type == vn) {
                count = sscanf(c, "%lf %lf %lf\n", &tv[0], &tv[1], &tv[2]);
                if(count != 3) {
                    std::cout << "ERROR: normal should have 3 coordinates\n";
                    fclose(obj_file);
                    return false;
                }
                N.push_back(tv);
            }
            else if(type == vt) {
                count = sscanf(c, "%lf %lf %lf\n", &tv[0], &tv[1], &tv[2]);
                if(count == 3) {
                    TC.push_back(tv);
                }
                else if(count == 2) {
                    tv.resize(2);
                    TC.push_back(tv);
                }
                else {
                    std::cout << "ERROR:texture should have 3 or 2 coordinates\n";
                    fclose(obj_file);
                    return false;
                }
            }
            else if(type == f) {
                const auto & shift = [&V](const int i)->int
                {
                    return i<0 ? i+V.size() : i-1;
                };
                const auto & shift_t = [&TC](const int i)->int
                {
                    return i<0 ? i+TC.size() : i-1;
                };
                const auto & shift_n = [&N](const int i)->int
                {
                    return i<0 ? i+N.size() : i-1;
                };
                std::vector<Index > f;
                std::vector<Index > ftc;
                std::vector<Index > fn;

                // Read each "word" after type
                char word[MAX_LINE];
                int offset;
                while(sscanf(c,"%s%n",word,&offset) == 1)
                {
                    // adjust offset
                    c += offset;
                    // Process word
                    long int i,it,in;
                    if(sscanf(word,"%ld/%ld/%ld",&i,&it,&in) == 3)
                    {
                        f.push_back(shift(i));
                        ftc.push_back(shift_t(it));
                        fn.push_back(shift_n(in));
                    }else if(sscanf(word,"%ld/%ld",&i,&it) == 2)
                    {
                        f.push_back(shift(i));
                        ftc.push_back(shift_t(it));
                    }else if(sscanf(word,"%ld//%ld",&i,&in) == 2)
                    {
                        f.push_back(shift(i));
                        fn.push_back(shift_n(in));
                    }else if(sscanf(word,"%ld",&i) == 1)
                    {
                        f.push_back(shift(i));
                    }else
                    {
                        std::cout << "Error: readOBJ() face on line has invalid element format\n";
                        fclose(obj_file);
                        return false;
                    }
                }
                if ((f.size() > 0 && fn.size() == 0 && ftc.size() == 0) ||
                    (f.size() > 0 && fn.size() == f.size() &&
                     ftc.size() == 0) ||
                    (f.size() > 0 && fn.size() == 0 &&
                     ftc.size() == f.size()) ||
                    (f.size() > 0 && fn.size() == f.size() &&
                     ftc.size() == f.size())) {
                    // No matter what add each type to lists so that lists are the
                    // correct lengths
                    F.push_back(f);
                    FTC.push_back(ftc);
                    FN.push_back(fn);
                } else {
                    std::cout << "Error: readOBJ() face on line has invalid format\n";
                    fclose(obj_file);
                    return false;
                }

            }
            else if(strlen(type) >= 1 && (type[0] == '#' ||
                                          type[0] == 'g'  ||
                                          type[0] == 's'  ||
                                          strcmp("usemtl",type)==0 ||
                                          strcmp("mtllib",type)==0))
            {
                //ignore comments or other shit
            }
            else
            {
                std::cout <<"Warning: readOBJ() ignored non-comment line \n";
            }
        }
        else
        {
            // ignore empty line
        }
    }
    fclose(obj_file);

    assert(F.size() == FN.size());
    assert(F.size() == FTC.size());

    return true;
}



template <typename DerivedV, typename DerivedF>
bool ReadObj(
        const std::string& fileName,
        Eigen::PlainObjectBase <DerivedV>& V, 
        Eigen::PlainObjectBase <DerivedF>& F)
{    
    std::vector<std::vector<typename DerivedV::Scalar>> vV, vTC, VN;
    std::vector<std::vector<typename DerivedF::Scalar>> vF, vFTC, vFN;
    bool ret = ReadObj(fileName, vV, vTC, VN, vF, vFTC, vFN);
    if(!ret) return ret;

    V.resize(vV.size(), 3);
    F.resize(vF.size(), 3);

    for(Eigen::Index i = 0; i < V.rows(); i++)
        for(Eigen::Index j = 0; j < V.cols(); j++)
            V(i, j) = vV[i][j];

    for(Eigen::Index i = 0 ; i < F.rows(); i++)
        for(Eigen::Index j = 0; j < F.cols(); j++)
            F(i, j) = vF[i][j];

    return true;
}

template bool ReadObj<double, int>(const std::string& fileName,
                                   std::vector<std::vector<double>>& V,
                                   std::vector<std::vector<double>>& TC,
                                   std::vector<std::vector<double>>& N,
                                   std::vector<std::vector<int>>& F,
                                   std::vector<std::vector<int>>& FTC,
                                   std::vector<std::vector<int>>& FN);

template bool ReadObj<RowMatrixXd, RowMatrixXi>(
    const std::string& fileName, Eigen::PlainObjectBase<RowMatrixXd>& V,
    Eigen::PlainObjectBase<RowMatrixXi>& F);

#ifndef Morphing_H
#define Morphing_H

#include "ShapeMatch.h" 

class Morphing
{
public:
    Morphing();
    ~Morphing();

    void computeTransform();
    void addMeshVec(std::vector<double> &mesh_vec);

private:
    std::vector<std::vector<double>> meshes_vec;
    std::vector<Eigen::Matrix3d> A_mats;
    std::vector<Eigen::Vector3d> t_vecs;
};

#endif
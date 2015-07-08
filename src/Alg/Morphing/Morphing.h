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
    Eigen::Matrix3d *getAMat(size_t i_mesh) { return &A_mats[i_mesh]; };
    std::vector<double> &getMeshVec(size_t i_mesh) { return meshes_vec[i_mesh]; };
    void alignMeshes();

private:
    std::vector<std::vector<double>> meshes_vec;
    std::vector<Eigen::Matrix3d> A_mats;
    std::vector<Eigen::Vector3d> t_vecs;
    std::vector<Eigen::Matrix3d> R_mats;
};

#endif
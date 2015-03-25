#include "Morphing.h"

Morphing::Morphing()
{
    meshes_vec.clear();
    A_mats.clear();
    t_vecs.clear();
}

Morphing::~Morphing()
{

}

void Morphing::computeTransform()
{
    ShapeMatch *shape_match;
    shape_match = new ShapeMatch;
    A_mats.clear();
    t_vecs.clear();
    std::cout<<"Total meshes: "<<meshes_vec.size()<<"\n";

    // we set meshes_vec[0] as default template mesh
    shape_match->setTemplate(meshes_vec[0]);
    A_mats.push_back(Eigen::Matrix3d::Identity());
    t_vecs.push_back(Eigen::Vector3d::Zero());

    for (size_t i = 1; i < meshes_vec.size(); ++i)
    {
        shape_match->setTarget(meshes_vec[i]);
        A_mats.push_back(shape_match->computeA());
        t_vecs.push_back(shape_match->computeT());
    }

    for (size_t i = 0; i < A_mats.size(); ++i)
    {
        std::cout<<"A: \n"<<A_mats[i]<<"\n";
        std::cout<<"t: \n"<<t_vecs[i]<<"\n";
    }

    delete shape_match;
}

void Morphing::addMeshVec(std::vector<double> &mesh_vec)
{
    meshes_vec.push_back(mesh_vec);
}
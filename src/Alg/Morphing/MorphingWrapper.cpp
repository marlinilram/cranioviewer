#include "MorphingWrapper.h"

MorphingWrapper::MorphingWrapper()
{
    meshes.clear();
    morphing_hander = new Morphing;
}

MorphingWrapper::~MorphingWrapper()
{
    delete morphing_hander;
}

void MorphingWrapper::loadMesh(std::string fName, vtkSmartPointer<vtkRenderer> renderer)
{
    meshes.push_back(Mesh(fName, renderer));

    std::vector<Mesh>::iterator cur_mesh = meshes.end()-1;
    int meshNum = cur_mesh->getMeshData()->GetNumberOfPoints();
    std::vector<double> temp_mesh_vec;
    temp_mesh_vec.clear();

    for (size_t i = 0; i < meshNum; ++i) {
        double* temp = cur_mesh->getMeshData()->GetPoint(i);
        temp_mesh_vec.push_back(temp[0]);
        temp_mesh_vec.push_back(temp[1]);
        temp_mesh_vec.push_back(temp[2]);
    }

    morphing_hander->addMeshVec(temp_mesh_vec);
}
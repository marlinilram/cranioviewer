#include "MorphingWrapper.h"

MorphingWrapper::MorphingWrapper()
{
    meshes.clear();
    morphing_handler = new Morphing;
}

MorphingWrapper::~MorphingWrapper()
{
    delete morphing_handler;
    for (size_t i = 0; i < meshes.size(); ++i)
    {
        delete meshes[i];
    }
}

void MorphingWrapper::loadMesh(std::string fName, vtkSmartPointer<vtkRenderer> renderer)
{
    Mesh *ptr;
    ptr = new Mesh(fName, renderer);
    meshes.push_back(ptr);
    std::vector<Mesh*>::iterator cur_mesh = meshes.end()-1;
    int meshNum = (*cur_mesh)->getMeshData()->GetNumberOfPoints();
    (*cur_mesh)->setMeshCenter();
    double *mesh_center = (*cur_mesh)->getMeshCenter();
    std::vector<double> temp_mesh_vec;
    temp_mesh_vec.clear();

    vtkSmartPointer<vtkPolyData> mesh_data =  (*cur_mesh)->getMeshData();
    vtkSmartPointer<vtkPoints> new_points = mesh_data->GetPoints();

    for (size_t i = 0; i < meshNum; ++i) {
        double* temp = (*cur_mesh)->getMeshData()->GetPoint(i);
        temp_mesh_vec.push_back(temp[0]-mesh_center[0]);
        temp_mesh_vec.push_back(temp[1]-mesh_center[1]);
        temp_mesh_vec.push_back(temp[2]-mesh_center[2]);
        new_points->SetPoint(i, temp_mesh_vec[3*i], temp_mesh_vec[3*i+1], temp_mesh_vec[3*i+2]);
    }
    new_points->Modified();
    (*cur_mesh)->resetMesh(mesh_data);
    morphing_handler->addMeshVec(temp_mesh_vec);
}

void MorphingWrapper::doMorphing(double *paras, int n_paras)
{
    if (n_paras != meshes.size())
    {
        std::cout<<"Warning: num of control parameters isn't match with internal number of meshes\n";
        return;
    }

    vtkSmartPointer<vtkMatrix4x4> mat = vtkSmartPointer<vtkMatrix4x4>::New();
    mat->Zero();
    for (size_t i_mesh = 0; i_mesh < meshes.size(); ++i_mesh)
    {
        double A_element;
        for (size_t i = 0; i < 3; ++i)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                A_element = mat->GetElement(i, j);
                mat->SetElement(i, j, A_element + paras[i_mesh]*morphing_handler->getAMat(i_mesh)(i, j));
            }
        }
    }
    mat->SetElement(3, 3, 1.0);

    meshes[0]->getTransform()->Identity();
    meshes[0]->getTransform()->SetMatrix(mat);
    meshes[0]->reloadTransform();
}

void MorphingWrapper::setCenterMesh(size_t i_mesh)
{
    for (size_t i = 0; i < meshes.size(); ++i)
    {
        if (i != i_mesh)
            meshes[i]->getActor()->SetVisibility(0);
    }
}
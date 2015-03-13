#ifndef Mesh_H
#define Mesh_H

#include <vtkOBJReader.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkCenterOfMass.h>
#include <string>
#include <ctime>

class Mesh
{
public:
    Mesh();
    Mesh(std::string fName);
    Mesh(std::string fName, vtkSmartPointer<vtkRenderer> mainWin_renderer);
    ~Mesh();

    void setMeshData(std::string fName);
    void setRenderer(vtkSmartPointer<vtkRenderer> mainWin_renderer);
    void updateTransform(double *tr_vals);
    void reloadTransform();
    vtkSmartPointer<vtkPolyData> getMeshData() { return mapper->GetInput(); };
    vtkSmartPointer<vtkTransform> getTransform() { return widget_transform; };
    void resetMesh(vtkSmartPointer<vtkPolyData> new_mesh_data);
    double *getMeshCenter() { return mesh_center; };
    void setMeshCenter(double new_center[3] = nullptr);

private:
    void computeCenter(double bounds[6]);

private:
    vtkSmartPointer<vtkPolyDataMapper> mapper;
    vtkSmartPointer<vtkActor> actor;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkTransform> widget_transform;
    vtkSmartPointer<vtkMatrix4x4> widget_transform_mat;

    double mesh_center[3];
};

#endif
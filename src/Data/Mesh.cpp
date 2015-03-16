#include "Mesh.h"

Mesh::Mesh()
{
    mapper = nullptr;
    actor = nullptr;
    renderer = nullptr;
    widget_transform = nullptr;
    widget_transform_mat = nullptr;
    kdtree_picker = nullptr;
}

Mesh::~Mesh()
{
    if (actor)
    {
        renderer->RemoveActor(actor);
        std::cout<<"remove mesh actor\n";
    }
}

Mesh::Mesh(std::string fName)
{
    setMeshData(fName);
}

Mesh::Mesh(std::string fName, vtkSmartPointer<vtkRenderer> mainWin_renderer)
{
    setMeshData(fName);
    setRenderer(mainWin_renderer);
}

void Mesh::setMeshData(std::string fName)
{
    vtkSmartPointer<vtkOBJReader> obj_reader = vtkSmartPointer< vtkOBJReader>::New();

    obj_reader->SetFileName(fName.c_str());
    obj_reader->Update();

    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(obj_reader->GetOutput());

    actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    double bounds[6];
    actor->GetBounds(bounds);
    computeCenter(bounds);

    widget_transform = vtkSmartPointer<vtkTransform>::New();
    widget_transform->PostMultiply();

    widget_transform_mat = vtkSmartPointer<vtkMatrix4x4>::New();
    widget_transform_mat->Identity();
}

void Mesh::computeCenter(double bounds[6])
{
    mesh_center[0] = (bounds[0] + bounds[1]) / 2;
    mesh_center[1] = (bounds[2] + bounds[3]) / 2;
    mesh_center[2] = (bounds[4] + bounds[5]) / 2;
    std::cout<<mesh_center[0]<<"\t"<<mesh_center[1]<<"\t"<<mesh_center[2]<<"\n";
}

void Mesh::setRenderer(vtkSmartPointer<vtkRenderer> mainWin_renderer)
{
    renderer = mainWin_renderer;
    renderer->AddActor(actor);

    renderer->ResetCamera();
}

void Mesh::updateTransform(double *tr_vals)
{
    widget_transform->Identity();
    widget_transform->Translate(-mesh_center[0], -mesh_center[1], -mesh_center[2]);
    widget_transform->RotateX(tr_vals[3]);
    widget_transform->RotateY(tr_vals[4]);
    widget_transform->RotateZ(tr_vals[5]);
    widget_transform->Scale(tr_vals[6], tr_vals[6], tr_vals[6]);
    widget_transform->Translate(mesh_center[0], mesh_center[1], mesh_center[2]);
    widget_transform->Translate(tr_vals[0], tr_vals[1], tr_vals[2]);
    widget_transform->GetMatrix(widget_transform_mat);

    actor->SetUserTransform(widget_transform);
}

void Mesh::reloadTransform()
{
    actor->SetUserTransform(widget_transform);

    //vtkSmartPointer<vtkMatrix4x4> mat = vtkSmartPointer<vtkMatrix4x4>::New();
    //widget_transform->GetMatrix(mat);
    //mat->Print(std::cout);
}

void Mesh::resetMesh(vtkSmartPointer<vtkPolyData> new_mesh_data)
{
    renderer->RemoveActor(actor);
    mapper->SetInputData(new_mesh_data);
    actor->SetMapper(mapper);

    renderer->AddActor(actor);
}

void Mesh::setMeshCenter(double new_center[3])
{
    vtkSmartPointer<vtkCenterOfMass> c_mass = vtkSmartPointer<vtkCenterOfMass>::New();

    c_mass->SetInputData(mapper->GetInput());
    c_mass->SetUseScalarsAsWeights(false);
    c_mass->Update();
    c_mass->GetCenter(mesh_center);
    //mesh_center[0] = new_center[0];
    //mesh_center[1] = new_center[1];
    //mesh_center[2] = new_center[2];
}

void Mesh::setKDTreeLocator()
{
    kdtree_picker = vtkSmartPointer<vtkKdTreePointLocator>::New();
    kdtree_picker->SetDataSet(mapper->GetInput());
    kdtree_picker->BuildLocator();
}

int Mesh::getClosestPtID(double world_coord[3])
{
    if (kdtree_picker)
    {
        int vertex_id = -1;
        vertex_id = kdtree_picker->FindClosestPoint(world_coord);
        double *vertex_coord;
        vertex_coord = mapper->GetInput()->GetPoint(vertex_id);

        //Create a sphere
        vtkSmartPointer<vtkSphereSource> sphereSource =
            vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(vertex_coord[0], vertex_coord[1], vertex_coord[2]);
        sphereSource->SetRadius(1.0);

        //Create a mapper and actor
        vtkSmartPointer<vtkPolyDataMapper> pick_mapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
        pick_mapper->SetInputConnection(sphereSource->GetOutputPort());

        vtkSmartPointer<vtkActor> pick_actor =
            vtkSmartPointer<vtkActor>::New();
        pick_actor->SetMapper(pick_mapper);
        pick_actor->GetProperty()->SetDiffuseColor(0.0, 1.0, 1.0);

        renderer->AddActor(pick_actor);

        return vertex_id;
    }

    return -1;
}
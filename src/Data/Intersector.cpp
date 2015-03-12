#include "Intersector.h"

Intersector::Intersector()
{
    mapper = nullptr;
    actor = nullptr;
    renderer = nullptr;
    transform_filter = nullptr;
    cutter = nullptr;
    stripper = nullptr;
    mesh = nullptr;
    image_slice = nullptr;
}

Intersector::~Intersector()
{
    if (actor)
    {
        renderer->RemoveActor(actor);
        std::cout<<"remove intersection actor\n";
    }
}

void Intersector::setImageSlice(ImageSlice *slice)
{
    image_slice = slice;

    renderer = image_slice->getRenderer();
}

void Intersector::setMesh(Mesh *cur_mesh)
{
    mesh = cur_mesh;
}

void Intersector::setRenderer(vtkSmartPointer<vtkRenderer> mainWin_renderer)
{
    renderer = mainWin_renderer;
}

void Intersector::setCutter()
{
    transform_filter = vtkSmartPointer< vtkTransformPolyDataFilter >::New();
    transform_filter->SetInputData(mesh->getMeshData());
    transform_filter->SetTransform(mesh->getTransform());

    cutter = vtkSmartPointer< vtkCutter >::New();
    cutter->SetInputConnection(transform_filter->GetOutputPort());
    cutter->SetCutFunction(image_slice->getSliceMapper()->GetSlicePlane());

    stripper = vtkSmartPointer< vtkStripper >::New();
    stripper->SetInputConnection(cutter->GetOutputPort());
    stripper->Update();

    mapper = vtkSmartPointer< vtkPolyDataMapper >::New();
    mapper->SetInputConnection(stripper->GetOutputPort());

    actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetDiffuseColor(0.9, 0.447, 0.0745);

    renderer->AddActor(actor);
}
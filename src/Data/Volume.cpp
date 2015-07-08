#include "Volume.h"

Volume::Volume()
{
    mapper = nullptr;
    actor = nullptr;
    renderer = nullptr;
    polydata = nullptr;
    vertex_glyph_filter = nullptr;
    nii_img = nullptr;

    isovalue = 0;
    isowidth = 0;
}

Volume::~Volume()
{
    if (actor)
    {
        renderer->RemoveActor(actor);
        std::cout<<"remove image volume actor\n";
    }
}

void Volume::setVolumeImg(NiiLoader *img)
{
    nii_img = img;
}

void Volume::setRenderer(vtkSmartPointer<vtkRenderer> mainWin_renderer)
{
    renderer = mainWin_renderer;
}

void Volume::dispTargetVolume()
{
    int voxel_num = 0;
    std::vector<double> voxel_data;// = nii_img->extractSkullVertex(isovalue, isowidth, voxel_num);
    std::vector<double> voxel_norm;
    nii_img->mcSkullVertex(voxel_data, voxel_norm, isovalue);

    std::cout<<"isovalue: "<<isovalue<<"\tisowidth: "<<isowidth<<"\tvoxel num: "<<voxel_num<<"\n";

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (size_t i = 0; i < voxel_data.size()/3; ++i)
    {
        points->InsertNextPoint(voxel_data[3*i + 0], voxel_data[3*i + 1], voxel_data[3*i + 2]);
    }

    if (actor) renderer->RemoveActor(actor);

    polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);

    vertex_glyph_filter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    vertex_glyph_filter->SetInputData(polydata);
    vertex_glyph_filter->Update();

    mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(vertex_glyph_filter->GetOutput());

    actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(2);
    actor->GetProperty()->SetDiffuseColor(0,1,1);

    renderer->AddActor(actor);
}

void Volume::noDispTargetVolume()
{
    if (actor) 
   {
        renderer->RemoveActor(actor);
        actor = nullptr;
    }
}

void Volume::setISOVal(int value)
{
    isovalue = value;
}

void Volume::setISOWidth(int value)
{
    isowidth = value;
}
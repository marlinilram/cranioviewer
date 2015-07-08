#ifndef Intersector_H
#define Intersector_H

#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkCutter.h>
#include <vtkStripper.h>
#include <vtkImageSlice.h>
#include <vtkTransform.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkPlane.h>
#include <vtkImplicitFunction.h>
#include <vtkPropPicker.h>
#include <vtkWorldPointPicker.h>
#include <vtkSphereSource.h>

#include "ImageSlice.h"
#include "Mesh.h"

class Intersector
{
public:
    Intersector();
    ~Intersector();

    void setRenderer(vtkSmartPointer<vtkRenderer> mainWin_renderer);
    void setImageSlice(ImageSlice *slice);
    void setMesh(Mesh *cur_mesh);
    void setCutter();
    void getWorldCoord(double x, double y, double coord[3] = nullptr);
    vtkSmartPointer<vtkRenderer> getRenderer() { return renderer; };
    void setThickness(int value);

private:
    vtkSmartPointer<vtkPolyDataMapper> mapper;
    vtkSmartPointer<vtkActor> actor;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkTransformPolyDataFilter> transform_filter;
    vtkSmartPointer<vtkCutter> cutter;
    vtkSmartPointer<vtkStripper> stripper;
    Mesh *mesh;
    ImageSlice *image_slice;
};

#endif
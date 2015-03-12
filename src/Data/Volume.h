#ifndef Volume_H
#define Volume_H

#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPoints.h>

#include "NiiLoader.h"

class Volume
{
public:
    Volume();
    ~Volume();

    void setVolumeImg(NiiLoader *img);
    void setRenderer(vtkSmartPointer<vtkRenderer> mainWin_renderer);
    void dispTargetVolume();
    void noDispTargetVolume();
    void setISOVal(int value);
    void setISOWidth(int value);
    int getISOVal() { return isovalue; };
    int getISOWidth() { return isowidth; };

private:
    vtkSmartPointer<vtkPolyDataMapper> mapper;
    vtkSmartPointer<vtkActor> actor;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkPolyData> polydata;
    vtkSmartPointer<vtkVertexGlyphFilter> vertex_glyph_filter;

    NiiLoader *nii_img;
    int isovalue;
    int isowidth;
};

#endif
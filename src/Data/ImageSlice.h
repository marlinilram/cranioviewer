#ifndef ImageSlice_H
#define ImageSlice_H

#include <string>

#include <vtkSmartPointer.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkPropPicker.h>
#include <vtkMath.h>
#include <vtkImageProperty.h>
#include <vtkLookupTable.h>

#include "NiiLoader.h"

class ImageSlice
{
public:
    ImageSlice();
    ImageSlice(vtkSmartPointer<vtkImageData> data, std::string orient, vtkSmartPointer<vtkRenderer> mainWin_renderer, bool img);
    ~ImageSlice();
    vtkSmartPointer<vtkImageSlice> getActor();
    void setRenderer(vtkSmartPointer<vtkRenderer> mainWin_renderer);
    void setSliceNum(float slice_pos);
    void setParallel(bool on);
    void setImgCamera();
    void setImgColorWin(double win);
    void setImgColorLev(double lev);
    void setImgColorLUT(double range[2]);
    double getPixelVal(double x, double y);
    vtkSmartPointer<vtkRenderer> getRenderer() { return renderer; };
    vtkSmartPointer<vtkImageSliceMapper> getSliceMapper() { return mapper; };

private:
    void setImgData(vtkSmartPointer<vtkImageData> data);
    void setOrient(std::string orient);

private:
    vtkSmartPointer<vtkImageData> nii_data;

    vtkSmartPointer<vtkImageSliceMapper> mapper;
    vtkSmartPointer<vtkImageSlice> actor;
    vtkSmartPointer<vtkRenderer> renderer;

    int max_slice_num;
    std::string img_orient;
    bool parallel_on;
};

#endif
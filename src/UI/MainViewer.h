#ifndef MainViewer_H
#define MainViewer_H

#include <QObject>
#include <QSlider>
#include <vtkXMLImageDataWriter.h>
#include "ImageSlice.h"
#include "Mesh.h"
#include "Volume.h"


class MainViewer : public QObject
{
    Q_OBJECT

public:
    MainViewer();
    ~MainViewer();

    void updateSlice(int orient, float slice_pos);
    void setSliders(QSlider *mainWin_sliders[3]);
    void setImgData(NiiLoader *data);
    void setImgData(std::string fName);
    void setImgData(vtkSmartPointer<vtkImageData> data);
    void setMeshData(std::string fName);
    void setRenderer(vtkSmartPointer<vtkRenderer> mainWin_renderer);
    NiiLoader *getImgData() { return nii_img; };
    Mesh *getMeshData() { return temp_mesh; };
    Volume *getVolumeData() { return img_volume; };
    void clearImage();
    void clearData();
    void clearMesh();
    void saveMesh(std::string f_name);
    void saveImg(std::string f_name);

    void updateVolumeView();
    void runMC(std::string fName);

private slots:
    void updateYZPlaneView(int value);
    void updateXZPlaneView(int value);
    void updateXYPlaneView(int value);
    void updateTransform(double *tr_vals);
    void updateISOVal(int value);
    void updateISOWidth(int value);
    void showVolume(int state);

signals:
    void updateRenderers();

private:                  
    // data
    NiiLoader *nii_img;
    ImageSlice *image_slices[3];
    Volume *img_volume;
    Mesh *temp_mesh;

    // UI widget ptr
    QSlider *sliders[3];

    // viewer handler
    vtkSmartPointer<vtkRenderer> renderer;

    // bool controller
    bool show_volume;
    bool show_mesh;
    bool show_slices;
};

#endif
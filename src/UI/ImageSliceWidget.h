#ifndef ImageSliceWidget_H
#define ImageSliceWidget_H

#include <QObject>
#include <QSlider>

#include "ImageSlice.h"
#include "Intersector.h"

class ImageSliceWidget : public QObject
{
    Q_OBJECT

public:
    ImageSliceWidget();
    ImageSliceWidget(ImageSlice *slice, QSlider *slider);
    ~ImageSliceWidget();

    void setSlider(QSlider *slider);
    void setSlice(vtkSmartPointer<vtkImageData> img_data, std::string orient, vtkSmartPointer<vtkRenderer> plane_renderer);
    void addSlice(vtkSmartPointer<vtkImageData> img_data);
    void setMesh(Mesh *mesh);
    void setIntersector(Mesh *mesh);
    void clearSlice();
    void clearIntersector();
    ImageSlice *getImageSlice() { return image_slice; };
    Intersector *getIntersector() { return intersector; };
    void setSliceLUT(double range[2]);

public slots:
    void updatePlaneView(int value);
    void setSliceColorWin(double win);
    void setSliceColorWin(int win) { setSliceColorWin(static_cast<double>(win)); };
    void setSliceColorLev(double lev);
    void setSliceColorLev(int lev) { setSliceColorLev(static_cast<double>(lev)); };
    void setIntersectionThick(int value);
    void toggleImgDisp();

signals:
    void updateRenderers();

private:
    ImageSlice *image_slice;
    Intersector *intersector;

    QSlider *image_slider;
};

#endif
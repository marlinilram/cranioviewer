#include "ImageSliceWidget.h"

ImageSliceWidget::ImageSliceWidget()
{
    image_slice = nullptr;
    intersector = nullptr;
    image_slider = nullptr;
}

ImageSliceWidget::~ImageSliceWidget()
{
    clearSlice();
    clearIntersector();
    std::cout<<"delete image slice widget finished\n";
}

ImageSliceWidget::ImageSliceWidget(ImageSlice *slice, QSlider *slider)
{
    image_slice = slice;
    image_slider = slider;
}

void ImageSliceWidget::updatePlaneView(int value)
{
    float slice_pos = (float)value/image_slider->maximum();
    if (image_slice)
        image_slice->setSliceNum(slice_pos);
    emit(updateRenderers());
}

void ImageSliceWidget::setSlider(QSlider *slider)
{
    image_slider = slider;
}

void ImageSliceWidget::setSlice(vtkSmartPointer<vtkImageData> img_data, std::string orient, vtkSmartPointer<vtkRenderer> plane_renderer)
{
    image_slice = new ImageSlice(img_data, orient, plane_renderer, true);
}

void ImageSliceWidget::addSlice(vtkSmartPointer<vtkImageData> img_data)
{
  image_slice->addImgData(img_data);
}

void ImageSliceWidget::setMesh(Mesh *mesh)
{
    if (!intersector)
    {
        intersector = new Intersector;

        intersector->setImageSlice(image_slice);
        intersector->setMesh(mesh);
        intersector->setCutter();
    }
    else
    {
        std::cout<<"current intersector ptr isn't null\n";
    }
}

void ImageSliceWidget::clearSlice()
{
    if (image_slice)
    {
        delete image_slice;
        image_slice = nullptr;
        std::cout<<"delete image slice\n";
    }
    else
    {
        std::cout<<"no image slice\n";
    }
}

void ImageSliceWidget::clearIntersector()
{
    if (intersector)
    {
        delete intersector;
        intersector = nullptr;
        std::cout<<"delete intersector\n";
    }
    else
    {
        std::cout<<"no intersector\n";
    }
}

void ImageSliceWidget::setSliceColorWin(double win)
{
  if (image_slice)
  {
    image_slice->setImgColorWin(win);
    emit(updateRenderers());
  }
}

void ImageSliceWidget::setSliceColorLev(double lev)
{
  if (image_slice)
  {
    image_slice->setImgColorLev(lev);
    emit(updateRenderers());
  }
}

void ImageSliceWidget::setSliceLUT(double range[2])
{
    image_slice->setImgColorLUT(range);
}

void ImageSliceWidget::setIntersectionThick(int value)
{
  if (intersector)
  {
    intersector->setThickness(value);
    emit(updateRenderers());
  }
}

void ImageSliceWidget::toggleImgDisp()
{
  if (image_slice)
  {
    image_slice->togImgDisp();
  }
  if (intersector)
  {
    intersector->setCutter();
  }
  emit(updateRenderers());
}
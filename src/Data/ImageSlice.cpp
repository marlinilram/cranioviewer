#include "ImageSlice.h"

ImageSlice::ImageSlice()
{
    renderer = nullptr;

    parallel_on = false;

    cur_image = 0;
}

ImageSlice::~ImageSlice()
{
  if (!actors.empty())
  {
    for (size_t i = 0; i < actors.size(); ++i)
    {
      renderer->RemoveActor(actors[i]);
      std::cout<<"remove " << i << "th image slice actor\n";
    }
  }

}

ImageSlice::ImageSlice(vtkSmartPointer<vtkImageData> data, std::string orient, vtkSmartPointer<vtkRenderer> mainWin_renderer, bool img)
  : renderer(nullptr), parallel_on(false), cur_image(0)
{
    setImgData(data);
    setOrient(orient);
    setRenderer(mainWin_renderer);

    if (img) 
   {
        setImgCamera();
        setParallel(true);
    }
    else 
    {
        setParallel(false);
        renderer->ResetCamera();
    }
}
vtkSmartPointer<vtkImageSlice> ImageSlice::getActor()
{
    return actors[cur_image];
}

void ImageSlice::setRenderer(vtkSmartPointer<vtkRenderer> mainWin_renderer)
{
    renderer = mainWin_renderer;
    renderer->AddViewProp(actors.back());
}

void ImageSlice::setImgData(vtkSmartPointer<vtkImageData> data)
{
    images_data.push_back(data);

    vtkSmartPointer<vtkImageSliceMapper> mapper = vtkSmartPointer<vtkImageSliceMapper>::New();
    mapper->SetInputData(images_data.back());
    mappers.push_back(mapper);
    vtkSmartPointer<vtkImageSlice> actor = vtkSmartPointer<vtkImageSlice>::New();
    actor->SetMapper(mappers.back());
    actors.push_back(actor);
}

void ImageSlice::addImgData(vtkSmartPointer<vtkImageData> data)
{
  images_data.push_back(data);

  vtkSmartPointer<vtkImageSliceMapper> mapper = vtkSmartPointer<vtkImageSliceMapper>::New();
  mapper->SetInputData(images_data.back());
  mappers.push_back(mapper);
  vtkSmartPointer<vtkImageSlice> actor = vtkSmartPointer<vtkImageSlice>::New();
  actor->SetMapper(mappers.back());
  actors.push_back(actor);
  renderer->AddViewProp(actors.back());
  actor->VisibilityOff();

  int dims[3];
  images_data.back()->GetDimensions(dims);
  if (img_orient == "XY")
  {
    mappers.back()->SetOrientationToZ();
    mappers.back()->SetSliceNumber(dims[2]/2);
    max_slice_num = dims[2];
  }
  else if (img_orient == "XZ")
  {
    mappers.back()->SetOrientationToY();
    mappers.back()->SetSliceNumber(dims[1]/2);
    max_slice_num = dims[1];
  }
  else if (img_orient == "YZ")
  {
    mappers.back()->SetOrientationToX();
    mappers.back()->SetSliceNumber(dims[0]/2);
    max_slice_num = dims[0];
    //renderer->GetActiveCamera()->SetPosition(1,0,0);
    //renderer->GetActiveCamera()->SetViewUp(0,0,1);
    //renderer->GetActiveCamera()->SetFocalPoint(0,0,0);
  }
}

void ImageSlice::togImgDisp()
{
  int max_image = actors.size();
  actors[cur_image]->VisibilityOff();
  int cur_slice_num = mappers[cur_image]->GetSliceNumber();
  cur_image = (cur_image + 1) % max_image;
  mappers[cur_image]->SetSliceNumber(cur_slice_num);
  actors[cur_image]->VisibilityOn();
}

void ImageSlice::setOrient(std::string orient)
{
    img_orient = orient;
    int dims[3];
    images_data[cur_image]->GetDimensions(dims);

    //std::cout<<"dims: "<<dims[0]<<"\t"<<dims[1]<<"\t"<<dims[2]<<"\n";

    if (orient == "XY")
    {
        mappers[cur_image]->SetOrientationToZ();
        mappers[cur_image]->SetSliceNumber(dims[2]/2);
        max_slice_num = dims[2];
    }
    else if (orient == "XZ")
    {
        mappers[cur_image]->SetOrientationToY();
        mappers[cur_image]->SetSliceNumber(dims[1]/2);
        max_slice_num = dims[1];
    }
    else if (orient == "YZ")
    {
        mappers[cur_image]->SetOrientationToX();
        mappers[cur_image]->SetSliceNumber(dims[0]/2);
        max_slice_num = dims[0];
    }
}

void ImageSlice::setSliceNum(float slice_pos)
{
    mappers[cur_image]->SetSliceNumber(slice_pos*max_slice_num);
    //std::cout<< slice_pos << "\t" << max_slice_num << "\n";
    renderer->Render();
}

void ImageSlice::setParallel(bool on)
{
    parallel_on = on;
    if (parallel_on)
        renderer->GetActiveCamera()->SetParallelProjection(1);
    else
        renderer->GetActiveCamera()->SetParallelProjection(0);
}

void ImageSlice::setImgCamera()
{
  int dims[3];
  images_data[cur_image]->GetDimensions(dims);
  std::cout<<"dims: "<<dims[0]<<"\t"<<dims[1]<<"\t"<<dims[2]<<"\n";
  double spaces[3];
  images_data[cur_image]->GetSpacing(spaces);
  std::cout<<"spaces: "<<spaces[0]<<"\t"<<spaces[1]<<"\t"<<spaces[2]<<"\n";
  std::cout<<"max slice number: "<<max_slice_num<<"\n";

    if (img_orient == "YZ") 
    {
      renderer->GetActiveCamera()->SetPosition(spaces[0]*(max_slice_num+1), spaces[1]*dims[1]/2, spaces[2]*dims[2]/2);
      renderer->GetActiveCamera()->SetViewUp(0,0,1);
      renderer->GetActiveCamera()->SetFocalPoint(0, spaces[1]*dims[1]/2, spaces[2]*dims[2]/2);
      renderer->GetActiveCamera()->SetParallelScale(spaces[2]*dims[2]);
    }
    else if (img_orient == "XZ")
    {
      renderer->GetActiveCamera()->SetPosition(spaces[0]*dims[0]/2,spaces[1]*(max_slice_num+1), spaces[2]*dims[2]/2);
      renderer->GetActiveCamera()->SetViewUp(0,0,1);
      renderer->GetActiveCamera()->SetFocalPoint(spaces[0]*dims[0]/2,0, spaces[2]*dims[2]/2);
      renderer->GetActiveCamera()->SetParallelScale(spaces[2]*dims[2]);
    }
    else if (img_orient == "XY")
    {
        renderer->GetActiveCamera()->SetPosition(spaces[0]*dims[0]/2,spaces[1]*dims[1]/2,spaces[2]*(max_slice_num+1));
        renderer->GetActiveCamera()->SetViewUp(0,1,0);
        renderer->GetActiveCamera()->SetFocalPoint(spaces[0]*dims[0]/2,spaces[1]*dims[1]/2,0);
        renderer->GetActiveCamera()->SetParallelScale(spaces[1]*dims[1]);
    }
    renderer->GetActiveCamera()->SetParallelProjection(1);
    //renderer->ResetCamera();
}

double ImageSlice::getPixelVal(double x, double y)
{
    // x y is actor coord
    vtkSmartPointer<vtkPropPicker> picker = vtkSmartPointer<vtkPropPicker>::New();
    picker->Pick(x, y, 0.0, renderer);

    double pos[3];
    picker->GetPickPosition(pos);
    //std::cout<<"picked position: "<<pos[0]<<"\t"<<pos[1]<<"\t"<<pos[2]<<"\n";

    //double bounds[6];
    //actor->GetBounds(bounds);
    //std::cout<<"actor bounds: ";
    //for (size_t i = 0; i < 6; ++i)
    //{
    //    std::cout<<bounds[i]<<"\t";
    //}
    //std::cout<<"\n";

    // X 0 Y 1 Z 2
    int img_coord[3] = {0,0,0};
    int img_extent[6] = {0,0,0,0,0,0};
    int orient = mappers[cur_image]->GetOrientation();
    double img_spcing[3];
    images_data[cur_image]->GetSpacing(img_spcing);
    images_data[cur_image]->GetExtent(img_extent);
    if (orient == 0)
    {
        img_coord[0] = mappers[cur_image]->GetSliceNumber();
        img_coord[1] = vtkMath::Round(pos[1]/img_spcing[1]);
        img_coord[2] = vtkMath::Round(pos[2]/img_spcing[2]);
    }
    else if (orient == 1)
    {
        img_coord[0] = vtkMath::Round(pos[0]/img_spcing[0]);
        img_coord[1] = mappers[cur_image]->GetSliceNumber();
        img_coord[2] = vtkMath::Round(pos[2]/img_spcing[2]);
    }
    else if (orient == 2)
    {
        img_coord[0] = vtkMath::Round(pos[0]/img_spcing[0]);
        img_coord[1] = vtkMath::Round(pos[1]/img_spcing[1]);
        img_coord[2] = mappers[cur_image]->GetSliceNumber();
    }
#ifdef DEBUG_PRINT
    std::string data_type(nii_data->GetScalarTypeAsString());
    std::cout<<"image data type: "<<data_type<<"\n";
#endif
    if (img_coord[0] >= img_extent[0] && img_coord[0] <= img_extent[1] 
    && img_coord[1] >= img_extent[2] && img_coord[1] <= img_extent[3]
    && img_coord[2] >= img_extent[4] && img_coord[2] <= img_extent[5])
    {
            return images_data[cur_image]->GetScalarComponentAsDouble(img_coord[0], img_coord[1], img_coord[2], 0);
    }
    
    return 0.0;
}

void ImageSlice::setImgColorWin(double win)
{
  double range[2];
  images_data[cur_image]->GetScalarRange(range);
  double n_win = (range[1] - range[0]) * win / 100;
  actors[cur_image]->GetProperty()->SetColorWindow(n_win);
  std::cout<<n_win<<"\n";
}

void ImageSlice::setImgColorLev(double lev)
{
  double range[2];
  images_data[cur_image]->GetScalarRange(range);
  double n_lev = range[0] + (range[1] - range[0]) * lev / 100;
  actors[cur_image]->GetProperty()->SetColorLevel(n_lev);
  std::cout<<n_lev<<"\n";
}

void ImageSlice::setImgColorLUT(double range[2])
{
    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetNumberOfTableValues(256);
    lut->SetRange(range[0], range[1]);
    lut->Build();

    actors[cur_image]->GetProperty()->SetLookupTable(lut);
}

void ImageSlice::setVisible(int state)
{
  if (state)
  {
    actors[cur_image]->VisibilityOn();
  }
  else
  {
    actors[cur_image]->VisibilityOff();
  }
}
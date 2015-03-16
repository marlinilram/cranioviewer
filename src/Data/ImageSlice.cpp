#include "ImageSlice.h"

ImageSlice::ImageSlice()
{
    nii_data = nullptr;

    mapper = nullptr;
    actor = nullptr;
    renderer = nullptr;

    parallel_on = false;
}

ImageSlice::~ImageSlice()
{
     if (actor)
     {
         renderer->RemoveActor(actor);
         std::cout<<"remove image slice actor\n";
     }
}

ImageSlice::ImageSlice(vtkSmartPointer<vtkImageData> data, std::string orient, vtkSmartPointer<vtkRenderer> mainWin_renderer, bool img)
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
    return actor;
}

void ImageSlice::setRenderer(vtkSmartPointer<vtkRenderer> mainWin_renderer)
{
    renderer = mainWin_renderer;
    renderer->AddViewProp(actor);
}

void ImageSlice::setImgData(vtkSmartPointer<vtkImageData> data)
{
    nii_data = data;

    mapper = vtkSmartPointer<vtkImageSliceMapper>::New();
    mapper->SetInputData(nii_data);
    actor = vtkSmartPointer<vtkImageSlice>::New();
    actor->SetMapper(mapper);
}

void ImageSlice::setOrient(std::string orient)
{
    img_orient = orient;
    int dims[3];
    nii_data->GetDimensions(dims);

    //std::cout<<"dims: "<<dims[0]<<"\t"<<dims[1]<<"\t"<<dims[2]<<"\n";

    if (orient == "XY")
    {
        mapper->SetOrientationToZ();
        mapper->SetSliceNumber(dims[2]/2);
        max_slice_num = dims[2];
    }
    else if (orient == "XZ")
    {
        mapper->SetOrientationToY();
        mapper->SetSliceNumber(dims[1]/2);
        max_slice_num = dims[1];
    }
    else if (orient == "YZ")
    {
        mapper->SetOrientationToX();
        mapper->SetSliceNumber(dims[0]/2);
        max_slice_num = dims[0];
    }
}

void ImageSlice::setSliceNum(float slice_pos)
{
    mapper->SetSliceNumber(slice_pos*max_slice_num);
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
    if (img_orient == "YZ") 
    {
        renderer->GetActiveCamera()->SetPosition(1,0,0);
        renderer->GetActiveCamera()->SetViewUp(0, 0,1);
    }
    else if (img_orient == "XZ")
    {
        renderer->GetActiveCamera()->SetPosition(0,1,0);
        renderer->GetActiveCamera()->SetViewUp(0,0,1);
    }
    else if (img_orient == "XY")
    {
        renderer->GetActiveCamera()->SetPosition(0,0,1);
        renderer->GetActiveCamera()->SetViewUp(0,1,0);
    }

    renderer->GetActiveCamera()->SetParallelProjection(1);
    renderer->ResetCamera();
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
    int orient = mapper->GetOrientation();
    double img_spcing[3];
    nii_data->GetSpacing(img_spcing);
    nii_data->GetExtent(img_extent);
    if (orient == 0)
    {
        img_coord[0] = mapper->GetSliceNumber();
        img_coord[1] = vtkMath::Round(pos[1]/img_spcing[1]);
        img_coord[2] = vtkMath::Round(pos[2]/img_spcing[2]);
    }
    else if (orient == 1)
    {
        img_coord[0] = vtkMath::Round(pos[0]/img_spcing[0]);
        img_coord[1] = mapper->GetSliceNumber();
        img_coord[2] = vtkMath::Round(pos[2]/img_spcing[2]);
    }
    else if (orient == 2)
    {
        img_coord[0] = vtkMath::Round(pos[0]/img_spcing[0]);
        img_coord[1] = vtkMath::Round(pos[1]/img_spcing[1]);
        img_coord[2] = mapper->GetSliceNumber();
    }
#ifdef DEBUG_PRINT
    std::string data_type(nii_data->GetScalarTypeAsString());
    std::cout<<"image data type: "<<data_type<<"\n";
#endif
    if (img_coord[0] >= img_extent[0] && img_coord[0] <= img_extent[1] 
    && img_coord[1] >= img_extent[2] && img_coord[1] <= img_extent[3]
    && img_coord[2] >= img_extent[4] && img_coord[2] <= img_extent[5])
    {
            return nii_data->GetScalarComponentAsDouble(img_coord[0], img_coord[1], img_coord[2], 0);
    }
    
    return 0.0;
}

void ImageSlice::setImgColorWin(double win)
{
    actor->GetProperty()->SetColorWindow(win);
}

void ImageSlice::setImgColorLev(double lev)
{
    actor->GetProperty()->SetColorLevel(lev);
}

void ImageSlice::setImgColorLUT(double range[2])
{
    vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetNumberOfTableValues(256);
    lut->SetRange(range[0], range[1]);
    lut->Build();

    actor->GetProperty()->SetLookupTable(lut);
}
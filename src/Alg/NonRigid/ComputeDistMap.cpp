#include "ComputeDistMap.h"

ComputeDistMap::ComputeDistMap()
{
    bone_img = nullptr;
    dist_map = nullptr;
}

ComputeDistMap::~ComputeDistMap()
{

}

void ComputeDistMap::runFastMarching()
{

    FastMarchingType::Pointer marcher = FastMarchingType::New();
    marcher->SetStoppingCriterion( criterion );

    // specify the size of the output image
    int *dims = bone_img->GetDimensions();
    FloatImageType::SizeType size = {{dims[0],dims[1],dims[2]}};
    marcher->SetOutputSize( size );

    // setup a 'alive image'
    FloatImageType::Pointer AliveImage = FloatImageType::New();
    AliveImage->SetLargestPossibleRegion( region );
    AliveImage->SetBufferedRegion( region );
    AliveImage->Allocate();
    AliveImage->FillBuffer( 0.0 );
    //FloatImageType::IndexType index;
    //index[0] = 10;
    //index[1] = 10;
    //index[2] = 10;
    //AliveImage->SetPixel(index, 1.0);
    //index[0] = 50;
    //index[1] = 50;
    //index[2] = 50;
    //AliveImage->SetPixel(index, 1.0);
    setActiveImg(AliveImage, true);
    // setup a 'trial image'
    FloatImageType::Pointer TrialImage = FloatImageType::New();
    TrialImage->SetLargestPossibleRegion( region );
    TrialImage->SetBufferedRegion( region );
    TrialImage->Allocate();
    TrialImage->FillBuffer( 0.0 );
    setTrialImg(TrialImage, AliveImage);

    marcher->SetInput( ITK_NULLPTR );
    marcher->SetSpeedConstant(1.0);
    marcher->SetOutputRegion(region);
    double *ori_spacing = bone_img->GetSpacing();
    FloatImageType::SpacingType output_spacing;
    output_spacing[0] = ori_spacing[0];
    output_spacing[1] = ori_spacing[1];
    output_spacing[2] = ori_spacing[2];
    marcher->SetOutputSpacing(output_spacing);

    AdaptorType::Pointer adaptor = AdaptorType::New();
    adaptor->SetAliveImage( AliveImage.GetPointer() );
    adaptor->SetAliveValue( 1.0 );
    adaptor->SetTrialImage( TrialImage.GetPointer() );
    adaptor->SetTrialValue( 1.0 );
    adaptor->Update();
    marcher->SetAlivePoints( adaptor->GetAlivePoints() );
    marcher->SetTrialPoints( adaptor->GetTrialPoints() );

    marcher->Update();


    ConnectorType::Pointer connector = ConnectorType::New();
    connector->SetInput(marcher->GetOutput());
    connector->Update();
    dist_map = vtkSmartPointer<vtkImageData>::New();
    dist_map->DeepCopy(connector->GetOutput());
    double *img_spacing = bone_img->GetSpacing();
    //dist_map->SetSpacing(img_spacing[0], img_spacing[1], img_spacing[2]);

    std::cout<<"compute dist map finished\n";

#ifdef DEBUG_PRINT
    std::cout<<"input image dims: "<<dims[0]<<"\t"<<dims[1]<<"\t"<<dims[2]<<"\n";
#endif
    dims = dist_map->GetDimensions();
#ifdef DEBUG_PRINT
    std::cout<<"output image dims: "<<dims[0]<<"\t"<<dims[1]<<"\t"<<dims[2]<<"\n";
#endif             

    //vtkSmartPointer<vtkImageActor> itk_img_actor = vtkSmartPointer<vtkImageActor>::New();
    //vtkSmartPointer<vtkRenderer> itk_img_ren = vtkSmartPointer<vtkRenderer>::New();
    //itk_img_ren->AddActor(itk_img_actor);
    //vtkSmartPointer<vtkRenderWindow> itk_img_renWin = vtkSmartPointer<vtkRenderWindow>::New();
    //itk_img_renWin->AddRenderer(itk_img_ren);

    //vtkSmartPointer<vtkImageViewer> itk_img_viewer = vtkSmartPointer<vtkImageViewer>::New();
    //vtkSmartPointer<vtkRenderWindowInteractor> itk_img_iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    //itk_img_iren->SetRenderWindow(itk_img_renWin);
    //itk_img_viewer->SetInputData(dist_map);
    //itk_img_viewer->SetupInteractor(itk_img_iren);
    //itk_img_viewer->Render();
    //itk_img_iren->Initialize();
    //itk_img_iren->Start();
}

vtkSmartPointer<vtkImageData> ComputeDistMap::computeOuterDist()
{
    FastMarchingType::Pointer marcher = FastMarchingType::New();
    marcher->SetStoppingCriterion( criterion );
    marcher->SetOutputSize( size );
    marcher->SetInput( ITK_NULLPTR );
    marcher->SetSpeedConstant(1.0);
    marcher->SetOutputRegion(region);
    marcher->SetOutputSpacing(output_spacing);

    // setup a 'alive image'
    FloatImageType::Pointer AliveImage = FloatImageType::New();
    AliveImage->SetLargestPossibleRegion( region );
    AliveImage->SetBufferedRegion( region );
    AliveImage->SetSpacing(output_spacing);
    AliveImage->Allocate();
    AliveImage->FillBuffer( 0.0 );
    //FloatImageType::IndexType index;
    //index[0] = 5;
    //index[1] = 5;
    //index[2] = 5;
    //AliveImage->SetPixel(index, 1.0);
    //index[0] = 495;
    //index[1] = 495;
    //index[2] = 495;
    //AliveImage->SetPixel(index, 1.0);
    setActiveImg(AliveImage, true);
    // setup a 'trial image'
    FloatImageType::Pointer TrialImage = FloatImageType::New();
    TrialImage->SetLargestPossibleRegion( region );
    TrialImage->SetBufferedRegion( region );
    TrialImage->SetSpacing(output_spacing);
    TrialImage->Allocate();
    TrialImage->FillBuffer( 0.0 );
    setTrialImg(TrialImage, AliveImage);

    AdaptorType::Pointer adaptor = AdaptorType::New();
    adaptor->SetAliveImage( AliveImage.GetPointer() );
    adaptor->SetAliveValue( 1.0 );
    adaptor->SetTrialImage( TrialImage.GetPointer() );
    adaptor->SetTrialValue( 1.0 );
    adaptor->Update();
    marcher->SetAlivePoints( adaptor->GetAlivePoints() );
    marcher->SetTrialPoints( adaptor->GetTrialPoints() );

    std::clock_t begin = std::clock();
    marcher->Update();
    std::clock_t end = std::clock();
    std::cout<<"ITK fast marching filter finished. Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\n";

    ConnectorType::Pointer connector = ConnectorType::New();
    connector->SetInput(marcher->GetOutput());
    connector->Update();
    vtkSmartPointer<vtkImageData> out_disp_map = connector->GetOutput();
    return out_disp_map;
}

vtkSmartPointer<vtkImageData> ComputeDistMap::computeInnerDist()
{
    FastMarchingType::Pointer marcher = FastMarchingType::New();
    marcher->SetStoppingCriterion( criterion );
    marcher->SetOutputSize( size );
    marcher->SetInput( ITK_NULLPTR );
    marcher->SetSpeedConstant(1.0);
    marcher->SetOutputRegion(region);
    marcher->SetOutputSpacing(output_spacing);

    // setup a 'alive image'
    FloatImageType::Pointer AliveImage = FloatImageType::New();
    AliveImage->SetLargestPossibleRegion( region );
    AliveImage->SetBufferedRegion( region );
    AliveImage->SetSpacing(output_spacing);
    AliveImage->Allocate();
    AliveImage->FillBuffer( 1.0 );
    setActiveImg(AliveImage, false);
    // setup a 'trial image'
    FloatImageType::Pointer TrialImage = FloatImageType::New();
    TrialImage->SetLargestPossibleRegion( region );
    TrialImage->SetBufferedRegion( region );
    TrialImage->SetSpacing(output_spacing);
    TrialImage->Allocate();
    TrialImage->FillBuffer( 0.0 );
    setTrialImg(TrialImage, AliveImage);

    AdaptorType::Pointer adaptor = AdaptorType::New();
    adaptor->SetAliveImage( AliveImage.GetPointer() );
    adaptor->SetAliveValue( 1.0 );
    adaptor->SetTrialImage( TrialImage.GetPointer() );
    adaptor->SetTrialValue( 1.0 );
    adaptor->Update();
    marcher->SetAlivePoints( adaptor->GetAlivePoints() );
    marcher->SetTrialPoints( adaptor->GetTrialPoints() );

    marcher->Update();

    ConnectorType::Pointer connector = ConnectorType::New();
    connector->SetInput(marcher->GetOutput());
    connector->Update();
    vtkSmartPointer<vtkImageData> in_disp_map = connector->GetOutput();
    return in_disp_map;
}

void ComputeDistMap::computeFinalDistMap()
{
    std::cout<<"compute outward dist map\n";
    std::clock_t begin;
    std::clock_t end;
    begin = std::clock();
    vtkSmartPointer<vtkImageData> outside = vtkSmartPointer<vtkImageData>::New();
    outside->DeepCopy(computeOuterDist());
    end = std::clock();
    std::cout<<"Compute outward dist map finished. Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\n";

    std::cout<<"compute inward dist map\n";
    begin = std::clock();
    vtkSmartPointer<vtkImageData> inside = vtkSmartPointer<vtkImageData>::New();
    inside->DeepCopy(computeInnerDist());
    end = std::clock();
    std::cout<<"Compute outward dist map finished. Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\n";

    std::cout<<"combine inward dist map and outward dist map\n";
    begin = std::clock();
    int *dims = bone_img->GetDimensions();
    float *out_ptr = static_cast<float *>(outside->GetScalarPointer());
    float *in_ptr = static_cast<float *>(inside->GetScalarPointer());
    for (size_t k = 0; k < dims[2]; ++k)
    {
        for (size_t j = 0; j < dims[1]; ++j)
        {
            for(size_t i = 0; i < dims[0]; ++i)
            {
                if (*out_ptr == 1.0)
                {   
                    *out_ptr = 1 - (*in_ptr);
                }
                else
                    *out_ptr = (*out_ptr) - 1;

                ++out_ptr;
                ++in_ptr;
            }
        }
    }
    end = std::clock();
    std::cout<<"output dims: "<<outside->GetDimensions()[0]<<"\t"<<outside->GetDimensions()[1]<<"\t"<<outside->GetDimensions()[2]<<"\n";
    std::cout<<"compute final dist map finished. Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\n";
    dist_map = vtkSmartPointer<vtkImageData>::New();
    dist_map->DeepCopy(outside);
}

void ComputeDistMap::setVTKImg(vtkSmartPointer< vtkImageData> img_data, double mid_interval, double wid_interval)
{
    bone_img = img_data;
    if (mid_interval != 0 && mid_interval != 0)
    {
        max_interval = mid_interval + wid_interval;
        min_interval = mid_interval - wid_interval;
    }
    else
    {
        max_interval = 1000;
        min_interval = 100;
    }

    // stop criteria
    criterion = CriterionType::New();
    criterion->SetThreshold( 500 );

    // specify the size of the output image
    int *dims = bone_img->GetDimensions();
    size[0] = dims[0];
    size[1] = dims[1];
    size[2] = dims[2];
    //size[0] = 500;
    //size[1] = 500;
    //size[2] = 500;

    std::cout<<"dist map size: "<<size<<"\n";

    // region of the output image
    region.SetSize( size );

    double *ori_spacing = bone_img->GetSpacing();
    output_spacing[0] = ori_spacing[0];
    output_spacing[1] = ori_spacing[1];
    output_spacing[2] = ori_spacing[2];
    //output_spacing[0] = 1;
    //output_spacing[1] = 1;
    //output_spacing[2] = 1;
}

void ComputeDistMap::setActiveImg(FloatImageType::Pointer active_ptr, bool out_tag)
{
    // index[slice][row][col] namely index[k][j][i]
    std::cout<<"fill active image for fast marching\n";
    std::clock_t begin = clock();
    FloatImageType::IndexType index;
    for (size_t k = 0; k < size[2]; ++k)
    {
        for (size_t j = 0; j < size[1]; ++j)
        {
            for(size_t i = 0; i < size[0]; ++i)
            {
                index[2] = k;
                index[1] = j;
                index[0] = i;

                short cur_val = static_cast<short *>(bone_img->GetScalarPointer(i, j, k))[0];
                if(cur_val <= max_interval && cur_val >= min_interval)
                {
                    if (out_tag)
                        active_ptr->SetPixel(index, 1.0);
                    else
                        active_ptr->SetPixel(index, 0.0);
                }
            }
        }
    }
    
    std::clock_t end = clock();

    std::cout<<"filling finished. Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\n";
}

void ComputeDistMap::setTrialImg(FloatImageType::Pointer trial_ptr, FloatImageType::Pointer active_ptr)
{
    // index[slice][row][col] namely index[k][j][i]
    std::cout<<"filling trial image for fast marching\n";
    std::clock_t begin = std::clock();
    FloatImageType::IndexType index;
    for (size_t k = 0; k < size[2]; ++k)
    {
        for (size_t j = 0; j < size[1]; ++j)
        {
            for (size_t i = 0; i < size[0]; ++i)
            {
                index[2] = k;
                index[1] = j;
                index[0] = i;
                if (active_ptr->GetPixel(index) == 1.0)
                {
#ifdef DEBUG_PRINT
                    std::cout<<"found one active point\n";
#endif
                    if (index[0] - 1 >= 0)
                    {
                        --index[0];
                        if (active_ptr->GetPixel(index) == 0.0)
                        {
                            trial_ptr->SetPixel(index, 1.0);
#ifdef DEBUG_PRINT
                            std::cout<<"found one trial point\n";
#endif
                        }
                        ++index[0];
                    }
                    if (index[0] + 1 <size[0])
                    {
                        ++index[0];
                        if (active_ptr->GetPixel(index) == 0.0)
                        {
                            trial_ptr->SetPixel(index, 1.0);
#ifdef DEBUG_PRINT
                            std::cout<<"found one trial point\n";
#endif
                        }
                        --index[0];
                    }

                    if (index[1] - 1 >= 0)
                    {
                        --index[1];
                        if (active_ptr->GetPixel(index) == 0.0)
                        {
                            trial_ptr->SetPixel(index, 1.0);
#ifdef DEBUG_PRINT
                            std::cout<<"found one trial point\n";
#endif
                        }
                        ++index[1];
                    }
                    if (index[1] + 1 <size[1])
                    {
                        ++index[1];
                        if (active_ptr->GetPixel(index) == 0.0)
                        {
                            trial_ptr->SetPixel(index, 1.0);
#ifdef DEBUG_PRINT
                            std::cout<<"found one trial point\n";
#endif
                        }
                        --index[1];
                    }

                    if (index[2] - 1 >= 0)
                    {
                        --index[2];
                        if (active_ptr->GetPixel(index) == 0.0)
                        {    
                            trial_ptr->SetPixel(index, 1.0);
#ifdef DEBUG_PRINT
                            std::cout<<"found one trial point\n";
#endif
                        }
                        ++index[2];
                    }
                    if (index[2] + 1 <size[2])
                    {
                        ++index[2];
                        if (active_ptr->GetPixel(index) == 0.0)
                        {
                            trial_ptr->SetPixel(index, 1.0);
#ifdef DEBUG_PRINT
                            std::cout<<"found one trial point\n";
#endif
                        }
                        --index[2];
                    }
                }
            }
        }
    }
    std::clock_t end = std::clock();
    std::cout<<"filling finished. Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\n";
}

void ComputeDistMap::gaussianSmooth(vtkSmartPointer<vtkImageData> img)
{
    typedef itk::DiscreteGaussianImageFilter< FloatImageType, ShortImageType > GaussianFilterType;

    FloatImageType::Pointer AliveImage = FloatImageType::New();
    AliveImage->SetLargestPossibleRegion( region );
    AliveImage->SetBufferedRegion( region );
    AliveImage->Allocate();
    AliveImage->FillBuffer( 0.0 );
    AliveImage->SetSpacing(output_spacing);

    std::cout<<"fill itk image for gaussian filter\n";
    std::clock_t begin = clock();
    FloatImageType::IndexType index;
    for (size_t k = 0; k < size[2]; ++k)
    {
        for (size_t j = 0; j < size[1]; ++j)
        {
            for(size_t i = 0; i < size[0]; ++i)
            {
                index[2] = k;
                index[1] = j;
                index[0] = i;

                short cur_val = static_cast<short *>(img->GetScalarPointer(i, j, k))[0];
                if(cur_val <= max_interval && cur_val >= min_interval)
                {
                        AliveImage->SetPixel(index, 500);
                }
            }
        }
    }

    std::clock_t end = clock();

    std::cout<<"filling finished. Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\n";

    GaussianFilterType::Pointer gaussian_filter = GaussianFilterType::New();
    gaussian_filter->SetInput(AliveImage);
    gaussian_filter->SetVariance(1.0);
    gaussian_filter->SetUseImageSpacingOn();
    gaussian_filter->Update();

    ShortConnectorType::Pointer connector = ShortConnectorType::New();
    connector->SetInput(gaussian_filter->GetOutput());
    connector->Update();
    bone_img = vtkSmartPointer<vtkImageData>::New();
    bone_img->DeepCopy(connector->GetOutput());

    min_interval = 1.0;
    computeFinalDistMap();
}
#include "MainViewer.h"

MainViewer::MainViewer()
{
    for (size_t i = 0; i < 3; ++i)
    {
        image_slices[i] = nullptr;
        sliders[i] = nullptr;
    }
    
    nii_img = nullptr;
    temp_mesh = nullptr;
    img_volume = nullptr;

    show_volume = false;
    show_mesh = false;
    show_slices = false;

    axis = vtkSmartPointer<vtkAxesActor>::New();
    axis->SetTotalLength(20,20,20);
}

MainViewer::~MainViewer()
{
    clearImage();
    clearMesh();

    renderer->RemoveActor(axis);
    std::cout<<"delete main viewer finished\n";
}

void MainViewer::setImgData(NiiLoader *data)
{
    image_slices[0] = new ImageSlice(data->getData(), "YZ", renderer, false);
    image_slices[1] = new ImageSlice(data->getData(), "XZ", renderer, false);
    image_slices[2] = new ImageSlice(data->getData(), "XY", renderer, false);

    show_slices = true;

    img_volume = new Volume;
    img_volume->setVolumeImg(data);
    img_volume->setRenderer(renderer);
}

void MainViewer::setImgData(vtkSmartPointer<vtkImageData> data)
{
    image_slices[0] = new ImageSlice(data, "YZ", renderer, false);
    image_slices[1] = new ImageSlice(data, "XZ", renderer, false);
    image_slices[2] = new ImageSlice(data, "XY", renderer, false);

    show_slices = true;
}


void MainViewer::addImgData(vtkSmartPointer<vtkImageData> data)
{
  image_slices[0]->addImgData(data);
  image_slices[1]->addImgData(data);
  image_slices[2]->addImgData(data);

  show_slices = true;
}
void MainViewer::setImgData(std::string fName)
{
    // set img data
    if (!nii_img)
    {
        nii_img = new NiiLoader(fName.c_str());
        setImgData(nii_img);
    }
    else
    {
        std::cout<<"Warning: image data ptr isn't null\n";
    }
}

void MainViewer::setRenderer(vtkSmartPointer<vtkRenderer> mainWin_renderer)
{
    renderer = mainWin_renderer;
    renderer->AddActor(axis);
}

void MainViewer::setMeshData(std::string fName)
{
    if (!temp_mesh)
    {
        temp_mesh = new Mesh(fName, renderer);
    }
    else
    {
        std::cout<<"Warning: temp_mesh ptr isn't null\n";
    }
}

void MainViewer::updateSlice(int orient, float slice_pos)
{
    if (image_slices[orient])
        image_slices[orient]->setSliceNum(slice_pos);
    emit(updateRenderers());
}

void MainViewer::updateYZPlaneView(int value)
{
    float slice_pos = (float)value/sliders[0]->maximum();
    updateSlice(0, slice_pos);
}

void MainViewer::updateXZPlaneView(int value)
{
    float slice_pos = (float)value/sliders[1]->maximum();
    updateSlice(1, slice_pos);
}

void MainViewer::updateXYPlaneView(int value)
{
    float slice_pos = (float)value/sliders[2]->maximum();
    updateSlice(2, slice_pos);
}

void MainViewer::setSliders(QSlider *mainWin_sliders[3])
{
    for (size_t i = 0; i < 3; ++i)
    {
        sliders[i] = mainWin_sliders[i];
    }
}

void MainViewer::updateTransform(double *tr_vals)
{
    if (temp_mesh)
    {
        temp_mesh->updateTransform(tr_vals);
    }

    emit(updateRenderers());
}

void MainViewer::clearImage()
{
    if (img_volume)
    {
        delete img_volume;
        img_volume = nullptr;
        std::cout<<"delete img_volume\n";
    }
    else
    {
        std::cout<<"no img_volume\n";
    }
    
    for (size_t i = 0; i < 3; ++i)
    {
        if (image_slices[i])
        {
            delete image_slices[i];
            image_slices[i] = nullptr;
            std::cout<<"delete image slice: "<<i<<"\n";
        }
        else
        {
            std::cout<<"no image slice: "<<i<<"\n";
        }
    }

    if (nii_img)
    {
        delete nii_img;
        nii_img = nullptr;
        std::cout<<"delete nii_image\n";
    }
    else
    {
        std::cout<<"no nii_image\n";
    }
}

void MainViewer::clearMesh()
{
    if (temp_mesh)
    {
        delete temp_mesh;
        temp_mesh = nullptr;
        std::cout<<"delete temp_mesh\n";
    }
    else
    {
        std::cout<<"no temp_mesh\n";
    }
}

void MainViewer::clearData()
{
    clearImage();
    clearMesh();
}

void MainViewer::updateISOVal(int value)
{
    img_volume->setISOVal(value);
    updateVolumeView();
}

void MainViewer::updateISOWidth(int value)
{
    img_volume->setISOWidth(value);
    updateVolumeView();
}

void MainViewer::showVolume(int state)
{
    switch(state)
    {
    case 0:
        show_volume = false;
        break;
    case 2:
        show_volume = true;
        break;
    default:
            break;
    }

    updateVolumeView();
}

void MainViewer::updateVolumeView()
{
    if (img_volume)
    {
        if (show_volume)
            img_volume->dispTargetVolume();
        else
            img_volume->noDispTargetVolume();

        emit(updateRenderers());
    }
}

void MainViewer::showMesh(int state)
{
  if (temp_mesh)
  {
    temp_mesh->setVisible(state);
    emit(updateRenderers());
  }
}

void MainViewer::showImage(int state)
{
  if (nii_img)
  {
    for (int i = 0; i < 3; ++i)
    {
      image_slices[i]->setVisible(state);
    }
    emit(updateRenderers());
  }
}

void MainViewer::saveMesh(std::string f_name)
{
    vtkSmartPointer<vtkOBJWriter> writer = vtkSmartPointer<vtkOBJWriter>::New();
    writer->SetInputData(temp_mesh->getMeshData());
    writer->SetFileName(f_name.c_str());
    writer->Update();
}

void MainViewer::saveImg(std::string f_name)
{
    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(f_name.c_str());
    writer->SetInputData(image_slices[0]->getSliceMapper()->GetInput());
    writer->Write();
}

void MainViewer::runMC(std::string fName)
{
    if (!nii_img)
    {
        return;
    }

    vtkSmartPointer<vtkMarchingCubes> surface = vtkSmartPointer<vtkMarchingCubes>::New();
    surface->SetInputData(nii_img->getData());
    surface->ComputeGradientsOn();
    surface->SetNumberOfContours(1);
    surface->SetValue(0, 125);
    surface->Update();

    vtkSmartPointer<vtkPolyDataConnectivityFilter> cfilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    cfilter->SetInputConnection(surface->GetOutputPort());
    cfilter->SetExtractionModeToLargestRegion();
    //cfilter->SetExtractionModeToPointSeededRegions();
    //cfilter->SetClosestPoint(x,y,z);
    cfilter->Update();

    vtkSmartPointer<vtkOBJWriter> writer = vtkSmartPointer<vtkOBJWriter>::New();
    writer->SetInputData(cfilter->GetOutput());
    writer->SetFileName(fName.c_str());
    writer->Update();
}
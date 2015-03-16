#include "CranioViewer.h"

CranioViewer::CranioViewer()
{
    setupUi(this);

    m_3DViewerRenderer = nullptr;
    m_YZViewerRenderer = nullptr;
    m_XZViewerRenderer = nullptr;
    m_XYViewerRenderer = nullptr;

    for (size_t i = 0; i < 3; ++i)
    {
        image_slice_widgets[i] = new ImageSliceWidget;
    }
    main_viewer = new MainViewer;
    tr_widget = new TrWidget;

    icp = new ICPWrapper;
    non_rigid = new NonRigidWrapper;

    event_handler = new VTKEventHandler;

    initVTKWin();

    connect( m_OpenAction, SIGNAL( triggered() ), this, SLOT( onOpenSlot() ) );
    connect( m_LoadTemplateMesh, SIGNAL( triggered() ), this, SLOT( loadMesh() ) );
    connect( actionTestITK, SIGNAL( triggered() ), this, SLOT( testITK() ) );
    connect( actionSaveMesh, SIGNAL( triggered() ), this, SLOT( saveMesh() ) );
    connect( actionSaveImg, SIGNAL( triggered() ), this, SLOT( saveImg() ) );

    connect( m_PushButtonICP, SIGNAL( clicked() ), this, SLOT( runICP() ) );
    connect( m_PushButtonNonRigidIter, SIGNAL( clicked() ), this, SLOT( nonRigidIter() ) );
    connect( m_PushButtonNonRigid, SIGNAL( clicked() ), this, SLOT( loadOutDistMap() ) );

    setSliceUIWidget();
    setTransUIWidget();
    setVTKEventHandler();
}

CranioViewer::~CranioViewer()
{
}

void CranioViewer::initVTKWin()
{
    if (m_3DViewerRenderer)
        m_3DViewer->GetRenderWindow()->RemoveRenderer(m_3DViewerRenderer);
    m_3DViewerRenderer = vtkSmartPointer<vtkRenderer>::New();
    m_3DViewerRenderer->SetBackground(0.3, 0.3, 0.3);
    m_3DViewer->GetRenderWindow()->AddRenderer(m_3DViewerRenderer);

    if (m_YZViewerRenderer)
        m_YZViewer->GetRenderWindow()->RemoveRenderer(m_YZViewerRenderer);
    m_YZViewerRenderer = vtkSmartPointer<vtkRenderer>::New();
    m_YZViewerRenderer->SetBackground(0.3, 0.3, 0.3);
    m_YZViewer->GetRenderWindow()->AddRenderer(m_YZViewerRenderer);
    m_YZViewer->GetRenderWindow()->GetInteractor()->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleImage>::New());

    if (m_XZViewerRenderer)
        m_XZViewer->GetRenderWindow()->RemoveRenderer(m_XZViewerRenderer);
    m_XZViewerRenderer = vtkSmartPointer<vtkRenderer>::New();
    m_XZViewerRenderer->SetBackground(0.3, 0.3, 0.3);
    m_XZViewer->GetRenderWindow()->AddRenderer(m_XZViewerRenderer);
    m_XZViewer->GetRenderWindow()->GetInteractor()->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleImage>::New());

    if (m_XYViewerRenderer)
        m_XYViewer->GetRenderWindow()->RemoveRenderer(m_XYViewerRenderer);
    m_XYViewerRenderer = vtkSmartPointer<vtkRenderer>::New();
    m_XYViewerRenderer->SetBackground(0.3, 0.3, 0.3);
    m_XYViewer->GetRenderWindow()->AddRenderer(m_XYViewerRenderer);
    m_XYViewer->GetRenderWindow()->GetInteractor()->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleImage>::New());
}
void CranioViewer::onOpenSlot()
{
    QString filter;
    filter = "nii image file (*.nii)";

    QDir dir;
    QString fileName = QFileDialog::getOpenFileName( this, QString(tr("Open Image")), dir.absolutePath() , filter );
    if ( fileName.isEmpty() == true ) return;

    // 支持带中文路径的读取
    QByteArray ba = fileName.toLocal8Bit();
    const char *fileName_str = ba.data();
    std::string temp = fileName_str;
    //output_filename = temp.substr(temp.find_last_of('/')+1);

    // clear data
    main_viewer->clearImage();
    for (size_t i = 0; i < 3; ++i)
    {
        image_slice_widgets[i]->clearSlice();
        image_slice_widgets[i]->clearIntersector();
    }
    // reinit
    initVTKWin();
    main_viewer->setRenderer(m_3DViewerRenderer);
    main_viewer->setImgData(temp);

    std::string orients[3] = {"YZ", "XZ", "XY"};
    vtkSmartPointer<vtkRenderer> plane_renderers[3]=
    {m_YZViewerRenderer, m_XZViewerRenderer, m_XYViewerRenderer};
    for (size_t i = 0; i < 3; ++i)
    {
        image_slice_widgets[i]->setSlice(main_viewer->getImgData()->getData(), orients[i], plane_renderers[i]);
    }

    if (main_viewer->getMeshData())
    {
        for (size_t i = 0; i < 3; ++i)
        {
            image_slice_widgets[i]->setMesh(main_viewer->getMeshData());
        }
        main_viewer->getMeshData()->setRenderer(m_3DViewerRenderer);
    }


    //if (temp_mesh)
    //{
    //    temp_mesh->setRenderer(m_3DViewerRenderer);

    //    for (size_t i = 0; i < 3; ++i)
    //    {
    //        intersectors[i] = new Intersector;
    //        
    //        intersectors[i]->setImageSlice(image_slices[i]);
    //        intersectors[i]->setMesh(temp_mesh);
    //        intersectors[i]->setCutter();
    //    }

    //    
    //}

    updateRenderers();
}

void CranioViewer::updateRenderers()
{
    m_3DViewer->GetRenderWindow()->Render();
    m_XYViewer->GetRenderWindow()->Render();
    m_XZViewer->GetRenderWindow()->Render();
    m_YZViewer->GetRenderWindow()->Render();
}

void CranioViewer::setSliceUIWidget()
{
    QSlider *sliders[3];
    sliders[0] = m_YZSlider;
    sliders[1] = m_XZSlider;
    sliders[2] = m_XYSlider;

    image_slice_widgets[0]->setSlider(m_YZSlider);
    image_slice_widgets[1]->setSlider(m_XZSlider);
    image_slice_widgets[2]->setSlider(m_XYSlider);
    main_viewer->setSliders(sliders);

    m_CheckBoxShowCTData->setCheckable(true);
    m_CheckBoxShowCTData->setChecked(true);
    m_PushButtonClrCTData->setEnabled(true);
    
    connect( m_YZSlider, SIGNAL( sliderMoved( int ) ), image_slice_widgets[0], SLOT( updatePlaneView( int ) ) );
    connect( m_YZSlider, SIGNAL( valueChanged( int ) ), image_slice_widgets[0], SLOT( updatePlaneView( int ) ) );
    connect( m_YZSlider, SIGNAL( sliderMoved( int ) ), main_viewer, SLOT( updateYZPlaneView( int ) ) );
    connect( m_YZSlider, SIGNAL( valueChanged( int ) ), main_viewer, SLOT( updateYZPlaneView( int ) ) );
    connect( m_XZSlider, SIGNAL( sliderMoved( int ) ), image_slice_widgets[1], SLOT( updatePlaneView( int ) ) );
    connect( m_XZSlider, SIGNAL( valueChanged( int ) ), image_slice_widgets[1], SLOT( updatePlaneView( int ) ) );
    connect( m_XZSlider, SIGNAL( sliderMoved( int ) ), main_viewer, SLOT( updateXZPlaneView( int ) ) );
    connect( m_XZSlider, SIGNAL( valueChanged( int ) ), main_viewer, SLOT( updateXZPlaneView( int ) ) );
    connect( m_XYSlider, SIGNAL( sliderMoved( int ) ), image_slice_widgets[2], SLOT( updatePlaneView( int ) ) );
    connect( m_XYSlider, SIGNAL( valueChanged( int ) ), image_slice_widgets[2], SLOT( updatePlaneView( int ) ) );
    connect( m_XYSlider, SIGNAL( sliderMoved( int ) ), main_viewer, SLOT( updateXYPlaneView( int ) ) );
    connect( m_XYSlider, SIGNAL( valueChanged( int ) ), main_viewer, SLOT( updateXYPlaneView( int ) ) );

    for (size_t i = 0; i < 3; ++i)
    {       
        connect( image_slice_widgets[i], SIGNAL( updateRenderers() ), this, SLOT( updateRenderers() ) );
    }
    connect( main_viewer, SIGNAL( updateRenderers() ), this, SLOT( updateRenderers() ) );
    
    // set volume display widget
    connect( m_showVolumeCheckBox, SIGNAL( stateChanged( int ) ), main_viewer, SLOT( showVolume( int ) ) );
    connect( m_ISOValSlider, SIGNAL( valueChanged( int ) ), main_viewer, SLOT( updateISOVal( int ) ) );
    connect( m_ISOValSlider, SIGNAL( sliderMoved( int ) ), main_viewer, SLOT( updateISOVal( int ) ) );
    connect( m_ISOValSpinBox, SIGNAL( valueChanged( int ) ), m_ISOValSlider, SLOT( setValue( int ) ) );
    connect( m_ISOValSlider, SIGNAL( valueChanged( int ) ), m_ISOValSpinBox, SLOT( setValue( int ) ) );
    connect( m_ISOWidthSpinBox, SIGNAL( valueChanged( int ) ), main_viewer, SLOT( updateISOWidth( int ) ) );

}

void CranioViewer::loadMesh()
{
    QString filter;
    filter = "obj file (*.obj)";

    QDir dir;
    QString fileName = QFileDialog::getOpenFileName( this, QString(tr("Open Obj File")), dir.absolutePath() , filter );
    if ( fileName.isEmpty() == true ) return;

    // 支持带中文路径的读取
    QByteArray ba = fileName.toLocal8Bit();
    const char *fileName_str = ba.data();
    std::string temp = fileName_str;

    // clear data
    main_viewer->clearMesh();
    for (size_t i = 0; i < 3; ++i)
    {
        image_slice_widgets[i]->clearIntersector();
    }

    // reinit
    main_viewer->setRenderer(m_3DViewerRenderer);
    main_viewer->setMeshData(temp);

    if (main_viewer->getImgData())
    {
        for (size_t i = 0; i < 3; ++i)
        {
            image_slice_widgets[i]->setMesh(main_viewer->getMeshData());
        }
    }

    updateRenderers();

}

void CranioViewer::setTransUIWidget()
{
    QSlider *sliders[7];
    QDoubleSpinBox *spinBox[5];

    sliders[0] = m_SliderTransLR;
    sliders[1] = m_SliderTransPA;
    sliders[2] = m_SliderTransIS;
    sliders[3] = m_SliderRotationLR;
    sliders[4] = m_SliderRotationPA;
    sliders[5] = m_SliderRotationIS;
    sliders[6] = m_SliderScale;

    spinBox[0] = m_SpinBoxTransLR;
    spinBox[1] = m_SpinBoxTransPA;
    spinBox[2] = m_SpinBoxTransIS;
    spinBox[3] = m_SpinBoxTransMin;
    spinBox[4] = m_SpinBoxTransMax;

    tr_widget->setWidgets(sliders, spinBox);

    m_CheckBoxShowSkull->setCheckable(true);
    m_CheckBoxShowSkull->setChecked(true);

    connect( m_SliderTransLR, SIGNAL( sliderMoved( int ) ), tr_widget, SLOT( updateSliderTransLR( int ) ) );
    connect( m_SliderTransLR, SIGNAL( valueChanged( int ) ), tr_widget, SLOT( updateSliderTransLR( int ) ) );
    connect( m_SpinBoxTransLR, SIGNAL( valueChanged( double ) ), tr_widget, SLOT( updateSpinBoxTransLR( double ) ) );
    connect( m_SliderTransPA, SIGNAL( sliderMoved( int ) ), tr_widget, SLOT( updateSliderTransPA( int ) ) );
    connect( m_SliderTransPA, SIGNAL( valueChanged( int ) ), tr_widget, SLOT( updateSliderTransPA( int ) ) );
    connect( m_SpinBoxTransPA, SIGNAL( valueChanged( double ) ), tr_widget, SLOT( updateSpinBoxTransPA( double ) ) );
    connect( m_SliderTransIS, SIGNAL( sliderMoved( int ) ), tr_widget, SLOT( updateSliderTransIS( int ) ) );
    connect( m_SliderTransIS, SIGNAL( valueChanged( int ) ), tr_widget, SLOT( updateSliderTransIS( int ) ) );
    connect( m_SpinBoxTransIS, SIGNAL( valueChanged( double ) ), tr_widget, SLOT( updateSpinBoxTransIS( double ) ) );
    connect( m_SpinBoxTransMax, SIGNAL( valueChanged( double ) ), tr_widget, SLOT( updateTransMax( double ) ) );
    connect( m_SpinBoxTransMin, SIGNAL( valueChanged( double ) ), tr_widget, SLOT( updateTransMin( double ) ) );
    connect( m_SliderRotationLR, SIGNAL( sliderMoved( int ) ), tr_widget, SLOT( updateSliderRotationLR( int ) ) );
    connect( m_SliderRotationLR, SIGNAL( valueChanged( int ) ), tr_widget, SLOT( updateSliderRotationLR( int ) ) );
    connect( m_SliderRotationPA, SIGNAL( sliderMoved( int ) ), tr_widget, SLOT( updateSliderRotationPA( int ) ) );
    connect( m_SliderRotationPA, SIGNAL( valueChanged( int ) ), tr_widget, SLOT( updateSliderRotationPA( int ) ) );
    connect( m_SliderRotationIS, SIGNAL( sliderMoved( int ) ), tr_widget, SLOT( updateSliderRotationIS( int ) ) );
    connect( m_SliderRotationIS, SIGNAL( valueChanged( int ) ), tr_widget, SLOT( updateSliderRotationIS( int ) ) );
    connect( m_SliderScale, SIGNAL( sliderMoved( int ) ), tr_widget, SLOT( updateSliderScale( int ) ) );
    connect( m_SliderScale, SIGNAL( valueChanged( int ) ), tr_widget, SLOT( updateSliderScale( int ) ) );

    connect( main_viewer, SIGNAL( updateRenderers() ), this, SLOT( updateRenderers() ) );
    connect( tr_widget, SIGNAL( updateTransform( double * ) ), main_viewer, SLOT( updateTransform( double * ) ) );
}

void CranioViewer::runICP()
{
    if (main_viewer->getImgData() && main_viewer->getMeshData())
    {
        connect( icp, SIGNAL( updateRenderers() ), this, SLOT( updateRenderers() ) );
        connect( icp, SIGNAL( resetTrans() ), this, SLOT( resetTrans() ) );

        icp->setTempData(main_viewer->getMeshData());
        icp->setImage(main_viewer->getImgData());
        icp->runICP();

        disconnect( icp, SIGNAL( resetTrans() ), this, SLOT( resetTrans() ) );
        disconnect( icp, SIGNAL( updateRenderers() ), this, SLOT( updateRenderers() ) );
    }
    //if (!icp) delete icp;
    //icp = new ICPWrapper;        
}

void CranioViewer::testITK()
{
    ComputeDistMap distmap;
    distmap.setVTKImg(main_viewer->getImgData()->getData(), main_viewer->getVolumeData()->getISOVal(), main_viewer->getVolumeData()->getISOWidth());
    //distmap.computeFinalDistMap();
    distmap.gaussianSmooth(main_viewer->getImgData()->getData());

    main_viewer->clearImage();
    for (size_t i = 0; i < 3; ++i)
    {
        image_slice_widgets[i]->clearSlice();
        image_slice_widgets[i]->clearIntersector();
    }
    // reinit
    initVTKWin();
    main_viewer->setRenderer(m_3DViewerRenderer);
    main_viewer->setImgData(distmap.getDistMap());

    std::string orients[3] = {"YZ", "XZ", "XY"};
    vtkSmartPointer<vtkRenderer> plane_renderers[3]=
    {m_YZViewerRenderer, m_XZViewerRenderer, m_XYViewerRenderer};
    for (size_t i = 0; i < 3; ++i)
    {
        image_slice_widgets[i]->setSlice(distmap.getDistMap(), orients[i], plane_renderers[i]);
        image_slice_widgets[i]->setSliceColorWin(10.0);
        image_slice_widgets[i]->setSliceColorLev(0.0);
    }

    if (main_viewer->getMeshData())
    {
        for (size_t i = 0; i < 3; ++i)
        {
            image_slice_widgets[i]->setMesh(main_viewer->getMeshData());
        }
        main_viewer->getMeshData()->setRenderer(m_3DViewerRenderer);
    }

    updateRenderers();

    //non_rigid->getNonRigid()->setDistMap(distmap.getDistMap());
    //non_rigid->getNonRigid()->setMesh(main_viewer->getMeshData());
}

void CranioViewer::setVTKEventHandler()
{
    QVTKWidget *qvtk_widgets[4] = { m_3DViewer, m_YZViewer, m_XZViewer, m_XYViewer };
    event_handler->setQVtkWidgets(qvtk_widgets);
    event_handler->setStatusBar(m_StatusBar);
    event_handler->setImageSliceWidgets(image_slice_widgets);
    event_handler->setMainViewer(main_viewer);

    connect( m_PushButtonAddCrsp, SIGNAL( clicked() ), event_handler, SLOT( initUserCrsp() ) );
}

void CranioViewer::resetTrans()
{
    m_SliderTransLR->setValue(20000);
    m_SliderTransPA->setValue(20000);
    m_SliderTransIS->setValue(20000);
    m_SliderRotationLR->setValue(0);
    m_SliderRotationPA->setValue(0);
    m_SliderRotationIS->setValue(0);
    m_SliderScale->setValue(100);
    tr_widget->resetTrans();
}

void CranioViewer::nonRigidIter()
{
    size_t n_iter = spinBoxOutIter->value();

    non_rigid->getNonRigid()->setGradMaxIter(spinBoxOutIter->value());

    non_rigid->getNonRigid()->setLamdArap(doubleSpinBoxLamdArap->value());

    non_rigid->getNonRigid()->setGradStep(doubleSpinBoxGradStep->value());

    non_rigid->getNonRigid()->setLamdUserCrsp(doubleSpinBoxLamdUserCrsp->value());

    non_rigid->getNonRigid()->setUserCrsp(event_handler->getPickedIds(), event_handler->getCrspPos());


    for (size_t i = 0; i < n_iter; ++i)
    {
    non_rigid->getNonRigid()->optStep();
    updateRenderers();
    }
}

void CranioViewer::loadOutDistMap()
{
    QString filter;
    filter = "vti image file (*.vti)";

    QDir dir;
    QString fileName = QFileDialog::getOpenFileName( this, QString(tr("Open Image")), dir.absolutePath() , filter );
    if ( fileName.isEmpty() == true ) return;

    // 支持带中文路径的读取
    QByteArray ba = fileName.toLocal8Bit();
    const char *fileName_str = ba.data();

    vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
    reader->SetFileName(fileName_str);
    reader->Update();

    vtkSmartPointer<vtkImageData> dist_map = vtkSmartPointer<vtkImageData>::New();
    dist_map->DeepCopy(reader->GetOutput());

    main_viewer->clearImage();
    for (size_t i = 0; i < 3; ++i)
    {
        image_slice_widgets[i]->clearSlice();
        image_slice_widgets[i]->clearIntersector();
    }
    // reinit
    initVTKWin();
    main_viewer->setRenderer(m_3DViewerRenderer);
    main_viewer->setImgData(dist_map);

    std::string orients[3] = {"YZ", "XZ", "XY"};
    vtkSmartPointer<vtkRenderer> plane_renderers[3]=
    {m_YZViewerRenderer, m_XZViewerRenderer, m_XYViewerRenderer};
    for (size_t i = 0; i < 3; ++i)
    {
        image_slice_widgets[i]->setSlice(dist_map, orients[i], plane_renderers[i]);
        image_slice_widgets[i]->setSliceColorWin(10.0);
        image_slice_widgets[i]->setSliceColorLev(0.0);
        //double lut_range[2] = {-5.0, 5.0};
        //image_slice_widgets[i]->setSliceLUT(lut_range);
    }

    if (main_viewer->getMeshData())
    {
        for (size_t i = 0; i < 3; ++i)
        {
            image_slice_widgets[i]->setMesh(main_viewer->getMeshData());
        }
        main_viewer->getMeshData()->setRenderer(m_3DViewerRenderer);
    }

    updateRenderers();

    non_rigid->getNonRigid()->setDistMap(dist_map);
    non_rigid->getNonRigid()->setMesh(main_viewer->getMeshData());
}

void CranioViewer::saveMesh()
{
    main_viewer->saveMesh("output.obj");
}

void CranioViewer::saveImg()
{
    main_viewer->saveImg("output.vti");
}
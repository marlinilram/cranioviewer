#include "VTKEventHandler.h"

VTKEventHandler::VTKEventHandler()
{
    event_connector = vtkSmartPointer<vtkEventQtSlotConnect>::New();

    for (size_t i = 0; i < 4; ++i)
    {
        qvtk_widgets[i] = nullptr;
    }

    for (size_t i = 0; i < 3; ++i)
    {
        imge_slice_widgets[i] = nullptr;
    }

    status_bar = nullptr;
    main_viewer = nullptr;

    picked_ids.clear();
    crsp_pos.clear();
}

VTKEventHandler::~VTKEventHandler()
{

}

void VTKEventHandler::setQVtkWidgets(QVTKWidget *mainWin_widgets[4])
{
    for (size_t i = 0; i < 4; ++i)
    {
        qvtk_widgets[i] = mainWin_widgets[i];
        event_connector->Connect(qvtk_widgets[i]->GetRenderWindow()->GetInteractor(), vtkCommand::MouseMoveEvent, this, SLOT(displayCoord(vtkObject *)));
        event_connector->Connect(qvtk_widgets[i]->GetRenderWindow()->GetInteractor(), vtkCommand::LeftButtonPressEvent, this, SLOT(getWorldCoord(vtkObject *)));
        event_connector->Connect(qvtk_widgets[i]->GetRenderWindow()->GetInteractor(), vtkCommand::KeyReleaseEvent, this, SLOT(setUserCrsp(vtkObject *)));
    }
}

void VTKEventHandler::setImageSliceWidgets(ImageSliceWidget *widgets[3])
{
    for (size_t i = 0; i < 3; ++i)
    {
        imge_slice_widgets[i] = widgets[i];
    }
}

void VTKEventHandler::displayCoord(vtkObject *obj)
{
    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::SafeDownCast(obj);

    int event_pos[2];
    iren->GetEventPosition(event_pos);

    vtkSmartPointer<vtkRenderer> cur_renderer = iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
    double val = 0.0;
    for (size_t i = 0; i < 3; ++i)
    {
        //std::cout<<cur_renderer.GetPointer()<<"\t"<<imge_slice_widgets[i]->getImageSlice()->getRenderer().GetPointer()<<"\n";
        if (imge_slice_widgets[i]->getImageSlice() && cur_renderer.GetPointer() == imge_slice_widgets[i]->getImageSlice()->getRenderer().GetPointer())
        {
            val = imge_slice_widgets[i]->getImageSlice()->getPixelVal(event_pos[0], event_pos[1]);
        }
        if (imge_slice_widgets[i]->getIntersector() && cur_renderer.GetPointer() == imge_slice_widgets[i]->getIntersector()->getRenderer().GetPointer())
        {
            imge_slice_widgets[i]->getIntersector()->getWorldCoord(event_pos[0], event_pos[1]);
        }
    }


    QString str;
    str.sprintf("x=%d : y=%d : val=%f", event_pos[0], event_pos[1], val);
    
    status_bar->showMessage(str);
}

void VTKEventHandler::setStatusBar(QStatusBar *mainWin_statusBar)
{
    status_bar = mainWin_statusBar;
}

void VTKEventHandler::getWorldCoord(vtkObject *obj)
{
    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::SafeDownCast(obj);

    int event_pos[2];
    iren->GetEventPosition(event_pos);

    double coord[3] = {0,0,0};

    vtkSmartPointer<vtkRenderer> cur_renderer = iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
    for (size_t i = 0; i < 3; ++i)
    {
        if (imge_slice_widgets[i]->getIntersector() && cur_renderer.GetPointer() == imge_slice_widgets[i]->getIntersector()->getRenderer().GetPointer())
        {
            imge_slice_widgets[i]->getIntersector()->getWorldCoord(event_pos[0], event_pos[1], coord);
            //main_viewer->getMeshData()->getClosestPtID(coord);
        }
    }
}

void VTKEventHandler::setMainViewer(MainViewer *viewer)
{
    main_viewer = viewer;
}

void VTKEventHandler::initUserCrsp()
{
    main_viewer->getMeshData()->setKDTreeLocator();

    picked_ids.clear();
    crsp_pos.clear();
}

void VTKEventHandler::setUserCrsp(vtkObject *obj)
{
    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::SafeDownCast(obj);
    std::string key = iren->GetKeySym();

    int event_pos[2];
    iren->GetEventPosition(event_pos);

    vtkSmartPointer<vtkRenderer> cur_renderer = iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
    for (size_t i = 0; i < 3; ++i)
    {
        //std::cout<<cur_renderer.GetPointer()<<"\t"<<imge_slice_widgets[i]->getImageSlice()->getRenderer().GetPointer()<<"\n";
        if (imge_slice_widgets[i]->getIntersector() && cur_renderer.GetPointer() == imge_slice_widgets[i]->getIntersector()->getRenderer().GetPointer())
        {
            double coord[3] = {0,0,0};
            imge_slice_widgets[i]->getIntersector()->getWorldCoord(event_pos[0], event_pos[1], coord);
            if (key == "1")
            {
                int picked_id = main_viewer->getMeshData()->getClosestPtID(coord);
                picked_ids.push_back(picked_id);
                std::cout<<"picked id: "<<picked_id<<"\tCurrent Crsps: "<<picked_ids.size()<<"\n";
            }
            else if (key == "2")
            {
                std::cout<<"crsp pos: "<<coord[0]<<"\t"<<coord[1]<<"\t"<<coord[2]<<"\tCurrent Crsps: "<<crsp_pos.size()/3<<"\n";
                crsp_pos.push_back(coord[0]);
                crsp_pos.push_back(coord[1]);
                crsp_pos.push_back(coord[2]);
            }
        }
    }
}
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
    }


    QString str;
    str.sprintf("x=%d : y=%d : val=%f", event_pos[0], event_pos[1], val);
    
    status_bar->showMessage(str);
}

void VTKEventHandler::setStatusBar(QStatusBar *mainWin_statusBar)
{
    status_bar = mainWin_statusBar;
}
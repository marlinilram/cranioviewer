#ifndef VTKEventHandler_H
#define VTKEventHandler_H

#include <QObject>
#include <QStatusBar>

#include <vtkObject.h>
#include <vtkEvent.h>
#include <vtkEventQtSlotConnect.h>
#include <vtkSmartPointer.h>
#include <QVTKWidget.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRendererCollection.h>

#include "ImageSliceWidget.h"


class VTKEventHandler : public QObject
{
    Q_OBJECT

public:
    VTKEventHandler();
    ~VTKEventHandler();

    void setQVtkWidgets(QVTKWidget *mainWin_widgets[4]);
    void setStatusBar(QStatusBar *mainWin_statusBar);
    void setImageSliceWidgets(ImageSliceWidget *widgets[3]);

private slots:
    void displayCoord(vtkObject *obj);

private:
    vtkSmartPointer<vtkEventQtSlotConnect> event_connector;
    QVTKWidget *qvtk_widgets[4]; 
    QStatusBar *status_bar;

    ImageSliceWidget *imge_slice_widgets[3];
};

#endif
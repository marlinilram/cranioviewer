#ifndef VTKEventHandler_H
#define VTKEventHandler_H

#include <vector>

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
#include "MainViewer.h"


class VTKEventHandler : public QObject
{
    Q_OBJECT

public:
    VTKEventHandler();
    ~VTKEventHandler();

    void setQVtkWidgets(QVTKWidget *mainWin_widgets[4]);
    void setStatusBar(QStatusBar *mainWin_statusBar);
    void setImageSliceWidgets(ImageSliceWidget *widgets[3]);
    void setMainViewer(MainViewer *viewer);
    std::vector<int> &getPickedIds() { return picked_ids; };
    std::vector<double> &getCrspPos() { return crsp_pos; };

private slots:
    void displayCoord(vtkObject *obj);
    void getWorldCoord(vtkObject *obj);
    void initUserCrsp();
    void setUserCrsp(vtkObject *obj);

private:
    vtkSmartPointer<vtkEventQtSlotConnect> event_connector;
    QVTKWidget *qvtk_widgets[4]; 
    QStatusBar *status_bar;

    ImageSliceWidget *imge_slice_widgets[3];
    MainViewer *main_viewer;

    std::vector<int> picked_ids;
    std::vector<double> crsp_pos;
};

#endif
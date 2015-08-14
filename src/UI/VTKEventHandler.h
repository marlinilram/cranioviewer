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
#include <vtkAxesActor.h>

#include "ImageSliceWidget.h"
#include "MainViewer.h"
#include "NonRigidWrapper.h"

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
    void setNonRigid(NonRigidWrapper* nonRigid);
    std::vector<int> &getPickedIds() { return picked_ids; };
    std::vector<double> &getCrspPos() { return crsp_pos; };
    std::vector<double> &getCenterPos() { return centerPos; };

private slots:
    void displayCoord(vtkObject *obj);
    void getWorldCoord(vtkObject *obj);
    void initUserCrsp();
    void setUserCrsp(vtkObject *obj);
    void selectMode();
    void getSelectCenter(vtkObject *obj);
    void getSelectRad(vtkObject *obj);
    void regionMoveMode();
    void initRegionMoveMode(vtkObject *obj);
    void regionMoveCall(vtkObject *obj);
    void setLocalSmooth(int state);

private:
    vtkSmartPointer<vtkEventQtSlotConnect> event_connector;
    QVTKWidget *qvtk_widgets[4]; 
    QStatusBar *status_bar;

    ImageSliceWidget *imge_slice_widgets[3];
    MainViewer *main_viewer;

    std::vector<int> picked_ids;
    std::vector<double> crsp_pos;
    int cur_crsp_count;

    std::vector<vtkSmartPointer<vtkActor>> picked_vertices_actors;

    bool select_mode;
    bool regionMove_mode;
    std::vector<double> centerPos;
    int centerTransformPtId;
    std::vector<int> centerPtIds;
    std::vector<double> centerOldPos;
    NonRigidWrapper* nonRigid;
    vtkSmartPointer<vtkAxesActor> centerTrackball;
};

#endif
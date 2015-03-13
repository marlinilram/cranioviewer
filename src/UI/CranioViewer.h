#ifndef CranioViewer_H
#define CranioViewer_H

#include <QMainWindow>
#include "ui_CranioViewer.h"

#include <QFileDialog>
#include <QDir>

#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>

#include "ImageSliceWidget.h"
#include "MainViewer.h"
#include "TrWidget.h"

#include "ICPWrapper.h"
#include "ComputeDistMap.h"
#include "NonRigidWrapper.h"

#include "VTKEventHandler.h"

class CranioViewer : public QMainWindow, public Ui::CranioViewer
{
    Q_OBJECT

public:
    CranioViewer();
    ~CranioViewer();

    void initVTKWin();
    void setSliceUIWidget();
    void setTransUIWidget();
    void setVTKEventHandler();

private slots:
    void onOpenSlot();
    void loadMesh();
    void updateRenderers();
    void runICP();
    void testITK();
    void resetTrans();
    void nonRigidIter();
    void loadOutDistMap();
    void saveMesh();
    void saveImg();

private:
    vtkSmartPointer<vtkRenderer> m_3DViewerRenderer;
    vtkSmartPointer<vtkRenderer> m_YZViewerRenderer;
    vtkSmartPointer<vtkRenderer> m_XZViewerRenderer;
    vtkSmartPointer<vtkRenderer> m_XYViewerRenderer;

    // TODO: 尽量把真正的data放到mainviewer和imageslicewidget里
    ImageSliceWidget *image_slice_widgets[3];
    MainViewer *main_viewer;
    TrWidget *tr_widget; 

    // alg
    ICPWrapper *icp;
    NonRigidWrapper *non_rigid;

    // vtk event
    VTKEventHandler *event_handler;
};
#endif
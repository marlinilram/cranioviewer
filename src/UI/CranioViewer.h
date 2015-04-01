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
#include "MorphingViewer.h"
#include "DistMapConfigDialog.h"

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
    void computeDistMap();
    void resetTrans();
    void nonRigidIter();
    void inflateIter();
    void loadOutDistMap();
    void saveMesh();
    void saveImg();
    void showControlPanel();
    void showMorphingViewer();
    void testMC();

private:
    vtkSmartPointer<vtkRenderer> m_3DViewerRenderer;
    vtkSmartPointer<vtkRenderer> m_YZViewerRenderer;
    vtkSmartPointer<vtkRenderer> m_XZViewerRenderer;
    vtkSmartPointer<vtkRenderer> m_XYViewerRenderer;

    ImageSliceWidget *image_slice_widgets[3];
    MainViewer *main_viewer;
    TrWidget *tr_widget; 
    MorphingViewer *morphing_viewer;
    DistMapConfigDialog *dist_map_config;

    // alg
    ICPWrapper *icp;
    NonRigidWrapper *non_rigid;

    // vtk event
    VTKEventHandler *event_handler;
};
#endif
#ifndef NonRigidWrapper_H
#define NonRigidWrapper_H

#include <QObject>
#include "NonRigid.h"
#include "MainViewer.h"

class NonRigidWrapper : public QObject
{
    Q_OBJECT

public:
    NonRigidWrapper();
    ~NonRigidWrapper();

    NonRigid *getNonRigid() { return non_rigid; };
    void NonRigidIter();
    void drawCrspLines();
    void clearCrspLines();
    void setRenderer(vtkSmartPointer<vtkRenderer> renderer) { this->renderer = renderer; };
    void setMainViewer(MainViewer *main_viewer) { this->main_viewer = main_viewer; };
    void initMCKDTree(int isoval);

signals:
    void updateRenderers();

private:
    NonRigid *non_rigid;
    MainViewer *main_viewer;
    std::vector<vtkSmartPointer<vtkActor>> crsp_lines_actors;
    vtkSmartPointer<vtkRenderer> renderer;
};

#endif
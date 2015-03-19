#ifndef MorphingViewer_H
#define MorphingViewer_H

#include <QWidget>
#include "ui_MorphingViewer.h"

#include <QFileDialog>
#include <QDir>

#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>

#include "MorphingWrapper.h"

class MorphingViewer : public QWidget, public Ui::MorphingViewer
{
    Q_OBJECT

public:
    MorphingViewer();
    ~MorphingViewer();

public slots:
    void show();
    void addMesh();
    void updateRenderer();

private:
    vtkSmartPointer<vtkRenderer> morphing_renderer;
    MorphingWrapper *morphing_wrapper;
};

#endif
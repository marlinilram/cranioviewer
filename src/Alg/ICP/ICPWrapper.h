#ifndef ICPWrapper_H
#define ICPWrapper_H

#include <iostream>
#include <vector>

#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>

#include <QObject>

#include "NiiLoader.h"
#include "icpPointToPlane.h"
#include "icpPointToPoint.h"
#include "Mesh.h"
#include "Eigen\Eigen"



class ICPWrapper : public QObject
{
    Q_OBJECT

public:
    ICPWrapper();
    ~ICPWrapper();

   void setTempData(Mesh *mesh);
   void setTargetData(std::vector<double> &data);
   void setUserTrans(vtkSmartPointer<vtkMatrix4x4> mat);
   void setImage(NiiLoader *img);
   void runICPStep();
   void runICP();
   void setIter(int n_iter) { max_iter = n_iter; };

signals:
   void updateRenderers();
   void resetTrans();

private:
    int max_iter;
    NiiLoader *nii_img;
    Mesh *temp_mesh;
    IcpPointToPlane *icp_pt_plane;
    IcpPointToPoint *icp_pt_pt;
    std::vector<double> temp_mesh_vec;
    vtkSmartPointer<vtkMatrix4x4> icp_trans_mat;
    vtkSmartPointer<vtkTransform> icp_trans;
};

#endif
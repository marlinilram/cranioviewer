#ifndef NonRigid_H
#define NonRigid_H

#include "Deform.h"
#include "Mesh.h"
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkImageGradient.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <limits>

using namespace Eigen;

class NonRigid : public Deform
{
public:
    NonRigid();
    ~NonRigid();

    double computeArap(VectorXf &p_vec, VectorXf &g_vec);
    double computeDistEnergy(VectorXf &p_vec, VectorXf &g_vec);
    double computeGradEnergy(VectorXf &p_vec, VectorXf &g_vec);
    void setdvec();
    void setDistMap(vtkSmartPointer<vtkImageData> img); 
    void setMesh(Mesh *mesh_data);
    double computeVertexGrad(Vector3f &v, Vector3f &grad);
    double trilinearInterpolation(double d[3], double corner_vals[2][2][2]);
    void optStep();
    void buildAdjList();
    void setUserTrans(vtkSmartPointer<vtkMatrix4x4> mat);
    void initOpt();
    void setLamdArap(double lamd) { lamd_arap = lamd; };
    void setGradStep(double lamd) { grad_step = lamd; };
    void setGradMaxIter(int num) { grad_max_iter = num; };
    void setUserCrsp(std::vector<int> &v_ids, std::vector<double> &crsp_pts);
    void setLamdUserCrsp(double lamd) { lamd_userCrsp = lamd; };
    void computeUserCrsp(VectorXf &p_vec, VectorXf &g_vec);


protected:
    MatrixX3f d_cur;
    MatrixX3f user_crsp;
    std::vector<int> user_v_ids;
    vtkSmartPointer<vtkImageData> dist_map;
    vtkSmartPointer<vtkImageData> dist_gradient;
    std::vector<double> temp_mesh_vec;

    double lamd_arap;
    double grad_step;
    double lamd_userCrsp;
    int grad_max_iter;

    Mesh *mesh_data;
};

#endif
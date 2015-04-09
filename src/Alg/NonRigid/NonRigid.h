#ifndef NonRigid_H
#define NonRigid_H

#include "Deform.h"
#include "Mesh.h"
#include "kdtree.h"

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
    void setLamdInflate(double lamd) { lamd_inflate = lamd; };
    void setLamdDist(double lamd) { lamd_dist = lamd; };
    void computeUserCrsp(VectorXf &p_vec, VectorXf &g_vec);

    void computePPrimeNormal();
    void computeInflateDir(VectorXf &p_vec, VectorXf &g_vec);
    void inflateOptStep();

    void buildKDTree(std::vector<double> &T_Pts);
    void optStepNRICP();

protected:
    MatrixX3f d_cur;
    MatrixX3f user_crsp;
    Matrix3Xf P_Prime_N;
    vector<Eigen::Vector3i> face_list_ori;
    std::vector<int> user_v_ids;
    vtkSmartPointer<vtkImageData> dist_map;
    vtkSmartPointer<vtkImageData> dist_gradient;
    std::vector<double> temp_mesh_vec;

    double lamd_arap;
    double grad_step;
    double lamd_userCrsp;
    double lamd_inflate;
    double lamd_dist;
    int grad_max_iter;

    Mesh *mesh_data;

    kdtree::KDTree *tree_target;
    kdtree::KDTreeArray tree_data_target;
};

#endif
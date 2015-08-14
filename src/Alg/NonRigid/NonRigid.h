#ifndef NonRigid_H
#define NonRigid_H

#include "Deform.h"
#include "Mesh.h"
#include "kdtree.h"

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
    void setScalarImg(vtkSmartPointer<vtkImageData> img);
    void setScalarThreshold(double threshold) { this->scalar_threshold = threshold; };
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
    void computeBoundCrsp(VectorXf &p_vec, VectorXf &g_vec);
    int buildBoundCrsp(std::vector<int> &v_ids, MatrixX3f &bound_crsp, std::vector<double> &bound_crsp_w);
    bool searchBound(int v_id, Vector3f &v_crsp, bool out = true);
    void refineOptStep();

    int buildBoundCrsp2(MatrixX3f &bound_crsp, SparseMatrix<float>& L_cur);
    void refineOptStep2();
    bool searchBound(int v_id, float v_crsp[3], double& bound_dist, double bound_thresh);

    void buildKDTree(std::vector<double> &T_Pts);
    std::vector<double>& getTreeNormal() { return tree_normal_target; };
    void optStepNRICP();
    std::vector<double>& getCrspLines() { return crsp_lines; };
    void setVisCrspLines(bool state) { vis_crsp_lines = state; };
    bool getVisCrspLines() { return this->vis_crsp_lines; };

    // submesh registration
    void setSubMesh(std::map<vtkIdType, int>& regionMap, std::vector<int>& regionBound);
    void refineSubMesh(std::vector<double>& centerPos);
    int buildSubBoundCrsp(Matrix3Xf &bound_crsp, SparseMatrix<float>& L_cur);
    void update_SubRi();
    void setSubdvec(Matrix3Xf& d_sub);
    void setLocalSmoothMode(int state) { local_smooth_mode = state; };

protected:
    MatrixX3f d_cur;
    MatrixX3f user_crsp;
    Matrix3Xf P_Prime_N;
    vector<Eigen::Vector3i> face_list_ori;
    std::vector<int> user_v_ids;
    vtkSmartPointer<vtkImageData> dist_map;
    vtkSmartPointer<vtkImageData> dist_gradient;
    vtkSmartPointer<vtkImageData> scalar_img;
    double scalar_threshold;
    std::vector<double> temp_mesh_vec;

    double lamd_arap;
    double grad_step;
    double lamd_userCrsp;
    double lamd_inflate;
    double lamd_dist;
    int grad_max_iter;
    bool vis_crsp_lines;
    bool has_dist_map;
    bool has_mesh;
    bool has_scalar_img;

    Mesh *mesh_data;

    kdtree::KDTree *tree_target;
    kdtree::KDTreeArray tree_data_target;
    std::vector<double> tree_normal_target;
    std::vector<double> crsp_lines;

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<float>> cholSubMesh;
    SparseMatrix<float> L_Submesh;
    SparseMatrix<float> L_SubmeshAux;
    std::map<vtkIdType, int> regionMap;
    std::vector<int> NewPtOrder;
    std::vector<int> OldToNewPtOrder;
    int NewPtNum;
    int NewCenterNum;
    std::vector<int> RUpdateList;
    bool local_smooth_mode;
};

#endif
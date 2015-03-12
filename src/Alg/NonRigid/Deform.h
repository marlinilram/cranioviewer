#ifndef _Deform_H
#define _Deform_H

#include <iostream>
#include <vector>
#include <math.h>
#include "Eigen\Eigen"
#include "Eigen\Sparse"
#include "WunderSVD3x3.h"

using namespace std;

class Deform
{
public:
	Deform(){};
	~Deform(){};
	Deform(double *P, int P_Num, vector<vector<int>> &adj_list, vector<Eigen::Vector3i> &face_list);

	float *do_Deform(vector<double> &T, vector<int> &idx_T, vector<double> &N);
	float *get_P_Prime();
	float *do_Deform_Iter(double &delta);
	void set_linear_sys(vector<double> &T, vector<int> &idx_T, vector<double> &N);
	void set_linear_sys(vector<double> &T, vector<int> &idx_T, vector<double> &N, vector<double> &U, vector<long long> &idx_U);

protected:
	void update_Ri();
	double update_P_Prime();
	float compute_wij(double *p1, double *p2, double *p3, double *p4 = nullptr);
	void find_share_Vertex(int i, int j, vector<vector<int>> &adj_list, vector<Eigen::Vector3i> &face_list, vector<int> &share_Vertex);

protected:
	Eigen::Matrix3Xf P_Prime;
	Eigen::Matrix3Xf P;
	vector<vector<int>> adj_list;
    vector<Eigen::Vector3i> face_list;
	Eigen::SparseMatrix<float> Weight;
	Eigen::SparseMatrix<float> L;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<float>> chol;
	Eigen::MatrixX3f d;
	vector<Eigen::Matrix3f> R;
	int max_iter;
	double min_delta;
	float lamd_deform;
	float lamd_norm_dist;
	float lamd_user_crsp;
};

#endif
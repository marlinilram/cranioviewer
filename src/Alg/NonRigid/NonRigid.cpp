#include "NonRigid.h"

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkWeakPointer.h>
#include <vtkImageCast.h>
#include <vtkImageGradient.h>
#include <vtkImageGradientMagnitude.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <limits>

NonRigid::NonRigid()
{
  tree_target = nullptr;
  dist_map = nullptr;
  scalar_img = nullptr;
  has_dist_map = false;
  has_mesh = false;
  has_scalar_img = false;
  scalar_threshold = 150.0;
  local_smooth_mode = false;
}

NonRigid::~NonRigid()
{
  if (tree_target)
  {
    delete tree_target;
    std::cout<<"Deleted KD Tree.\n";
  }
}

double NonRigid::computeArap(VectorXf &p_vec, VectorXf &g_vec)
{
  VectorXf d_vec(VectorXf::Map(d_cur.data(), d_cur.cols()*d_cur.rows()));

  g_vec = L*p_vec - d_vec;

  return (0.5*p_vec.transpose()*L*p_vec - d_vec.transpose()*p_vec)(0, 0);
}

double NonRigid::computeDistEnergy(VectorXf &p_vec, VectorXf &g_vec)
{
  double dist_e = 0;
  Vector3f cur_v(0,0,0);
  Vector3f cur_g(0,0,0);
  size_t P_Num = p_vec.size()/3;
  for (size_t i = 0; i < P_Num; ++i)
  {
    cur_v << p_vec(i + 0*P_Num), p_vec(i + 1*P_Num), p_vec(i + 2*P_Num);
    double cur_dist_val = computeVertexGrad(cur_v, cur_g);
    //g_vec(i + 0*P_Num) = cur_g(0);
    //g_vec(i + 1*P_Num) = cur_g(1);
    //g_vec(i + 2*P_Num) = cur_g(2);
    if (cur_dist_val > 0.0)
    {
      g_vec(i + 0*P_Num) = cur_g(0);
      g_vec(i + 1*P_Num) = cur_g(1);
      g_vec(i + 2*P_Num) = cur_g(2);
      dist_e += cur_dist_val * cur_dist_val;
    }
    else
    {
      g_vec(i + 0*P_Num) = -P_Prime_N(0, i);
      g_vec(i + 1*P_Num) = -P_Prime_N(1, i);
      g_vec(i + 2*P_Num) = -P_Prime_N(2, i);
    }
  }
  return dist_e;
}

void NonRigid::setdvec()
{
  d_cur = MatrixX3f::Zero(adj_list.size(), 3);
  for (size_t i = 0; i < adj_list.size(); ++i)
  {
    for (size_t j = 0; j < adj_list[i].size(); ++j)
    {
      d_cur.row(i) += ((Weight.coeffRef(i, adj_list[i][j])/2)*(R[i]+R[adj_list[i][j]])*(P.col(i) - P.col(adj_list[i][j]))).transpose();
    }
  }
}

void NonRigid::setScalarImg(vtkSmartPointer<vtkImageData> img)
{
  if (scalar_img != nullptr) scalar_img->Print(std::cout);
  scalar_img = img;
  has_scalar_img = true;
}

void NonRigid::setDistMap(vtkSmartPointer<vtkImageData> img)
{
  if (dist_map != nullptr)  dist_map->Print(std::cout);
  //dist_map = vtkSmartPointer<vtkImageData>::New();
  //dist_map->DeepCopy(img);
  dist_map = img;

  std::cout<<"dist map scalar size: "<<dist_map->GetNumberOfScalarComponents()<<"\n";

  vtkSmartPointer<vtkImageGradient> g_filter = vtkSmartPointer<vtkImageGradient>::New();

  g_filter->SetInputData(dist_map);
  g_filter->SetDimensionality(3);
  g_filter->Update();

  vtkSmartPointer<vtkImageCast> cast_filter = vtkSmartPointer<vtkImageCast>::New();
  cast_filter->SetInputConnection(g_filter->GetOutputPort());
  cast_filter->SetOutputScalarTypeToFloat();
  cast_filter->Update();

  dist_gradient = (cast_filter->GetOutput());

  has_dist_map = true;

  //std::cout<<"dist gradient scalar size: "<<dist_gradient->GetNumberOfScalarComponents()<<"\n";
}

double NonRigid::computeVertexGrad(Vector3f &v, Vector3f &grad)
{
  // 1. get the cell contains vertex and its 8 corner points
  // 2. compute gradients for this 8 points
  // 3. compute trilinear interpolation

  double x[3] = {v[0], v[1], v[2]};
  int idx[3];
  double pcoords[3];
  if (!dist_map->ComputeStructuredCoordinates(x, idx, pcoords))
  {
    grad = Vector3f(0,0,0);
    return 0;
  }


  // idx indicates the origin point x0 y0 z0 of found cell
  // pcoords is parametric xd yd zd for trilinear interpolation
  // get the left 7 points
  double grad_x[2][2][2];
  double grad_y[2][2][2];
  double grad_z[2][2][2];
  double f_dist[2][2][2];
  double cur_grad[3];
  for (size_t i = 0; i < 2; ++i)
  {
    for (size_t j = 0; j < 2; ++j)
    {
      for (size_t k = 0; k < 2; ++k)
      {
        // or get gradient here
        //corner_vals[i][j][k] = dist_map->GetScalarComponentAsDouble(idx[0]+i, idx[1]+j, idx[2]+k, 0);
        grad_x[i][j][k] = dist_gradient->GetScalarComponentAsDouble(idx[0]+i, idx[1]+j, idx[2]+k, 0);
        grad_y[i][j][k] = dist_gradient->GetScalarComponentAsDouble(idx[0]+i, idx[1]+j, idx[2]+k, 1);
        grad_z[i][j][k] = dist_gradient->GetScalarComponentAsDouble(idx[0]+i, idx[1]+j, idx[2]+k, 2);
        f_dist[i][j][k] = dist_map->GetScalarComponentAsDouble(idx[0]+i, idx[1]+j, idx[2]+k, 0);
      }
    }
  }

  // interpolation
  grad = Vector3f(0,0,0);
  grad(0) = this->trilinearInterpolation(pcoords, grad_x);
  grad(1) = this->trilinearInterpolation(pcoords, grad_y);
  grad(2) = this->trilinearInterpolation(pcoords, grad_z);
  return this->trilinearInterpolation(pcoords, f_dist);
}

double NonRigid::trilinearInterpolation(double d[3], double corner_vals[2][2][2])
{
  double c00 = corner_vals[0][0][0]*(1.0-d[0]) + corner_vals[1][0][0]*d[0];
  double c10 = corner_vals[0][1][0]*(1.0-d[0]) + corner_vals[1][1][0]*d[0];
  double c01 = corner_vals[0][0][1]*(1.0-d[0]) + corner_vals[1][0][1]*d[0];
  double c11 = corner_vals[0][1][1]*(1.0-d[0]) + corner_vals[1][1][1]*d[0];

  double c0 = c00*(1.0-d[1]) + c10*d[1];
  double c1 = c01*(1.0-d[1]) + c11*d[1];

  return c0*(1.0-d[2]) + c1*d[2];
}

double NonRigid::computeGradEnergy(VectorXf &p_vec, VectorXf &g_vec)
{
  VectorXf arap_g = g_vec;
  VectorXf dist_g = g_vec;
  VectorXf userCrsp_g = g_vec;
  double f_e = 0;
  computeArap(p_vec, arap_g);

  f_e += computeDistEnergy(p_vec, dist_g);

  computeUserCrsp(p_vec, userCrsp_g);

  g_vec = lamd_dist*dist_g + lamd_arap*arap_g + lamd_userCrsp*userCrsp_g;
  g_vec.normalize();
  return f_e;
}

void NonRigid::optStep()
{
  if (!has_dist_map || !has_mesh) return;
  size_t P_Num = adj_list.size();
  VectorXf p_vec(3*P_Num);
  p_vec.segment(0, P_Num) = P_Prime.row(0);
  p_vec.segment(0+P_Num, P_Num) = P_Prime.row(1);
  p_vec.segment(0+2*P_Num, P_Num) = P_Prime.row(2);
  VectorXf g_vec = VectorXf::Zero(p_vec.size());
  double f_val = numeric_limits<double>::max();
  double f_cur_val = numeric_limits<double>::max();

  vtkSmartPointer<vtkPolyData> mesh = mesh_data->getMeshData();
  vtkSmartPointer<vtkPoints> new_points = mesh->GetPoints();

  std::clock_t begin = std::clock();
  update_Ri();
  setdvec();
  for (size_t i = 0; i < 1; ++i)
  {
    //begin = std::clock();
    computePPrimeNormal();
    f_cur_val = computeGradEnergy(p_vec, g_vec);
    //end = std::clock();
    //std::cout<<"Computed Gradient finished. Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\nfunction value: "<<f_val<<"delta: "<<f_val-f_cur_val<<"\n";
    //if (f_val - f_cur_val < 1e-3) break;
    std::cout << "Dist Energy: " << f_cur_val << "\n";
    f_val = f_cur_val;
    p_vec = p_vec - grad_step*g_vec;
  }

  P_Prime.row(0) = p_vec.segment(0, P_Num);
  P_Prime.row(1) = p_vec.segment(0+P_Num, P_Num);
  P_Prime.row(2) = p_vec.segment(0+2*P_Num, P_Num);

  std::clock_t end = std::clock();
  std::cout<<"One opt step finished. Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\n";

  if (new_points->GetNumberOfPoints() == p_vec.size()/3)
  {
    for (size_t i = 0; i != new_points->GetNumberOfPoints(); ++i) {
      new_points->SetPoint(i, p_vec[i], p_vec[i+P_Num], p_vec[i+2*P_Num]);
    }
    new_points->Modified();
    mesh_data->resetMesh(mesh);
  }
  else std::cout<<"Error: number of points doesn't match\n";
}

void NonRigid::setMesh(Mesh *mesh_data)
{
  this->mesh_data = mesh_data;
  lamd_arap = 10;
  grad_max_iter = 25;
  grad_step = 0.01;

  std::clock_t begin;
  std::clock_t end;
  begin = std::clock();

  buildAdjList();

  end = std::clock();
  std::cout<<"Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\n";

  // init mesh data vector to current position
  vtkSmartPointer<vtkMatrix4x4> trans_mat = vtkSmartPointer<vtkMatrix4x4>::New();
  mesh_data->getTransform()->GetMatrix(trans_mat);
  int meshNum = mesh_data->getMeshData()->GetNumberOfPoints();
  temp_mesh_vec.clear();
  for (size_t i = 0; i < meshNum; ++i) {
    double* temp = mesh_data->getMeshData()->GetPoint(i);
    temp_mesh_vec.push_back(temp[0]);
    temp_mesh_vec.push_back(temp[1]);
    temp_mesh_vec.push_back(temp[2]);
  }

  setUserTrans(trans_mat);

  begin = std::clock();

  initOpt();

  end = std::clock();
  std::cout<<"Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\n";

  has_mesh = true;
}

void NonRigid::buildAdjList()
{
  // build adjacent list
  std::cout<<"building adjacent list\n";
  adj_list.clear();
  vtkSmartPointer<vtkPolyData> mesh = mesh_data->getMeshData();
  for (vtkIdType i = 0; i != mesh->GetNumberOfPoints(); ++i) {
    vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(i, cellIdList);
    vector<int> v_adj_list;
    for (vtkIdType j = 0; j != cellIdList->GetNumberOfIds(); ++j) {
      vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
      mesh->GetCellPoints(cellIdList->GetId(j), pointIdList);
      if (pointIdList->GetId(0) == i) {
        v_adj_list.push_back(pointIdList->GetId(1));
        v_adj_list.push_back(pointIdList->GetId(2));
      }
      else if (pointIdList->GetId(1) == i){
        v_adj_list.push_back(pointIdList->GetId(0));
        v_adj_list.push_back(pointIdList->GetId(2));
      }
      else {
        v_adj_list.push_back(pointIdList->GetId(0));
        v_adj_list.push_back(pointIdList->GetId(1));
      }
    }
    sort(v_adj_list.begin(), v_adj_list.end());
    vector<int>::iterator iter = unique(v_adj_list.begin(), v_adj_list.end());
    v_adj_list.erase(iter, v_adj_list.end());
    v_adj_list.shrink_to_fit();
    adj_list.push_back(v_adj_list);
  }


  // build face list
  face_list.clear();
  face_list_ori.clear();
  for (vtkIdType i = 0; i != mesh->GetNumberOfCells(); ++i) {
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
    mesh->GetCellPoints(i, pointIdList);
    vector<int> f;
    f.push_back(pointIdList->GetId(0));
    f.push_back(pointIdList->GetId(1));
    f.push_back(pointIdList->GetId(2));
    face_list_ori.push_back((Eigen::Map<Eigen::Vector3i>(&f[0])));
    sort(f.begin(), f.end());
    face_list.push_back(Eigen::Map<Eigen::Vector3i>(&f[0]));
  }
  std::cout<<"building finished\n";
}

void NonRigid::setUserTrans(vtkSmartPointer<vtkMatrix4x4> mat)
{
  if (!temp_mesh_vec.empty())
  {
    // convert vertices to eigen 4*n matrix
    Eigen::MatrixXd v_mat(4, temp_mesh_vec.size()/3);
    v_mat << Eigen::Map<Eigen::MatrixXd>(&temp_mesh_vec[0], 3, temp_mesh_vec.size()/3), Eigen::RowVectorXd::Ones(temp_mesh_vec.size()/3);

    // set transform
    Eigen::Matrix4d trans_mat;
    for (size_t i = 0; i < 4; ++i)
    {
      for (size_t j = 0; j < 4; ++j)
      {
        trans_mat(i, j) = mat->GetElement(i, j);
      }
    }

    Eigen::MatrixXd new_v_mat = trans_mat * v_mat;
    temp_mesh_vec.clear();
    for (size_t i = 0; i < new_v_mat.cols(); ++i)
    {
      temp_mesh_vec.push_back(new_v_mat(0, i)/new_v_mat(3, i));
      temp_mesh_vec.push_back(new_v_mat(1, i)/new_v_mat(3, i));
      temp_mesh_vec.push_back(new_v_mat(2, i)/new_v_mat(3, i));
    }
  }
}

void NonRigid::initOpt()
{
  std::cout<<"init optimization\n";

  size_t P_Num = temp_mesh_vec.size()/3;
  P.resize(3, P_Num);
  for (int i = 0; i != P_Num; ++i) {P.col(i) << temp_mesh_vec[3*i], temp_mesh_vec[3*i+1], temp_mesh_vec[3*i+2];}
  P_Prime = P;
  P_Prime_N.resize(3, P_Num);
  P_Prime_N.setOnes();
  R = vector<Matrix3f>(P_Num, Matrix3f::Identity());

  vector<Triplet<float>> weight_list;
  weight_list.reserve(3*7*P_Num); // each vertex may have about 7 vertices connected
  vector<Triplet<float>> weight_sum;
  weight_sum.reserve(3*P_Num);
  for (decltype(adj_list.size()) i = 0; i != adj_list.size(); ++i) {
    float wi = 0;
    for (decltype(adj_list[i].size()) j = 0; j != adj_list[i].size(); ++j) {
      int id_j = adj_list[i][j];
      //if (i < id_j) {
      vector<int> share_Vertex;

      find_share_Vertex(i, id_j, adj_list, face_list, share_Vertex);			
      float wij = 0;
      if (share_Vertex.size()==2) wij = compute_wij(&temp_mesh_vec[3*i], &temp_mesh_vec[3*id_j], &temp_mesh_vec[3*share_Vertex[0]], &temp_mesh_vec[3*share_Vertex[1]]);
      else wij = compute_wij(&temp_mesh_vec[3*i], &temp_mesh_vec[3*id_j], &temp_mesh_vec[3*share_Vertex[0]]);

      wi += wij;
      weight_list.push_back(Triplet<float>(i, id_j, wij));
      weight_list.push_back(Triplet<float>(i+P_Num, id_j+P_Num, wij));
      weight_list.push_back(Triplet<float>(i+2*P_Num, id_j+2*P_Num, wij));
      //weight_list.push_back(Triplet<float>(id_j, i, wij));
      //}
    }
    if (wi < 0.1) cout << "Edge Weight Sum Warning: " << wi << "\tvertex: "<<i<<"\n";
    weight_sum.push_back(Triplet<float>(i, i, wi));
    weight_sum.push_back(Triplet<float>(i+P_Num, i+P_Num, wi));
    weight_sum.push_back(Triplet<float>(i+2*P_Num, i+2*P_Num, wi));
    for (int id_r = 0; id_r != 3; ++id_r)
      for (int id_c = 0; id_c != 3; ++id_c)
        weight_sum.push_back(Triplet<float>(i+id_r*P_Num, i+id_c*P_Num, 0));
    //cout << R[i] << endl;
  }
  SparseMatrix<float> Weight_sum(3*P_Num, 3*P_Num);
  Weight_sum.setFromTriplets(weight_sum.begin(), weight_sum.end());
  Weight.resize(3*P_Num, 3*P_Num);
  Weight.setFromTriplets(weight_list.begin(), weight_list.end());
  L =  Weight_sum - Weight;
  chol.analyzePattern(L);
  std::cout<<"init opt finished\n";

  //std::ofstream f_debug("arap.txt");
  //if (f_debug)
  //{
  //    f_debug<<P<<"\n";
  //    for (size_t i = 0; i < R.size(); ++i)
  //    {
  //        f_debug<<R[i]<<"\n";
  //    }
  //    //f_debug<<L<<"\n";
  //    f_debug.close();
  //}
}

void NonRigid::setUserCrsp(std::vector<int> &v_ids, std::vector<double> &crsp_pts)
{
  user_crsp = MatrixX3f::Zero(adj_list.size(), 3);
  user_v_ids = v_ids;

  for (size_t i = 0; i < v_ids.size(); ++i)
  {
    user_crsp.row(v_ids[i]) = RowVector3f(crsp_pts[3*i+0], crsp_pts[3*i+1], crsp_pts[3*i+2]);
  }
}

void NonRigid::computeUserCrsp(VectorXf &p_vec, VectorXf &g_vec)
{
  size_t P_Num = adj_list.size();
  g_vec = VectorXf::Zero(3*P_Num);

  for (size_t i = 0; i < user_v_ids.size(); ++i)
  {
    g_vec(user_v_ids[i] + 0*P_Num) = p_vec(user_v_ids[i] + 0*P_Num);
    g_vec(user_v_ids[i] + 1*P_Num) = p_vec(user_v_ids[i] + 1*P_Num);
    g_vec(user_v_ids[i] + 2*P_Num) = p_vec(user_v_ids[i] + 2*P_Num);
  }

  g_vec = g_vec - Map<VectorXf>(user_crsp.data(), 3*P_Num, 1);
}

void NonRigid::computePPrimeNormal()
{
  Eigen::Vector3f edge_0;
  Eigen::Vector3f edge_1;
  Eigen::Vector3f cur_f_normal;
  size_t P_Num = adj_list.size();

  P_Prime_N.setZero();

  for (size_t i = 0; i < face_list_ori.size(); ++i)
  {
    edge_0 = P_Prime.col(face_list_ori[i](1)) - P_Prime.col(face_list_ori[i](0));
    edge_1 = P_Prime.col(face_list_ori[i](2)) - P_Prime.col(face_list_ori[i](1));

    cur_f_normal = edge_0.cross(edge_1);
    cur_f_normal.normalize();

    P_Prime_N.col(face_list_ori[i](0)) += cur_f_normal;
    P_Prime_N.col(face_list_ori[i](1)) += cur_f_normal;
    P_Prime_N.col(face_list_ori[i](2)) += cur_f_normal;
  }

  for (size_t i = 0; i < P_Num; ++i)
  {
    P_Prime_N.col(i).normalize();
  }
}

void NonRigid::inflateOptStep()
{
  // if D(vi) < 0, move towards normal direction
  // else, move towards gradient direction

  size_t P_Num = adj_list.size();
  VectorXf p_vec(3*P_Num);
  p_vec.segment(0, P_Num) = P_Prime.row(0);
  p_vec.segment(0+P_Num, P_Num) = P_Prime.row(1);
  p_vec.segment(0+2*P_Num, P_Num) = P_Prime.row(2);
  VectorXf g_vec = VectorXf::Zero(p_vec.size());
  double f_val = numeric_limits<double>::max();
  double f_cur_val = numeric_limits<double>::max();

  vtkSmartPointer<vtkPolyData> mesh = mesh_data->getMeshData();
  vtkSmartPointer<vtkPoints> new_points = mesh->GetPoints();

  std::clock_t begin = std::clock();
  update_Ri();
  setdvec();
  for (size_t i = 0; i < grad_max_iter; ++i)
  {
    //begin = std::clock();
    computePPrimeNormal();
    computeInflateDir(p_vec, g_vec);
    //end = std::clock();
    //std::cout<<"Computed Gradient finished. Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\nfunction value: "<<f_val<<"delta: "<<f_val-f_cur_val<<"\n";

    p_vec = p_vec - grad_step*g_vec;
    P_Prime.row(0) = p_vec.segment(0, P_Num);
    P_Prime.row(1) = p_vec.segment(0+P_Num, P_Num);
    P_Prime.row(2) = p_vec.segment(0+2*P_Num, P_Num);
  }

  std::clock_t end = std::clock();
  std::cout<<"One opt step finished. Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\n";

  if (new_points->GetNumberOfPoints() == p_vec.size()/3)
  {
    for (size_t i = 0; i != new_points->GetNumberOfPoints(); ++i) {
      new_points->SetPoint(i, p_vec[i], p_vec[i+P_Num], p_vec[i+2*P_Num]);
    }
    new_points->Modified();
    mesh_data->resetMesh(mesh);
  }
  else std::cout<<"Error: number of points doesn't match\n";
}

void NonRigid::computeInflateDir(VectorXf &p_vec, VectorXf &g_vec)
{
  VectorXf arap_g = g_vec;
  VectorXf inflate_g = g_vec;
  VectorXf userCrsp_g = g_vec;
  VectorXf boundCrsp_g = g_vec;
  double f_e = 0;
  computeArap(p_vec, arap_g);
  computeBoundCrsp(p_vec, boundCrsp_g);
  computeUserCrsp(p_vec, userCrsp_g);

  Vector3f cur_v(0,0,0);
  Vector3f cur_g(0,0,0);
  size_t P_Num = p_vec.size()/3;
  for (size_t i = 0; i < P_Num; ++i)
  {
    cur_v << p_vec(i + 0*P_Num), p_vec(i + 1*P_Num), p_vec(i + 2*P_Num);
    if (computeVertexGrad(cur_v, cur_g) > 0.0)
    {
      inflate_g(i + 0*P_Num) = cur_g(0);
      inflate_g(i + 1*P_Num) = cur_g(1);
      inflate_g(i + 2*P_Num) = cur_g(2);
    }
    else
    {
      inflate_g(i + 0*P_Num) = -P_Prime_N(0, i);
      inflate_g(i + 1*P_Num) = -P_Prime_N(1, i);
      inflate_g(i + 2*P_Num) = -P_Prime_N(2, i);
    }
  }

  g_vec = lamd_inflate*inflate_g + lamd_arap*arap_g + lamd_userCrsp*userCrsp_g + lamd_userCrsp*boundCrsp_g;
  g_vec.normalize();
}

void NonRigid::buildKDTree(std::vector<double> &T_Pts)
{
  std::cout<<"\nBuilding KD Tree\n";
  tree_data_target.resize(boost::extents[T_Pts.size()/3][3]);

  for (size_t m = 0; m < T_Pts.size()/3; ++m)
  {
    for(size_t n = 0; n < 3; ++n)
    {
      tree_data_target[m][n] = (float)T_Pts[m*3+n];
    }
  }

  if (tree_target)
  {
    delete tree_target;
    std::cout<<"Deleted KD Tree.\n";
  }

  tree_target = new kdtree::KDTree(tree_data_target);
  std::cout<<"Building KD Tree finished\n";
}

void NonRigid::optStepNRICP()
{
  std::cout<<"current dist map ref count: "<<dist_map->GetReferenceCount()<<"\n";
}

bool NonRigid::searchBound(int v_id, Vector3f &v_crsp, bool out)
{
  Vector3d x;
  int idx[3];
  double pcoords[3];
  double f_dist[2][2][2];
  double step = 0.0;
  double step_length = 0.05;
  int search_iter = 0;

  //x = (P_Prime.col(v_id) + step * P_Prime_N.col(v_id)).cast<double>();

  //if (!dist_map->ComputeStructuredCoordinates(x.data(), idx, pcoords))
  //{
  //  return false;
  //}

  //std::vector<float> query(3);
  //query[0] = x[0];
  //query[1] = x[1];
  //query[2] = x[2];
  //kdtree::KDTreeResultVector result;
  //tree_target->n_nearest(query, 1, result);
  //v_crsp << 
  //  tree_target->the_data[result[0].idx][0],
  //  tree_target->the_data[result[0].idx][1],
  //  tree_target->the_data[result[0].idx][2];
  //return true;

  do 
  {
    step = step + step_length;//(out ? 0.05 : -0.05);
    x = (P_Prime.col(v_id) + step * P_Prime_N.col(v_id)).cast<double>();

    if (!dist_map->ComputeStructuredCoordinates(x.data(), idx, pcoords))
    {
      return false;
    }

    for (size_t i = 0; i < 2; ++i)
    {
      for (size_t j = 0; j < 2; ++j)
      {
        for (size_t k = 0; k < 2; ++k)
        {
          f_dist[i][j][k] = dist_map->GetScalarComponentAsDouble(idx[0]+i, idx[1]+j, idx[2]+k, 0);
        }
      }
    }

    ++search_iter;
  } while ((this->trilinearInterpolation(pcoords, f_dist)) <= 0);

  //if (search_iter <= 50)
  {
    v_crsp = (P_Prime.col(v_id) + (step - (out ? step_length/2 : -step_length/2)) * P_Prime_N.col(v_id));
    return true;
  }
  //else return false;
}

int NonRigid::buildBoundCrsp(std::vector<int> &v_ids, MatrixX3f &bound_crsp, std::vector<double> &bound_crsp_w)
{
  size_t P_Num = P_Prime.cols();
  bound_crsp = MatrixX3f::Zero(P_Num, 3);
  crsp_lines.clear();

  std::vector<bool> is_in_user_crsp(P_Num, false);
  for (size_t i = 0; i < user_v_ids.size(); ++i)
  {
    is_in_user_crsp[user_v_ids[i]] = true;
  }

  Vector3f v_crsp;
  Vector3d x;
  int idx[3];
  double pcoords[3];
  double f_dist[2][2][2];
  int num_outBoundCrsp = 0;
  int num_inBoundCrsp = 0;
  for (size_t v_id = 0; v_id < P_Num; ++v_id)
  {
    if (is_in_user_crsp[v_id] == true) continue;

    x = (P_Prime.col(v_id)).cast<double>();

    if (!dist_map->ComputeStructuredCoordinates(x.data(), idx, pcoords))
    {
      continue;
    }

    for (size_t i = 0; i < 2; ++i)
    {
      for (size_t j = 0; j < 2; ++j)
      {
        for (size_t k = 0; k < 2; ++k)
        {
          f_dist[i][j][k] = dist_map->GetScalarComponentAsDouble(idx[0]+i, idx[1]+j, idx[2]+k, 0);
        }
      }
    }

    double cur_dist_val = this->trilinearInterpolation(pcoords, f_dist);
    double cur_crsp_dist = 0.0;
    if (cur_dist_val >= 0) 
    {
      bound_crsp.row(v_id) = 10 * lamd_arap * RowVector3f(x[0], x[1], x[2]);

      v_ids.push_back(v_id);
      bound_crsp_w.push_back(10 * lamd_arap);
      //continue;
    }
    else if (cur_dist_val < 0)
    {
      if(searchBound(v_id, v_crsp, true)) 
      {
        cur_crsp_dist = (x-v_crsp.cast<double>()).norm();
        v_ids.push_back(v_id);
        bound_crsp.row(v_id) = RowVector3f(v_crsp[0], v_crsp[1], v_crsp[2]) / cur_crsp_dist;
        bound_crsp_w.push_back(1 / cur_crsp_dist);
        ++num_outBoundCrsp;
        crsp_lines.push_back(x[0]);
        crsp_lines.push_back(x[1]);
        crsp_lines.push_back(x[2]);
        crsp_lines.push_back(v_crsp[0]);
        crsp_lines.push_back(v_crsp[1]);
        crsp_lines.push_back(v_crsp[2]);
      }
    }
    //else
    //{
    //  //if(searchBound(v_id, v_crsp, false) && (x-v_crsp.cast<double>()).norm() <= 3.0) 
    //  //{
    //  //  v_ids.push_back(v_id);
    //  //  bound_crsp.row(v_id) = RowVector3f(v_crsp[0], v_crsp[1], v_crsp[2]);
    //  //  ++num_inBoundCrsp;
    //  //}
    //  bound_crsp.row(v_id) = RowVector3f(x[0], x[1], x[2]);

    //  v_ids.push_back(v_id);
    //  bound_crsp_w.push_back(10*lamd_arap);
    //}
  }

  std::cout<<"search in: "<<num_inBoundCrsp<<"\tsearch out: "<<num_outBoundCrsp<<"\n";
  return num_inBoundCrsp + num_outBoundCrsp;
}

void NonRigid::computeBoundCrsp(VectorXf &p_vec, VectorXf &g_vec)
{
  size_t P_Num = P_Prime.cols();
  std::vector<int> v_ids;
  MatrixX3f bound_crsp = MatrixX3f::Zero(P_Num, 3);
  g_vec = VectorXf::Zero(3*P_Num);

  //buildBoundCrsp(v_ids, bound_crsp);

  for (size_t i = 0; i < v_ids.size(); ++i)
  {
    g_vec(v_ids[i] + 0*P_Num) = p_vec(v_ids[i] + 0*P_Num);
    g_vec(v_ids[i] + 1*P_Num) = p_vec(v_ids[i] + 1*P_Num);
    g_vec(v_ids[i] + 2*P_Num) = p_vec(v_ids[i] + 2*P_Num);
  }

  g_vec = g_vec - Map<VectorXf>(bound_crsp.data(), 3*P_Num, 1);

  std::cout<<"Num of bound crsp: "<<v_ids.size()<<"\n";
}

void NonRigid::refineOptStep()
{
  if (!has_dist_map || !has_mesh) return;
  size_t P_Num = adj_list.size();
  std::vector<int> v_ids;
  MatrixX3f bound_crsp = MatrixX3f::Zero(P_Num, 3);
  std::vector<double> bound_crsp_w;
  VectorXf P_Prime_vec;

  vtkSmartPointer<vtkPolyData> mesh = mesh_data->getMeshData();
  vtkSmartPointer<vtkPoints> new_points = mesh->GetPoints();

  std::clock_t begin = std::clock();

  computePPrimeNormal();
  int num_boundCrsp = buildBoundCrsp(v_ids, bound_crsp, bound_crsp_w);
  std::cout<<"Num of bound crsp: "<<num_boundCrsp<<"\n";
  if (vis_crsp_lines) return;

  SparseMatrix<float> L_cur = lamd_arap * L;

  for (size_t i = 0; i < v_ids.size(); ++i)
  {
    L_cur.coeffRef(v_ids[i], v_ids[i]) += lamd_userCrsp * bound_crsp_w[i];
    L_cur.coeffRef(v_ids[i] + P_Num, v_ids[i] + P_Num) += lamd_userCrsp * bound_crsp_w[i];
    L_cur.coeffRef(v_ids[i] + 2 * P_Num, v_ids[i] + 2 * P_Num) += lamd_userCrsp * bound_crsp_w[i];
  }
  for (size_t i = 0; i < user_v_ids.size(); ++i)
  {
    L_cur.coeffRef(user_v_ids[i], user_v_ids[i]) += lamd_userCrsp;
    L_cur.coeffRef(user_v_ids[i] + P_Num, user_v_ids[i] + P_Num) += lamd_userCrsp;
    L_cur.coeffRef(user_v_ids[i] + 2 * P_Num, user_v_ids[i] + 2 * P_Num) += lamd_userCrsp;
  }
  chol.factorize(L_cur);

  for (size_t i = 0; i < grad_max_iter; ++i)
  {
    update_Ri();
    setdvec();

    P_Prime_vec = chol.solve(
      lamd_arap * VectorXf::Map(d_cur.data(), d_cur.cols()*d_cur.rows())
      + lamd_userCrsp * VectorXf::Map(bound_crsp.data(), bound_crsp.cols()*bound_crsp.rows())
      + lamd_userCrsp * VectorXf::Map(user_crsp.data(), user_crsp.cols()*user_crsp.rows()));

    P_Prime.row(0) = P_Prime_vec.segment(0, P_Num);
    P_Prime.row(1) = P_Prime_vec.segment(0+P_Num, P_Num);
    P_Prime.row(2) = P_Prime_vec.segment(0+2*P_Num, P_Num);
  }

  std::clock_t end = std::clock();
  std::cout<<"One opt step finished. Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\n";

  for (size_t i = 0; i != new_points->GetNumberOfPoints(); ++i) 
  {
    new_points->SetPoint(i, P_Prime_vec[i], P_Prime_vec[i+P_Num], P_Prime_vec[i+2*P_Num]);
  }
  new_points->Modified();
  mesh_data->resetMesh(mesh);
}

bool NonRigid::searchBound(int v_id, float v_crsp[3], double& bound_dist, double bound_thresh)
{
  Vector3d x;
  int idx[3];
  double pcoords[3];
  double f_dist[2][2][2];
  double step = 0.0;
  double step_length = 0.1;
  int search_iter = 0;
  int max_search_iter = 100;

  do 
  {
    step = step + step_length;//(out ? 0.05 : -0.05);
    x = (P_Prime.col(v_id) + step * P_Prime_N.col(v_id)).cast<double>();

    if (!scalar_img->ComputeStructuredCoordinates(x.data(), idx, pcoords))
    {
      return false;
    }

    for (size_t i = 0; i < 2; ++i)
    {
      for (size_t j = 0; j < 2; ++j)
      {
        for (size_t k = 0; k < 2; ++k)
        {
          f_dist[i][j][k] = scalar_img->GetScalarComponentAsDouble(idx[0]+i, idx[1]+j, idx[2]+k, 0);
        }
      }
    }

    ++search_iter;
  } while ((this->trilinearInterpolation(pcoords, f_dist)) >= bound_thresh && search_iter < max_search_iter);

  if (search_iter < max_search_iter)
  {
    // also means bound_dist <= 5
    bound_dist = step;
    v_crsp[0] = x[0];
    v_crsp[1] = x[1];
    v_crsp[2] = x[2];
    return true;
  }
  else return false;
}

int NonRigid::buildBoundCrsp2(MatrixX3f &bound_crsp, SparseMatrix<float>& L_cur)
{
  size_t P_Num = P_Prime.cols();
  bound_crsp = MatrixX3f::Zero(P_Num, 3);
  crsp_lines.clear();

  std::vector<bool> is_in_user_crsp(P_Num, false);
  for (size_t i = 0; i < user_v_ids.size(); ++i)
  {
    is_in_user_crsp[user_v_ids[i]] = true;
  }

  Vector3d x;
  std::vector<float> query(3);
  kdtree::KDTreeResultVector result;
  int idx[3];
  double pcoords[3];
  double f_dist[2][2][2];
  int num_outBoundCrsp = 0;
  int num_inBoundCrsp = 0;
  double cur_dist_val;
  double cur_crsp_dist;
  int T_idx;
  for (size_t v_id = 0; v_id < P_Num; ++v_id)
  {
    if (is_in_user_crsp[v_id] == true) continue;

    x = (P_Prime.col(v_id)).cast<double>();

    if (!scalar_img->ComputeStructuredCoordinates(x.data(), idx, pcoords))
    {
      continue;
    }

    for (size_t i = 0; i < 2; ++i)
    {
      for (size_t j = 0; j < 2; ++j)
      {
        for (size_t k = 0; k < 2; ++k)
        {
          f_dist[i][j][k] = scalar_img->GetScalarComponentAsDouble(idx[0]+i, idx[1]+j, idx[2]+k, 0);
        }
      }
    }

    cur_dist_val = this->trilinearInterpolation(pcoords, f_dist);
    T_idx = -1;

    if (cur_dist_val >= scalar_threshold) 
    {
      // it is inside of the bone
      // search closest point in same half plane
      //std::cout<<"begin search crsp inside the bone.\n";
      if (searchBound(v_id, &query[0], cur_crsp_dist, scalar_threshold)) 
      {
        tree_target->n_nearest(query, 1, result);
        T_idx = result[0].idx;
        //query[0] = x[0];
        //query[1] = x[1];
        //query[2] = x[2];
        //tree_target->r_nearest(query, cur_crsp_dist*cur_crsp_dist, result);

        //Vector3f cur_v_normal = P_Prime_N.col(v_id);
        //int temp_T_idx;
        //double cur_min_dist = numeric_limits<double>::max();
        //for (int result_id = 0; result_id < result.size(); ++result_id)
        //{
        //  temp_T_idx = result[result_id].idx;
        //  if (
        //    (cur_v_normal[0]*tree_normal_target[3*temp_T_idx + 0]
        //  + cur_v_normal[1]*tree_normal_target[3*temp_T_idx + 1]
        //  + cur_v_normal[2]*tree_normal_target[3*temp_T_idx + 2]) >= 0)
        //  {

        //    Vector3d cur_crsp(
        //      tree_target->the_data[temp_T_idx][0],
        //      tree_target->the_data[temp_T_idx][1],
        //      tree_target->the_data[temp_T_idx][2]);
        //    if ((cur_crsp - x).norm() < cur_min_dist)
        //    {
        //      cur_min_dist = (cur_crsp - x).norm();
        //      T_idx = temp_T_idx;
        //    }
        //  }
        //}
        if (T_idx != -1) ++num_inBoundCrsp;
        //else
        //{
          //std::cout<<"range: "<<cur_crsp_dist<<"\tpoints in range: "<<result.size()<<"\tcrsp dist: "<<cur_min_dist<<"\n";
        //}
      }
      //std::cout<<"found a feasible crsp (T_idx: "<<T_idx<< ") to exocranium.\n";
    }
    else
    {
      // it is outside of the bone
      // search closest point

      //std::cout<<"begin search closest crsp.\n";
      if (searchBound(v_id, &query[0], cur_crsp_dist, scalar_threshold)) 
      {
        tree_target->n_nearest(query, 1, result);
        T_idx = result[0].idx;
      } // deal with possible holes
      else
      {
        query[0] = x[0];
        query[1] = x[1];
        query[2] = x[2];
        tree_target->n_nearest(query, 1, result);
        T_idx = result[0].idx;
      }
      ++num_outBoundCrsp;
      //std::cout<<"found a feasible crsp (T_idx: "<<T_idx<< ") to closest boundary.\n";
    }

    if (T_idx != -1)
    {
      bound_crsp.row(v_id) = (
        tree_target->the_data[T_idx][0]*tree_normal_target[3*T_idx + 0]
      + tree_target->the_data[T_idx][1]*tree_normal_target[3*T_idx + 1]
      + tree_target->the_data[T_idx][2]*tree_normal_target[3*T_idx + 2])
        * RowVector3f(
        tree_normal_target[3*T_idx + 0],
        tree_normal_target[3*T_idx + 1],
        tree_normal_target[3*T_idx + 2]);

      for (int id_r = 0; id_r != 3; ++id_r)
      {
        for (int id_c = 0; id_c != 3; ++id_c)
        {
          L_cur.coeffRef(v_id + id_r*P_Num, v_id + id_c*P_Num) += 
            lamd_arap
            * tree_normal_target[3*T_idx + id_r]
          * tree_normal_target[3*T_idx + id_c];
        }
      }

      crsp_lines.push_back(x[0]);
      crsp_lines.push_back(x[1]);
      crsp_lines.push_back(x[2]);
      crsp_lines.push_back(tree_target->the_data[T_idx][0]);
      crsp_lines.push_back(tree_target->the_data[T_idx][1]);
      crsp_lines.push_back(tree_target->the_data[T_idx][2]);
    }
  }

  std::cout<<"search in: "<<num_inBoundCrsp<<"\tsearch out: "<<num_outBoundCrsp<<"\n";
  return num_inBoundCrsp + num_outBoundCrsp;
}

void NonRigid::refineOptStep2()
{
  if (!has_scalar_img || !has_mesh) return;
  size_t P_Num = adj_list.size();
  std::vector<int> v_ids;
  MatrixX3f bound_crsp = MatrixX3f::Zero(P_Num, 3);
  VectorXf P_Prime_vec;

  vtkSmartPointer<vtkPolyData> mesh = mesh_data->getMeshData();
  vtkSmartPointer<vtkPoints> new_points = mesh->GetPoints();


  computePPrimeNormal();
  SparseMatrix<float> L_cur = lamd_arap * L;
  int num_boundCrsp = buildBoundCrsp2(bound_crsp, L_cur);
  std::cout<<"Num of bound crsp: "<<num_boundCrsp<<"\n";
  if (vis_crsp_lines) return;
  std::clock_t begin = std::clock();

  for (size_t i = 0; i < user_v_ids.size(); ++i)
  {
    L_cur.coeffRef(user_v_ids[i], user_v_ids[i]) += lamd_userCrsp;
    L_cur.coeffRef(user_v_ids[i] + P_Num, user_v_ids[i] + P_Num) += lamd_userCrsp;
    L_cur.coeffRef(user_v_ids[i] + 2 * P_Num, user_v_ids[i] + 2 * P_Num) += lamd_userCrsp;
  }

  // dif from last shape
  //for (int i = 0; i < P_Num; ++i)
  //{
  //  L_cur.coeffRef(i, i) += lamd_arap / 3;
  //  L_cur.coeffRef(i + P_Num, i + P_Num) += lamd_arap / 3;
  //  L_cur.coeffRef(i + 2 * P_Num, i + 2 * P_Num) += lamd_arap / 3;
  //}
  //VectorXf P_Prime_last_vec(3*P_Num);
  //P_Prime_last_vec.segment(0, P_Num) = P_Prime.row(0);
  //P_Prime_last_vec.segment(0+P_Num, P_Num) = P_Prime.row(1);
  //P_Prime_last_vec.segment(0+2*P_Num, P_Num) = P_Prime.row(2);

  std::clock_t end = std::clock();
  double updateL = double(end-begin)/CLOCKS_PER_SEC;

  begin = std::clock();
  chol.factorize(L_cur);
  end = std::clock();
  double factorizeL = double(end-begin)/CLOCKS_PER_SEC;

  begin = std::clock();

  for (size_t i = 0; i < grad_max_iter; ++i)
  {
    update_Ri();
    setdvec();

    P_Prime_vec = chol.solve(
      lamd_arap * VectorXf::Map(d_cur.data(), d_cur.cols()*d_cur.rows())
      + lamd_arap * VectorXf::Map(bound_crsp.data(), bound_crsp.cols()*bound_crsp.rows())
      + lamd_userCrsp * VectorXf::Map(user_crsp.data(), user_crsp.cols()*user_crsp.rows())
      );//+ lamd_arap / 3 * P_Prime_last_vec);

    P_Prime.row(0) = P_Prime_vec.segment(0, P_Num);
    P_Prime.row(1) = P_Prime_vec.segment(0+P_Num, P_Num);
    P_Prime.row(2) = P_Prime_vec.segment(0+2*P_Num, P_Num);
  }

  end = std::clock();
  double solveShape = double(end-begin)/CLOCKS_PER_SEC;
  std::cout<<"One opt step finished. Update L: "<<updateL<<" Factorize L: "<<factorizeL<<" Solve Shape: "<<solveShape<<"\n";

  for (size_t i = 0; i != new_points->GetNumberOfPoints(); ++i) 
  {
    new_points->SetPoint(i, P_Prime_vec[i], P_Prime_vec[i+P_Num], P_Prime_vec[i+2*P_Num]);
  }
  new_points->Modified();
  mesh_data->resetMesh(mesh);
}

void NonRigid::setSubMesh(std::map<vtkIdType, int>& regionMap, std::vector<int>& regionBound)
{
  // the region is classified to 3 types
  // 1. center vertices
  // 2. vertices we want to update
  // 3. boundary
  // all vertices need to recompute its Ri every iteration
  // but only need to recompute 2's d matrix

  // A11X1 + A12X2 = D1
  // A11 is the submatrix of L contains 2
  // A12

  this->regionMap = regionMap;
  NewPtOrder.clear();
  NewPtNum = 0;
  std::map<vtkIdType, int>::iterator it_regionMap;
  std::vector<int> CenterPtIds;
  for (it_regionMap = regionMap.begin(); it_regionMap != regionMap.end(); ++it_regionMap)
  {
    if (it_regionMap->second == 0)
    {
      CenterPtIds.push_back(it_regionMap->first);
    }
    else
    {
      NewPtOrder.push_back(it_regionMap->first);
    }
  }
  NewPtNum = NewPtOrder.size();
  NewCenterNum = CenterPtIds.size();
  NewPtOrder.insert(NewPtOrder.end(), CenterPtIds.begin(), CenterPtIds.end());
  RUpdateList = regionBound;
  RUpdateList.insert(RUpdateList.end(), NewPtOrder.begin(), NewPtOrder.end());
  std::cout<<"New Pt Num: "<<NewPtNum<<"\n";
  std::cout<<"New Pt Num with center: "<<NewPtOrder.size()<<"\n";
  std::cout<<"Pt Num in region map: "<<regionMap.size()<<"\n";

  int P_Num = adj_list.size();
  for (int i = 0; i < P_Num; ++i)
  {
    it_regionMap = regionMap.find(i);
    if (it_regionMap == regionMap.end())
    {
      NewPtOrder.push_back(i);
    }
  }
  std::cout<<"Total P_Num: "<<P_Num<<" Size of New Pt Order: "<<NewPtOrder.size()<<"\n";
  OldToNewPtOrder.resize(P_Num, 0);
  for (int i = 0; i < P_Num; ++i)
  {
    OldToNewPtOrder[NewPtOrder[i]] = i;
  }

  SparseMatrix<float> L_Rearrange(3*P_Num, 3*P_Num);
  vector<Triplet<float>> elem_L_SubMesh;
  for (int k = 0; k < L.outerSize(); ++k)
  {
    for (SparseMatrix<float>::InnerIterator it(L, k); it; ++it)
    {
      int NewRowId = OldToNewPtOrder[it.row() % P_Num];
      int NewColId = OldToNewPtOrder[it.col() % P_Num];

      elem_L_SubMesh.push_back(Triplet<float>(
        3*NewRowId + it.row()/P_Num, 
        3*NewColId + it.col()/P_Num,
        it.value()));
    }
  }
  L_Rearrange.setFromTriplets(elem_L_SubMesh.begin(), elem_L_SubMesh.end());

  L_Submesh = L_Rearrange.block(0,0, 3*NewPtNum, 3*NewPtNum);
  L_SubmeshAux = L_Rearrange.block(0, 3*NewPtNum, 3*NewPtNum, 3*(P_Num-NewPtNum));

  cholSubMesh.analyzePattern(L_Submesh);
}

void NonRigid::refineSubMesh(std::vector<double>& centerPos)
{
  if (!has_scalar_img || !has_mesh) return;

  vtkSmartPointer<vtkPolyData> mesh = mesh_data->getMeshData();
  vtkSmartPointer<vtkPoints> new_points = mesh->GetPoints();


  computePPrimeNormal();
  SparseMatrix<float> L_cur = lamd_arap * L_Submesh;
  Matrix3Xf bound_crsp = Matrix3Xf::Zero(3, NewPtNum);
  Matrix3Xf user_crsp_sub = Matrix3Xf::Zero(3, NewPtNum);
  int num_boundCrsp = 0;
  if (!local_smooth_mode)
  {
    num_boundCrsp = buildSubBoundCrsp(bound_crsp, L_cur);
  }
  std::cout<<"Num of bound crsp: "<<num_boundCrsp<<"\n";
  if (vis_crsp_lines) return;

  std::clock_t begin = std::clock();

  std::map<vtkIdType, int>::iterator it_regionMap;
  for (size_t i = 0; i < user_v_ids.size(); ++i)
  {
    int NewPtIdx = OldToNewPtOrder[user_v_ids[i]];
    if (NewPtIdx < NewPtNum)
    {
      L_cur.coeffRef(3*NewPtIdx + 0, 3*NewPtIdx + 0) += lamd_userCrsp;
      L_cur.coeffRef(3*NewPtIdx + 1, 3*NewPtIdx + 1) += lamd_userCrsp;
      L_cur.coeffRef(3*NewPtIdx + 2, 3*NewPtIdx + 2) += lamd_userCrsp;
      user_crsp_sub.col(NewPtIdx) = user_crsp.row(user_v_ids[i]).transpose();
    }
  }

  // dif from last shape
  for (int i = 0; i < NewPtNum; ++i)
  {
    L_cur.coeffRef(3*i + 0, 3*i + 0) += lamd_arap / 3;
    L_cur.coeffRef(3*i + 1, 3*i + 1) += lamd_arap / 3;
    L_cur.coeffRef(3*i + 2, 3*i + 2) += lamd_arap / 3;
  }

  std::clock_t end = std::clock();
  double updateL = double(end-begin)/CLOCKS_PER_SEC;

  begin = std::clock();
  cholSubMesh.factorize(L_cur);
  end = std::clock();
  double factorizeL = double(end-begin)/CLOCKS_PER_SEC;

  begin = std::clock();

  VectorXf P_Prime_vec;
  VectorXf P_Prime_vec_full(3*NewPtOrder.size());
  for (int i = 0; i < NewPtOrder.size(); ++i)
  {
    P_Prime_vec_full[3*i + 0] = P_Prime(0, NewPtOrder[i]);
    P_Prime_vec_full[3*i + 1] = P_Prime(1, NewPtOrder[i]);
    P_Prime_vec_full[3*i + 2] = P_Prime(2, NewPtOrder[i]);
  }
  for (int i = 0; i < NewCenterNum; ++i)
  {
    P_Prime_vec_full[3*(i+NewPtNum) + 0] = centerPos[3*i + 0];
    P_Prime_vec_full[3*(i+NewPtNum) + 1] = centerPos[3*i + 1];
    P_Prime_vec_full[3*(i+NewPtNum) + 2] = centerPos[3*i + 2];
  }

  Matrix3Xf d_cur_sub(3, NewPtNum);
  VectorXf aux_submesh = L_SubmeshAux * P_Prime_vec_full.segment(3*NewPtNum, 3*(NewPtOrder.size()-NewPtNum));

  // dif from last shape
  VectorXf P_Prime_last_vec = P_Prime_vec_full.segment(0, 3*NewPtNum);

  for (size_t i = 0; i < grad_max_iter; ++i)
  {
    update_SubRi();
    setSubdvec(d_cur_sub);

    P_Prime_vec = cholSubMesh.solve(
      lamd_arap * (VectorXf::Map(d_cur_sub.data(), d_cur_sub.cols()*d_cur_sub.rows())
      - aux_submesh)
      + lamd_arap * VectorXf::Map(bound_crsp.data(), bound_crsp.cols()*bound_crsp.rows())
      + lamd_userCrsp * VectorXf::Map(user_crsp_sub.data(), user_crsp_sub.cols()*user_crsp_sub.rows())
      + lamd_arap / 3 * P_Prime_last_vec);

    P_Prime_vec_full.segment(0, 3*NewPtNum) = P_Prime_vec;

    for (int j = 0; j < NewPtNum + NewCenterNum; ++j)
    {
      P_Prime(0, NewPtOrder[j]) = P_Prime_vec_full[3*j + 0];
      P_Prime(1, NewPtOrder[j]) = P_Prime_vec_full[3*j + 1];
      P_Prime(2, NewPtOrder[j]) = P_Prime_vec_full[3*j + 2];
    }
  }

  end = std::clock();
  double solveShape = double(end-begin)/CLOCKS_PER_SEC;
  std::cout<<"One opt step finished. Update L: "<<updateL<<" Factorize L: "<<factorizeL<<" Solve Shape: "<<solveShape<<"\n";

  for (size_t i = 0; i != new_points->GetNumberOfPoints(); ++i) 
  {
    new_points->SetPoint(i, P_Prime(0,i), P_Prime(1,i), P_Prime(2,i));
  }
  mesh_data->getActor()->GetMapper()->Modified();
}

int NonRigid::buildSubBoundCrsp(Matrix3Xf &bound_crsp, SparseMatrix<float>& L_cur)
{
  bound_crsp = Matrix3Xf::Zero(3, NewPtNum);
  crsp_lines.clear();

  std::vector<bool> is_in_user_crsp(NewPtNum, false);
  for (size_t i = 0; i < user_v_ids.size(); ++i)
  {
    if (OldToNewPtOrder[user_v_ids[i]] < NewPtNum)
    {
      is_in_user_crsp[OldToNewPtOrder[user_v_ids[i]]] = true;
    }
  }

  Vector3d x;
  std::vector<float> query(3);
  kdtree::KDTreeResultVector result;
  int idx[3];
  double pcoords[3];
  double f_dist[2][2][2];
  int num_outBoundCrsp = 0;
  int num_inBoundCrsp = 0;
  double cur_dist_val;
  double cur_crsp_dist;
  int T_idx;
  for (size_t v_id = 0; v_id < NewPtNum; ++v_id)
  {
    if (is_in_user_crsp[v_id] == true) continue;

    x = (P_Prime.col(NewPtOrder[v_id])).cast<double>();

    if (!scalar_img->ComputeStructuredCoordinates(x.data(), idx, pcoords))
    {
      continue;
    }

    for (size_t i = 0; i < 2; ++i)
    {
      for (size_t j = 0; j < 2; ++j)
      {
        for (size_t k = 0; k < 2; ++k)
        {
          f_dist[i][j][k] = scalar_img->GetScalarComponentAsDouble(idx[0]+i, idx[1]+j, idx[2]+k, 0);
        }
      }
    }

    cur_dist_val = this->trilinearInterpolation(pcoords, f_dist);
    T_idx = -1;

    if (cur_dist_val >= scalar_threshold) 
    {
      if (searchBound(NewPtOrder[v_id], &query[0], cur_crsp_dist, scalar_threshold)) 
      {
        tree_target->n_nearest(query, 1, result);
        T_idx = result[0].idx;

        if (T_idx != -1) ++num_inBoundCrsp;

      }
    }
    else
    {
      // it is outside of the bone
      // search closest point

      //std::cout<<"begin search closest crsp.\n";
      if (searchBound(NewPtOrder[v_id], &query[0], cur_crsp_dist, scalar_threshold)) 
      {
        tree_target->n_nearest(query, 1, result);
        T_idx = result[0].idx;
      } // deal with possible holes
      else
      {
        query[0] = x[0];
        query[1] = x[1];
        query[2] = x[2];
        tree_target->n_nearest(query, 1, result);
        T_idx = result[0].idx;
      }
      ++num_outBoundCrsp;
      //std::cout<<"found a feasible crsp (T_idx: "<<T_idx<< ") to closest boundary.\n";
    }

    if (T_idx != -1)
    {
      bound_crsp.col(v_id) = (
        tree_target->the_data[T_idx][0]*tree_normal_target[3*T_idx + 0]
      + tree_target->the_data[T_idx][1]*tree_normal_target[3*T_idx + 1]
      + tree_target->the_data[T_idx][2]*tree_normal_target[3*T_idx + 2])
        * Vector3f(
        tree_normal_target[3*T_idx + 0],
        tree_normal_target[3*T_idx + 1],
        tree_normal_target[3*T_idx + 2]);

      for (int id_r = 0; id_r != 3; ++id_r)
      {
        for (int id_c = 0; id_c != 3; ++id_c)
        {
          L_cur.coeffRef(3*v_id + id_r, 3*v_id + id_c) += 
            lamd_arap
            * tree_normal_target[3*T_idx + id_r]
          * tree_normal_target[3*T_idx + id_c];
        }
      }

      crsp_lines.push_back(x[0]);
      crsp_lines.push_back(x[1]);
      crsp_lines.push_back(x[2]);
      crsp_lines.push_back(tree_target->the_data[T_idx][0]);
      crsp_lines.push_back(tree_target->the_data[T_idx][1]);
      crsp_lines.push_back(tree_target->the_data[T_idx][2]);
    }
  }

  std::cout<<"search in: "<<num_inBoundCrsp<<"\tsearch out: "<<num_outBoundCrsp<<"\n";
  return num_inBoundCrsp + num_outBoundCrsp;
}

void NonRigid::setSubdvec(Matrix3Xf& d_sub)
{
  d_sub = Matrix3Xf::Zero(3, NewPtNum);
  for (int i = 0; i < NewPtNum; ++i)
  {
    for (size_t j = 0; j < adj_list[NewPtOrder[i]].size(); ++j)
    {
      d_sub.col(i) += ((Weight.coeffRef(NewPtOrder[i], adj_list[NewPtOrder[i]][j])/2)*(R[NewPtOrder[i]]+R[adj_list[NewPtOrder[i]][j]])*(P.col(NewPtOrder[i]) - P.col(adj_list[NewPtOrder[i]][j])));
    }
  }
}

void NonRigid::update_SubRi()
{
  Matrix3f Si;
  MatrixXf Di;
  Matrix3Xf Pi_Prime;
  Matrix3Xf Pi;
  for (int i = 0; i < RUpdateList.size(); ++i) 
  {
    Di = MatrixXf::Zero(adj_list[RUpdateList[i]].size(), adj_list[RUpdateList[i]].size());
    Pi_Prime.resize(3, adj_list[RUpdateList[i]].size());
    Pi.resize(3, adj_list[RUpdateList[i]].size());
    // if there is not any single unconnected point this for loop can have a more efficient representation
    for (decltype(adj_list[RUpdateList[i]].size()) j = 0; j != adj_list[RUpdateList[i]].size(); ++j) {
      Di(j, j) = Weight.coeffRef(RUpdateList[i], adj_list[RUpdateList[i]][j]);
      Pi.col(j) = P.col(RUpdateList[i]) - P.col(adj_list[RUpdateList[i]][j]);
      Pi_Prime.col(j) = P_Prime.col(RUpdateList[i]) - P_Prime.col(adj_list[RUpdateList[i]][j]);
    }
    Si = Pi * Di * Pi_Prime.transpose();
    Matrix3f Ui;
    Vector3f Wi;
    Matrix3f Vi;
    wunderSVD3x3(Si, Ui, Wi, Vi);
    R[RUpdateList[i]] = Vi * Ui.transpose();
  }
}
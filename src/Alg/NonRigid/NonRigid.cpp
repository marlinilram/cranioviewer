#include "NonRigid.h"

NonRigid::NonRigid()
{
  tree_target = nullptr;
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
        dist_e += computeVertexGrad(cur_v, cur_g);
        g_vec(i + 0*P_Num) = cur_g(0);
        g_vec(i + 1*P_Num) = cur_g(1);
        g_vec(i + 2*P_Num) = cur_g(2);
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

void NonRigid::setDistMap(vtkSmartPointer<vtkImageData> img)
{
    dist_map = img;
    std::cout<<"dist map scalar size: "<<dist_map->GetNumberOfScalarComponents()<<"\n";

    vtkSmartPointer<vtkImageGradient> g_filter = vtkSmartPointer<vtkImageGradient>::New();

    g_filter->SetInputData(dist_map);
    g_filter->SetDimensionality(3);
    g_filter->Update();

    dist_gradient =(g_filter->GetOutput());

    std::cout<<"dist gradient scalar size: "<<dist_gradient->GetNumberOfScalarComponents()<<"\n";
}

double NonRigid::computeVertexGrad(Vector3f &v, Vector3f &grad)
{
    // 1. get the cell contains vertex and its 8 corner points
    // 2. compute gradients for this 8 points
    // 3. compute trilinear interpolation

    double x[3] = {v[0], v[1], v[2]};
    int idx[3];
    double pcoords[3];
    dist_map->ComputeStructuredCoordinates(x, idx, pcoords);

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
        f_cur_val = computeGradEnergy(p_vec, g_vec);
        //end = std::clock();
        //std::cout<<"Computed Gradient finished. Elapsed time: "<<double(end-begin)/CLOCKS_PER_SEC<<"\nfunction value: "<<f_val<<"delta: "<<f_val-f_cur_val<<"\n";
        if (f_val - f_cur_val < 1e-3) break;
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
        if (wi < 0.1) cout << "Edge Weight Sum Warning: " << wi << endl;
        weight_sum.push_back(Triplet<float>(i, i, wi));
        weight_sum.push_back(Triplet<float>(i+P_Num, i+P_Num, wi));
        weight_sum.push_back(Triplet<float>(i+2*P_Num, i+2*P_Num, wi));
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

    P_Prime_N.setOnes();

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

void NonRigid::computeInflateDir(VectorXf &p_vec, VectorXf &g_vec)
{
    VectorXf arap_g = g_vec;
    VectorXf inflate_g = g_vec;
    VectorXf userCrsp_g = g_vec;
    double f_e = 0;
    computeArap(p_vec, arap_g);
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

    g_vec = lamd_inflate*inflate_g + lamd_arap*arap_g + lamd_userCrsp*userCrsp_g;
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

}
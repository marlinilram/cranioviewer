#include "ICPWrapper.h"

ICPWrapper::ICPWrapper()
{
    icp_pt_plane = nullptr;
    icp_pt_pt = nullptr;
    temp_mesh = nullptr;
    nii_img = nullptr;
    max_iter = 100;

    icp_trans_mat = vtkSmartPointer<vtkMatrix4x4>::New();
    icp_trans = vtkSmartPointer<vtkTransform>::New();
    icp_trans->PostMultiply();
}

ICPWrapper::~ICPWrapper()
{
    delete icp_pt_pt;
    delete icp_pt_plane;
}

void ICPWrapper::setTempData(Mesh *mesh)
{
    temp_mesh = mesh;

    // init mesh data vector to current position
    vtkSmartPointer<vtkMatrix4x4> trans_mat = vtkSmartPointer<vtkMatrix4x4>::New();
    temp_mesh->getTransform()->GetMatrix(trans_mat);
    int meshNum = temp_mesh->getMeshData()->GetNumberOfPoints();
    temp_mesh_vec.clear();
    for (size_t i = 0; i < meshNum; ++i) {
        double* temp = temp_mesh->getMeshData()->GetPoint(i);
        temp_mesh_vec.push_back(temp[0]);
        temp_mesh_vec.push_back(temp[1]);
        temp_mesh_vec.push_back(temp[2]);
    }
    
    setUserTrans(trans_mat);
}

void ICPWrapper::setUserTrans(vtkSmartPointer<vtkMatrix4x4> mat)
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

void ICPWrapper::setTargetData(std::vector<double> &data)
{
    if (icp_pt_pt)
    {
        delete icp_pt_pt;
    }

    icp_pt_pt = new IcpPointToPoint(&data[0], data.size()/3, 3);
}

void ICPWrapper::setImage(NiiLoader *img)
{
    nii_img = img;
}

void ICPWrapper::runICP()
{
    int iso_val = 120;
    int iso_width = 20 ;
    int voxel_num = 0;
    setTargetData(nii_img->extractSkullVertex(iso_val, iso_width, voxel_num));

    if (icp_pt_pt->checkWorkable()) {

        vector<int32_t> active;
        for (int32_t i=0; i<temp_mesh_vec.size()/3; i++) active.push_back(i);

        double last_error = std::numeric_limits<double>::max();
        double cur_error = -1;
        //int max_iter = 100;
        double min_delta = 1e-3;

        Matrix R = Matrix::eye(3);
        Matrix t(3,1);
        icp_trans_mat->Identity();
        temp_mesh->getTransform()->Concatenate(icp_trans);

        for (int32_t iter=0; iter<max_iter; iter++) {

            cur_error = icp_pt_pt->fitStep4Public(&temp_mesh_vec[0], temp_mesh_vec.size()/3, R, t, active);
            cout << "cur delta = " << cur_error << "\t" << "the " << iter << "th iteration" << endl;
            //cout<<"R: "<<R<<"\n"<<"t: "<<t<<"\n";

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) icp_trans_mat->SetElement(i, j, R.val[i][j]);
            }
            for (int i = 0; i < 3; ++i) icp_trans_mat->SetElement(i, 3, t.val[i][0]);
            icp_trans->Identity();
            icp_trans->Concatenate(icp_trans_mat);
            temp_mesh->reloadTransform();
            emit(updateRenderers());
            //system("pause");
            if (abs(last_error - cur_error) < min_delta)
                break;
            last_error = cur_error;
        }

        // give the transform to model and clear up all transform in UI
        setUserTrans(icp_trans_mat);

        vtkSmartPointer<vtkPolyData> mesh_data = temp_mesh->getMeshData();
        vtkSmartPointer<vtkPoints> new_points = mesh_data->GetPoints();
        if (new_points->GetNumberOfPoints() == temp_mesh_vec.size()/3)
        {
            for (size_t i = 0; i != new_points->GetNumberOfPoints(); ++i) {
                new_points->SetPoint(i, temp_mesh_vec[3*i], temp_mesh_vec[3*i+1], temp_mesh_vec[3*i+2]);
            }
            new_points->Modified();
            temp_mesh->resetMesh(mesh_data);
            temp_mesh->setMeshCenter();
        }
        else std::cout<<"Error: number of points doesn't match\n";

        emit(resetTrans());
    }

    delete icp_pt_pt;
    icp_pt_pt = nullptr;

}

void ICPWrapper::runICPStep()
{

}
#include "ShapeMatch.h"

ShapeMatch::ShapeMatch()
{

}

ShapeMatch::~ShapeMatch()
{

}

void ShapeMatch::setTarget(std::vector<double> &target_vec)
{
    Eigen::Map<Eigen::Matrix3Xd> target_mat(&target_vec[0], 3, target_vec.size()/3);

    target_center = Eigen::Vector3d(target_mat.row(0).mean(), target_mat.row(0).mean(), target_mat.row(0).mean());
    target_mat.row(0) = (target_mat.row(0).array() - target_center(0)).matrix();
    target_mat.row(1) = (target_mat.row(1).array() - target_center(1)).matrix();
    target_mat.row(2) = (target_mat.row(2).array() - target_center(2)).matrix();
}

void ShapeMatch::setTemplate(std::vector<double> &template_vec)
{
    Eigen::Map<Eigen::Matrix3Xd> template_mat(&template_vec[0], 3, template_vec.size()/3);

    template_center = Eigen::Vector3d(template_mat.row(0).mean(), template_mat.row(1).mean(), template_mat.row(2).mean());
    template_mat.row(0) = (template_mat.row(0).array() - template_center(0)).matrix();
    template_mat.row(1) = (template_mat.row(1).array() - template_center(1)).matrix();
    template_mat.row(2) = (template_mat.row(2).array() - template_center(2)).matrix();
}

Eigen::Matrix3d ShapeMatch::computeA()
{
    Eigen::MatrixXd C_n = Eigen::MatrixXd::Identity(target_mat.cols(), target_mat.cols()) - Eigen::MatrixXd::Ones(target_mat.cols(), target_mat.cols()) / target_mat.cols();
    Eigen::Matrix3d A_pq = target_mat*C_n*template_mat.transpose();
    Eigen::Matrix3d A_qq = template_mat*C_n*template_mat.transpose();

    return A_pq*A_qq;
}

Eigen::Vector3d ShapeMatch::computeT()
{
    return target_center - template_center;
}
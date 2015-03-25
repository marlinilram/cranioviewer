#ifndef ShapeMatch_H
#define ShapeMatch_H

#include "Eigen\Eigen"
#include <vector>
#include <iostream>
#include <fstream>

#include "WunderSVD3x3.h"

class ShapeMatch
{
public:
    ShapeMatch();
    ~ShapeMatch();

    void setTemplate(std::vector<double> &template_vec);
    void setTarget(std::vector<double> &target_vec);

    Eigen::Matrix3d computeA();
    Eigen::Vector3d computeT();

private:
    Eigen::Matrix3Xd template_mat;
    Eigen::Matrix3Xd target_mat;
    Eigen::Vector3d template_center;
    Eigen::Vector3d target_center;
};

#endif
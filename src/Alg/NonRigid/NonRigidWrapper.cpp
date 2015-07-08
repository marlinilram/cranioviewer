#include "NonRigidWrapper.h"
#include <vtkLineSource.h>

NonRigidWrapper::NonRigidWrapper()
{
    non_rigid = new NonRigid;
}

NonRigidWrapper::~NonRigidWrapper()
{
    delete non_rigid;
}

void NonRigidWrapper::drawCrspLines()
{
  clearCrspLines();

  std::vector<double>& crsp_lines = non_rigid->getCrspLines();
  std::cout<<"crsp_lines: "<<crsp_lines.size()/6<<"\n";
  for (size_t i = 0; i < crsp_lines.size()/6; ++i)
  {
    vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
    lineSource->SetPoint1(&crsp_lines[i*6]);
    lineSource->SetPoint2(&crsp_lines[i*6+3]);
    lineSource->Update();

    // Visualize
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(lineSource->GetOutputPort());
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetLineWidth(2);
    actor->GetProperty()->SetDiffuseColor(0,1,0);
    crsp_lines_actors.push_back(actor);
    renderer->AddActor(crsp_lines_actors.back());
  }
}

void NonRigidWrapper::clearCrspLines()
{
  for (size_t i = 0; i < crsp_lines_actors.size(); ++i)
  {
    renderer->RemoveActor(crsp_lines_actors[i]);
  }
  crsp_lines_actors.clear();
}

void NonRigidWrapper::initMCKDTree(int isoval)
{
  isoval = (isoval == 0) ? 150 : isoval;
  std::vector<double> voxel_data;// = nii_img->extractSkullVertex(isovalue, isowidth, voxel_num);
  std::vector<double>& voxel_norm = non_rigid->getTreeNormal();
  main_viewer->getImgData()->mcSkullVertex(voxel_data, voxel_norm, isoval, 1);

  non_rigid->buildKDTree(voxel_data);
}
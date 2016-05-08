#include "RegionMoveStyle.h"

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkPicker.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>

#include "Mesh.h"

vtkStandardNewMacro(RegionMoveStyle);

RegionMoveStyle::RegionMoveStyle()
{
  this->StartPosition[0] = this->StartPosition[1] = 0;
  this->EndPosition[0] = this->EndPosition[1] = 0;
  this->Moving = 0;
  this->MotionVec[0] = this->MotionVec[1] = this->MotionVec[2] = 0;
  this->LastMotionVec[0] = this->LastMotionVec[1] = this->LastMotionVec[2] = 0;
  this->StartDPosition[0] = this->StartDPosition[1] = this->StartDPosition[2] = 0; // in display
  this->StartWPosition[0] = this->StartWPosition[1] = this->StartWPosition[2] = 0;
}

RegionMoveStyle::~RegionMoveStyle()
{

}

void RegionMoveStyle::OnMouseMove()
{
  if (!this->Interactor || !this->Moving)
  {
    return;
  }

  // pressed and move
  int NewEndPosition[2];
  NewEndPosition[0] = this->Interactor->GetEventPosition()[0];
  NewEndPosition[1] = this->Interactor->GetEventPosition()[1];
  int *size = this->Interactor->GetRenderWindow()->GetSize();
  if (NewEndPosition[0] > (size[0]-1))
  {
    NewEndPosition[0] = size[0]-1;
  }
  if (NewEndPosition[0] < 0)
  {
    NewEndPosition[0] = 0;
  }
  if (NewEndPosition[1] > (size[1]-1))
  {
    NewEndPosition[1] = size[1]-1;
  }
  if (NewEndPosition[1] < 0)
  {
    NewEndPosition[1] = 0;
  }

  if (sqrt((NewEndPosition[0]-this->EndPosition[0])*(NewEndPosition[0]-this->EndPosition[0])
    + (NewEndPosition[1]-this->EndPosition[1])*(NewEndPosition[1]-this->EndPosition[1])) > 3)
  {
    this->EndPosition[0] = NewEndPosition[0];
    this->EndPosition[1] = NewEndPosition[1];

    double NewWPosition[4];
    vtkRenderer* ren = this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
    this->ComputeDisplayToWorld(ren,
      this->StartDPosition[0] + this->EndPosition[0] - this->StartPosition[0],
      this->StartDPosition[1] + this->EndPosition[1] - this->StartPosition[1],
      this->StartDPosition[2], 
      NewWPosition);

    MotionVec[0] = NewWPosition[0] - StartWPosition[0];
    MotionVec[1] = NewWPosition[1] - StartWPosition[1];
    MotionVec[2] = NewWPosition[2] - StartWPosition[2];
  }
}

void RegionMoveStyle::OnLeftButtonDown()
{
  if (!this->Interactor)
  {
    return;
  }

  this->Moving = 1;



  this->StartPosition[0] = this->Interactor->GetEventPosition()[0];
  this->StartPosition[1] = this->Interactor->GetEventPosition()[1];
  this->EndPosition[0] = this->StartPosition[0];
  this->EndPosition[1] = this->StartPosition[1];

  this->InvokeEvent(vtkCommand::StartInteractionEvent);
}

void RegionMoveStyle::OnLeftButtonUp()
{
  if (!this->Interactor || !this->Moving)
  {
    return;
  }

  this->Moving = 0;
  this->InvokeEvent(vtkCommand::EndInteractionEvent);
}

void RegionMoveStyle::SetStartWPosition(double StartWPosition[3])
{
  this->StartWPosition[0] = StartWPosition[0];
  this->StartWPosition[1] = StartWPosition[1];
  this->StartWPosition[2] = StartWPosition[2];

  vtkRenderWindow *renWin = this->Interactor->GetRenderWindow();
  vtkRenderer* ren = renWin->GetRenderers()->GetFirstRenderer();

  this->ComputeWorldToDisplay(ren,
    StartWPosition[0],
    StartWPosition[1],
    StartWPosition[2],
    this->StartDPosition);
}

double* RegionMoveStyle::GetLastMotionVec()
{
  return LastMotionVec;
}

double* RegionMoveStyle::GetMotionVec()
{
  return MotionVec;
}

int RegionMoveStyle::GetMovingTag()
{
   return this->Moving;
}
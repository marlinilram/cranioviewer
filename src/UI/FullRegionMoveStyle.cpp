#include "FullRegionMoveStyle.h"

#include <vtkObjectFactory.h>
#include "vtkCamera.h"
#include "vtkCellPicker.h"
#include "vtkCallbackCommand.h"
#include <vtkCommand.h>
#include "vtkMath.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkProp3D.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include <vtkSmartPointer.h>
#include "vtkTransform.h"
#include <vtkAxesActor.h>

vtkStandardNewMacro(FullRegionMoveStyle);

FullRegionMoveStyle::FullRegionMoveStyle()
{

}

FullRegionMoveStyle::~FullRegionMoveStyle()
{

}

//----------------------------------------------------------------------------
void FullRegionMoveStyle::OnMouseMove()
{
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  switch (this->State)
  {
  case VTKIS_ROTATE:
    this->FindPokedRenderer(x, y);
    this->Rotate();
    this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
    break;

  case VTKIS_PAN:
    this->FindPokedRenderer(x, y);
    this->Pan();
    this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
    break;

  case VTKIS_DOLLY:
    this->FindPokedRenderer(x, y);
    this->Dolly();
    this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
    break;

  case VTKIS_SPIN:
    this->FindPokedRenderer(x, y);
    this->Spin();
    this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
    break;

  case VTKIS_USCALE:
    this->FindPokedRenderer(x, y);
    this->UniformScale();
    this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
    break;
  }
}

//----------------------------------------------------------------------------
void FullRegionMoveStyle::OnLeftButtonDown()
{
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  this->FindPokedRenderer(x, y);
  this->FindPickedActor(x, y);
  if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
  {
    return;
  }

  //this->GrabFocus(this->EventCallbackCommand);
  if (this->Interactor->GetShiftKey())
  {
    this->StartPan();
  }
  else if (this->Interactor->GetControlKey())
  {
    this->StartSpin();
  }
  else
  {
    this->StartRotate();
  }
}

//----------------------------------------------------------------------------
void FullRegionMoveStyle::OnLeftButtonUp()
{
  switch (this->State)
  {
  case VTKIS_PAN:
    this->EndPan();
    break;

  case VTKIS_SPIN:
    this->EndSpin();
    break;

  case VTKIS_ROTATE:
    this->EndRotate();
    break;
  }

  if ( this->Interactor )
  {
    this->ReleaseFocus();
  }
}

//----------------------------------------------------------------------------
void FullRegionMoveStyle::OnRightButtonDown()
{
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  this->FindPokedRenderer(x, y);
  this->FindPickedActor(x, y);
  if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
  {
    return;
  }

  //this->GrabFocus(this->EventCallbackCommand);
  if (this->Interactor->GetControlKey())
  {
    this->StartDolly();
  }
  else
  {
    this->StartPan();
  }
}

//----------------------------------------------------------------------------
void FullRegionMoveStyle::OnRightButtonUp()
{
  switch (this->State)
  {
  case VTKIS_DOLLY:
    this->EndDolly();
    break;

  case VTKIS_PAN:
    this->EndPan();
    break;
  }

  if ( this->Interactor )
  {
    this->ReleaseFocus();
  }
}

//----------------------------------------------------------------------------
void FullRegionMoveStyle::OnMiddleButtonDown()
{
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  this->FindPokedRenderer(x, y);
  this->FindPickedActor(x, y);
  if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
  {
    return;
  }

  this->GrabFocus(this->EventCallbackCommand);
  this->StartUniformScale();
}

//----------------------------------------------------------------------------
void FullRegionMoveStyle::OnMiddleButtonUp()
{
  switch (this->State)
  {
  case VTKIS_USCALE:
    this->EndUniformScale();
    break;
  }

  if ( this->Interactor )
  {
    this->ReleaseFocus();
  }
}

//----------------------------------------------------------------------------
void FullRegionMoveStyle::Rotate()
{
  if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
  {
    return;
  }

  vtkRenderWindowInteractor *rwi = this->Interactor;
  vtkCamera *cam = this->CurrentRenderer->GetActiveCamera();

  // First get the origin of the assembly
  double *obj_center = vtkTransform::SafeDownCast(this->InteractionProp->GetUserTransform())->GetPosition();

  // GetLength gets the length of the diagonal of the bounding box
  double boundRadius = this->InteractionProp->GetLength() * 0.1;

  // Get the view up and view right vectors
  double view_up[3], view_look[3], view_right[3];

  cam->OrthogonalizeViewUp();
  cam->ComputeViewPlaneNormal();
  cam->GetViewUp(view_up);
  vtkMath::Normalize(view_up);
  cam->GetViewPlaneNormal(view_look);
  vtkMath::Cross(view_up, view_look, view_right);
  vtkMath::Normalize(view_right);

  // Get the furtherest point from object position+origin
  double outsidept[3];

  outsidept[0] = obj_center[0] + view_right[0] * boundRadius;
  outsidept[1] = obj_center[1] + view_right[1] * boundRadius;
  outsidept[2] = obj_center[2] + view_right[2] * boundRadius;

  // Convert them to display coord
  double disp_obj_center[3];

  this->ComputeWorldToDisplay(obj_center[0], obj_center[1], obj_center[2],
    disp_obj_center);

  this->ComputeWorldToDisplay(outsidept[0], outsidept[1], outsidept[2],
    outsidept);

  double radius = sqrt(vtkMath::Distance2BetweenPoints(disp_obj_center,
    outsidept));
  double nxf = (rwi->GetEventPosition()[0] - disp_obj_center[0]) / radius;

  double nyf = (rwi->GetEventPosition()[1] - disp_obj_center[1]) / radius;

  double oxf = (rwi->GetLastEventPosition()[0] - disp_obj_center[0]) / radius;

  double oyf = (rwi->GetLastEventPosition()[1] - disp_obj_center[1]) / radius;

  if (((nxf * nxf + nyf * nyf) <= 1.0) &&
    ((oxf * oxf + oyf * oyf) <= 1.0))
  {
    double newXAngle = vtkMath::DegreesFromRadians( asin( nxf ) );
    double newYAngle = vtkMath::DegreesFromRadians( asin( nyf ) );
    double oldXAngle = vtkMath::DegreesFromRadians( asin( oxf ) );
    double oldYAngle = vtkMath::DegreesFromRadians( asin( oyf ) );

    double scale[3];
    scale[0] = scale[1] = scale[2] = 1.0;

    double **rotate = new double*[2];

    rotate[0] = new double[4];
    rotate[1] = new double[4];

    rotate[0][0] = newXAngle - oldXAngle;
    rotate[0][1] = view_up[0];
    rotate[0][2] = view_up[1];
    rotate[0][3] = view_up[2];

    rotate[1][0] = oldYAngle - newYAngle;
    rotate[1][1] = view_right[0];
    rotate[1][2] = view_right[1];
    rotate[1][3] = view_right[2];


    this->Prop3DTransform(this->InteractionProp,
      obj_center,
      2,
      rotate,
      scale);

    delete [] rotate[0];
    delete [] rotate[1];
    delete [] rotate;

    if (this->AutoAdjustCameraClippingRange)
    {
      this->CurrentRenderer->ResetCameraClippingRange();
    }

    rwi->Render();
  }
}

//----------------------------------------------------------------------------
void FullRegionMoveStyle::Spin()
{
  if ( this->CurrentRenderer == NULL || this->InteractionProp == NULL )
  {
    return;
  }

  vtkRenderWindowInteractor *rwi = this->Interactor;
  vtkCamera *cam = this->CurrentRenderer->GetActiveCamera();

  // Get the axis to rotate around = vector from eye to origin

  double *obj_center = this->InteractionProp->GetCenter();

  double motion_vector[3];
  double view_point[3];

  if (cam->GetParallelProjection())
  {
    // If parallel projection, want to get the view plane normal...
    cam->ComputeViewPlaneNormal();
    cam->GetViewPlaneNormal( motion_vector );
  }
  else
  {
    // Perspective projection, get vector from eye to center of actor
    cam->GetPosition( view_point );
    motion_vector[0] = view_point[0] - obj_center[0];
    motion_vector[1] = view_point[1] - obj_center[1];
    motion_vector[2] = view_point[2] - obj_center[2];
    vtkMath::Normalize(motion_vector);
  }

  double disp_obj_center[3];

  this->ComputeWorldToDisplay(obj_center[0], obj_center[1], obj_center[2],
    disp_obj_center);

  double newAngle =
    vtkMath::DegreesFromRadians( atan2( rwi->GetEventPosition()[1] - disp_obj_center[1],
    rwi->GetEventPosition()[0] - disp_obj_center[0] ) );

  double oldAngle =
    vtkMath::DegreesFromRadians( atan2( rwi->GetLastEventPosition()[1] - disp_obj_center[1],
    rwi->GetLastEventPosition()[0] - disp_obj_center[0] ) );

  double scale[3];
  scale[0] = scale[1] = scale[2] = 1.0;

  double **rotate = new double*[1];
  rotate[0] = new double[4];

  rotate[0][0] = newAngle - oldAngle;
  rotate[0][1] = motion_vector[0];
  rotate[0][2] = motion_vector[1];
  rotate[0][3] = motion_vector[2];

  this->Prop3DTransform( this->InteractionProp,
    obj_center,
    1,
    rotate,
    scale );

  delete [] rotate[0];
  delete [] rotate;

  if ( this->AutoAdjustCameraClippingRange )
  {
    this->CurrentRenderer->ResetCameraClippingRange();
  }

  rwi->Render();
}

//----------------------------------------------------------------------------
void FullRegionMoveStyle::Pan()
{
  if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
  {
    return;
  }

  vtkRenderWindowInteractor *rwi = this->Interactor;

  // Use initial center as the origin from which to pan

  double *obj_center = this->InteractionProp->GetCenter();

  double disp_obj_center[3], new_pick_point[4];
  double old_pick_point[4], motion_vector[3];

  this->ComputeWorldToDisplay(obj_center[0], obj_center[1], obj_center[2],
    disp_obj_center);

  this->ComputeDisplayToWorld(rwi->GetEventPosition()[0],
    rwi->GetEventPosition()[1],
    disp_obj_center[2],
    new_pick_point);

  this->ComputeDisplayToWorld(rwi->GetLastEventPosition()[0],
    rwi->GetLastEventPosition()[1],
    disp_obj_center[2],
    old_pick_point);

  motion_vector[0] = new_pick_point[0] - old_pick_point[0];
  motion_vector[1] = new_pick_point[1] - old_pick_point[1];
  motion_vector[2] = new_pick_point[2] - old_pick_point[2];

  if (this->InteractionProp->GetUserMatrix() != NULL)
  {
    vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
    t->Identity();
    t->PostMultiply();
    t->Concatenate(this->InteractionProp->GetUserTransform());
    t->Translate(motion_vector[0], motion_vector[1], motion_vector[2]);
    this->InteractionProp->SetUserTransform(t);
  }
  else
  {
    this->InteractionProp->AddPosition(motion_vector[0],
      motion_vector[1],
      motion_vector[2]);
  }

  if (this->AutoAdjustCameraClippingRange)
  {
    this->CurrentRenderer->ResetCameraClippingRange();
  }

  rwi->Render();
}

//----------------------------------------------------------------------------
void FullRegionMoveStyle::Dolly()
{
  if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
  {
    return;
  }

  vtkRenderWindowInteractor *rwi = this->Interactor;
  vtkCamera *cam = this->CurrentRenderer->GetActiveCamera();

  double view_point[3], view_focus[3];
  double motion_vector[3];

  cam->GetPosition(view_point);
  cam->GetFocalPoint(view_focus);

  double *center = this->CurrentRenderer->GetCenter();

  int dy = rwi->GetEventPosition()[1] - rwi->GetLastEventPosition()[1];
  double yf = dy / center[1] * this->MotionFactor;
  double dollyFactor = pow(1.1, yf);

  dollyFactor -= 1.0;
  motion_vector[0] = (view_point[0] - view_focus[0]) * dollyFactor;
  motion_vector[1] = (view_point[1] - view_focus[1]) * dollyFactor;
  motion_vector[2] = (view_point[2] - view_focus[2]) * dollyFactor;

  if (this->InteractionProp->GetUserMatrix() != NULL)
  {
    vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
    t->Identity();
    t->PostMultiply();
    t->Concatenate(this->InteractionProp->GetUserTransform());
    t->Translate(motion_vector[0], motion_vector[1], motion_vector[2]);
    this->InteractionProp->SetUserTransform(t);
  }
  else
  {
    this->InteractionProp->AddPosition(motion_vector);
  }

  if (this->AutoAdjustCameraClippingRange)
  {
    this->CurrentRenderer->ResetCameraClippingRange();
  }

  rwi->Render();
}

//----------------------------------------------------------------------------
void FullRegionMoveStyle::UniformScale()
{
  if (this->CurrentRenderer == NULL || this->InteractionProp == NULL)
  {
    return;
  }

  vtkRenderWindowInteractor *rwi = this->Interactor;

  int dy = rwi->GetEventPosition()[1] - rwi->GetLastEventPosition()[1];

  double *obj_center = this->InteractionProp->GetCenter();
  double *center = this->CurrentRenderer->GetCenter();

  double yf = dy / center[1] * this->MotionFactor;
  double scaleFactor = pow(1.1, yf);

  double **rotate = NULL;

  double scale[3];
  scale[0] = scale[1] = scale[2] = scaleFactor;

  this->Prop3DTransform(this->InteractionProp,
    obj_center,
    0,
    rotate,
    scale);

  if (this->AutoAdjustCameraClippingRange)
  {
    this->CurrentRenderer->ResetCameraClippingRange();
  }

  rwi->Render();
}


void FullRegionMoveStyle::FindPickedActor(int x, int y)
{
  // do nothing
}

void FullRegionMoveStyle::SetInteractionProp(vtkProp3D* InteractionProp)
{
  this->InteractionProp = InteractionProp;
}

void FullRegionMoveStyle::Prop3DTransform(vtkProp3D *prop3D,
                                                       double *boxCenter,
                                                       int numRotation,
                                                       double **rotate,
                                                       double *scale)
{
  vtkMatrix4x4 *oldMatrix = vtkMatrix4x4::New();
  prop3D->GetMatrix(oldMatrix);

  double orig[3];
  prop3D->GetOrigin(orig);

  vtkSmartPointer<vtkTransform> newTransform = vtkSmartPointer<vtkTransform>::New();
  newTransform->PostMultiply();
  if (prop3D->GetUserMatrix() != NULL)
  {
    newTransform->Concatenate(prop3D->GetUserTransform());
  }
  else
  {
    newTransform->SetMatrix(oldMatrix);
  }

  newTransform->Translate(-(boxCenter[0]), -(boxCenter[1]), -(boxCenter[2]));

  for (int i = 0; i < numRotation; i++)
  {
    newTransform->RotateWXYZ(rotate[i][0], rotate[i][1],
      rotate[i][2], rotate[i][3]);
  }

  if ((scale[0] * scale[1] * scale[2]) != 0.0)
  {
    newTransform->Scale(scale[0], scale[1], scale[2]);
  }

  newTransform->Translate(boxCenter[0], boxCenter[1], boxCenter[2]);

  // now try to get the composit of translate, rotate, and scale
  newTransform->Translate(-(orig[0]), -(orig[1]), -(orig[2]));
  newTransform->PreMultiply();
  newTransform->Translate(orig[0], orig[1], orig[2]);

  if (prop3D->GetUserMatrix() != NULL)
  {
    prop3D->SetUserTransform(newTransform);
  }
  else
  {
    prop3D->SetPosition(newTransform->GetPosition());
    prop3D->SetScale(newTransform->GetScale());
    prop3D->SetOrientation(newTransform->GetOrientation());
  }
  oldMatrix->Delete();
}

int FullRegionMoveStyle::GetEventState()
{
  return this->State;
}
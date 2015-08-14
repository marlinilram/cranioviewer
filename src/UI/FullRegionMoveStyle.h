#ifndef FullRegionMoveStyle_H
#define FullRegionMoveStyle_H

#include <vtkInteractorStyleTrackballActor.h>

class FullRegionMoveStyle : public vtkInteractorStyleTrackballActor
{
public:
  static FullRegionMoveStyle* New();
  vtkTypeMacro(FullRegionMoveStyle, vtkInteractorStyleTrackballActor);

  // Description:
  // Event bindings controlling the effects of pressing mouse buttons
  // or moving the mouse.
  void OnMouseMove();
  void OnLeftButtonDown();
  void OnLeftButtonUp();
  void OnMiddleButtonDown();
  void OnMiddleButtonUp();
  void OnRightButtonDown();
  void OnRightButtonUp();

  // These methods for the different interactions in different modes
  // are overridden in subclasses to perform the correct motion. Since
  // they might be called from OnTimer, they do not have mouse coord parameters
  // (use interactor's GetEventPosition and GetLastEventPosition)
  void Rotate();
  void Spin();
  void Pan();
  void Dolly();
  void UniformScale();

  void SetInteractionProp(vtkProp3D* InteractionProp);

  int GetEventState();

protected:
  FullRegionMoveStyle();
  ~FullRegionMoveStyle();

  void FindPickedActor(int x, int y);

  void Prop3DTransform(vtkProp3D *prop3D,
    double *boxCenter,
    int NumRotation,
    double **rotate,
    double *scale);

private:
  FullRegionMoveStyle(const FullRegionMoveStyle&);
  void operator=(const FullRegionMoveStyle&);
};

#endif
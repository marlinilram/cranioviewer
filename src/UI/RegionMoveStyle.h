#ifndef RegionMoveStyle_H
#define RegionMoveStyle_H

#include <vtkInteractorStyle.h>

class RegionMoveStyle : public vtkInteractorStyle
{
public:
  static RegionMoveStyle* New();
  vtkTypeMacro(RegionMoveStyle, vtkInteractorStyle);

  virtual void OnMouseMove();
  virtual void OnLeftButtonDown();
  virtual void OnLeftButtonUp();

  void SetStartWPosition(double StartWPosition[3]);
  double* GetLastMotionVec();
  double* GetMotionVec();

  int GetMovingTag();

protected:
  RegionMoveStyle();
  ~RegionMoveStyle();

  int StartPosition[2];
  int EndPosition[2];
  int Moving;
  double MotionVec[3];
  double LastMotionVec[3];
  double StartDPosition[3]; // in display
  double StartWPosition[3];

private:
  RegionMoveStyle(const RegionMoveStyle&);
  void operator=(const RegionMoveStyle&);
};

#endif
#ifndef SelectModeStyle_H
#define SelectModeStyle_H

#include <vtkInteractorStyleDrawPolygon.h>

class SelectModeStyle : public vtkInteractorStyleDrawPolygon
{
public:
  static SelectModeStyle* New();
  vtkTypeMacro(SelectModeStyle, vtkInteractorStyleDrawPolygon);

  virtual void OnMouseMove();
  virtual void OnLeftButtonDown();
  virtual void OnLeftButtonUp();

  void GetStartPos(int StartPos[2]);
  void GetEndPos(int EndPos[2]);
  double GetRadius();
  int GetMovingTag();
  void SetStartPtId(vtkIdType id);
  vtkIdType GetStartPtId();

protected:
  SelectModeStyle();
  ~SelectModeStyle();

  void DrawPolygon();
  void DrawPixels(int StartPos[2],
    int EndPos[2], unsigned char *pixels, int *size);

  vtkIdType StartPtId;

private:
  SelectModeStyle(const SelectModeStyle&);
  void operator=(const SelectModeStyle&);
};

#endif